#define NO_IMPORT_ARRAY

#include <pybindings.h>

/* #include <numeric> */
/* #include <algorithm> */
#include <iostream>
#include <boost/python.hpp>

#include <container_pybindings.h>
#include <G3SuperTimestream.h>
#include <cereal/types/utility.hpp>

#include <FLAC/stream_encoder.h>
#include <bzlib.h>

// Debugging variables for compressors.
// Must range from 0 (silent) to 4.
#define SO3G_BZ2_VERBOSITY 0

static
std::string get_bz2_error_string(int err) {
	std::ostringstream s;
	switch(err) {
	case BZ_CONFIG_ERROR:
		s << "BZ_CONFIG_ERROR (library compilation issue)";
		break;
	case BZ_PARAM_ERROR:
		s <<  "BZ_PARAM_ERROR (bad blocksize, verbosity, etc)";
		break;
	case BZ_MEM_ERROR:
		s << "BZ_MEM_ERROR (not enough memory is available)";
		break;
	case BZ_OK:
		s << "BZ_OK (no problem)";
	default:
		s << "Unknown error code " << err;
	}
	return s.str();
}

static
std::string get_algo_error_string(std::string var_name, int algo_code)
{
	std::ostringstream s;
	s << "No support for compression algorithm " << var_name << "=" << algo_code;
	return s.str();
}

// Split each datum in x into:
//
//    x[i] = r[i] + y[i]
//
// where r is a multiple of snap_to, -step_at <= y[i] < step_at.  The
// algorithm tries to have r[i] change slowly, meaning that r[i+1]
// will not differ from r[i] unless necessary to satisfy the
// range restriction on y[i].
//
// This is only implemented for snap_to and step_at both powers of 2!
//
// Let's allow r and x to point to same memory.
template <typename T>
static
int rebranch(int32_t *y, T *r, T *x, int n_samps, T snap_to, T step_at)
{
	T branch = 0;
	int fails = 0;
	for (int i=0; i < n_samps; i++) {
		T new_branch = (x[i] + snap_to / 2) & ~(snap_to - 1);
		if (new_branch - branch >= step_at ||
		    new_branch - branch <  -step_at)
			branch = new_branch;
		y[i] = x[i] - branch;
		// Don't update r until you use x -- they are allowed to overlap.
		r[i] = branch;
		if (y[i] < -step_at || y[i] >= step_at)
			fails++;
	}
	return fails;
}

FLAC__StreamEncoderWriteStatus flac_encoder_write_cb(
	const FLAC__StreamEncoder *encoder,
	const FLAC__byte buffer[],
	size_t bytes,
	unsigned samples,
	unsigned current_frame,
	void *client_data)
{
	auto fb = (struct G3SuperTimestream::flac_block *)client_data;
	assert(bytes + fb->count <= fb->size);
	memcpy(fb->buf + fb->count, buffer, bytes);
	fb->count += bytes;
	return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

static
struct G3SuperTimestream::flac_block encode_flac_bz2(
	int8_t data_algo, PyArrayObject *array, float precision)
{
	struct G3SuperTimestream::flac_block fb;
	int n = PyArray_NBYTES(array);

	fb.buf = new char[n];
	fb.size = n;
	fb.count = 0;
	fb.offsets.push_back(0);
	fb.precision = precision;

	int n_samps = PyArray_DIMS(array)[1];
	int32_t d[n_samps];
	const int32_t *chan_ptrs[1] = {d};

	char r[n_samps * sizeof(int64_t)];
	auto r32 = (int32_t*)r;
	auto r64 = (int64_t*)r;

	npy_intp type_num = PyArray_TYPE(array);
	char *src = (char*)(PyArray_DATA(array));
	FLAC__StreamEncoder *encoder = nullptr;
	if (data_algo & G3SuperTimestream::ALGO_DO_FLAC)
		encoder = FLAC__stream_encoder_new();

	int32_t M = (1 << 24);
	for (int i=0; i< PyArray_DIMS(array)[0]; i++, src+=PyArray_STRIDES(array)[0]) {
		if (encoder != nullptr) {
			int fails = 0;
			if (type_num == NPY_INT32) {
				fails = rebranch<int32_t>(d, r32, (int32_t*)src, n_samps, M, M/2);
			} else if (type_num == NPY_INT64) {
				fails = rebranch<int64_t>(d, r64, (int64_t*)src, n_samps, M, M/2);
			} else if (type_num == NPY_FLOAT32) {
				for (int j=0; j < n_samps; j++)
					r32[j] = roundf(((float*)src)[j] / precision);
				fails = rebranch<int32_t>(d, r32, r32, n_samps, M, M/2);
			} else if (type_num == NPY_FLOAT64) {
				for (int j=0; j < n_samps; j++)
					r64[j] = round(((double*)src)[j] / precision);
				fails = rebranch<int64_t>(d, r64, r64, n_samps, M, M/2);
			} else
				throw g3supertimestream_exception("Invalid array type encountered.");
			if (fails > 0)
				throw g3supertimestream_exception("Data fail.");

			FLAC__stream_encoder_set_channels(encoder, 1);
			FLAC__stream_encoder_set_bits_per_sample(encoder, 24);
			FLAC__stream_encoder_set_compression_level(encoder, 1);
			FLAC__stream_encoder_init_stream(
				encoder, flac_encoder_write_cb, NULL, NULL, NULL, (void*)(&fb));
			FLAC__stream_encoder_process(encoder, chan_ptrs, n_samps);
			FLAC__stream_encoder_finish(encoder);
		} else {
			memcpy(r, src, n_samps * PyArray_ITEMSIZE(array));
			if (type_num == NPY_FLOAT32) {
				for (int j=0; j<n_samps; j++)
					((int32_t*)r)[j] = roundf(((float*)r)[j] / precision);
			} else if (type_num == NPY_FLOAT64) {
				for (int j=0; j<n_samps; j++)
					((int64_t*)r)[j] = round(((double*)r)[j] / precision);
			}
		}
		fb.offsets.push_back(fb.count);

		// And the bz2
		if (data_algo & G3SuperTimestream::ALGO_DO_BZ) {
			assert(fb.size > fb.count);
			unsigned int n_write = fb.size - fb.count;
			int err = BZ2_bzBuffToBuffCompress(
				fb.buf + fb.count, &n_write, r, n_samps * PyArray_ITEMSIZE(array),
				5, SO3G_BZ2_VERBOSITY, 1);
			if (err != BZ_OK)
				throw g3supertimestream_exception(get_bz2_error_string(err));
			fb.count += n_write;
		}
		fb.offsets.push_back(fb.count);
	}
	if (encoder != nullptr)
		FLAC__stream_encoder_delete(encoder);

	return fb;
}


/* G3SuperTimestream */

std::string G3SuperTimestream::Description() const
{
	std::ostringstream s;
	s << "G3SuperTimestream("
	  << names.size() << ", " << times.size() << ")";
	return s.str();
}

std::string G3SuperTimestream::Summary() const
{
	return Description();
}

template <class A> void G3SuperTimestream::load(A &ar, unsigned v)
{
	using namespace cereal;
	ar & make_nvp("parent", base_class<G3FrameObject>(this));

	ar & make_nvp("times_algo", options.times_algo);

	if (options.times_algo == ALGO_NONE) {
		ar & make_nvp("times", times);
	} else if (options.times_algo == ALGO_DO_BZ) {
		int n_samps;
		ar & make_nvp("n_samps", n_samps);
		times.resize(n_samps);

		int max_bytes = 0;
		ar & make_nvp("comp_bytes", max_bytes);

		char *buf = (char*)malloc(max_bytes);
		ar & make_nvp("times_data", binary_data(buf, max_bytes));
		unsigned int n_decomp = n_samps * sizeof(times[0]);
		int err = BZ2_bzBuffToBuffDecompress(
			     (char*)&times[0], &n_decomp, buf, max_bytes,
			     1, SO3G_BZ2_VERBOSITY);
		if (err != BZ_OK)
			throw g3supertimestream_exception(get_bz2_error_string(err));
		free(buf);
	} else
		throw g3supertimestream_exception(
			get_algo_error_string("times_algo", options.times_algo));

	ar & make_nvp("names", names);

	// Read the desc.
        ar & make_nvp("type_num", desc.type_num);
        ar & make_nvp("ndim", desc.ndim);
	ar & make_nvp("shape", desc.shape);
        ar & make_nvp("nbytes", desc.nbytes);

	ar & make_nvp("data_algo", options.data_algo);
	if (options.data_algo == ALGO_NONE) {
		assert(PyArray_EquivByteorders(NPY_NATIVE, NPY_LITTLE)); // Check it's not 1997
		array = (PyArrayObject*)
			PyArray_SimpleNew(desc.ndim, desc.shape, desc.type_num);
		ar & make_nvp("data_raw", binary_data((char*)PyArray_DATA(array),
						      PyArray_NBYTES(array)));
	} else {
		// Read the flacblock
		flac = new struct flac_block;
		ar & make_nvp("precision", flac->precision);
		ar & make_nvp("offsets", flac->offsets);
		ar & make_nvp("payload_bytes", flac->count);
		flac->buf = new char[flac->count];
		ar & make_nvp("payload", binary_data(flac->buf, flac->count));
	}
}

template <class A> void G3SuperTimestream::save(A &ar, unsigned v) const
{
	using namespace cereal;
	ar & make_nvp("parent", base_class<G3FrameObject>(this));

	ar & make_nvp("times_algo", options.times_algo);

	if (options.times_algo == ALGO_NONE) {
		// No compression.
		ar & make_nvp("times", times);
	} else if (options.times_algo == ALGO_DO_BZ) {
		int n_samps = times.size();
		ar & make_nvp("n_samps", n_samps);

		unsigned int max_bytes = n_samps * sizeof(times[0]);
		char *buf = (char*)malloc(max_bytes);
		int err = BZ2_bzBuffToBuffCompress(
			buf, &max_bytes, (char*)&times[0], max_bytes,
			     5, SO3G_BZ2_VERBOSITY, 1);
		if (err != BZ_OK)
			throw g3supertimestream_exception(get_bz2_error_string(err));
		ar & make_nvp("comp_bytes", max_bytes);
		ar & make_nvp("times_data", binary_data(buf, max_bytes));
		free(buf);
	} else
		throw g3supertimestream_exception(
			get_algo_error_string("times_algo", options.times_algo));

	ar & make_nvp("names", names);

	// Write the desc.
        ar & make_nvp("type_num", desc.type_num);
        ar & make_nvp("ndim", desc.ndim);
	ar & make_nvp("shape", desc.shape);
        ar & make_nvp("nbytes", desc.nbytes);

	ar & make_nvp("data_algo", options.data_algo);
	if (options.data_algo == ALGO_NONE) {
		assert(array != nullptr);
		// Check the endianness
		auto this_descr = PyArray_DESCR(array);
		assert(PyArray_EquivByteorders(this_descr.byte_order, NPY_LITTLE));

		// Might as well use numpy to repack it properly...
		PyArrayObject *contig = PyArray_GETCONTIGUOUS(array);
		ar & make_nvp("data_raw", binary_data((char*)PyArray_DATA(contig),
						      PyArray_NBYTES(contig)));
		Py_DECREF((PyObject*)contig);
	} else {
		struct flac_block *_flac = flac;
		if (_flac == nullptr) {
			// Encode to a copy.
			_flac = new struct flac_block;
			*_flac = encode_flac_bz2(options.data_algo, array, options.precision);
		}

		// Write the flacblock
		ar & make_nvp("precision", _flac->precision);
		ar & make_nvp("offsets", _flac->offsets);
		ar & make_nvp("payload_bytes", _flac->count);
		ar & make_nvp("payload", binary_data(_flac->buf, _flac->count));

		if (_flac != flac) {
			delete _flac->buf;
			delete _flac;
		}
	}
}

bool G3SuperTimestream::Encode() {
	if (array == nullptr)
		return false;

	// Compress the array data.
	if (options.data_algo == ALGO_NONE)
		return false;
	else {
		flac = new struct flac_block;
		*flac = encode_flac_bz2(options.data_algo, array, options.precision);
		Py_XDECREF(array);
		array = nullptr;
	}
	return true;
}

struct flac_helper {
	int bytes_remaining;
	char *src;
	char *dest;
	float precision;
};

FLAC__StreamDecoderReadStatus read_callback(
	const FLAC__StreamDecoder *decoder, FLAC__byte buffer[], size_t *bytes, void *client_data)
{
	auto fh = (struct flac_helper *)client_data;
	/* printf(" ... read %i (remaining: %i)\n", *bytes, fh->bytes_remaining); */
	if (fh->bytes_remaining == 0) {
		*bytes = 0;
		return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM;
	}
	if (fh->bytes_remaining < *bytes)
		*bytes = fh->bytes_remaining;
	memcpy(buffer, fh->src, *bytes);
	fh->bytes_remaining -= *bytes;
	fh->src += *bytes;
	return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
}

template <typename T>
FLAC__StreamDecoderWriteStatus write_callback_int(
	const FLAC__StreamDecoder *decoder, const FLAC__Frame *frame, const FLAC__int32 *const buffer[], void *client_data)
{
	auto fh = (struct flac_helper *)client_data;
	int n = frame->header.blocksize;
	for (int i=0; i<n; i++)
		((T*)fh->dest)[i] = buffer[0][i];
	fh->dest += n * sizeof(T);
	return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}

static void flac_decoder_error_cb(const FLAC__StreamDecoder *decoder,
				  FLAC__StreamDecoderErrorStatus status, void *client_data)
{

	switch (status) {
	case FLAC__STREAM_DECODER_ERROR_STATUS_LOST_SYNC:
		printf("FLAC decoding error (lost sync)");
	case FLAC__STREAM_DECODER_ERROR_STATUS_BAD_HEADER:
		printf("FLAC decoding error (bad header)");
	case FLAC__STREAM_DECODER_ERROR_STATUS_FRAME_CRC_MISMATCH:
		printf("FLAC decoding error (CRC mismatch)");
	case FLAC__STREAM_DECODER_ERROR_STATUS_UNPARSEABLE_STREAM:
		printf("FLAC decoding error (unparseable stream)");
	default:
		printf("FLAC decoding error (%d)", status);
	}
}

template <typename T>
void expand_branch(struct flac_helper *fh, int n_bytes, int nsamps,
	char *temp)
{
	bool own_temp = (temp == nullptr);
	unsigned int temp_size = nsamps * sizeof(T);
	if (own_temp)
		temp = new char[temp_size];

	int err = BZ2_bzBuffToBuffDecompress(
		temp, &temp_size, fh->src, n_bytes,
		1, SO3G_BZ2_VERBOSITY);
	if (err != BZ_OK)
		throw g3supertimestream_exception(get_bz2_error_string(err));
	// Add it in ...
	for (int i=0; i<nsamps; i++)
		((T*)fh->dest)[i] += ((T*)temp)[i];

	if (own_temp)
		delete temp;
}

bool G3SuperTimestream::Decode()
{
	if (flac == nullptr)
		return false;

	assert (options.data_algo != 0);

	FLAC__StreamDecoder *decoder = nullptr;
	if (options.data_algo & ALGO_DO_FLAC)
		decoder = FLAC__stream_decoder_new();
	array = (PyArrayObject*)
		PyArray_ZEROS(desc.ndim, desc.shape, desc.type_num, 0);

	struct flac_helper helper;
	helper.precision = flac->precision;

	FLAC__StreamDecoderWriteCallback this_write_callback;
	void (*expand_func)(struct flac_helper *, int, int, char*) = nullptr;
	int elsize = 0;

	switch (desc.type_num) {
	case NPY_INT32:
	case NPY_FLOAT32:
		this_write_callback = &write_callback_int<int32_t>;
		expand_func = expand_branch<int32_t>;
		elsize = sizeof(int32_t);
		break;
	case NPY_INT64:
	case NPY_FLOAT64:
		this_write_callback = &write_callback_int<int64_t>;
		expand_func = expand_branch<int64_t>;
		elsize = sizeof(int64_t);
		break;
	default:
		throw g3supertimestream_exception("Invalid array type encountered.");
	}

	char *temp = new char[PyArray_SHAPE(array)[1] * elsize];

	for (int i=0; i<desc.shape[0]; i++) {
		char* this_data = (char*)PyArray_DATA(array) + PyArray_STRIDES(array)[0]*i;
		if (decoder != nullptr) {
			helper.dest = this_data;
			helper.src = flac->buf + flac->offsets[2*i];
			helper.bytes_remaining = flac->offsets[2*i+1] - flac->offsets[2*i];

			FLAC__stream_decoder_init_stream(
				decoder, read_callback, NULL, NULL, NULL, NULL,
				*this_write_callback, NULL, flac_decoder_error_cb,
				(void*)&helper);

			FLAC__stream_decoder_process_until_end_of_stream(decoder);
			FLAC__stream_decoder_finish(decoder);
		}

		// And bz2 bit
		if (options.data_algo & ALGO_DO_BZ) {
			helper.src = flac->buf + flac->offsets[2*i+1];
			helper.dest = this_data;
			expand_func(&helper, flac->offsets[2*i+2] - flac->offsets[2*i+1],
				    PyArray_SHAPE(array)[1], temp);
		}

		// Now convert for precision.
		if (desc.type_num == NPY_FLOAT32) {
			auto src = (int32_t*)this_data;
			auto dest = (float*)this_data;
			for (int j=0; j<PyArray_SHAPE(array)[1]; j++)
				dest[j] = (float)src[j] * flac->precision;
		} else if (desc.type_num == NPY_FLOAT64) {
			auto src = (int64_t*)this_data;
			auto dest = (double*)this_data;
			for (int j=0; j<PyArray_SHAPE(array)[1]; j++)
				dest[j] = src[j] * flac->precision;
		}
	}
	delete temp;
	if (decoder != nullptr)
		FLAC__stream_decoder_delete(decoder);

	// Destroy the flac bundle.
	delete flac->buf;
	delete flac;
	flac = nullptr;

	return true;
}

int G3SuperTimestream::Options(int data_algo, int times_algo, float precision)
{
	if (data_algo >= 0)
		options.data_algo = data_algo;
	if (times_algo >= 0)
		options.times_algo = times_algo;
	if (precision >= 0)
		options.precision = precision;
	return 0;
}


G3SuperTimestream::G3SuperTimestream() {
	options.times_algo = ALGO_DO_BZ;
	options.data_algo = ALGO_DO_FLAC | ALGO_DO_BZ;
	array = nullptr;
	flac = nullptr;
}

G3SuperTimestream::~G3SuperTimestream()
{
	if (array != nullptr)
		Py_XDECREF(array);
	if (flac != nullptr) {
		delete flac->buf;
		delete flac;
	}
}

static
void safe_set_times(G3SuperTimestream &self, G3VectorTime _times)
{
	// Only allow this if it doesn't upset consistency.  We will
	// assume that, coming in, we're internally consistent.
	if (_times.size() != self.times.size() && self.times.size() != 0) {
		std::ostringstream s;
		s << "Cannot set .times because it conflicts with "
		  << "the established number of samples (" << self.times.size()
		  << ").";
		throw g3supertimestream_exception(s.str());
	}
	self.times = _times;
}

static
void safe_set_names(G3SuperTimestream &self, G3VectorString _names)
{
	// Only allow this if it doesn't upset consistency.  We will
	// assume that, coming in, we're internally consistent.
	if (_names.size() != self.names.size() && self.names.size() != 0) {
		std::ostringstream s;
		s << "Cannot set .names because it conflicts with "
		  << "the established number of channels (" << self.names.size()
		  << ").";
		throw g3supertimestream_exception(s.str());
	}
	self.names = _names;
}

static
void safe_set_data(G3SuperTimestream &self, const bp::object object_in)
{
	// Note this function, as invoked here, might return a
	// reference or create a new array.
	PyObject *ob = PyArray_FromAny(object_in.ptr(), NULL, 0, 0, 0, NULL);
	if (ob == NULL)
		throw g3supertimestream_exception("Could not decode array.");

	PyArrayObject *_array = reinterpret_cast<PyArrayObject*>(ob);

	if (PyArray_NDIM(_array) != 2) {
		Py_XDECREF(ob);
		throw g3supertimestream_exception("Bad ndim.");
	}
	if (PyArray_DIMS(_array)[0] != self.names.size()) {
		Py_XDECREF(ob);
		throw g3supertimestream_exception("Bad shape[0].");
	}
	if (PyArray_DIMS(_array)[1] != self.times.size()) {
		Py_XDECREF(ob);
		throw g3supertimestream_exception("Bad shape[1].");
	}

	// Clear cached array or compressed data.
	if (self.array) {
		Py_XDECREF((PyObject*)self.array);
		self.array = nullptr;
	}
	if (self.flac) {
		delete self.flac->buf;
		delete self.flac;
		self.flac = nullptr;
	}

	self.desc.ndim = PyArray_NDIM(_array);
	self.desc.type_num = PyArray_TYPE(_array);

	self.desc.nbytes = PyArray_NBYTES(_array);
	for (int i=0; i<self.desc.ndim; i++)
		self.desc.shape[i] = PyArray_DIMS(_array)[i];

	self.array = _array;
}

static
bp::object safe_get_data(G3SuperTimestream &self)
{
	if (self.array == nullptr)
		self.Decode();
	return bp::object(bp::handle<>(bp::borrowed(reinterpret_cast<PyObject*>(self.array))));
}

static
bp::object safe_get_dtype(G3SuperTimestream &self)
{
	if (self.array == nullptr)
		return bp::object(); // not good enough...
	return bp::object(bp::handle<>(bp::borrowed(reinterpret_cast<PyObject*>(
							    PyArray_DESCR(self.array)->typeobj))));
}


G3_SPLIT_SERIALIZABLE_CODE(G3SuperTimestream);

static void translate_ValueError(g3supertimestream_exception const& e)
{
	PyErr_SetString(PyExc_ValueError, e.msg_for_python().c_str());
}


PYBINDINGS("so3g")
{
	EXPORT_FRAMEOBJECT(G3SuperTimestream, init<>(), "G3SuperTimestream()")
		.add_property("times", &G3SuperTimestream::times, &safe_set_times,
			      "Times vector.  Setting this stores a copy, but getting returns a reference.")
		.add_property("names", &G3SuperTimestream::names, &safe_set_names,
			      "Names vector.  Setting this stores a copy, but getting returns a reference.")
		.add_property("data", &safe_get_data, &safe_set_data,
			      "Data array.")
		.add_property("dtype", &safe_get_dtype, "Numpy dtype of underlying array.")
		.def("encode", &G3SuperTimestream::Encode, "Compress.")
		.def("decode", &G3SuperTimestream::Decode, "De-compress.")
		.def("options", &G3SuperTimestream::Options,
		     (bp::arg("data_algo")=-1, bp::arg("times_algo")=-1, bp::arg("precision")=-1.),
		     "Get/set compression options.")
		;
	register_pointer_conversions<G3SuperTimestream>();

	bp::register_exception_translator<g3supertimestream_exception>(&translate_ValueError);
}
