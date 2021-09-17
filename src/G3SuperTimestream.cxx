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
struct G3SuperTimestream::flac_block encode_flac(PyArrayObject *array)
{
	struct G3SuperTimestream::flac_block fb;
	int n = PyArray_NBYTES(array);

	fb.buf = new char[n];
	fb.size = n;
	fb.count = 0;
	fb.offsets.push_back(0);

	int n_samps = PyArray_DIMS(array)[1];
	int32_t d[n_samps];

	char *src = (char*)(PyArray_DATA(array));
	const int32_t *chan_ptrs[1];
	FLAC__StreamEncoder *encoder = FLAC__stream_encoder_new();

	int32_t M = (1 << 24);
	for (int i=0; i< PyArray_DIMS(array)[0]; i++, src+=PyArray_STRIDES(array)[0]) {
		FLAC__stream_encoder_set_channels(encoder, 1);
		FLAC__stream_encoder_set_bits_per_sample(encoder, 24);
		FLAC__stream_encoder_set_compression_level(encoder, 1);
		FLAC__stream_encoder_init_stream(
			encoder, flac_encoder_write_cb, NULL, NULL, NULL, (void*)(&fb));
		chan_ptrs[0] = (int32_t*)d;
		int32_t pivot = round(double(((int32_t*)src)[0]) / M) * M;
		int warnings = 0;
		for (int i=0; i < n_samps; i++) {
			d[i] = ((int32_t*)src)[i] + pivot;
			if (d[i] >= M/2 || d[i] < -M/2)
				warnings++;
		}
		FLAC__stream_encoder_process(encoder, chan_ptrs, n_samps);
		FLAC__stream_encoder_finish(encoder);
		fb.offsets.push_back(fb.count);
		fb.pivots.push_back(pivot);
		fb.warnings.push_back(warnings);
	}
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
	ar & make_nvp("times", times);
	ar & make_nvp("names", names);

	// Read the desc.
        ar & make_nvp("type_num", desc.type_num);
        ar & make_nvp("ndim", desc.ndim);
	ar & make_nvp("shape", desc.shape);
        ar & make_nvp("item_size", desc.item_size);
        ar & make_nvp("nbytes", desc.nbytes);

	// Read the flacblock
	flac = new struct flac_block;
	ar & make_nvp("precision", flac->precision);
	ar & make_nvp("payload_bytes", flac->count);
	ar & make_nvp("offsets", flac->offsets);
	ar & make_nvp("pivots", flac->pivots);
	flac->buf = new char[flac->count];
	ar & make_nvp("payload", binary_data(flac->buf, flac->count));
}

template <class A> void G3SuperTimestream::save(A &ar, unsigned v) const
{
	using namespace cereal;
	ar & make_nvp("parent", base_class<G3FrameObject>(this));
	ar & make_nvp("times", times);
	ar & make_nvp("names", names);

	struct flac_block *_flac = flac;
	if (_flac == nullptr) {
		// Encode to a copy.
		_flac = new struct flac_block;
		*_flac = encode_flac(array);
        }

	// Write the desc.
        ar & make_nvp("type_num", desc.type_num);
        ar & make_nvp("ndim", desc.ndim);
	ar & make_nvp("shape", desc.shape);
        ar & make_nvp("item_size", desc.item_size);
        ar & make_nvp("nbytes", desc.nbytes);

	// Write the flacblock
	//ar & make_nvp("codec", 1);
	ar & make_nvp("precision", _flac->precision);
	ar & make_nvp("payload_bytes", _flac->count);
	ar & make_nvp("offsets", _flac->offsets);
	ar & make_nvp("pivots", _flac->pivots);
	ar & make_nvp("payload", binary_data(_flac->buf, _flac->count));
}

bool G3SuperTimestream::Encode(int codec) {
	if (array == nullptr)
		return false;

	// Compress the array data.
	flac = new struct flac_block;
	*flac = encode_flac(array);

	Py_XDECREF(array);
	array = nullptr;

	return true;
}

struct flac_helper {
	int bytes_remaining;
	char *src;
	int32_t *dest;
	int32_t pivot;
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

FLAC__StreamDecoderWriteStatus write_callback(
	const FLAC__StreamDecoder *decoder, const FLAC__Frame *frame, const FLAC__int32 *const buffer[], void *client_data)
{
	auto fh = (struct flac_helper *)client_data;
	int n = frame->header.blocksize;
	for (int i=0; i<n; i++)
		*(fh->dest++) = buffer[0][i] + fh->pivot;
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

bool G3SuperTimestream::Decode() {
	/* printf("decode %p %p\n", array, flac); */

	if (flac == nullptr)
		return false;

	array = (PyArrayObject*)
		PyArray_SimpleNew(desc.ndim, desc.shape, desc.type_num);

	auto decoder = FLAC__stream_decoder_new();

	for (int i=0; i<desc.shape[0]; i++) {
		struct flac_helper helper;
		helper.src = flac->buf + flac->offsets[i];
		helper.bytes_remaining = flac->offsets[i+1] - flac->offsets[i];
		helper.dest = (int32_t*)PyArray_DATA(array) + desc.shape[1] * i;
		helper.pivot = flac->pivots[i];

		FLAC__stream_decoder_init_stream(
			decoder, read_callback, NULL, NULL, NULL, NULL,
			write_callback, NULL, flac_decoder_error_cb, (void*)&helper);

		FLAC__stream_decoder_process_until_end_of_stream(decoder);
		FLAC__stream_decoder_finish(decoder);
	}

        FLAC__stream_decoder_delete(decoder);

	// Kill the flac bundle.
	delete flac->buf;
	delete flac;
	flac = nullptr;

	return true;
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

	if (self.flac) {
		delete self.flac->buf;
		delete self.flac;
		self.flac = nullptr;
	}

	self.desc.ndim = PyArray_NDIM(_array);
	self.desc.type_num = PyArray_TYPE(_array);

	self.desc.item_size = PyArray_ITEMSIZE(_array);
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
		;
	register_pointer_conversions<G3SuperTimestream>();

	bp::register_exception_translator<g3supertimestream_exception>(&translate_ValueError);
}
