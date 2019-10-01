 #include <pybindings.h>
#include <serialization.h>
#include <G3VectorFixed.h>

#include <cereal/types/vector.hpp>

#include <FLAC/stream_encoder.h>
#include <cmath>

//
// FLAC stuff. This is all ripped out of G3Timestream, probably
// without any modification.
//

template<typename A>
struct FlacDecoderCallbackArgs {
    A *inbuf;
    std::vector<double> *outbuf;
    size_t pos;
    size_t nbytes;
};

enum FLACNaNFlag {
    NoNan = 0,
    AllNan = 1,
    SomeNan = 2
};

static FLAC__StreamEncoderWriteStatus flac_encoder_write_cb(
    const FLAC__StreamEncoder *encoder, const FLAC__byte buffer[], size_t bytes,
    unsigned samples, unsigned current_frame, void *client_data)
{
    std::vector<uint8_t> *outbuf = (std::vector<uint8_t> *)(client_data);

    outbuf->insert(outbuf->end(), buffer, buffer + bytes);
    return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

template<typename A>
static FLAC__StreamDecoderReadStatus flac_decoder_read_cb(
    const FLAC__StreamDecoder *decoder, FLAC__byte buffer[], size_t *bytes,
    void *client_data)
{
    FlacDecoderCallbackArgs<A> *args =
        (FlacDecoderCallbackArgs<A> *)(client_data);

    ssize_t bytes_left = ssize_t(args->nbytes) - args->pos;

    if (bytes_left <= 0 || *bytes == 0) {
        *bytes = 0;
        return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM;
    } else if (*bytes >= size_t(bytes_left)) {
        *bytes = bytes_left;
        args->inbuf->template loadBinary<1>(buffer, bytes_left);
        args->pos += bytes_left;
        return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
    } else {
        args->inbuf->template loadBinary<1>(buffer, *bytes);
        args->pos += *bytes;
        return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
    }
}
  
template<typename A>
static FLAC__StreamDecoderWriteStatus flac_decoder_write_cb(
    const FLAC__StreamDecoder *decoder, const FLAC__Frame *frame,
    const FLAC__int32 *const buffer[], void *client_data)
{
    FlacDecoderCallbackArgs<A> *args =
        (FlacDecoderCallbackArgs<A> *)(client_data);

    size_t oldsize = args->outbuf->size();
    args->outbuf->resize(oldsize + frame->header.blocksize);
    for (size_t i = 0; i < frame->header.blocksize; i++)
        (*args->outbuf)[oldsize + i] = buffer[0][i];
    return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}
  
static void flac_decoder_error_cb(const FLAC__StreamDecoder *decoder,
                                  FLAC__StreamDecoderErrorStatus status, void *client_data)
{

    switch (status) {
    case FLAC__STREAM_DECODER_ERROR_STATUS_LOST_SYNC:
        log_fatal("FLAC decoding error (lost sync)");
    case FLAC__STREAM_DECODER_ERROR_STATUS_BAD_HEADER:
        log_fatal("FLAC decoding error (bad header)");
    case FLAC__STREAM_DECODER_ERROR_STATUS_FRAME_CRC_MISMATCH:
        log_fatal("FLAC decoding error (CRC mismatch)");
    case FLAC__STREAM_DECODER_ERROR_STATUS_UNPARSEABLE_STREAM:
        log_fatal("FLAC decoding error (unparseable stream)");
    default:
        log_fatal("FLAC decoding error (%d)", status);
    }
}

//
// Serialization.  This is modified only slightly from
// G3Timestream... we need to store the precision, we do not need to
// store the timestamps.
//

template <class A> void G3VectorFixed::save(A &ar, unsigned v) const
{
    ar & cereal::make_nvp("G3FrameObject",
                          cereal::base_class<G3FrameObject>(this));
    ar & cereal::make_nvp("precision", precision);
    ar & cereal::make_nvp("flac_level", flac_level);
        
    if (flac_level == 0) {
        // No compression.
        ar & cereal::make_nvp("data",
                              cereal::base_class<std::vector<double> >(this));
    } else {
        std::vector<int32_t> inbuf;
        std::vector<uint8_t> outbuf;
        const int32_t *chanmap[1];
        uint8_t nanflag;
        size_t nans = 0;

        // Copy to 24-bit integers
        inbuf.resize(size());
        for (size_t i = 0; i < size(); i++)
            inbuf[i] = ((int32_t(round((*this)[i] / precision)) & 0x00ffffff) << 8)
                >> 8;
        chanmap[0] = &inbuf[0];

        // Mark bad samples using an out-of-band signal. Since we
        // convert to 24-bit integers going into FLAC, which have no
        // out-of-range values for signalling, this requires special
        // care. Usually a timestream is either all-valid or
        // all-invalid, so signal that with a single byte flag. In the
        // rare case that only some samples are valid, store a
        // validity mask.
        std::vector<bool> nanbuf(size(), false);
        for (size_t i = 0; i < size(); i++) {
            if (!std::isfinite((*this)[i])) {
                nans++;
                nanbuf[i] = true;
                inbuf[i] = (i >=0 ? inbuf[-1] : 0);
            }
        }
        nanflag = SomeNan;
        if (nans == 0)
            nanflag = NoNan;
        else if (nans == size())
            nanflag = AllNan;
        ar & cereal::make_nvp("nanflag", nanflag);
        if (nanflag == SomeNan)
            ar & cereal::make_nvp("nanmask", nanbuf);

        // Now do FLAC encoding
        FLAC__StreamEncoder *encoder = FLAC__stream_encoder_new();
        FLAC__stream_encoder_set_channels(encoder, 1);
        // XXX: should assert if high-order 8 bits are not clear
        FLAC__stream_encoder_set_bits_per_sample(encoder, 24);
        FLAC__stream_encoder_set_compression_level(encoder, flac_level);
        FLAC__stream_encoder_init_stream(
            encoder, flac_encoder_write_cb, NULL, NULL, NULL, (void*)(&outbuf));
        FLAC__stream_encoder_process (encoder, chanmap, inbuf.size());
        FLAC__stream_encoder_finish(encoder);
        FLAC__stream_encoder_delete(encoder);

        ar & cereal::make_nvp("data", outbuf);
    }
}

template <class A> void G3VectorFixed::load(A &ar, unsigned v)
{
    G3_CHECK_VERSION(v);

    ar & cereal::make_nvp("G3FrameObject",
                          cereal::base_class<G3FrameObject>(this));
    ar & cereal::make_nvp("precision", precision);
    ar & cereal::make_nvp("flac_level", flac_level);

    if (flac_level) {
        FlacDecoderCallbackArgs<A> callback;
        uint8_t nanflag;
        std::vector<bool> nanbuf;

        callback.inbuf = &ar;
        callback.outbuf = this;
        callback.pos = 0;

        // if (units != Counts)
        // 	log_fatal("Cannot use FLAC on non-counts timestreams");

        ar & cereal::make_nvp("nanflag", nanflag);
        if (nanflag == SomeNan)
            ar & cereal::make_nvp("nanmask", nanbuf);

        ar & cereal::make_size_tag(callback.nbytes);

        // Typical compression ratio: N bytes in input = N samples
        reserve(callback.nbytes);

        FLAC__StreamDecoder *decoder = FLAC__stream_decoder_new();
        FLAC__stream_decoder_init_stream(decoder,
                                         flac_decoder_read_cb<A>, NULL, NULL, NULL, NULL,
                                         flac_decoder_write_cb<A>, NULL, flac_decoder_error_cb,
                                         (void*)(&callback));
        FLAC__stream_decoder_process_until_end_of_stream(decoder);
        FLAC__stream_decoder_finish(decoder);
        FLAC__stream_decoder_delete(decoder);

        // Apply precision.
        for (size_t i=0; i < size(); i++)
            (*this)[i] *= precision;
                
        // Apply NaN mask
        if (nanflag == AllNan) {
            for (size_t i = 0; i < size(); i++)
                (*this)[i] = NAN;
        } else if (nanflag == SomeNan) {
            for (size_t i = 0; i < size(); i++)
                if (nanbuf[i])
                    (*this)[i] = NAN;
        }

    } else {
        ar & cereal::make_nvp("data",
                              cereal::base_class<std::vector<double> >(this));
    }
}

//
// Help the user with sanity checking on range and precision.
//

int G3VectorFixed::CheckPrecision()
{
    double tol = 1e-6;
    int fail_count = 0;
    for (size_t i=0; i<size(); i++) {
        double v = (*this)[i] / precision;
        if (fabs(v - round(v)) > tol)
            ++fail_count;
    }
    return fail_count;
}

int G3VectorFixed::CheckRange()
{
    const double max_i24 =  8388607;
    const double min_i24 = -8388608;
    
    int fail_count = 0;
    for (size_t i=0; i<size(); i++) {
        double v = round((*this)[i] / precision);
        if ((v < min_i24) || (v > max_i24))
            ++fail_count;
    }
    return fail_count;
}

std::string G3VectorFixed::Description() const
{
    std::ostringstream desc;
    desc << "G3VectorFixed(size=" << size() << ", precision=" << precision << ")";
    return desc.str();
}


G3_SPLIT_SERIALIZABLE_CODE(G3VectorFixed);

namespace {

    SET_LOGGER("G3VectorFixed");

    static int
    G3VectorFixed_getbuffer(PyObject *obj, Py_buffer *view, int flags)
    {
	if (view == NULL) {
            PyErr_SetString(PyExc_ValueError, "NULL view");
            return -1;
	}

	view->shape = NULL;

	boost::python::handle<> self(boost::python::borrowed(obj));
	boost::python::object selfobj(self);
	G3VectorFixedPtr ts = boost::python::extract<G3VectorFixedPtr>(selfobj)();
	view->obj = obj;
	view->buf = (void*)&(*ts)[0];
	view->len = ts->size() * sizeof(double);
	view->readonly = 0;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
            view->format = (char *)"d";
	else
            view->format = NULL;
	view->ndim = 1;

	// Abuse internal pointer in the absence of smalltable. This is safe
	// on all architectures except MIPS N32.
	view->internal = (void *)ts->size();
	view->shape = (Py_ssize_t *)(&view->internal);
	view->strides = &view->itemsize;
	view->suboffsets = NULL;

	// Try to hold onto our collective hats. This is still very dangerous if
	// the G3VectorFixed's underlying vector is resized.
	Py_INCREF(obj);

	return 0;
    }

    static int
    G3VectorFixed_nsamples(const G3VectorFixed &r)
    {
	return r.size();
    }

    static G3VectorFixedPtr
    G3VectorFixed_getslice(const G3VectorFixed &a, boost::python::slice slice)
    {
	using namespace boost::python;
	int start(0), stop(a.size()), step(1);

	// Normalize and check slice boundaries
	if (slice.start().ptr() != Py_None)
            start = extract<int>(slice.start())();
	if (slice.stop().ptr() != Py_None)
            stop = extract<int>(slice.stop())();
	if (slice.step().ptr() != Py_None)
            step = extract<int>(slice.step())();

	if (start < 0)
            start = a.size() + start;
	if (stop < 0)
            stop = a.size() + stop;

	if (start >= a.size() || start < 0)
            log_fatal("Start index %d out of range", start);
	if (stop > a.size() || stop < 0)
            log_fatal("Stop index %d out of range", stop);
	if (step >= a.size() || step <= 0)
            log_fatal("Step index %d out of range", step);
	if (start >= stop)
            log_fatal("Start index %d >= stop index %d", start, stop);

	// Get stop index corresponding to step parameter
        int count = ((stop - start) + (step - 1) / step) * step;

	// Build new TS
	G3VectorFixedPtr out(new G3VectorFixed(count, a.precision));
	for (int i = start, j = 0; i < stop; i += step, j++)
            (*out)[j] = a[i];
	
	return out;
    }
}

static PyBufferProcs timestream_bufferprocs;

PYBINDINGS("so3g") {
    namespace bp = boost::python;

    bp::object ts =
        EXPORT_FRAMEOBJECT(G3VectorFixed, init<>(),
                           "Vector of data encoded with a fixed-point scheme.  Data "
                           "can be serialized using FLAC compression, provided that "
                           "the dynamic range is less than 24 bits and the precision "
                           "(e.g. 0.01) has been specified.  Buffer protocol exposes "
                           "the internal storage for numpy access through "
                           "np.asarray(...).")
        .def(bp::init<G3Vector<double>, double>((bp::arg("data"), 
                                                 bp::arg("precision")=1.0)))
        .def_readwrite("precision", &G3VectorFixed::precision,
                       "Precision at which to store the data.")
        .def_readwrite("flac_level", &G3VectorFixed::flac_level,
                       "Level of flac compression to use.  0 to disable, 9 for max, "
                       "5 is usually enough.")
        .add_property("n_samples", &G3VectorFixed_nsamples,
                      "Number of samples in the timestream. Equivalent to len(ts)")
        .def("CheckPrecision", &G3VectorFixed::CheckPrecision,
             "Count the number of rounding failures, given current "
             "precision setting.")
        .def("CheckRange", &G3VectorFixed::CheckRange,
             "Count the number of range failures, given current "
             "precision setting.")
        .def("_cxxslice", G3VectorFixed_getslice, "Slice-only __getitem__")
	;
    register_pointer_conversions<G3VectorFixed>();

    // Add buffer protocol interface
    PyTypeObject *tsclass = (PyTypeObject *)ts.ptr();
    timestream_bufferprocs.bf_getbuffer = G3VectorFixed_getbuffer;
    tsclass->tp_as_buffer = &timestream_bufferprocs;

}

