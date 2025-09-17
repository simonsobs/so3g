#include <pybindings.h>
//#include <serialization.h>

#include <boost/python.hpp>
#include <boost/preprocessor.hpp>

// #include <G3Frame.h>
// #include <G3Data.h>
// #include <G3Module.h>
// #include <G3EventBuilder.h>
// #include <G3Pipeline.h>
// #include <G3Timestream.h>
// #include <G3SimpleLoggers.h>
// #include <G3Constants.h>

namespace bp = boost::python;

// The following implements the headerless module registration code
typedef std::map<std::string, std::vector<void (*)()> > module_reg_t;
static module_reg_t *modregs = NULL;

G3ModuleRegistrator::G3ModuleRegistrator(const char *mod, void (*def)())
{
	if (modregs == NULL)
		modregs = new module_reg_t;
	(*modregs)[mod].push_back(def);
}

void G3ModuleRegistrator::CallRegistrarsFor(const char *mod)
{
	for (auto i = (*modregs)[mod].begin(); i != (*modregs)[mod].end(); i++)
		(*i)();
}


// Nonsense boilerplate for POD vector numpy bindings
#define numpy_vector_infrastructure(T, name, conv) \
template <> \
boost::shared_ptr<std::vector<T> > \
container_from_object(boost::python::object v) \
{ \
	return numpy_container_from_object<std::vector<T> >(v); \
} \
static int \
vector_getbuffer_##name(PyObject *obj, Py_buffer *view, int flags) \
{ \
	return pyvector_getbuffer<T>(obj, view, flags, conv); \
} \
static PyBufferProcs vec_bufferprocs_##name; \
struct numpy_vector_from_python_##name { \
	numpy_vector_from_python_##name() { \
		boost::python::converter::registry::push_back( \
		    &convertible, &construct, \
		    boost::python::type_id<std::vector<T> >()); \
	} \
	static void *convertible(PyObject* obj_ptr) { \
		Py_buffer view; \
		if (PyObject_GetBuffer(obj_ptr, &view, \
		    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) { \
			PyErr_Clear(); \
			return NULL; \
		} \
		if (view.ndim == 0) { \
			PyBuffer_Release(&view); \
			return NULL; \
		} \
		PyBuffer_Release(&view); \
		return obj_ptr; \
	} \
	static void construct(PyObject* obj_ptr, \
	    boost::python::converter::rvalue_from_python_stage1_data* data) { \
		void* storage = ( \
		    (boost::python::converter::rvalue_from_python_storage<std::vector<T> >*)data)->storage.bytes; \
		new (storage) std::vector<T>; \
		boost::shared_ptr<std::vector<T> > swap_storage = numpy_container_from_object<std::vector<T> >(boost::python::object(boost::python::handle<>(boost::python::borrowed(obj_ptr)))); \
		((std::vector<T> *)(storage))->swap(*swap_storage); \
		data->convertible = storage; \
	} \
};

#if PY_MAJOR_VERSION < 3
#define numpy_vector_of(T, name, desc) \
{ \
	numpy_vector_from_python_##name(); \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##name.bf_getbuffer = vector_getbuffer_##name; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##name; \
	vdclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER; \
}
#else
#define numpy_vector_of(T, name, desc) \
{ \
	numpy_vector_from_python_##name(); \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##name.bf_getbuffer = vector_getbuffer_##name; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##name; \
}
#endif

numpy_vector_infrastructure(int64_t, int64_t, "q")
numpy_vector_infrastructure(uint64_t, uint64_t, "Q")
numpy_vector_infrastructure(int32_t, int32_t, "i")
numpy_vector_infrastructure(uint32_t, uint32_t, "I")
numpy_vector_infrastructure(double, double, "d")
numpy_vector_infrastructure(float, float, "f")

// Apple, for their own insane reasons, defines uint64_t as
// "unsigned long long" even on LP64 systems where longs are
// 64-bit. Because "long long" (not a standard C type!) is not
// actually the same type as "long", even when both are 64-bit
// integers, the uint64_t definition above does not do the right
// thing for size_t on 64-bit Apple systems.
//
// Thanks, Apple. "Think Different!"
#if defined(__APPLE__) && defined(__LP64__)
numpy_vector_infrastructure(size_t, size_t, "L")
numpy_vector_infrastructure(ssize_t, ssize_t, "l")
struct apple_size
{
	static PyObject* convert(const std::vector<size_t> &arg) {
		return boost::python::to_python_value<std::vector<uint64_t> >()(*(std::vector<uint64_t> *)(uintptr_t)(&arg));
	}
};

struct apple_ssize
{
	static PyObject* convert(const std::vector<ssize_t> &arg) {
		return boost::python::to_python_value<std::vector<int64_t> >()(*(std::vector<int64_t> *)(intptr_t)(&arg));
	}
};
#endif

template <> boost::shared_ptr<std::vector<std::complex<float> > >
numpy_container_from_object(boost::python::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<float> > >(v);
}
template <> boost::shared_ptr<std::vector<std::complex<double> > >
numpy_container_from_object(boost::python::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<double> > >(v);
}

numpy_vector_infrastructure(std::complex<double>, cxdouble, "Zd");
numpy_vector_infrastructure(std::complex<float>, cxfloat, "Zf");

PYBINDINGS("so3g") //SPT3G_PYTHON_MODULE(core)
{
	bp::docstring_options docopts(true, true, false);

	// Some POD types
	register_vector_of<bool>("Bool");
	numpy_vector_of(int64_t, int64_t, "Int64");
	numpy_vector_of(uint64_t, uint64_t, "UInt64");
	numpy_vector_of(int32_t, int32_t, "Int");
	numpy_vector_of(uint32_t, uint32_t, "UInt");

#if defined(__APPLE__) && defined(__LP64__)
	numpy_vector_from_python_size_t();
	numpy_vector_from_python_ssize_t();
	bp::to_python_converter<std::vector<size_t>, apple_size, false>();
	bp::to_python_converter<std::vector<ssize_t>, apple_ssize, false>();
#endif

	numpy_vector_of(double, double, "Double");
	numpy_vector_of(std::complex<double>, cxdouble, "ComplexDouble");
	numpy_vector_of(float, float, "Float");
	numpy_vector_of(std::complex<float>, cxfloat, "ComplexFloat");
	register_vector_of<std::string>("String");
	register_vector_of<unsigned char>("UnsignedChar");

	// Do everything else
	//G3ModuleRegistrator::CallRegistrarsFor("core");
}

