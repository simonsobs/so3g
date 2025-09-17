#pragma once

#include <boost/python.hpp>
namespace bp = boost::python;

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/facilities/overload.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/python.hpp>

#include <cstdint>

class G3ModuleRegistrator {
public:
	G3ModuleRegistrator(const char *mod, void (*def)());
	static void CallRegistrarsFor(const char *mod);
};

#define EXPORT_G3MODULE_AND(mod, T, init, docstring, other_defs)   \
	static void registerfunc##T() { \
		using namespace boost::python; \
		class_<T, bases<G3Module>, boost::shared_ptr<T>, \
		  boost::noncopyable>(#T, docstring, init) \
		    .def_readonly("__g3module__", true) \
                other_defs \
		; \
	} \
	static G3ModuleRegistrator register##T(mod, registerfunc##T);

#define EXPORT_G3MODULE(mod, T, init, docstring) \
    EXPORT_G3MODULE_AND(mod, T, init, docstring, )


#define PYBINDINGS(mod) \
	static void ___pybindings_registerfunc(); \
	static G3ModuleRegistrator ___pybindings_register(mod, ___pybindings_registerfunc); \
	static void ___pybindings_registerfunc() 

#define EXPORT_FRAMEOBJECT(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, boost::shared_ptr<T> >(#T, docstring, boost::python::initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

#define EXPORT_FRAMEOBJECT_NOINITNAMESPACE(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, boost::shared_ptr<T> >(#T, docstring, initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

// Declare a python module with a name and the name of its enclosing package scope.
// name should be be a bare token, while pkg should be a string literal, e.g.:
//     SPT3G_PYTHON_MODULE_2(foo, "spt3g.bar")
// for a package whose fully qualified name will be spt3g.bar.foo
#define SPT3G_PYTHON_MODULE_2(name, pkg) \
BOOST_PYTHON_MODULE(name) { \
	namespace bp = boost::python; \
	auto mod = bp::scope(); \
	std::string package_prefix = pkg; \
	std::string full_name = package_prefix + "." + bp::extract<std::string>(mod.attr("__name__"))(); \
	mod.attr("__name__") = full_name; \
	mod.attr("__package__") = package_prefix; \
	void BOOST_PP_CAT(spt3g_init_module_, name)(); \
	BOOST_PP_CAT(spt3g_init_module_, name)(); \
	if(PY_MAJOR_VERSION < 3){ \
		Py_INCREF(mod.ptr()); \
		PyDict_SetItemString(PyImport_GetModuleDict(),full_name.c_str(),mod.ptr()); \
	} \
} \
void BOOST_PP_CAT(spt3g_init_module_, name)()

// Declare a python module with the given name, assuming that the enclosing package 
// is the default "spt3g".
#define SPT3G_PYTHON_MODULE_1(name) SPT3G_PYTHON_MODULE_2(name, "spt3g")

// Declare a python module with a name and optionally the name of its enclosing package scope.
// name should be be a bare token, while if provided the enclosing package name should be a
// string literal.
// If the enclosing package name is not specified, it will default to "spt3g".
#define SPT3G_PYTHON_MODULE(...) BOOST_PP_OVERLOAD(SPT3G_PYTHON_MODULE_,__VA_ARGS__)(__VA_ARGS__)
