#pragma once

#include <boost/python.hpp>
#include <exception>

// Wrap a Py_buffer view with a destructor that releases the buffer.
// This allows us to do RAII; without it we'd have to explicitly
// cleanup any successful buffer views before throwing exceptions
// related to the content of the buffer.

class BufferWrapper {
public:
    Py_buffer view;
    BufferWrapper() {
        view.obj = NULL;
    }
    ~BufferWrapper() {
        PyBuffer_Release(&view);
    }
};

// so3g_exception is our internal base class, which defines the
// interface we use for converting C++ exceptions to python.

class so3g_exception : std::exception
{
public:
    virtual std::string msg_for_python() const throw() = 0;
};


// The exceptions below should be used when processing objects with
// the buffer protocol (probably numpy arrays).

class buffer_exception : public so3g_exception
{
public:
    std::string var_name;
    buffer_exception(std::string var_name) : var_name{var_name} {}

    std::string msg_for_python() const throw() {
	std::ostringstream s;
        s << "Argument '" << var_name << "' does not expose buffer protocol, "
            "is not contiguous, or does not export a format.";
	return s.str();
    }
};

class shape_exception : public so3g_exception
{
public:
    std::string var_name;
    std::string detail;
    shape_exception(std::string var_name, std::string detail) :
        var_name{var_name}, detail(detail) {}

    std::string msg_for_python() const throw() {
	std::ostringstream s;
        s << "Buffer '" << var_name << "' has incompatible shape: "
          << detail << ".";
	return s.str();
    }
};

class dtype_exception : public so3g_exception
{
public:
    std::string var_name;
    std::string type_str;
    dtype_exception(std::string var_name, std::string type_str) :
        var_name{var_name}, type_str{type_str} {}

    std::string msg_for_python() const throw() {
	std::ostringstream s;
        s << "Expected buffer '" << var_name << "' to contain items of type "
          << type_str << ".";
	return s.str();
    }
};

class agreement_exception : public so3g_exception
{
public:
    std::string var1, var2, prop;
    agreement_exception(std::string var1, std::string var2, std::string prop) :
        var1{var1}, var2{var2}, prop{prop} {}

    std::string msg_for_python() const throw() {
	std::ostringstream s;
        s << "Expected buffers '" << var1 << "' and '" << var2 << "' to have "
          << "the same " << prop << ".";
	return s.str();
    }
};

class general_agreement_exception : public so3g_exception
{
public:
    std::string text;
    general_agreement_exception(std::string text) :
        text{text} {}

    std::string msg_for_python() const throw() {
	return text;
    }
};

class general_type_exception : public so3g_exception
{
public:
    std::string text;
    general_type_exception(std::string text) :
        text{text} {}

    std::string msg_for_python() const throw() {
	return text;
    }
};
