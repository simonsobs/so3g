#pragma once

int get_dtype(const bp::object &);

template <typename T>
T _calculate_median(const T*, const int);
