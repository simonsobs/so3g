#pragma once

#include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>

#include <nanobind/nanobind.h>

namespace nb = nanobind;


class LookupTable
{
public:
    int n;
    double dx;
    std::vector<double> vals;

    LookupTable() {};
    LookupTable(int _n) : n(_n) {
        vals.resize(n);
    }

    double get_raw(double x) {
        double d = x / dx;
        int idx = (int)d;
        if (x < 0)
            return vals[0];
        if (idx >= n - 1)
            return vals[n - 1];
        d -= idx;
        return vals[idx]*(1-d) + vals[idx+1]*d;
    }
};

class asinTable : public LookupTable
{
public:
    asinTable(int _n) : LookupTable(_n) {
        dx = 1. / (n - 1);
        for (int i = 0; i < n; i++)
            vals[i] = asin(i*dx);
    }

    double get(double x) {
        if (x < 0)
            return -get_raw(-x);
        return get_raw(x);
    }
};

class atan2Table : public LookupTable
{
public:
    atan2Table(int _n) : LookupTable(_n) {
        vals.resize(n);
        dx = 1. / (n - 1);
        for (int i = 0; i < n; i++)
            vals[i] = atan2(i*dx, 1.);
    }

    double get(double y, double x) {
        // Pedantic note: this doesn't return exactly the same thing
        // as atan2... it doesn't recognize -0.
        if (y < 0)
            return -get(-y, x);
        if (x < 0)
            return M_PI - get(y, -x);
        if (y == 0)
            return vals[0];
        if (x < y)
            return M_PI/2 - get_raw(x/y);
        return get_raw(y/x);
    }
};

void register_so_linterp(nb::module_ & m);
