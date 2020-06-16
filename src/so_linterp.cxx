#include <so_linterp.h>

#include <pybindings.h>
#include <boost/python.hpp>
namespace bp = boost::python;

double test_trig(int table_size, int verbose)
{
    // Report maximimum absolute discrepancy in angle.
    double worst = 0.;
    double lo = -1., hi = 1., step = .01;

    asinTable asin_lt(table_size);
    for (double x = lo; x < hi; x += step) {
        double y0 = asin(x);
        double y1 = asin_lt.get(x);
        worst = std::max(worst, abs(y1 - y0));
        if (verbose)
            std::cout << "asin(" << x << ")" << " "
                      << y0 << " " << y1 << " " << y1 - y0 << "\n";
    }

    atan2Table atan2_lt(table_size);
    for (double _x = -3; _x < 3.1; _x += 0.5) {
        for (double _y = lo; _y < hi; _y += step) {
            double y0 = atan2(_y, _x);
            double y1 = atan2_lt.get(_y, _x);
            worst = std::max(worst, abs(y1 - y0));
            if (verbose)
                std::cout << "atan2(" << _y << ", " << _x << ") "
                          << y0 << " " << y1 << " " << y1 - y0 << "\n";
        }
    }

    return worst;
}

PYBINDINGS("so3g")
{
    bp::def("test_trig", test_trig,
        "For use in test suite -- determines worst arctrig discrepancy, in radians.");
}
