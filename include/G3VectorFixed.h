// -*- mode: c++ -*-

#ifndef _G3_VECTOR_FIXED_H
#define _G3_VECTOR_FIXED_H

#include <G3Frame.h>
#include <G3TimeStamp.h>
#include <G3Vector.h>
#include <vector>
#include <map>

const int DEFAULT_FLAC_LEVEL = 5;
const double DEFAULT_PRECISION = 1.0;

class G3VectorFixed : public G3Vector<double> {
public:
    // Default constructor
    G3VectorFixed(double precision_=DEFAULT_PRECISION) :
        G3Vector<double>(), precision(precision_), flac_level(DEFAULT_FLAC_LEVEL) {}

    // Pre-allocated.
    G3VectorFixed(std::vector<double>::size_type s, double precision_=DEFAULT_PRECISION) :
        G3Vector<double>(s), precision(precision_), flac_level(DEFAULT_FLAC_LEVEL) {}

    // Copy construct, from floating vector type, with explicit precision.
    G3VectorFixed(const G3Vector<double> &data_in, double precision_) :
        G3Vector<double>(data_in), precision(precision_) {}

    // Copy constructor.
    G3VectorFixed(const G3VectorFixed &r) :
        G3Vector<double>(r), precision(r.precision),
        flac_level(r.flac_level) {}

    int CheckPrecision();
    int CheckRange();

    // Cereal.
    template <class A> void load(A &ar, unsigned v);
    template <class A> void save(A &ar, unsigned v) const;

    // Repr.
    std::string Description() const;
    std::string Summary() const { return Description(); };

    //! Precision multiplier.  E.g. 0.01 to handle values in
    //! [-2**23*0.01, (2**23-1)*0.01]
    double precision;
    //! FLAC compression level (1-9, 0 disables).
    int32_t flac_level;
private:
    SET_LOGGER("G3VectorFixed");
};

G3_POINTERS(G3VectorFixed);

namespace cereal {
    template <class A> struct specialize<A, G3VectorFixed,
                                         cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(G3VectorFixed, 0);

#endif
