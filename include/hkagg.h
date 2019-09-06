#pragma once

#include <G3Frame.h>
#include <G3Map.h>

#include <stdint.h>

using namespace std;

enum HKFrameType {
     session = 0,
     status = 1,
     data = 2,
};


class IrregBlock : public G3MapFrameObject {
    // Storage for a set of co-sampled data vectors, and the single
    // vector of associated timestamps.  The vectors can have
    // different types, but must be one of the explicitly handled
    // types.  Run Check() to confirm your object has been populated
    // properly.
public:
    G3VectorTime times;
    const IrregBlock Concatenate(const IrregBlock other);
    bool Check();

    string Description() const;
    string Summary() const;

    template <class A> void save(A &ar, unsigned v) const;
    template <class A> void load(A &ar, unsigned v);
};


class IrregBlockDouble : public G3FrameObject {
    // Stores a block of timestamped data.  This consists of named
    // vectors in .data, and a vector of timestamps in .t.  The user
    // should assure that all these vectors are the same length.
    //
    // In the present version, all data vectors as well as the
    // timestamps are doubles.
public:
    string prefix;
    G3MapVectorDouble data;
    G3VectorDouble t;

    string Description() const;
    string Summary() const;
    template <class A> void serialize(A &ar, unsigned v);
};

namespace cereal {
    template <class A> struct specialize<
        A, IrregBlock, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(IrregBlock, 0);
G3_SERIALIZABLE(IrregBlockDouble, 0);
