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

class IrregBlockDouble : public G3FrameObject {
public:
    string prefix;
    G3MapVectorDouble data;
    G3VectorDouble t;

    string Description() const;
    string Summary() const;
    template <class A> void serialize(A &ar, unsigned v);
};


G3_SERIALIZABLE(IrregBlockDouble, 0);
