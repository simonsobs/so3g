#pragma once

#include <G3Frame.h>
#include <G3Map.h>

#include <stdint.h>

using namespace std;

class SampleFlagSet: public G3MapVectorInt {
public:
    /* Representation of sample cuts. */
    std::map<string, G3MapVectorInt> flags;
    std::string name;
    
    /* Required. */
    string Description() const;
    string Summary() const;
    template <class A> void serialize(A &ar, unsigned v);
};

// class SampleFlags : public G3FrameObject {
// public:
//     /* Representation of sample cuts. */
//     std::vector<std::string> flag_names;
//     std::vector<G3
    
//     /* Required. */
//     string Description() const;
//     string Summary() const;
//     template <class A> void serialize(A &ar, unsigned v);
// };

G3_SERIALIZABLE(SampleFlagSet, 0);
