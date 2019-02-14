#pragma once

#include <G3Frame.h>
#include <stdint.h>
#include <string>

using namespace std;

// This is a C++ representation of FITS world coordinate system parameters.
// It should ideally have full wcslib support so one can query coordinates
// etc. But for now, it just stores and serializes the wcs header information.
// This can be changed later withou the python side needing to worry.
//
// We use the header string representation to communicate between python and C++.
// This is the standard serialization format, so it should not lose any important
// information. But it has a performance cost. The raw conversions from wcslib
// take 20 us from wcs to string and 60 us from string to wcs, but the python
// wrapper adds a lot of overhead, bringing it to 1 ms and 1.3 ms respectively.
//
// At the current level this whole class could be replaced by just a string, but
// the idea is that it will be expanded later.

class G3WCS : public G3FrameObject {
public:
    string header;

    G3WCS();
    G3WCS(const G3WCS &);
    G3WCS(const string & hdr);
    // Required for G3FrameObjects.
    string Description() const;
    template <class A> void serialize(A &ar, unsigned v);
};

G3_SERIALIZABLE(G3WCS, 0);
