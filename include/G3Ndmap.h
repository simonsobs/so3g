#pragma once

#include <G3Frame.h>
#include <G3Ndarray.h>
#include <G3WCS.h>

using namespace std;

// This is a C++ representation of ndmaps (enmaps) from pixell.
// Like G3Ndarray and G3WCS it should get some more useful functions
// from the C++ side. For now, its purpose is simply allowing python
// to store these objects in frames, and serialization.

class G3Ndmap : public G3FrameObject {
public:
    G3Ndarray data;
    G3WCS wcs;

    G3Ndmap();
    G3Ndmap(const G3Ndmap &);
    G3Ndmap(const G3Ndarray &, const G3WCS &);
    G3Ndmap(const bp::numpy::ndarray &, const string &);
    // Required for G3FrameObjects.
    string Description() const;
    template <class A> void serialize(A &ar, unsigned v);
};

G3_SERIALIZABLE(G3Ndmap, 0);
G3_POINTERS(G3Ndmap);
