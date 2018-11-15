#pragma once

#include <G3Frame.h>
#include <G3Map.h>

#include <stdint.h>

using namespace std;

class HKInfo : public G3FrameObject {
public:
    int32_t session_id;
    string hk_source;
    
    string Description() const;
    string Summary() const;
    template <class A> void serialize(A &ar, unsigned v);
};

G3_SERIALIZABLE(HKInfo, 0);
