#pragma once

#include <G3Frame.h>
#include <G3Map.h>
#include <G3TimeStamp.h>
#include <G3Timestream.h>
#include <Intervals.h>

#include <deque>
#include <stdint.h>

using namespace std;
namespace bp = boost::python;

template <typename T>
class Rebundler {
public:
    Rebundler() {};
    bp::object Process(G3FrameObjectPtr p);
    bp::object ExtractIntervalTime(G3Time start, G3Time end, bool flush);
    bp::object ExtractInterval(int length, bool flush);

private:
    deque<G3FrameObjectPtr> buffer;
    G3Time consumed_to;
};


typedef Rebundler<G3TimestreamMap> RebundlerPrimaryMap;
