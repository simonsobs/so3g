#define NO_IMPORT_ARRAY

#include <pybindings.h>

#include <iostream>
#include <limits>
#include <type_traits>

#include <boost/python.hpp>

#include <Rebundler.h>

template <>
bp::object Rebundler<G3TimestreamMap>::Process(
    G3FrameObjectPtr p)
{
    bp::list out;
    auto m = dynamic_cast<G3TimestreamMap *>(p.get());
    buffer.push_back(p);
    return out;
}

template <>
bp::object Rebundler<G3TimestreamMap>::ExtractIntervalTime(
    G3Time start_time, G3Time end_time, bool flush)
{
    // Drop all frames that are to the left of the domain.
    while (buffer.size()) {
        auto tm0 = dynamic_cast<G3TimestreamMap *>(buffer.front().get());
        G3Time stop = tm0->GetStopTime();
        if (stop < start_time)
            buffer.pop_front();
        else
            break;
    }

    // Do we cover the interval completely?
    if (!flush) {
        // That we've passed the end of requested interval.
        if (buffer.size() <= 0)
            return bp::object();
        auto tm0 = dynamic_cast<G3TimestreamMap *>(buffer.back().get());
        G3Time stop = tm0->GetStopTime();
        if (stop < end_time)
            return bp::object();
    }
    // Create the new object.
    bp::list out;
    if (buffer.size() <= 0)
        return out;

    vector<pair<int,int>> slices;
    int n = 0;
    double rate;
    pair<G3Time,G3Time> times;
    for (auto tm_ptr: buffer) {
        auto tm = dynamic_cast<G3TimestreamMap *>(tm_ptr.get());
        auto t0 = tm->GetStartTime();
        auto t1 = tm->GetStopTime();
        double _i0 = (double(start_time) - double(t0)) * tm->GetSampleRate();
        double _i1 = (double(end_time)   - double(t0)) * tm->GetSampleRate();
        int i0 = _i0 > 0 ? int(_i0) : 0;
        int i1 = _i1 > 0 ? int(_i1) : 0;
        if (i1 > tm->NSamples())
            i1 = tm->NSamples();
        slices.push_back(make_pair(i0,i1));
        if (n == 0) {
            rate = tm->GetSampleRate();
            times.first = G3Time(t0 + int64_t(i0 / rate));
        }
        n += i1 - i0;
    }
    times.second = G3Time(times.first + int64_t((n-1) / rate));
    consumed_to = G3Time(times.second + int64_t(0.5 / rate));
    
    G3TimestreamMap tm_out;
    auto tmr_ptr = buffer.front();
    auto tmr_ref = dynamic_cast<G3TimestreamMap *>(tmr_ptr.get());
    for (auto key_ref: *tmr_ref) {
        auto pv = new G3Timestream(n);
        tm_out[key_ref.first] = G3TimestreamPtr(pv);
        auto vi = pv->begin();
        pv->start = times.first;
        pv->stop = times.second;
        auto sli = slices.begin();
        auto tm_ptr = buffer.begin();
        while (tm_ptr != buffer.end()) {
            auto tm = dynamic_cast<G3TimestreamMap *>(tm_ptr->get());
            auto p = tm->at(key_ref.first);
            int count = sli->second - sli->first;
            if (count > 0) {
                memcpy(&(*vi), &(p->at(sli->first)), count*sizeof(*vi));
                vi += count;
            }
            sli++;
            tm_ptr++;
        }
    }
    out.append(tm_out);
    return out;
}


template <>
bp::object Rebundler<G3TimestreamMap>::ExtractInterval(
    int length, bool flush)
{
    // Do we have enough information to convert "length" to timestamps?
    for (auto tmr_ptr: buffer) {
        auto tmr = dynamic_cast<G3TimestreamMap *>(tmr_ptr.get());
        if (tmr->NSamples() == 0)
            continue;
        if (consumed_to > tmr->GetStopTime())
            continue;
        G3Time start = tmr->GetStartTime();
        if (start < consumed_to)
            start = consumed_to;
        G3Time stop = G3Time(start + int64_t(length / tmr->GetSampleRate()));
        return ExtractIntervalTime(start, stop, flush);
    }
    return bp::object();
}


using namespace boost::python;

#define EXPORT_REBUNDLER(DOMAIN_TYPE, CLASSNAME) \
    class_<CLASSNAME>(#CLASSNAME) \
    .def("Process", &CLASSNAME::Process, \
         "Add element.") \
    .def("ExtractIntervalTime", &CLASSNAME::ExtractIntervalTime, \
         "Rebundle into interval.")\
    .def("ExtractInterval", &CLASSNAME::ExtractInterval, \
         "Rebundle into interval.")


PYBINDINGS("so3g")
{
    EXPORT_REBUNDLER(G3TimestreamMap,  RebundlerPrimaryMap);
}
