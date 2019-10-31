#define NO_IMPORT_ARRAY
#include "so3g_numpy.h"

#include <pybindings.h>
#include <dataio.h>
#include <G3Map.h>
#include <G3Timestream.h>

#include <boost/python.hpp>

#include <chrono>
#include <iostream>

#include "exceptions.h"
#include "G3ThreadedExpander.h"

// Each thread is started in the constructor and sits idle until its
// particular slot of the jobs_ vector is filled with work to do.
// This is identified by a change of jobs_ state from "idle" to
// "working".  When the thread has completed the extraction of objects
// from the frame, it sets state to "done" and signals back to main
// thread.
//
// All threads will wait on the same condition variable, and will
// notify whenever a job state changes or the exit_now_ trigger has
// been set.

void G3ThreadedExpander::expander_thread_func(int index, G3ThreadedExpander *parent)
{
    while (true) {
        {
            std::unique_lock<std::mutex> lk(parent->lock_);
            parent->cv_.wait(lk, [&parent, index]{ return
                        parent->jobs_[index].state == FrameJob::working ||
                        parent->exit_now_; });
            if (parent->exit_now_)
                return;
            if (parent->jobs_[index].state != FrameJob::working)
                continue;
        }
        auto job = &parent->jobs_[index];
        
        auto uf = &job->uframe;
        uf->clear();

        auto f = job->frame;
        for (auto key: f->Keys()) {
            // Yes, this is the slow part.
            uf->operator[](key) =  f->operator[](key);
            f->Delete(key);
        }

        {
            std::unique_lock<std::mutex> lk(parent->lock_);
            parent->jobs_[index].state = FrameJob::done;
            parent->cv_.notify_all();
        }
    }
}


G3ThreadedExpander::G3ThreadedExpander(int n_threads) :
    n_threads_(n_threads)
{
    std::lock_guard<std::mutex> lk(lock_);
    exit_now_ = 0;
    index_in = 0;
    index_out = 0;

    FrameJob idle_job = {FrameJob::idle};
    bp::object none;
    for (int i=0; i<n_threads_; i++) {
        threads_.push_back(std::thread(expander_thread_func, i, this));
        jobs_.push_back(idle_job);
    }
}

G3ThreadedExpander::~G3ThreadedExpander()
{
    {
        std::lock_guard<std::mutex> lk(lock_);
        exit_now_ = 1;
        cv_.notify_all();
    }

    for (int i=0; i<n_threads_; i++)
        threads_[i].join();
}

void G3ThreadedExpander::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
    // In absence of threads, pass through.
    if (n_threads_ == 0) {
	out.push_back(frame);
        return;
    }

    bp::object none;
    while (true) {
        // Check our threads and pull out processed frames.
        std::unique_lock<std::mutex> lk(lock_);
        while (index_out < index_in) {
            int slot = index_out % n_threads_;
            auto job = &jobs_[slot];
            if (job->state == FrameJob::done) {
                out.push_back(job->frame);
                unpacked_frames.push_back(std::move(job->uframe));

                job->state = FrameJob::idle;
                job->frame = nullptr;
                ++index_out;
            } else if (frame->type == G3Frame::EndProcessing) {
                // Fine, we'll wait.
                cv_.wait(lk, [this,slot]{ return this->jobs_[slot].state == FrameJob::done; });
            } else {
                // No problem, move on.
                break;
            }
        }

        if (frame->type == G3Frame::EndProcessing) {
            out.push_back(frame);
            break;
        }

        // We must cache this frame somehow.
        int slot = index_in % n_threads_;
        auto job = &jobs_[slot];
        if (job->state == 0) {
            job->state = FrameJob::working;
            job->frame = frame;
            cv_.notify_all();
            frame = nullptr;
            ++index_in;
            break;
        } else {
            // Wait for the slot to free up, then loop around.
            cv_.wait(lk, [this,slot]{ return this->jobs_[slot].state == FrameJob::done; });
        }
    }
}

// Start of post-read extraction methods.
//
// Below this point are only function used after stream Processing is
// complete, i.e. after an EndProcessing has been pushed through.

bp::object G3ThreadedExpander::Extract() {
    // Transfer all stored frame objects into python dictionaries and
    // return.  Clear the stored objects.
    bp::list out;
    for (auto unpack: unpacked_frames) {
        bp::dict repack;
        for (auto item: unpack)
            repack[item.first] =
                bp::object(boost::const_pointer_cast<G3FrameObject>(item.second));
        out.append(repack);
    }
    unpacked_frames.clear();
    return out;
}

//! descend_to
//
// Examine each of the map(str,G3FrameObjectPtr) elements in
// unpacked_frames and attempt to descend to the item identified by
// path_, a list of strings.  Returns a list of successfully
// identified objects, at most one per input frame.  If a path cannot
// be descended because of type conversion issues, an exception is
// raised.  If a path cannot be descended because a key doesn't exist,
// the frame is skipped.

static
std::vector<G3FrameObjectConstPtr> descend_to(
    const std::vector<G3ThreadedExpander::UFrame> unpacked_frames,
    bp::list path_)
{
    // Convert "path" to native vector<string>.
    std::vector<std::string> path;
    for (int i=0; i<bp::len(path_); i++)
        path.push_back(bp::extract<std::string>(path_[i])());
    assert(path.size() > 0);

    std::vector<G3FrameObjectConstPtr> endpoints;
    for (auto unpack: unpacked_frames) {
        // Each unpacked frame is C++ map from string to *FrameObject.
        // We have to decode, then descend the rest of the tree.
        auto p = path.begin();
        auto it = unpack.find(*p);
        if (it == unpack.end())
            continue;
        auto target = it->second;
        for (++p; p != path.end(); ++p) {
            // Promote to a Map, search it, get new result.
            auto tmap = boost::dynamic_pointer_cast
                <const G3MapFrameObject>(target);
            if (tmap == nullptr) {
                std::ostringstream s;
                s << "Requested path element could not be decoded as a Map: ";
                for (auto q=path.begin(); q!=p; ++q)
                    s << *q << " -> ";
                s << *p << "!\n";
                throw key_crawling_exception(s.str());
            }
            auto it = tmap->find(*p);
            assert(it != tmap->end());
            target = it->second;
        }
        endpoints.push_back(target);
    }
    return endpoints;
}

struct ndarray_bundle {
    void *dest;
    PyObject *po;
    bp::object bpo;
};

struct ndarray_bundle create_ndarray_1d(int dtype, int n_samples)
{
    ndarray_bundle out;
    npy_intp dims[32];
    dims[0] = n_samples;
    out.po = PyArray_EMPTY(1, dims, dtype, 0);
    out.dest = (void*)PyArray_DATA((PyArrayObject *)out.po);
    out.bpo = bp::object(bp::handle<>(out.po));
    return out;
}

//! Restream(path)
//
// Concatenate a field from cached frame data into numpy array.
//
// The "path" must be a list of strings that address the target (a
// G3VectorDouble) in each unpacked object cache.
//
// Returns a numpy array, or None if array would be of 0 length.

bp::object G3ThreadedExpander::Restream(bp::list path)
{
    auto series = descend_to(unpacked_frames, path);

    // Scan for total length.
    int n_samples = 0;
    std::vector<G3VectorDoubleConstPtr> vo;
    for (auto fo: series) {
        auto v = boost::dynamic_pointer_cast<const G3VectorDouble>(fo);
        assert(v != nullptr);
        n_samples += v->size();
        vo.push_back(v);
    }

    if (n_samples == 0) {
        bp::object none;
        return none;
    }

    auto arr = create_ndarray_1d(NPY_FLOAT64, n_samples);
    auto dest = (double*)arr.dest;

     //Fill it.
    int offset = 0;
    for (auto v: vo) {
        int n = v->size();
        std::copy(v->begin(), v->end(), dest+offset);
        offset += n;
    }

    return arr.bpo;
}


//! RestreamTimestream(path, keys)
//
// Concatenates selected fields from cached frame data into numpy
// arrays.
//
// The "path" must be a list of strings that address the target (a
// G3TimestreamMap) in each unpacked object cache.  If keys is None,
// then all elements are copied out.  Otherwise, only the elements
// listed in keys are copied out.  If a field requested in keys does
// not exist in the first data frame, it will be ignored.
//
// comp_spec must either be None or an instructive tuple.
//
// Returns (keys, data, start, stop) where data is a dictionary of numpy
// arrays.
//

bp::object G3ThreadedExpander::RestreamTimestream(bp::list path, bp::object keys,
                                              bp::object comp_spec)
{
    auto series = descend_to(unpacked_frames, path);

    std::vector<G3FrameObjectConstPtr> comp_series0;
    std::vector<G3FrameObjectConstPtr> comp_series1;
    
    int compression_scheme = 0;
    if (comp_spec.ptr() != Py_None) {
        bp::tuple t = bp::extract<bp::tuple>(comp_spec)();
        assert(bp::extract<int>(t[0]) == 1);
        compression_scheme = 1;
        // Then assume it's a path to standardized compression dictionaries.
        comp_series0 = descend_to(unpacked_frames, bp::extract<bp::list>(t[1]));
        comp_series1 = descend_to(unpacked_frames, bp::extract<bp::list>(t[2]));
        assert(comp_series0.size() == comp_series1.size());
        assert(comp_series0.size() == series.size());
    }

    // Convert "keys" to non-python for safe use.
    bool get_all = (keys.ptr() == Py_None);
    std::vector<std::string> keys_;
    if (!get_all) {
        bp::list t = bp::extract<bp::list>(keys)();
        for (int i=0; i<bp::len(t); i++)
            keys_.push_back(bp::extract<std::string>(t[i])());
    }

    // Determine total length, save timestamps.
    std::vector<G3TimestreamMapConstPtr> ts_maps;
    int n_samples = 0;
    G3Time start, stop;
    for (auto target: series) {
        auto vm = boost::dynamic_pointer_cast<const G3TimestreamMap>(target);
        if (vm == nullptr) {
            std::ostringstream s;
            s << "A terminal path element could not be decoded as a "
              << "G3TimestreamMap.";
            throw value_decoding_exception(s.str());
        }
        if (n_samples==0)
            start = vm->GetStartTime();
        stop = vm->GetStopTime();
        ts_maps.push_back(vm);
        if (keys_.size() == 0) {
            // Populate it.
            for (auto item: *vm)
                keys_.push_back(item.first);
        }
        auto it = vm->find(keys_[0]);
        assert(it != vm->end());
        auto t = (*it).second;
        n_samples += t->size();
    }

    // Extract compression information, too.
    std::vector<std::vector<std::pair<double,double>>> comp_data;
    if (compression_scheme == 1) {
        // Find each key in compression frames and store the decoded
        // double in comp_data, which will have shape (n_frame, n_det, 2).
        for (int i=0; i<comp_series0.size(); i++) {
            std::vector<std::pair<double,double>> f;
            auto md0 = boost::dynamic_pointer_cast<const G3MapDouble>
                (comp_series0[i]);
            auto md1 = boost::dynamic_pointer_cast<const G3MapDouble>
                (comp_series1[i]);
            assert(md0 != nullptr && md1 != nullptr);
            for (auto k: keys_) {
                auto it0 = md0->find(k);
                auto it1 = md1->find(k);
                assert((it0 != md0->end()) && (it1 != md1->end()));
                f.push_back(std::make_pair(1./(*it0).second, (*it1).second));
            }
            comp_data.push_back(std::move(f));
        }
    }

    // Allocate output -- python again...
    std::vector<double*> dests_d;
    std::vector<float*> dests_f;
    bp::list ndarrays;
    struct ndarray_bundle arr;    
    for (int i=0; i<keys_.size(); i++) {
        switch(compression_scheme) {
        case 0:
            arr = create_ndarray_1d(NPY_FLOAT64, n_samples);
            dests_d.push_back((double*)arr.dest);
            ndarrays.append(arr.bpo);
            break;
        case 1:
            arr = create_ndarray_1d(NPY_FLOAT32, n_samples);
            dests_f.push_back((float*)arr.dest);
            ndarrays.append(arr.bpo);
        }
    }

    //Fill them.
#pragma omp parallel for
    for (int j=0; j<keys_.size(); j++) {
        int offset = 0;
        for (int i=0; i<ts_maps.size(); i++) {
            auto vm = ts_maps[i];
            auto it = vm->find(keys_[j]);
            assert(it != vm->end());
            auto t = (*it).second;
            int n = t->size();
            switch(compression_scheme) {
            case 0:
                std::copy(t->begin(), t->end(), dests_d[j] + offset);
                // for (int k=0; k<n; k++)
                //     dests_d[j][offset+k] = t->operator[](k);
                break;
            case 1:
                auto rescale = comp_data[i][j];
                for (int k=0; k<n; k++)
                    dests_f[j][offset+k] = t->operator[](k) * rescale.first + rescale.second;
                break;
            }
            offset += n;
        }
    }

    // What fields were loaded?
    bp::list keys_out;
    if (get_all) {
        for (auto t: keys_)
            keys_out.append(t);
    } else
        keys_out = bp::extract<bp::list>(keys)();

    return bp::make_tuple(keys_out, ndarrays, start, stop);
}

// End of post-read, single thread extraction methods.

PYBINDINGS("so3g") {
    using namespace boost::python;

    // Instead of EXPORT_G3MODULE since there are two constructors
    class_<G3ThreadedExpander, bases<G3Module>, boost::shared_ptr<G3ThreadedExpander>,
           boost::noncopyable>(
               "G3ThreadedExpander",
               "Extract objects from each frame and store them.  This "
               "extraction is done by a thread pool.  The extracted objects "
               "may then be expanded directly into numpy arrays, taking "
               "advantage of OpenMP.",
               init<int>((arg("n_threads")=0)))
        .def("Extract", &G3ThreadedExpander::Extract,
            "Returns the full unpacked_frames structure, as a list of dicts.")
        .def("Restream", &G3ThreadedExpander::Restream,
             "Extract a G3VectorDouble from the given path.")
        .def("RestreamTimestream", &G3ThreadedExpander::RestreamTimestream,
             "Extract vectors from the series of G3TimestreamMap at the "
             "given path.")
        .def_readonly("__g3module__", true)
	;
}

