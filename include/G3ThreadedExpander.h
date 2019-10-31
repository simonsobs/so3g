#ifndef _G3_EXPOOL_H
#define _G3_EXPOOL_H

#include <string>
#include <vector>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <boost/iostreams/filtering_stream.hpp>

#include <G3Module.h>

class G3ThreadedExpander : public G3Module {
public:
    G3ThreadedExpander(int n_threads);
    ~G3ThreadedExpander();

    // This is the storage area for unpacked objects.
    typedef std::map<std::string,G3FrameObjectConstPtr> UFrame;
    std::vector<UFrame> unpacked_frames;

    // Use of the Expander has two phases; first it is a G3Module that
    // processes frames through Process().  During that phase,
    // multiple pthreads may be in operation and the state and data
    // variables (including unpacked_frames) are protected by a mutex.
    // The threaded phase ends once EndProcessing has been passed
    // through Process.

    void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

    // In the second phase, contiguous data can be extracted from the
    // unpacked_frames using the methods below.

    bp::object Extract();
    bp::object RestreamTimestream(bp::list path,
                                  bp::object terms,
                                  bp::object comp_info);
    bp::object Restream(bp::list path);

private:

    int n_threads_;
    std::vector<std::thread> threads_;
    static void expander_thread_func(int index, G3ThreadedExpander *parent);

    // This lock protects all private variables.
    std::mutex lock_;
    std::condition_variable cv_;

    int exit_now_;
    int index_in;
    int index_out;

    // Work for each thread.
    typedef struct {
        enum {idle, working, done} state;
        G3FramePtr frame;
        UFrame uframe;
    } FrameJob;
    std::vector<FrameJob> jobs_;

    SET_LOGGER("G3ThreadedExpander");
};

G3_POINTERS(G3ThreadedExpander);

#endif
