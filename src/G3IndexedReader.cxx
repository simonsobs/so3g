#include <pybindings.h>
#include <dataio.h>
#include <G3IndexedReader.h>
#include "exceptions.h"

#include <boost/filesystem.hpp>

G3IndexedReader::G3IndexedReader(std::string filename, int n_frames_to_read) : 
    prefix_file_(false), n_frames_to_read_(n_frames_to_read),
    n_frames_read_(0)
{
	boost::filesystem::path fpath(filename);
	if (filename.find("://") == std::string::npos &&
	   (!boost::filesystem::exists(fpath) ||
	    !boost::filesystem::is_regular_file(fpath)))
		log_fatal("Could not find file %s", filename.c_str());
	StartFile(filename);
}

G3IndexedReader::G3IndexedReader(std::vector<std::string> filename, int n_frames_to_read) :
    prefix_file_(false), n_frames_to_read_(n_frames_to_read),
    n_frames_read_(0)
{
	if (filename.size() == 0)
		log_fatal("Empty file list provided to G3IndexedReader");

	for (auto i = filename.begin(); i != filename.end(); i++){
		boost::filesystem::path fpath(*i);
		if (i->find("://") == std::string::npos &&
		   (!boost::filesystem::exists(fpath) ||
		    !boost::filesystem::is_regular_file(fpath)))
			log_fatal("Could not find file %s", i->c_str());
		filename_.push_back(*i);
	}
	StartFile(filename_.front());
	filename_.pop_front();
}

void G3IndexedReader::StartFile(std::string path)
{
	log_info("Starting file %s\n", path.c_str());
	cur_file_ = path;
	g3_istream_from_path(stream_, path);
}

void G3IndexedReader::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	// Bail (ends processing) if too many frames have passed by
	if (!frame && n_frames_to_read_ > 0 && 
	    n_frames_read_ >= n_frames_to_read_)
		return;

	// If this is not a driving module and we haven't already, prepend the
	// file to the input frame
	if (frame) {
		if (!prefix_file_) {
			prefix_file_ = true;
			std::deque<G3FramePtr> subdeque;
			while (true) {
				Process(G3FramePtr(), subdeque);
				if (subdeque.size() == 0)
					break;
				for (auto i = subdeque.begin();
				    i != subdeque.end(); i++)
					out.push_back(*i);
				subdeque.clear();
			}
		}
		out.push_back(frame);
	}

	if (stream_.peek() == EOF) {
		if (filename_.size() > 0) {
			StartFile(filename_.front());
			filename_.pop_front();
		} else {
			// Stop processing
			return;
		}
	}
	frame = G3FramePtr(new G3Frame);
	try {
		frame->load(stream_);
	} catch (...) {
		log_error("Exception raised while reading file %s",
		    cur_file_.c_str());
		throw;
	}
	out.push_back(frame);
	n_frames_read_++;
}

int G3IndexedReader::Seek(int offset) {
    if (stream_.peek() == EOF && offset != Tell()) {
        log_error("Cannot seek; stream closed at EOF.");
        throw RuntimeError_exception("Cannot seek; stream closed at EOF.");
    }
    return boost::iostreams::seek(stream_, offset, std::ios_base::beg);
}

int G3IndexedReader::Tell() {
    return boost::iostreams::seek(stream_, 0, std::ios_base::cur);
}

PYBINDINGS("so3g") {
	using namespace boost::python;

	// Instead of EXPORT_G3MODULE since there are two constructors
	class_<G3IndexedReader, bases<G3Module>, boost::shared_ptr<G3IndexedReader>,
	    boost::noncopyable>("G3IndexedReader",
	      "Read frames from disk. Takes either the path to a file to read "
	      "or an iterable of files to be read in sequence. If "
	      "n_frames_to_read is greater than zero, will stop after "
	      "n_frames_to_read frames rather than at the end of the file[s].",
	    init<std::string, int>((arg("filename"),
	      arg("n_frames_to_read")=0)))
		.def(init<std::vector<std::string>, int>((arg("filename"), 
                                                          arg("n_frames_to_read")=0)))
            .def("Tell", &G3IndexedReader::Tell,
                 "Return the current byte offset from start of stream.")
            .def("Seek", &G3IndexedReader::Seek,
                 "Position the stream read pointer at specific byte offset. "
                 "Note that once EOF is reached, Seek does not work any more.")
		.def_readonly("__g3module__", true)
	;
}

