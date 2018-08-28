#ifndef _G3_INDREADER_H
#define _G3_INDREADER_H

#include <string>
#include <vector>
#include <boost/iostreams/filtering_stream.hpp>

#include <G3Module.h>

class G3IndexedReader : public G3Module {
public:
	G3IndexedReader(std::string filename, int n_frames_to_read = -1);
	G3IndexedReader(std::vector<std::string> filenames, int n_frames_to_read = -1);

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
        int Seek(int offset);
        int Tell();
private:
	void StartFile(std::string path);
	bool prefix_file_;
	std::string cur_file_;
	std::deque<std::string> filename_;
	boost::iostreams::filtering_istream stream_;
	int n_frames_to_read_;
	int n_frames_read_;

	SET_LOGGER("G3IndexedReader");
};

G3_POINTERS(G3IndexedReader);

#endif
