#pragma once

#include "so3g_numpy.h"
#include <G3Frame.h>
#include <G3Map.h>

#include <exception>
#include <stdint.h>

using namespace std;

class G3SuperTimestream : public G3FrameObject {
	// Storage for a 2d array with shape (n_dets, n_samps) along
	// with the vector of associated timestamps (n_samps), and the
	// vector of channel names (n_dets).  Serializes with lossless
	// compression of int32 and int64, and controllable
	// quantization and compression of float32 and lfoat64.
public:
	G3SuperTimestream();
	~G3SuperTimestream();

	G3VectorTime times;
	G3VectorString names;

	string Description() const;
	string Summary() const;

	bool Encode();
	bool Decode();
	int Options(int data_algo=-1, int times_algo=-1, float precision=-1.);

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;

	struct array_desc {
		npy_intp type_num;
		npy_intp ndim;
		npy_intp shape[32];
		npy_intp nbytes;
	};

	struct flac_block {
		float precision;
		int size;
		char *buf;
		int count;
		vector<int> offsets;
		vector<int32_t> pivots;
		vector<int> warnings;
	};

	enum algos {
		ALGO_NONE = 0,
		ALGO_DO_FLAC = (1 << 0),
		ALGO_DO_BZ = (1 << 1)
	};

	struct {
		int8_t times_algo;
		int8_t data_algo;
		float precision;
	} options;

	int is_encoded;
	struct array_desc desc;
	PyArrayObject *array;
	struct flac_block *flac;
};

// This specialization tells cereal to use G3SuperTimestream::load/save
// and not the base class' load/save.
namespace cereal {
	template <class A> struct specialize<
		A, G3SuperTimestream, cereal::specialization::member_load_save> {};
}

G3_POINTERS(G3SuperTimestream);
G3_SERIALIZABLE(G3SuperTimestream, 0);

class g3supertimestream_exception : std::exception
{
	// Exception raised when internal validity checks fail.  This will
	// also be mapped to some particular Python exception type.
public:
	std::string text;
	g3supertimestream_exception(std::string text) :	text{text} {};

	std::string msg_for_python() const throw() {
		return text;
	}
};
