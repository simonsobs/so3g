# This file is generated automatically by scanning a compiled
# C++ boost-python module, to facilitate documentation builds.
# It should not be present in the master branch!
# Edits will be lost.

# No handler for "__name__"

# No handler for "__doc__"

# No handler for "__package__"

# No handler for "__loader__"

# No handler for "__spec__"

# No handler for "__path__"

# No handler for "__file__"

# No handler for "__cached__"

# No handler for "__builtins__"

# No handler for "os"

# No handler for "np"

# No handler for "npversion"

# No handler for "npversplit"

# No handler for "spt3g"

# No handler for "spt3g_core"

# No handler for "load_pybindings"

def version():
  return "0.1.17"
# No handler for "useful_info"

# No handler for "greet"

class TestClass:
  """None"""
  pass
  def runme():
    """runme() -> None

    C++ signature :
        void runme(TestClass {lvalue})"""
    pass

class TestFrame:
  """TestFrame for demonstration."""
  pass

class IrregBlockDouble:
  """Data block for irregularly sampled data."""
  pass

class HKFrameType:
  """Identifier for generic HK streams."""
  pass

class G3Ndarray:
  """G3Ndarray default constructor"""
  pass
  def to_array():
    """to_array() -> object
    Get the wrapped numpy array

    C++ signature :
        boost::python::api::object to_array(G3Ndarray {lvalue})"""
    pass

class G3WCS:
  """G3WCS default constructor"""
  pass

class G3Ndmap:
  """G3Ndmap default constructor"""
  pass

class G3IndexedReader:
  """Read frames from disk. Takes either the path to a file to read or an iterable of files to be read in sequence. If n_frames_to_read is greater than zero, will stop after n_frames_to_read frames rather than at the end of the file[s]."""
  pass
  def Tell():
    """Tell() -> int
    Return the current byte offset from start of stream.

    C++ signature :
        int Tell(G3IndexedReader {lvalue})"""
    pass
  def Seek():
    """Seek(arg2) -> int
    Position the stream read pointer at specific byte offset. Note that once EOF is reached, Seek does not work any more.

    C++ signature :
        int Seek(G3IndexedReader {lvalue},int)"""
    pass

class IntervalsDouble:
  """A finite series of non-overlapping semi-open intervals on a domain of type: double."""
  pass
  def add_interval():
    """add_interval(start, end) -> IntervalsDouble
    Merge an interval into the set."""
    pass
  def append_interval_no_check():
    """append_interval_no_check(start, end) -> IntervalsDouble
    Append an interval to the set without checking for overlap or sequence."""
    pass
  def merge():
    """merge(arg2) -> IntervalsDouble
    Merge an Intervals into the set."""
    pass
  def intersect():
    """intersect(source) -> IntervalsDouble
    Intersect another doublewith this one."""
    pass
  def complement():
    """complement() -> IntervalsDouble
    Return the complement (over domain)."""
    pass
  def array():
    """array() -> object
    Return the intervals as a 2-d numpy array."""
    pass
  @staticmethod
  def from_array():
    """from_array(input_array) -> IntervalsDouble
    Return a IntervalsDouble based on an (n,2) ndarray."""
    pass
  @staticmethod
  def from_mask():
    """from_mask(input_array, n_bits) -> object
    Return a list of IntervalsDouble, extracted from the first 
    n_bits bits of input_array (a 1-d array of integer type)."""
    pass
  @staticmethod
  def mask():
    """mask(intervals_list, n_bits) -> object
    Return an ndarray bitmask from a list of IntervalsDouble.
    The dtype will be the smallest available to hold n_bits."""
    pass
  def copy():
    """copy() -> IntervalsDouble
    Get a new object with a copy of the data."""
    pass

class std_map_indexing_suite__MapIntervalsDoubleBaseMap_entry:
  """None"""
  pass
  def data():
    """data() -> IntervalsDouble
    K.data() -> the value associated with this pair.
    """
    pass
  def key():
    """key() -> str
    K.key() -> the key associated with this pair.
    """
    pass
  def first():
    """first() -> str
    K.first() -> the first item in this pair.
    """
    pass
  def second():
    """second() -> IntervalsDouble
    K.second() -> the second item in this pair.
    """
    pass

class MapIntervalsDouble:
  """Mapping from strings to Intervals over double."""
  pass
  def keys():
    """keys() -> list
    D.keys() -> list of D's keys
    """
    pass
  def has_key():
    """has_key(arg2) -> bool
    D.has_key(k) -> True if D has a key k, else False
    """
    pass
  def values():
    """values() -> list
    D.values() -> list of D's values
    """
    pass
  def items():
    """items() -> list
    D.items() -> list of D's (key, value) pairs, as 2-tuples
    """
    pass
  def clear():
    """clear() -> None
    D.clear() -> None.  Remove all items from D.
    """
    pass
  def copy():
    """copy() -> MapIntervalsDouble
    D.copy() -> a shallow copy of D
    """
    pass
  def get():
    """get(arg2, default_val) -> object
    D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
    """
    pass
  def pop():
    """pop( (MapIntervalsDouble)arg1, (str)arg2) -> object

pop( (MapIntervalsDouble)arg1, (str)arg2, (object)arg3) -> object :
    D.pop(k[,d]) -> v, remove specified key and return the corresponding value
    If key is not found, d is returned if given, otherwise KeyError is raised
    """
    pass
  def popitem():
    """popitem() -> object
    D.popitem() -> (k, v), remove and return some (key, value) pair as a
    2-tuple; but raise KeyError if D is empty
    """
    pass
  @staticmethod
  def fromkeys():
    """fromkeys(arg1, arg2) -> object
    MapIntervalsDouble.fromkeys(S,v) -> New MapIntervalsDouble with keys from S and values equal to v.
    """
    pass
  def update():
    """update(arg2) -> None
    D.update(E) -> None.  Update D from E: for k in E: D[k] = E[k]
    """
    pass
  def iteritems():
    """iteritems() -> object
    D.iteritems() -> an iterator over the (key, value) items of D
    """
    pass
  def iterkeys():
    """iterkeys() -> object
    D.iterkeys() -> an iterator over the keys of D
    """
    pass
  def itervalues():
    """itervalues() -> object
    D.itervalues() -> an iterator over the values of D
    """
    pass

class IntervalsInt:
  """A finite series of non-overlapping semi-open intervals on a domain of type: int64_t."""
  pass
  def add_interval():
    """add_interval(start, end) -> IntervalsInt
    Merge an interval into the set."""
    pass
  def append_interval_no_check():
    """append_interval_no_check(start, end) -> IntervalsInt
    Append an interval to the set without checking for overlap or sequence."""
    pass
  def merge():
    """merge(arg2) -> IntervalsInt
    Merge an Intervals into the set."""
    pass
  def intersect():
    """intersect(source) -> IntervalsInt
    Intersect another int64_twith this one."""
    pass
  def complement():
    """complement() -> IntervalsInt
    Return the complement (over domain)."""
    pass
  def array():
    """array() -> object
    Return the intervals as a 2-d numpy array."""
    pass
  @staticmethod
  def from_array():
    """from_array(input_array) -> IntervalsInt
    Return a IntervalsInt based on an (n,2) ndarray."""
    pass
  @staticmethod
  def from_mask():
    """from_mask(input_array, n_bits) -> object
    Return a list of IntervalsInt, extracted from the first 
    n_bits bits of input_array (a 1-d array of integer type)."""
    pass
  @staticmethod
  def mask():
    """mask(intervals_list, n_bits) -> object
    Return an ndarray bitmask from a list of IntervalsInt.
    The dtype will be the smallest available to hold n_bits."""
    pass
  def copy():
    """copy() -> IntervalsInt
    Get a new object with a copy of the data."""
    pass

class std_map_indexing_suite__MapIntervalsIntBaseMap_entry:
  """None"""
  pass
  def data():
    """data() -> IntervalsInt
    K.data() -> the value associated with this pair.
    """
    pass
  def key():
    """key() -> str
    K.key() -> the key associated with this pair.
    """
    pass
  def first():
    """first() -> str
    K.first() -> the first item in this pair.
    """
    pass
  def second():
    """second() -> IntervalsInt
    K.second() -> the second item in this pair.
    """
    pass

class MapIntervalsInt:
  """Mapping from strings to Intervals over int64_t."""
  pass
  def keys():
    """keys() -> list
    D.keys() -> list of D's keys
    """
    pass
  def has_key():
    """has_key(arg2) -> bool
    D.has_key(k) -> True if D has a key k, else False
    """
    pass
  def values():
    """values() -> list
    D.values() -> list of D's values
    """
    pass
  def items():
    """items() -> list
    D.items() -> list of D's (key, value) pairs, as 2-tuples
    """
    pass
  def clear():
    """clear() -> None
    D.clear() -> None.  Remove all items from D.
    """
    pass
  def copy():
    """copy() -> MapIntervalsInt
    D.copy() -> a shallow copy of D
    """
    pass
  def get():
    """get(arg2, default_val) -> object
    D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
    """
    pass
  def pop():
    """pop( (MapIntervalsInt)arg1, (str)arg2) -> object

pop( (MapIntervalsInt)arg1, (str)arg2, (object)arg3) -> object :
    D.pop(k[,d]) -> v, remove specified key and return the corresponding value
    If key is not found, d is returned if given, otherwise KeyError is raised
    """
    pass
  def popitem():
    """popitem() -> object
    D.popitem() -> (k, v), remove and return some (key, value) pair as a
    2-tuple; but raise KeyError if D is empty
    """
    pass
  @staticmethod
  def fromkeys():
    """fromkeys(arg1, arg2) -> object
    MapIntervalsInt.fromkeys(S,v) -> New MapIntervalsInt with keys from S and values equal to v.
    """
    pass
  def update():
    """update(arg2) -> None
    D.update(E) -> None.  Update D from E: for k in E: D[k] = E[k]
    """
    pass
  def iteritems():
    """iteritems() -> object
    D.iteritems() -> an iterator over the (key, value) items of D
    """
    pass
  def iterkeys():
    """iterkeys() -> object
    D.iterkeys() -> an iterator over the keys of D
    """
    pass
  def itervalues():
    """itervalues() -> object
    D.itervalues() -> an iterator over the values of D
    """
    pass

class IntervalsInt32:
  """A finite series of non-overlapping semi-open intervals on a domain of type: int32_t."""
  pass
  def add_interval():
    """add_interval(start, end) -> IntervalsInt32
    Merge an interval into the set."""
    pass
  def append_interval_no_check():
    """append_interval_no_check(start, end) -> IntervalsInt32
    Append an interval to the set without checking for overlap or sequence."""
    pass
  def merge():
    """merge(arg2) -> IntervalsInt32
    Merge an Intervals into the set."""
    pass
  def intersect():
    """intersect(source) -> IntervalsInt32
    Intersect another int32_twith this one."""
    pass
  def complement():
    """complement() -> IntervalsInt32
    Return the complement (over domain)."""
    pass
  def array():
    """array() -> object
    Return the intervals as a 2-d numpy array."""
    pass
  @staticmethod
  def from_array():
    """from_array(input_array) -> IntervalsInt32
    Return a IntervalsInt32 based on an (n,2) ndarray."""
    pass
  @staticmethod
  def from_mask():
    """from_mask(input_array, n_bits) -> object
    Return a list of IntervalsInt32, extracted from the first 
    n_bits bits of input_array (a 1-d array of integer type)."""
    pass
  @staticmethod
  def mask():
    """mask(intervals_list, n_bits) -> object
    Return an ndarray bitmask from a list of IntervalsInt32.
    The dtype will be the smallest available to hold n_bits."""
    pass
  def copy():
    """copy() -> IntervalsInt32
    Get a new object with a copy of the data."""
    pass

class std_map_indexing_suite__MapIntervalsInt32BaseMap_entry:
  """None"""
  pass
  def data():
    """data() -> IntervalsInt32
    K.data() -> the value associated with this pair.
    """
    pass
  def key():
    """key() -> str
    K.key() -> the key associated with this pair.
    """
    pass
  def first():
    """first() -> str
    K.first() -> the first item in this pair.
    """
    pass
  def second():
    """second() -> IntervalsInt32
    K.second() -> the second item in this pair.
    """
    pass

class MapIntervalsInt32:
  """Mapping from strings to Intervals over int32_t."""
  pass
  def keys():
    """keys() -> list
    D.keys() -> list of D's keys
    """
    pass
  def has_key():
    """has_key(arg2) -> bool
    D.has_key(k) -> True if D has a key k, else False
    """
    pass
  def values():
    """values() -> list
    D.values() -> list of D's values
    """
    pass
  def items():
    """items() -> list
    D.items() -> list of D's (key, value) pairs, as 2-tuples
    """
    pass
  def clear():
    """clear() -> None
    D.clear() -> None.  Remove all items from D.
    """
    pass
  def copy():
    """copy() -> MapIntervalsInt32
    D.copy() -> a shallow copy of D
    """
    pass
  def get():
    """get(arg2, default_val) -> object
    D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
    """
    pass
  def pop():
    """pop( (MapIntervalsInt32)arg1, (str)arg2) -> object

pop( (MapIntervalsInt32)arg1, (str)arg2, (object)arg3) -> object :
    D.pop(k[,d]) -> v, remove specified key and return the corresponding value
    If key is not found, d is returned if given, otherwise KeyError is raised
    """
    pass
  def popitem():
    """popitem() -> object
    D.popitem() -> (k, v), remove and return some (key, value) pair as a
    2-tuple; but raise KeyError if D is empty
    """
    pass
  @staticmethod
  def fromkeys():
    """fromkeys(arg1, arg2) -> object
    MapIntervalsInt32.fromkeys(S,v) -> New MapIntervalsInt32 with keys from S and values equal to v.
    """
    pass
  def update():
    """update(arg2) -> None
    D.update(E) -> None.  Update D from E: for k in E: D[k] = E[k]
    """
    pass
  def iteritems():
    """iteritems() -> object
    D.iteritems() -> an iterator over the (key, value) items of D
    """
    pass
  def iterkeys():
    """iterkeys() -> object
    D.iterkeys() -> an iterator over the keys of D
    """
    pass
  def itervalues():
    """itervalues() -> object
    D.itervalues() -> an iterator over the values of D
    """
    pass

class IntervalsTime:
  """A finite series of non-overlapping semi-open intervals on a domain of type: G3Time."""
  pass
  def add_interval():
    """add_interval(start, end) -> IntervalsTime
    Merge an interval into the set."""
    pass
  def append_interval_no_check():
    """append_interval_no_check(start, end) -> IntervalsTime
    Append an interval to the set without checking for overlap or sequence."""
    pass
  def merge():
    """merge(arg2) -> IntervalsTime
    Merge an Intervals into the set."""
    pass
  def intersect():
    """intersect(source) -> IntervalsTime
    Intersect another G3Timewith this one."""
    pass
  def complement():
    """complement() -> IntervalsTime
    Return the complement (over domain)."""
    pass
  def array():
    """array() -> object
    Return the intervals as a 2-d numpy array."""
    pass
  @staticmethod
  def from_array():
    """from_array(input_array) -> IntervalsTime
    Return a IntervalsTime based on an (n,2) ndarray."""
    pass
  @staticmethod
  def from_mask():
    """from_mask(input_array, n_bits) -> object
    Return a list of IntervalsTime, extracted from the first 
    n_bits bits of input_array (a 1-d array of integer type)."""
    pass
  @staticmethod
  def mask():
    """mask(intervals_list, n_bits) -> object
    Return an ndarray bitmask from a list of IntervalsTime.
    The dtype will be the smallest available to hold n_bits."""
    pass
  def copy():
    """copy() -> IntervalsTime
    Get a new object with a copy of the data."""
    pass

class std_map_indexing_suite__MapIntervalsTimeBaseMap_entry:
  """None"""
  pass
  def data():
    """data() -> IntervalsTime
    K.data() -> the value associated with this pair.
    """
    pass
  def key():
    """key() -> str
    K.key() -> the key associated with this pair.
    """
    pass
  def first():
    """first() -> str
    K.first() -> the first item in this pair.
    """
    pass
  def second():
    """second() -> IntervalsTime
    K.second() -> the second item in this pair.
    """
    pass

class MapIntervalsTime:
  """Mapping from strings to Intervals over G3Time."""
  pass
  def keys():
    """keys() -> list
    D.keys() -> list of D's keys
    """
    pass
  def has_key():
    """has_key(arg2) -> bool
    D.has_key(k) -> True if D has a key k, else False
    """
    pass
  def values():
    """values() -> list
    D.values() -> list of D's values
    """
    pass
  def items():
    """items() -> list
    D.items() -> list of D's (key, value) pairs, as 2-tuples
    """
    pass
  def clear():
    """clear() -> None
    D.clear() -> None.  Remove all items from D.
    """
    pass
  def copy():
    """copy() -> MapIntervalsTime
    D.copy() -> a shallow copy of D
    """
    pass
  def get():
    """get(arg2, default_val) -> object
    D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
    """
    pass
  def pop():
    """pop( (MapIntervalsTime)arg1, (str)arg2) -> object

pop( (MapIntervalsTime)arg1, (str)arg2, (object)arg3) -> object :
    D.pop(k[,d]) -> v, remove specified key and return the corresponding value
    If key is not found, d is returned if given, otherwise KeyError is raised
    """
    pass
  def popitem():
    """popitem() -> object
    D.popitem() -> (k, v), remove and return some (key, value) pair as a
    2-tuple; but raise KeyError if D is empty
    """
    pass
  @staticmethod
  def fromkeys():
    """fromkeys(arg1, arg2) -> object
    MapIntervalsTime.fromkeys(S,v) -> New MapIntervalsTime with keys from S and values equal to v.
    """
    pass
  def update():
    """update(arg2) -> None
    D.update(E) -> None.  Update D from E: for k in E: D[k] = E[k]
    """
    pass
  def iteritems():
    """iteritems() -> object
    D.iteritems() -> an iterator over the (key, value) items of D
    """
    pass
  def iterkeys():
    """iterkeys() -> object
    D.iterkeys() -> an iterator over the keys of D
    """
    pass
  def itervalues():
    """itervalues() -> object
    D.itervalues() -> an iterator over the values of D
    """
    pass

class BFilterParams:
  """None"""
  pass

class BFilterBank:
  """None"""
  pass
  def add():
    """add(arg2) -> BFilterBank

    C++ signature :
        BFilterBank {lvalue} add(BFilterBank {lvalue},BFilterParams)"""
    pass
  def init():
    """init(arg2) -> BFilterBank

    C++ signature :
        BFilterBank {lvalue} init(BFilterBank {lvalue},int)"""
    pass
  def apply():
    """apply(arg2, arg3) -> None

    C++ signature :
        void apply(BFilterBank {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass

class RangesInt32:
  """A finite series of non-overlapping semi-open intervals on a domain
of type: int32_t.

To create an empty object, instantiate with just a sample count:
``RangesInt32(count)``.

Alternately, consider convenience methods such as ``from_mask``,
``from_array``, and ``from_bitmask``; see below.

In addition to the methods explained below, note the that following
operators have been defined and perform as follows (where ``r1`` and
``r2`` are objects of this class:

- ``~r1`` is equivalent to ``r1.complement()``
- ``r1 *= r2`` is equivalent to ``r1.intersect(r2)``
- ``r1 += r2`` is equivalent to ``r1.merge(r2)``
- ``r1 * r2`` and ``r1 + r2`` behave as you might expect, returning a
  new object and leaving ``r1`` and ``r2`` unmodified.

The object also supports slicing.  For example, if ``r1`` has
count = 100 then r1[10:-5] returns a new object (not a view)
that has count = 85.  A data member ``reference`` keeps track
of the history of shifts in the first sample; for example if
r1.reference = 0 then r1[10:-5].reference will be -10.  This
variable can be interpreted as giving the logical index, in
the new index system, of where index=0 of the original object
would be found.  This is useful for bookkeeping in some cases.
"""
  pass
  def add_interval():
    """add_interval(start, end) -> RangesInt32
    Merge an interval into the set."""
    pass
  def append_interval_no_check():
    """append_interval_no_check(start, end) -> RangesInt32
    Append an interval to the set without checking for overlap or sequence."""
    pass
  def merge():
    """merge(src) -> RangesInt32
    Merge ranges from another RangesInt32 into this one."""
    pass
  def buffer():
    """buffer(buff) -> RangesInt32
    Buffer each interval by an amount specified by buff"""
    pass
  def buffered():
    """buffered(buff) -> RangesInt32
    Return an interval buffered by buff"""
    pass
  def close_gaps():
    """close_gaps(gap) -> RangesInt32
    Remove gaps between ranges less than gap"""
    pass
  def intersect():
    """intersect(src) -> RangesInt32
    Intersect another RangesInt32 with this one."""
    pass
  def complement():
    """complement() -> RangesInt32
    Return the complement (over domain)."""
    pass
  def zeros_like():
    """zeros_like() -> RangesInt32
    Return range of same length but no intervals"""
    pass
  def ones_like():
    """ones_like() -> RangesInt32
    Return range of same length and interval spanning count"""
    pass
  def ranges():
    """ranges() -> object
    Return the intervals as a 2-d numpy array of ranges."""
    pass
  @staticmethod
  def from_array():
    """from_array(data, count) -> RangesInt32
    The input data must be an (n,2) shape ndarray of int32. The integer count sets the domain of the object."""
    pass
  @staticmethod
  def from_bitmask():
    """from_bitmask(bitmask_array) -> object
    Return a list of RangesInt32 extracted from an ndarray encoding a bitmask."""
    pass
  @staticmethod
  def from_mask():
    """from_mask(bool_array) -> object
    Return a list of RangesInt32 extracted from an ndarray of bool."""
    pass
  @staticmethod
  def bitmask():
    """bitmask(ranges_list, n_bits) -> object
    Return an ndarray bitmask from a list ofRangesInt32.
    n_bits determines the output integer type.  Bits are assigned from 
    LSB onwards; use None in the list to skip a bit."""
    pass
  def mask():
    """mask() -> object
    Return a boolean mask from this Ranges object."""
    pass
  def copy():
    """copy() -> RangesInt32
    Get a new object with a copy of the data."""
    pass

class std_map_indexing_suite__MapRangesInt32BaseMap_entry:
  """None"""
  pass
  def data():
    """data() -> RangesInt32
    K.data() -> the value associated with this pair.
    """
    pass
  def key():
    """key() -> str
    K.key() -> the key associated with this pair.
    """
    pass
  def first():
    """first() -> str
    K.first() -> the first item in this pair.
    """
    pass
  def second():
    """second() -> RangesInt32
    K.second() -> the second item in this pair.
    """
    pass

class MapRangesInt32:
  """Mapping from strings to Ranges over int32_t."""
  pass
  def keys():
    """keys() -> list
    D.keys() -> list of D's keys
    """
    pass
  def has_key():
    """has_key(arg2) -> bool
    D.has_key(k) -> True if D has a key k, else False
    """
    pass
  def values():
    """values() -> list
    D.values() -> list of D's values
    """
    pass
  def items():
    """items() -> list
    D.items() -> list of D's (key, value) pairs, as 2-tuples
    """
    pass
  def clear():
    """clear() -> None
    D.clear() -> None.  Remove all items from D.
    """
    pass
  def copy():
    """copy() -> MapRangesInt32
    D.copy() -> a shallow copy of D
    """
    pass
  def get():
    """get(arg2, default_val) -> object
    D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
    """
    pass
  def pop():
    """pop( (MapRangesInt32)arg1, (str)arg2) -> object

pop( (MapRangesInt32)arg1, (str)arg2, (object)arg3) -> object :
    D.pop(k[,d]) -> v, remove specified key and return the corresponding value
    If key is not found, d is returned if given, otherwise KeyError is raised
    """
    pass
  def popitem():
    """popitem() -> object
    D.popitem() -> (k, v), remove and return some (key, value) pair as a
    2-tuple; but raise KeyError if D is empty
    """
    pass
  @staticmethod
  def fromkeys():
    """fromkeys(arg1, arg2) -> object
    MapRangesInt32.fromkeys(S,v) -> New MapRangesInt32 with keys from S and values equal to v.
    """
    pass
  def update():
    """update(arg2) -> None
    D.update(E) -> None.  Update D from E: for k in E: D[k] = E[k]
    """
    pass
  def iteritems():
    """iteritems() -> object
    D.iteritems() -> an iterator over the (key, value) items of D
    """
    pass
  def iterkeys():
    """iterkeys() -> object
    D.iterkeys() -> an iterator over the keys of D
    """
    pass
  def itervalues():
    """itervalues() -> object
    D.itervalues() -> an iterator over the values of D
    """
    pass

class RebundlerPrimaryMap:
  """None"""
  pass
  def Process():
    """Process(arg2) -> object
    Add element.

    C++ signature :
        boost::python::api::object Process(Rebundler<G3TimestreamMap> {lvalue},boost::shared_ptr<G3FrameObject>)"""
    pass
  def ExtractIntervalTime():
    """ExtractIntervalTime(arg2, arg3, arg4) -> object
    Rebundle into interval.

    C++ signature :
        boost::python::api::object ExtractIntervalTime(Rebundler<G3TimestreamMap> {lvalue},G3Time,G3Time,bool)"""
    pass
  def ExtractInterval():
    """ExtractInterval(arg2, arg3) -> object
    Rebundle into interval.

    C++ signature :
        boost::python::api::object ExtractInterval(Rebundler<G3TimestreamMap> {lvalue},int,bool)"""
    pass

class ProjEng_Flat_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_T_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_T_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_QU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_QU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_TQU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Flat_TQU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjFlat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_T_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_T_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_QU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_QU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_TQU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Quat_TQU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjQuat, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_T_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_T_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_QU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_QU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_TQU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CAR_TQU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_T_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_T_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_QU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_QU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_TQU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_CEA_TQU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_T_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_T_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_QU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_QU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_TQU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ARC_TQU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjARC, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_T_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_T_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_QU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_QU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_TQU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_TAN_TQU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjTAN, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_T_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_T_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_QU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_QU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_TQU_Tiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<Tiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, NearestNeighbor>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_ZEA_TQU_NonTiled_Bilinear:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjZEA, Pixelizor2_Flat<NonTiled, Bilinear>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Precomp_NonTiled:
  """None"""
  pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object from_map(ProjEng_Precomp<NonTiled> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_map(ProjEng_Precomp<NonTiled> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjEng_Precomp<NonTiled> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_Precomp_Tiled:
  """None"""
  pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object from_map(ProjEng_Precomp<Tiled> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_map(ProjEng_Precomp<Tiled> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjEng_Precomp<Tiled> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_HP_T_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_HP_QU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_HP_TQU_NonTiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<NonTiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_HP_T_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinT> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_HP_QU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class ProjEng_HP_TQU_Tiled:
  """None"""
  pass
  def coords():
    """coords(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object coords(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixels():
    """pixels(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object pixels(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_hits():
    """tile_hits(arg2, arg3) -> IntVector

    C++ signature :
        std::vector<int, std::allocator<int> > tile_hits(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object)"""
    pass
  def tile_ranges():
    """tile_ranges(arg2, arg3, arg4) -> object

    C++ signature :
        boost::python::api::object tile_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pointing_matrix():
    """pointing_matrix(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object pointing_matrix(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def pixel_ranges():
    """pixel_ranges(arg2, arg3, arg4, arg5) -> object

    C++ signature :
        boost::python::api::object pixel_ranges(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,int)"""
    pass
  def zeros():
    """zeros(arg2) -> object

    C++ signature :
        boost::python::api::object zeros(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object)"""
    pass
  def from_map():
    """from_map(arg2, arg3, arg4, arg5, arg6) -> object

    C++ signature :
        boost::python::api::object from_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_map():
    """to_map(arg2, arg3, arg4, arg5, arg6, arg7, arg8) -> object

    C++ signature :
        boost::python::api::object to_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass
  def to_weight_map():
    """to_weight_map(arg2, arg3, arg4, arg5, arg6, arg7) -> object

    C++ signature :
        boost::python::api::object to_weight_map(ProjectionEngine<ProjCAR, Pixelizor_Healpix<Tiled>, SpinTQU> {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object)"""
    pass

class G3SuperTimestream:
  """G3SuperTimestream()

  Construct with names and times uninitialized."""
  pass
  def encode():
    """Compress."""
    pass
  def decode():
    """Decompress."""
    pass
  def extract():
    """Decompress data subset into an array."""
    pass
  def calibrate():
    """calibrate(cal_factors)

Apply per-channel scale factors.  Note this puts you in float mode
(if you're not already) and modifies both .quanta and .data."""
    pass
  def options():
    """Set compression options."""
    pass

# No handler for "test_g3super"

# No handler for "test_trig"

# No handler for "nmat_detvecs_apply"

# No handler for "process_cuts"

# No handler for "translate_cuts"

# No handler for "get_gap_fill_poly"

# No handler for "get_gap_fill_poly64"

# No handler for "test_buffer_wrapper"

# No handler for "block_moment"

# No handler for "block_moment64"

# No handler for "block_minmax"

# No handler for "block_minmax64"

# No handler for "matched_jumps"

# No handler for "matched_jumps64"

# No handler for "find_quantized_jumps"

# No handler for "find_quantized_jumps64"

# No handler for "subtract_jump_heights"

# No handler for "subtract_jump_heights64"

# No handler for "clean_flag"

# No handler for "interp1d_linear"

# No handler for "detrend"

# No handler for "__version__"

# No handler for "config"

# No handler for "instance_config"

# No handler for "soframe"

# No handler for "hk"

# No handler for "proj"

