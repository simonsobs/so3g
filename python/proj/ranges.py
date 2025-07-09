import so3g
import numpy as np

"""Objects will self report as being of type "RangesInt32" rather than
Ranges.  But let's try to use so3g.proj.Ranges when testing types and
making new ones and stuff."""

Ranges = so3g.RangesInt32


class RangesMatrix():
    """This is a wrapper for multi-dimensional grids of Ranges objects.
    This can be used to store Ranges objects per-detector (the 2-d
    case) or per-thread and per-detector (a 3-d case, required if
    using OpenMP).  The right-most axis always corresponds to the
    sample axis carried by the Ranges objects.

    In addition to .shape and multi-dimensional slicing, it supports
    inversion (the ~ operator) and multiplication against other
    RangesMatrix or Ranges objects, with broadcasting rules similar to
    standard numpy array rules.

    """
    def __init__(self, items=[], child_shape=None, skip_shape_check=False):
        self.ranges = [x for x in items]
        if len(items):
            child_shape = items[0].shape
            if not skip_shape_check:
                assert all([(item.shape == child_shape) for item in items])
        elif child_shape is None:
            child_shape = ()
        self._child_shape = child_shape

    def __repr__(self):
        return 'RangesMatrix(' + ','.join(map(str, self.shape)) + ')'

    def __len__(self):
        return self.shape[0]

    def copy(self):
        return RangesMatrix([x.copy() for x in self.ranges],
                            child_shape=self.shape[1:])

    def zeros_like(self):
        return RangesMatrix([x.zeros_like() for x in self.ranges],
                            child_shape=self.shape[1:])
    
    def ones_like(self):
        return RangesMatrix([x.ones_like() for x in self.ranges],
                            child_shape=self.shape[1:])

    def buffer(self, buff):
        [x.buffer(buff) for x in self.ranges]
        ## just to make this work like Ranges.buffer()
        return self
    
    def buffered(self, buff):
        out = self.copy()
        [x.buffer(buff) for x in out.ranges]
        return out

    @property
    def shape(self):
        if len(self.ranges) == 0:
            return (0,) + self._child_shape
        return (len(self.ranges),) + self.ranges[0].shape

    def __getitem__(self, index):
        """RangesMatrix supports multi-dimensional indexing, slicing, and
        numpy-style advanced indexing (with some restrictions).

        The right-most axis of RangesMatrix has special restrictions,
        since it corresponds to Ranges objects: it can only accept
        slice-based indexing, not integer or advanced indexing.

        To guarantee that you get copies, rather than references to,
        the lowest level Ranges objects, make sure the index tuple
        includes a slice along the final (right-most) axis.  For
        example::

          # Starting point
          rm = RangesMatrix.zeros(10, 10000)

          # Get and modify an element in rm ...
          r = rm[0]
          r.add_interval(0, 10)  # This _does_ affect rm!

          # Get a copy of an element from rm.
          r = rm[0,:]
          r.add_interval(0, 10)  # This will not affect rm.

          # More generally, this is equivalent to rm1 = rm.copy():
          new_rm = rm[..., :]

        """
        if not isinstance(index, tuple):
            index = (index,)

        eidx = [i for i, a in enumerate(index) if a is Ellipsis]
        if len(eidx) == 1:
            # Fill in missing slices.
            eidx = eidx[0]
            n_free = len(self.shape) - sum([e is not None for e in index]) + 1
            index = index[:eidx] + tuple([slice(None)] * n_free) + index[eidx+1:]
        elif len(eidx) > 1:
            raise IndexError("An index can only have a single ellipsis ('...')")

        return _gibasic(self, index)

    def __add__(self, x):
        if isinstance(x, Ranges):
            return self.__class__([d + x for d in self.ranges])
        elif isinstance(x, RangesMatrix):
            # Check for shape compatibility.
            nd_a = len(self.shape)
            nd_b = len(x.shape)
            ndim = min(nd_a, nd_b)
            ok = [(a==1 or b==1 or a == b) for a,b in zip(self.shape[-ndim:], x.shape[-ndim:])]
            if not all(ok) or nd_a < nd_b:
                raise ValueError('Operands have incompatible shapes: %s %s' %
                                 self.shape, x.shape)
            if nd_a == nd_b:
                # Broadcast if either has shape 1...
                if x.shape[0] == 1:
                    return self.__class__([r + x[0] for r in self.ranges])
                elif self.shape[0] == 1:
                    return self.__class__([self.ranges[0] + _x for _x in x])
                elif self.shape[0] == x.shape[0]:
                    return self.__class__([r + d for r, d in zip(self.ranges, x)])
            return self.__class__([r + x for r in self.ranges])
        
    def __mul__(self, x):
        if isinstance(x, Ranges):
            return self.__class__([d * x for d in self.ranges])
        elif isinstance(x, RangesMatrix):
            # Check for shape compatibility.
            nd_a = len(self.shape)
            nd_b = len(x.shape)
            ndim = min(nd_a, nd_b)
            ok = [(a==1 or b==1 or a == b) for a,b in zip(self.shape[-ndim:], x.shape[-ndim:])]
            if not all(ok) or nd_a < nd_b:
                raise ValueError('Operands have incompatible shapes: %s %s' %
                                 self.shape, x.shape)
            if nd_a == nd_b:
                # Broadcast if either has shape 1...
                if x.shape[0] == 1:
                    return self.__class__([r * x[0] for r in self.ranges])
                elif self.shape[0] == 1:
                    return self.__class__([self.ranges[0] * _x for _x in x])
                elif self.shape[0] == x.shape[0]:
                    return self.__class__([r * d for r, d in zip(self.ranges, x)])
            return self.__class__([r * x for r in self.ranges])

    def __invert__(self):
        return self.__class__([~x for x in self.ranges])

    @staticmethod
    def concatenate(items, axis=0):
        """Static method to concatenate multiple RangesMatrix (or Ranges)
        objects along the specified axis.

        """
        # Check shape compatibility...
        s = list(items[0].shape)
        s[axis] = -1
        for item in items[1:]:
            s1 = list(item.shape)
            s1[axis] = -1
            if s != s1:
                raise ValueError('Contributed items must have same shape on non-cat axis.')
        def collect(items, join_depth):
            # Recurse through items in order to concatenate along axis join_depth.
            ranges = []
            if join_depth > 0:
                for co_items in zip(*items):
                    ranges.append(collect(co_items, join_depth - 1))
            else:
                if isinstance(items[0], RangesMatrix):
                    for item in items:
                        ranges.extend(item.ranges)
                else:
                    # The items are Ranges, then.
                    n, r = 0, []
                    for item in items:
                        r.extend(item.ranges() + n)
                        n += item.count
                    r = Ranges.from_array(
                        np.array(r, dtype='int32').reshape((-1, 2)), n)
                    return r
            return RangesMatrix(ranges, child_shape=items[0].shape[1:])
        return collect(items, axis)

    def close_gaps(self, gap_size=0):
        """Call close_gaps(gap_size) on all children.  Any ranges that abutt
        each other within the gap_size are merged into a single entry.
        Usually a gap_size of 0 is not possible, but if a caller is
        carefully using append_interval_no_check, then it can happen.

        """
        for r in self.ranges:
            r.close_gaps(gap_size)

    def get_stats(self):
        samples, intervals = [], []
        for r in self.ranges:
            if isinstance(r, Ranges):
                ra = r.ranges()
                samples.append(np.dot(ra, [-1, 1]).sum())
                intervals.append(ra.shape[0])
            else:
                sub = r.get_stats()
                samples.append(sum(sub['samples']))
                intervals.append(sum(sub['intervals']))
        return {
            'samples': samples,
            'intervals': intervals}

    @staticmethod
    def full(shape, fill_value):
        """Construct a RangesMatrix with the specified shape, initialized to
        fill_value.

        Args:
          shape (tuple of int): The shape.  Must have at least one
            element.  If there is only one element, a Ranges object is
            returned, not a RangesMatrix.
          fill_value (bool): True or False.  If False, the wrapped
            Ranges objects will have no intervals.  If True, the
            wrapped Ranges objects will all have a single interval
            spanning their entire domain.

        Returns:
          Ranges object (if shape has a single element) or a
          RangesMatrix object (len(shape) > 1).

        See also: zeros, ones.

        """
        assert fill_value in [True, False]
        if isinstance(shape, int):
            shape = (shape,)
        assert(len(shape) > 0)
        if len(shape) == 1:
            r = Ranges(shape[0])
            if fill_value:
                r = ~r
            return r
        return RangesMatrix([RangesMatrix.full(shape[1:], fill_value)
                             for i in range(shape[0])],
                            child_shape=shape[1:])

    @classmethod
    def zeros(cls, shape):
        """Equivalent to full(shape, False)."""
        return cls.full(shape, False)

    @classmethod
    def ones(cls, shape):
        """Equivalent to full(shape, True)."""
        return cls.full(shape, True)

    @classmethod
    def from_mask(cls, mask):
        """Create a RangesMatrix from a boolean mask.  The input mask can have
        any dimensionality greater than 0 but be aware that if ndim==1
        then a Ranges object, and not a RangesMatrix, is returned.

        Args:
          mask (ndarray): Must be boolean array with at least 1 dimension.

        Returns:
          RangesMatrix with the same shape as mask, with ranges
          corresponding to the intervals where mask was True.

        """
        assert(mask.ndim > 0)
        if mask.ndim == 1:
           return Ranges.from_mask(mask)
        if len(mask) == 0:
            return cls(child_shape=mask.shape[1:])
        # Recurse.
        return cls([cls.from_mask(m) for m in mask])

    def mask(self, dest=None):
        """Return the boolean mask equivalent of this object."""
        if dest is None:
            dest = np.empty(self.shape, bool)
        if len(self.ranges) and isinstance(self.ranges[0], Ranges):
            for d, r in zip(dest, self.ranges):
                d[:] = r.mask()
        else:
            # Recurse
            for d, rm in zip(dest, self.ranges):
                rm.mask(dest=d)
        return dest


# Support functions for RangesMatrix.__getitem__.  It's helpful to
# take this logic out of the class to handle some differences between
# Ranges and RangesMatrix.
#
# In _gibasic and _giadv, the target must be a Ranges or RangesMatrix,
# and index must be a tuple with no Ellipsis in it (but None are ok).
# The entry point is _gibasic, which will call _giadv if/when it
# encounteres an advanced index.

def _gibasic(target, index):
    if len(index) == 0:
        return target
    if index[0] is None:
        return RangesMatrix([_gibasic(target, index[1:])], skip_shape_check=True)
    is_rm = isinstance(target, RangesMatrix)
    if not is_rm and len(index) > 1:
        raise IndexError(f'Too many indices (extras: {index[1:]}).')
    if isinstance(index[0], (np.ndarray, list, tuple)):
        if is_rm:
            return _giadv(target, index)
        raise IndexError('Ranges (or last axis of RangesMatrix) '
                         'cannot use advanced indexing.')
    if isinstance(index[0], (int, np.int32, np.int64)):
        if is_rm:
            return _gibasic(target.ranges[index[0]], index[1:])
        raise IndexError('Cannot apply integer index to Ranges object.')
    if isinstance(index[0], slice):
        if is_rm:
            rm = RangesMatrix([_gibasic(r, index[1:]) for r in target.ranges[index[0]]],
                              child_shape=target.shape[1:], skip_shape_check=True)
            if rm.shape[0] == 0:
                # If your output doesn't have any .ranges, you need to
                # fake one in order to figure out how the dimensions
                # play out.
                fake_child = RangesMatrix.zeros(target.shape[1:])
                rm._child_shape = _gibasic(fake_child, index[1:]).shape
            return rm
        return target[index[0]]
    raise IndexError(f'Unexpected target[index]: {target}[{index}]')

def _giadv(target, index):
    adv_index, adv_axis, basic_index = [], [], []
    adr_axis = 0
    for axis, ind in enumerate(index):
        if isinstance(ind, (list, tuple, np.ndarray)):
            ind = np.asarray(ind)
            if ind.dtype == bool:
                if ind.ndim != 1 or len(ind) != target.shape[adr_axis]:
                    raise IndexError('index mask with shape '
                                     f'{ind.shape} is not compatible with '
                                     f'{target.shape[adr_axis]}')
                ind = ind.nonzero()[0]
            adv_index.append(ind)
            adv_axis.append(axis)
            basic_index.append(None)
        else:
            basic_index.append(ind)
        if ind is not None:
            adr_axis += 1
    assert(adv_axis[0] == 0)  # Don't call me until you need to.

    br = np.broadcast(*adv_index)

    # If zero size, short circuit.
    if br.size == 0:
        for _axis in adv_axis:
            basic_index[_axis] = 0
        child_thing = _gibasic(target, basic_index)
        return RangesMatrix.zeros(br.shape + child_thing.shape)

    # Note super_list will not be empty.
    super_list = []
    for idx_tuple in br:
        for _axis, _basic in zip(adv_axis, idx_tuple):
            basic_index[_axis] = _basic
        sub = _gibasic(target, basic_index)
        super_list.append(sub)

    return _giassemble(super_list, br.shape)

def _giassemble(items, shape):
    if len(shape) > 1:
        stride = len(items) // shape[0]
        groups = [items[i*stride:(i+1)*stride] for i in range(shape[0])]
        items = [_giassemble(g, shape[1:]) for g in groups]
    return RangesMatrix(items)
