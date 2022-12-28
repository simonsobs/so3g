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
    def __init__(self, items=[], child_shape=None):
        self.ranges = [x for x in items]
        if len(items):
            child_shape = items[0].shape
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
        if index is None:
            index = (None,)
        if isinstance(index, slice):
            index = (index,)
        if isinstance(index, tuple) and len(index) == 0:
            return self
        elif isinstance(index, tuple):
            if index[0] is None:
                #Insert a dimension.
                if 0 in self.shape:
                    return self.__class__([], self.shape)
                return self.__class__([self[index[1:]]], self.shape)

            if len(self.shape) == 2 and len(index) > 1 and index[1] is None:
                # This case corresponds to trying to inject a new
                # dimension before the last one, e.g. if shape is
                # (100, 10000) and index is [:,None].  This requires
                # special treatment because a simple Ranges object
                # can't self-promote like Ranges(10000)[None,:] ->
                # RangesMatrix(1,10000).
                new_index = tuple([index[0], slice(0,1)] + list(index[2:]))
                new_shape = (self.shape[0], 1, self.shape[1])
                new_self = RangesMatrix([RangesMatrix([r], new_shape[2:])
                                         for r in self.ranges], new_shape[1:])
                return new_self[new_index]

            if len(self.shape) == 2 and len(index) > 2:
                raise IndexError("Too many indices to RangesMatrix.")

            if isinstance(index[0], np.ndarray):
                if index[0].dtype is bool:
                    return RangesMatrix([self.ranges[i][index[1:]]
                                           for i in index[0].nonzero()[0]])
                return RangesMatrix([self.ranges[i][index[1:]] for i in index[0]],
                                    self.shape[1:])
            elif isinstance(index[0], slice):
                return RangesMatrix([d[index[1:]] for d in self.ranges[index[0]]],
                                    self.shape[1:])
            else:
                return self.ranges[index[0]][index[1:]]
        return self.ranges[index]
    
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
