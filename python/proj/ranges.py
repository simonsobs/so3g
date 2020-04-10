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
            elif isinstance(index[0], np.ndarray):
                if index[0].dtype is np.bool:
                    return RangesMatrix([self.ranges[i][index[1:]]
                                           for i in index[0].nonzero()[0]])
                return RangesMatrix([self.ranges[i][index[1:]] for i in index[0]],
                                    self.shape[1:])
            return RangesMatrix([d[index[1:]] for d in self.ranges[index[0]]],
                                self.shape[1:])
        return self.ranges[index]

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
            s1 = list(items[0].shape)
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
                    r = Ranges.from_array(np.array(r))
                    r.count = n
                    return r
            return RangesMatrix(ranges)
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
