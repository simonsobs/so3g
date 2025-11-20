import so3g
import spt3g.core as core
import numpy as np

def length_tests(iv, rows, indent_text='    '):
    for (lo, hi), len_exp, len_comp in rows:
        iv.add_interval(lo,hi)
        print(indent_text + 'add [', lo, ',', hi, ') ->', iv.array())
        assert(len(iv.array()) == len_exp)
        assert(len((-iv).array()) == len_comp)


print('Testing interface for all data types:')
for dtype in [
        so3g.IntervalsDouble,
        so3g.IntervalsInt,
]:
    print('    ', dtype)
    o = dtype()
    a = o.array()
    print('    ', o)
    print('    ', a, a.dtype, a.shape)


print()
print('Testing operations.')
iv = so3g.IntervalsDouble()
iv.add_interval(1., 2.)
iv.add_interval(3., 4.)
everything = (iv + ~iv)
assert(len(everything.array()) == 1)


print()
print('Testing simple interval merging:')
iv = so3g.IntervalsDouble()
length_tests(iv, [
    ((1., 2.), 1, 2),
    ((3., 4.), 2, 3),
    ((2., 3.), 1, 2)])


print()
print('Testing domain trimming:')
iv = so3g.IntervalsDouble(0., 1.)
print('   ', iv)
length_tests(iv, [
    (( 0.1,  0.2), 1, 2),
    ((-1.0,  0.0), 1, 2),
    (( 1.0,  2.0), 1, 2),
    ((-1.0,  0.1), 1, 1),
    ((-0.2,  1.1), 1, 0),
])

print()
print('Testing domain treatment on combination')
lo0, hi0 = 0., 10.
for lo1 in [-1, 5, 11]:
    for hi1 in [-1, 5, 11]:
        iv0 = so3g.IntervalsDouble(lo0, hi0)
        iv1 = so3g.IntervalsDouble(lo1, hi1)
        ivx = iv0 * iv1
        lo, hi = ivx.domain
        print('    ', iv0.domain, ' + ', iv1.domain, ' -> ',
              ivx.domain)
        assert( lo >= max(lo0, lo1) and
                (hi==lo or hi <= min(hi0, hi1)) )

print()
print('Testing import.')
iv0 = so3g.IntervalsDouble()\
          .add_interval(1, 2)\
          .add_interval(3, 4)


iv1 = iv0.copy()
print('iv0 = ', iv0)
print('iv1 = ', iv1, flush=True)
print('iv0.array() = ', iv0.array(), flush=True)
print('iv1.array() = ', iv1.array(), flush=True)
assert(np.all(iv0.array() == iv1.array()))

iv1 = so3g.IntervalsDouble.from_array(iv0.array())
print('iv1 = iv0.from_array(): ', iv1, flush=True)
assert(np.all(iv0.array() == iv1.array()))


print()
iv1 = so3g.IntervalsDouble()\
          .add_interval(0., 1.)\
          .add_interval(2., 3.)\
          .add_interval(4., 5.)

iv2 = so3g.IntervalsDouble()\
          .add_interval(1., 2.5)

assert(len((iv1 + iv2).array()) == 2)
assert(len((iv1 * iv2).array()) == 1)
assert(len((iv1 - iv2).array()) == 2)
assert(len((iv2 - iv1).array()) == 4)


print()
print('Interval <-> mask testing')
mask = np.zeros(20, np.uint16)
n_bit, target_bit = 16, 12
for ikill, nint in [(None, 0),
                    (19,   1),
                    ( 0,   2),
                    (10,   3),
                    (11,   3)]:
    if ikill is not None:
        mask[ikill] = (1<<target_bit)
    iv = so3g.IntervalsInt.from_mask(mask, n_bit)
    print(f"IntervalsInt.from_mask({mask},{n_bit}) -> {iv}", flush=True)
    assert(len(iv) == n_bit)
    for i in range(n_bit):
        if i == target_bit:
            assert(len(iv[i].array()) == nint)
        else:
            assert(len(iv[i].array()) == 0)

print(f"so3g.IntervalsInt.mask({iv}, -1)", flush=True)
mask1 = so3g.IntervalsInt.mask(iv,-1)
print(f"mask = {mask}, mask1 = {mask1}", flush=True)
assert(np.all(mask == mask1))

print('...bit-width checking works?')
try:
    mask3 = so3g.IntervalsInt.mask(iv,8)
except ValueError as e:
    mask3 = 'failed'
    print(' -> ', e)
assert(mask3 == 'failed')

print('...domain checking works?')
iv.append(so3g.IntervalsInt(-2,20))
try:
    mask3 = so3g.IntervalsInt.mask(iv,-1)
except Exception as e:
    mask3 = 'failed'
    print(' -> ', e)
assert(mask3 == 'failed')

# Type failing works?  Can't create mask from non-integer Intervals.
try:
    mask3 = so3g.IntervalsDouble.mask([], 8)
except ValueError:
    mask3 = 'failed'
assert(mask3 == 'failed')
