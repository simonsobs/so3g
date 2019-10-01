import so3g
from spt3g import core
import numpy as np
import os

# Get a copy of that ArgumentError...
ArgumentError = None
try:
    m0 = so3g.G3VectorFixed('blech')
except Exception as e:
    ArgumentError = e.__class__

# Test constructors.
print('Testing constructors...')
test_data = [0., 1., 2., 3., 4., 5.5]
prec = 0.1
for good_data in [
        test_data,
        [int(_d) for _d in test_data],
        np.array(test_data),
        core.G3VectorDouble(test_data),
        core.G3Timestream(test_data),
]:
    # Test with and without precision argument.
    m0 = so3g.G3VectorFixed(good_data, prec)
    assert(m0.precision == prec)
    m0 = so3g.G3VectorFixed(good_data)
    assert(m0.precision == 1.)


# Ok that these fail.
for bad_input in [
        np.array(test_data, int),
]:
    try:
        m0 = so3g.G3VectorFixed(bad_input)
        raise RuntimeError('Succeeded unexpectedly: %s' % bad_input)
    except ArgumentError:
        pass


# Test buffer protocol.
d = np.array(m0)
assert(len(d) == m0.n_samples)


# Compression test.

test_filename = 'test.g3'
n = 1000000
precision = 0.01
print('Test set has %i points, using %s.' % (n, test_filename))
test_data = np.sin(np.arange(n) * np.pi/180) + np.random.normal(size=n) * .1
test_data = (test_data / precision).round().astype(int) * precision

vd = core.G3VectorDouble(test_data)
vi = core.G3VectorInt([1,2,3,4])
vt = core.G3Timestream(vi)

m0 = so3g.G3VectorFixed(vd, .01)
m0.flac_level = 0

for level in [0, 1, 5, 9, 100]:
    print('Testing compression level %i' % level)
    m0.flac_level = level
    f = core.G3Frame()
    f['x'] = m0
    core.G3Writer(test_filename).Process(f)
    print('  Compressed file size is %i' % os.path.getsize(test_filename))
    f = core.G3Reader(test_filename).Process(None)[0]
    print('  Readback: ', np.array(f['x'])[:10])
    print('  Max error: ', max(abs(np.array(f['x']) - np.array(vd))))
    print()

# Test out of bounds:
limit = test_data.max() * .9
y = (test_data / limit * 2**23).round()
m = so3g.G3VectorFixed(test_data / limit * 2**23, 1.)
print('Range check: ', m.CheckRange(), ((y > 2**23-1) + (y < -2**23)).sum())

# Test precision.
m0 = so3g.G3VectorFixed(test_data)
m1 = so3g.G3VectorFixed(test_data, precision)
print('Precision check [0, lots]: ', m1.CheckPrecision(), m0.CheckPrecision())
