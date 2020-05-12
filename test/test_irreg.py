import unittest

import so3g
from spt3g import core
SEC = core.G3Units.sec

import numpy as np

def get_test_block(length, keys, offset=0):
    type_cycle = [core.G3VectorDouble, core.G3VectorInt]
    t0 = core.G3Time('2019-01-01T12:30:00') + offset*SEC
    m = core.G3TimesampleMap()
    m.times = core.G3VectorTime([t0 + t*SEC for t in np.arange(length)])
    for i,k in enumerate(keys):
        y = (np.random.uniform(size=length) * 100).astype(int)
        m[k] = type_cycle[i % len(type_cycle)](y)
    return m
    
class TestIrregBlock(unittest.TestCase):
    def test_00_internal_checks(self):
        # Valid block.
        m = get_test_block(100, ['x', 'y', 'z'])
        m.check()

        # Construct invalid blocks.
        m = core.G3TimesampleMap()
        t0 = core.G3Time('2019-01-01T12:30:00')
        m.times = core.G3VectorTime([t0, t0 + 10*SEC, t0 + 20*SEC])
        with self.assertRaises(ValueError):
            m['x'] = core.G3VectorDouble([1,2])

    def test_10_concat(self):
        # Test concatenation.
        key_list = ['x', 'y', 'z']
        m0 = get_test_block(100, key_list)
        m1 = get_test_block(200, key_list, 100)
        for fail_type, fail_vec in [
                (ValueError, get_test_block(200, key_list + ['extra'], 100)),
                (ValueError, get_test_block(200, key_list[:-1], 100)),
                ]:
            with self.assertRaises(fail_type):
                m0.concatenate(fail_vec)

    def test_20_serialization(self):
        m0 = get_test_block(100, ['x', 'y'])
        m1 = get_test_block(200, ['x', 'y'], 100)
        m2 = m0.concatenate(m1)
        m0.check()
        m1.check()
        m2.check()
        f = core.G3Frame()
        f['irreg0'] = m0
        f['irreg1'] = m1
        core.G3Writer('test.g3').Process(f)
        f = core.G3Reader('test.g3').Process(None)[0]
        f['irreg0'].check()
        f['irreg1'].check()
        f['irreg0'].concatenate(f['irreg1'])['x']

if __name__ == '__main__':
    unittest.main()
