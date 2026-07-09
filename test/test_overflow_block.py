import unittest

import so3g
import numpy as np


class TestOverflowBlock(unittest.TestCase):

    def test_overflow(self):
        a = np.zeros((4050, 603260), dtype='float32')
        so3g.block_minmax(a, a, 800, 2, 0)


if __name__ == '__main__':
    unittest.main()
