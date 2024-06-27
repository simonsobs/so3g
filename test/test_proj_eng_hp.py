import unittest
import itertools

from so3g import proj

import numpy as np

DEG = np.pi/180


def get_scan():
    dt = np.arange(0., 600, .1)
    T = 10.
    phase = (dt % (2*T)) / T
    phase[phase > 1] = 2 - phase[phase > 1]
    v = 2.
    az = v * (phase - .5) + 60.
    el = 50. + dt*0
    times = 1654041600 + dt
    return times, az*DEG, el*DEG


def get_basics():
    t, az, el = get_scan()
    csl = proj.CelestialSightLine.az_el(t, az, el, weather='vacuum', site='so')
    fp = proj.FocalPlane.from_xieta(['a', 'b'], [0., .1*DEG], [0, .1*DEG])
    asm = proj.Assembly.attach(csl, fp)
    return ((t, az, el), asm)


class TestProjEngHP(unittest.TestCase):
    """Test ProjectionistHealpix
       Based on TestProjEng
    """
    
    def test_00_basic(self):
        scan, asm = get_basics()
        nside = 128
        p = proj.ProjectionistHealpix.for_healpix(nside)
        sig = np.ones((2, len(scan[0])), 'float32')
        for comps in ['T', 'QU', 'TQU']:
            m = p.to_map(sig, asm, comps=comps)
            assert(np.any(m != 0))
            w = p.to_weights(asm, comps=comps)
            assert(np.any(w != 0))
            m = np.asarray(m, dtype=np.float64)
            s = p.from_map(m, asm, comps=comps)

        # Does det_weights seem to work?
        m = p.to_map(sig, asm, comps='T',
                     det_weights=np.array([0., 0.], dtype='float32'))[0]
        assert(np.all(m==0))

        # Raise if pointing invalid.
        asm.dets[1,2] = np.nan
        with self.assertRaises(ValueError):
           p.to_map(sig, asm, comps='T')
        with self.assertRaises(ValueError):
           p.to_weights(asm, comps='T')

    def test_10_tiled(self):
        scan, asm = get_basics()
        nside = 128
        for nside_tile in [8, 'auto']:
            p = proj.ProjectionistHealpix.for_healpix(nside, nside_tile)
            p.active_tiles = p.get_active_tiles(asm)['active_tiles']
            sig = np.ones((2, len(scan[0])), 'float32')
            for comps in ['T', 'QU', 'TQU']:
                m = p.to_map(sig, asm, comps=comps)
                m2 = [tile for tile in m if tile is not None]
                assert(np.any(m2))
                w = p.to_weights(asm, comps=comps)
                w2 = [tile for tile in w if tile is not None]
                assert(np.any(w2))
            # Identify active subtiles?
            print(p.active_tiles)
        
    def test_20_threads(self):
        for (tiled, interpol, method) in itertools.product(
                [False, True],
                ['nearest'],
                ['simple', 'tiles']):
            # For error messages ...
            detail = f'(method={method}, tiled={tiled}, interpol={interpol})'
            scan, asm = get_basics()
            nside = 128
            if tiled:
                nside_tile = 8
            else:
                nside_tile = None
                
            p = proj.ProjectionistHealpix.for_healpix(nside, nside_tile, interpol=interpol)
            sig = np.ones((2, len(scan[0])), 'float32')
            n_threads = 3

            if method in ['tiles'] and not tiled:
                with self.assertRaises(RuntimeError, msg=
                                       f'Expected assignment to fail ({detail})'):
                    threads = p.assign_threads(asm, method=method, n_threads=n_threads)
                continue
            else:
                threads = p.assign_threads(asm, method=method, n_threads=n_threads)
            # This may need to be generalized if we implement fancier threads schemes.
            self.assertIsInstance(threads, list,
                                  msg=f'a thread assignment routine did not return a list ({detail})')

            # Make sure the threads cover the TOD, or not,
            # depending on clipped.
            counts0 = threads[0].mask().sum(axis=0)
            counts1 = np.zeros(counts0.shape, int)

            self.assertEqual(threads[0].shape, (n_threads,) + sig.shape,
                             msg=f'a thread bunch has wrong shape ({detail})')

            for t in threads[1:]:
                counts1 += t.mask().sum(axis=0)
                self.assertEqual(t.shape[1:], sig.shape,
                                 msg=f'a thread bunch has unexpected shape ({detail})')

            target = set([1])
            self.assertEqual(set((counts0 + counts1).ravel()), target,
                             msg=f'threads does not cover TOD ({detail})')
            # Only the first segment should be non-empty, unless bilinear.
            if interpol == 'nearest':
                self.assertEqual(counts1.sum(), 0)


if __name__ == '__main__':
    unittest.main()
