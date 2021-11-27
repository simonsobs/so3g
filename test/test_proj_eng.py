import unittest

from so3g import proj

import numpy as np

# Don't require pixell for testing
try:
    from pixell import enmap
    pixell_found = True
except ModuleNotFoundError:
    pixell_found = False

requires_pixell = unittest.skipIf(pixell_found is False, "pixell not found")

DEG = np.pi/180


# Making these testing-support functions available in the library
# would be great.

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

    # And a map ... of where?
    ra, dec = csl.coords()[:, :2].T
    ra0, dec0 = ra.mean(), dec.mean()
    shape, wcs = enmap.geometry((dec0, ra0), res=(.01*DEG, -0.01*DEG),
                                shape=(300, 300), proj='tan', ref=(dec0, ra0))
    return ((t, az, el), asm, (shape, wcs))


class TestProjEng(unittest.TestCase):
    """Test the Projectionist and supporting structures.

    """
    @requires_pixell
    def test_00_basic(self):
        scan, asm, (shape, wcs) = get_basics()
        p = proj.Projectionist.for_geom(shape, wcs)
        sig = np.ones((2, len(scan[0])), 'float32')
        for comps in ['T', 'QU', 'TQU']:
            m = p.to_map(sig, asm, comps=comps)[0]
            assert(np.any(m != 0))
            w = p.to_weights(asm, comps=comps)[0, 0]
            assert(np.any(w != 0))

    @requires_pixell
    def test_10_tiled(self):
        scan, asm, (shape, wcs) = get_basics()
        p = proj.Projectionist.for_tiled(shape, wcs, (150, 150))
        sig = np.ones((2, len(scan[0])), 'float32')
        for comps in ['T', 'QU', 'TQU']:
            m = p.to_map(sig, asm, comps=comps)[0][0]
            assert(np.any(m != 0))
            w = p.to_weights(asm, comps=comps)[0][0, 0]
            assert(np.any(w != 0))
        # Identify active subtiles?
        p = proj.Projectionist.for_tiled(shape, wcs, (20, 20))
        print(p.active_tiles)
        p2 = p.get_active_tiles(asm, assign=2)
        print(p2)

    @requires_pixell
    def test_20_threads(self):
        scan, asm, (shape, wcs) = get_basics()
        p = proj.Projectionist.for_geom(shape, wcs)
        sig = np.ones((2, len(scan[0])), 'float32')
        n_threads = 3
        for method in proj.wcs.THREAD_ASSIGNMENT_METHODS:
            print(f'Assigning threads using {method}...')
            if method not in ['tiles']:
                threads = p.assign_threads(asm, method=method, n_threads=n_threads)
                self.assertEqual(threads.shape, (n_threads,) + sig.shape)
            else:
                with self.assertRaises(RuntimeError):
                    threads = p.assign_threads(asm, method=method, n_threads=n_threads)

    @requires_pixell
    def test_21_threads_tiled(self):
        scan, asm, (shape, wcs) = get_basics()
        sig = np.ones((len(asm.dets), len(scan[0])), 'float32')
        n_threads = 3
        comps = 'T'
        for method in proj.wcs.THREAD_ASSIGNMENT_METHODS:
            p = proj.Projectionist.for_tiled(shape, wcs, (150, 150), active_tiles=False)
            print(f'Assigning threads using {method}...')
            threads = p.assign_threads(asm, method=method, n_threads=n_threads)
            self.assertEqual(threads.shape, (n_threads,) + sig.shape)


if __name__ == '__main__':
    unittest.main()
