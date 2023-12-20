import unittest
import itertools

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


def get_basics(clipped=True):
    t, az, el = get_scan()
    csl = proj.CelestialSightLine.az_el(t, az, el, weather='vacuum', site='so')
    fp = proj.FocalPlane.from_xieta(['a', 'b'], [0., .1*DEG], [0, .1*DEG])
    asm = proj.Assembly.attach(csl, fp)

    # And a map ... of where?
    ra, dec = csl.coords()[:, :2].T
    ra0, dec0 = ra.mean(), dec.mean()
    shape = (250, 300)  # This will clip some samples
    if not clipped:
        shape = (200, 350)  # This will not.
    shape, wcs = enmap.geometry((dec0, ra0), res=(.01*DEG, -0.01*DEG),
                                shape=shape, proj='tan', ref=(dec0, ra0))
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
        for (clipped, tiled, interpol, method) in itertools.product(
                [False, True],
                [False, True],
                ['nearest', 'bilinear'],
                proj.wcs.THREAD_ASSIGNMENT_METHODS):
            # For error messages ...
            detail = f'(method={method}, tiled={tiled}, clipped={clipped}, interpol={interpol})'
            scan, asm, (shape, wcs) = get_basics(clipped=clipped)
            if tiled:
                p = proj.Projectionist.for_tiled(shape, wcs, (150, 150), active_tiles=False,
                                                 interpol=interpol)
            else:
                p = proj.Projectionist.for_geom(shape, wcs, interpol=interpol)
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

            target = set([0,1]) if clipped else set([1])
            self.assertEqual(set((counts0 + counts1).ravel()), target,
                             msg=f'threads does not cover TOD ({detail})')
            # Only the first segment should be non-empty, unless bilinear.
            if interpol == 'nearest':
                self.assertEqual(counts1.sum(), 0)


if __name__ == '__main__':
    unittest.main()
