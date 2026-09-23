import unittest
import itertools

import so3g
from so3g import proj

import numpy as np

# Don't require pixell for testing
try:
    from pixell import enmap, wcsutils
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


def get_basics(clippedness=0):
    t, az, el = get_scan()
    csl = proj.CelestialSightLine.az_el(t, az, el, weather='vacuum', site='so')
    fp = proj.FocalPlane.from_xieta([0., .1*DEG], [0, .1*DEG])
    asm = proj.Assembly.attach(csl, fp)

    # And a map ... of where?
    ra, dec = csl.coords()[:, :2].T
    ra0, dec0 = ra.mean(), dec.mean()
    shape = (200, 350)  # This will not clip samples...
    if clippedness == 1:
        shape = (250, 300)  # But this will
    elif clippedness == 2:
        ra0 = ra0 - 10*DEG   # And this should miss the TOD entirely.
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

        target = np.zeros((1, ) + m.shape[:-1])
        with self.assertRaises(RuntimeError):
            p.to_map(sig, asm, comps='T', output=target)

        # Does det_weights seem to work?
        m = p.to_map(sig, asm, comps='T',
                     det_weights=np.array([0., 0.], dtype='float32'))[0]
        assert(np.all(m==0))

        # Can't assign to quat fields, so do
        # it this way instead
        asm.fplane.quats[1] = asm.fplane.quats[1]*np.nan
        with self.assertRaises(ValueError):
           p.to_map(sig, asm, comps='T')
        with self.assertRaises(ValueError):
           p.to_weights(asm, comps='T')

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
        for (clippedness, tiled, interpol, method) in itertools.product(
                [0, 1, 2],
                [False, True],
                ['nearest', 'bilinear'],
                proj.wcs.THREAD_ASSIGNMENT_METHODS):
            # For error messages ...
            detail = f'(method={method}, tiled={tiled}, clippedness={clippedness}, interpol={interpol})'
            scan, asm, (shape, wcs) = get_basics(clippedness=clippedness)
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

            target = {1: set([0,1]),
                      0: set([1]),
                      2: set([0])}[clippedness]
            self.assertEqual(set((counts0 + counts1).ravel()), target,
                             msg=f'threads does not cover TOD ({detail})')
            # Only the first segment should be non-empty, unless bilinear.
            if interpol == 'nearest':
                self.assertEqual(counts1.sum(), 0)

    @requires_pixell
    def test_30_bilin(self):
        """1-sample single boresight-pointed detector with no coordinate transformation"""
        # Trivial geometry where pixel coordinate and sky coordinate are the same thing
        shape = (2,2)
        wcs   = wcsutils.explicit(crval=[0,0], crpix=[1,1], cdelt=[1,1], ctype=["RA---CAR","DEC--CAR"])
        # Cases we will consider
        cases = [
            # A single sample hitting (0.1, 0.7)
            np.array([[0.1],[0.7]]),
            # A single sample with integer coordinates
            np.array([[1.0],[0.0]]),
        ]
        for x, y in cases:
            whole = np.all(x==np.round(x)) and np.all(y==np.round(y))
            csl   = proj.CelestialSightLine.for_lonlat(x*DEG, y*DEG)
            fp    = proj.FocalPlane.from_xieta([0.0],[0.0])
            asm   = proj.Assembly.attach(csl, fp)
            dtype = np.float32
            tod   = np.ones((1,1), dtype)
            # Expected result for NN and bilin
            wy    = np.sum(np.round([1-y,y]),1)
            wx    = np.sum(np.round([1-x,x]),1)
            targ_nn = wy[:,None]*wx[None,:]
            wy    = np.sum([1-y,y],1)
            wx    = np.sum([1-x,x],1)
            targ_li = wy[:,None]*wx[None,:]
            # Nearest neighbor untiled
            p     = proj.Projectionist.for_geom(shape, wcs, interpol="nearest")
            m_nn  = p.to_map(tod, asm, comps="T")[0]
            assert np.allclose(m_nn, targ_nn)
            # Bilinear untiled
            p     = proj.Projectionist.for_geom(shape, wcs, interpol="bilinear")
            m_li  = p.to_map(tod, asm, comps="T")[0]
            assert np.allclose(m_li, targ_li)
            if whole: assert np.allclose(m_nn, m_li)
            # Nearest neighbor tiled
            p     = proj.Projectionist.for_tiled(shape, wcs, (10, 10), interpol="nearest")
            m_nn  = p.to_map(tod, asm, comps="T")[0][0]
            assert np.allclose(m_nn, targ_nn)
            p     = proj.Projectionist.for_tiled(shape, wcs, (10, 10), interpol="bilinear")
            m_li  = p.to_map(tod, asm, comps="T")[0][0]
            assert np.allclose(m_li, targ_li)
            if whole: assert np.allclose(m_nn, m_li)

    def test_40_det_left_engine(self):
        """ProjEng with det_left should compute q_ofs * q_bore."""
        t, az, el = get_scan()
        csl = proj.CelestialSightLine.naive_az_el(t[:1000], az[:1000], el[:1000])
        fp = proj.FocalPlane.from_xieta([0., .1*DEG, -.3*DEG], [0, .1*DEG, .2*DEG],
                                        [0., 10*DEG, 45*DEG])
        pe = so3g.ProjEng_Quat_TQU_NonTiled((1, 1, 1., 1., 1., 1.))
        self.assertFalse(pe.det_left)
        q_normal = pe.coords(csl.Q, fp.quats, None)
        pe.det_left = True
        self.assertTrue(pe.det_left)
        q_left = pe.coords(csl.Q, fp.quats, None)
        for i, q in enumerate(fp.quats):
            np.testing.assert_allclose(np.asarray(q_normal[i]),
                                       np.asarray(csl.Q * q), atol=1e-12)
            np.testing.assert_allclose(np.asarray(q_left[i]),
                                       np.asarray(q * csl.Q), atol=1e-12)

    @requires_pixell
    def test_41_det_left_projectionist(self):
        """A det_left Assembly should project like the equivalent
        explicit per-detector boresight."""
        scan, asm, (shape, wcs) = get_basics()
        Q = asm.Q
        qmid = Q[len(Q)//2]
        # Detector quats that, applied on the left, give pointing close
        # to the usual Q * q_det, so that samples land in the map.
        ofs = [qmid * q * ~qmid for q in asm.fplane.quats]
        fp_left = proj.FocalPlane(quats=np.array([[q.a, q.b, q.c, q.d] for q in ofs]))
        asm_left = proj.Assembly.attach(Q, fp_left, det_left=True)
        refs = [proj.Assembly.attach(q * Q, proj.FocalPlane.boresight())
                for q in ofs]

        for tiled in [False, True]:
            if tiled:
                p = proj.Projectionist.for_tiled(shape, wcs, (50, 50))
            else:
                p = proj.Projectionist.for_geom(shape, wcs)
            pix = p.get_pixels(asm_left)
            coords = p.get_coords(asm_left)
            for i, ref in enumerate(refs):
                pix_ref = p.get_pixels(ref)[0]
                assert np.any(np.asarray(pix_ref)[:, 0] >= 0)
                np.testing.assert_array_equal(np.asarray(pix[i]), np.asarray(pix_ref))
                np.testing.assert_allclose(np.asarray(coords[i]),
                                           np.asarray(p.get_coords(ref)[0]), atol=1e-10)

        # Round trip through to_map / from_map.
        p = proj.Projectionist.for_geom(shape, wcs)
        sig = np.ones((len(ofs), len(Q)), 'float32')
        m = p.to_map(sig, asm_left, comps='T')
        m_ref = sum(p.to_map(sig[i:i+1], ref, comps='T') for i, ref in enumerate(refs))
        np.testing.assert_allclose(m, m_ref)
        tod = p.from_map(m, asm_left, comps='T')
        for i, ref in enumerate(refs):
            np.testing.assert_allclose(tod[i], p.from_map(m, ref, comps='T')[0])

        # Thread assignment should work (domdir falls back to simple).
        for method in ['simple', 'domdir']:
            threads = p.assign_threads(asm_left, method=method, n_threads=2)
            m_thr = p.to_map(sig, asm_left, comps='T', threads=threads)
            np.testing.assert_allclose(m_thr, m)

if __name__ == '__main__':
    unittest.main()
