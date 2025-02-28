import unittest

import so3g
import numpy as np

from scipy.interpolate import interp1d
from scipy.signal import welch, windows


class TestPolyFill(unittest.TestCase):
    """Test the polynomial gap filling."""

    def test_00_basic(self):
        n = 1000
        p0 = [1.0, 0.2, -0.14]
        t = np.arange(1000) * 0.1
        sig0 = np.polyval(p0, t)[None, :]

        ix, nx = 200, 20
        r = so3g.RangesInt32(n)
        r.add_interval(ix, ix + nx)
        s = r.mask()

        BAD_DATA = -100

        for dtype, fillfunc, tolerance in [
            ("float32", so3g.get_gap_fill_poly, 1e-4),
            ("float64", so3g.get_gap_fill_poly64, 1e-10),
        ]:
            sig1 = sig0.astype(dtype)

            # Write the fill into ex array (inplace=False)
            sig = sig1.copy()
            sig[0, s] = BAD_DATA
            ex = np.zeros((s.sum()), dtype=sig.dtype)
            fillfunc([r], sig, 10, 2, False, ex)
            np.testing.assert_allclose(ex, sig1[0, s], rtol=tolerance)
            np.testing.assert_allclose(BAD_DATA, sig[0, s], rtol=tolerance)

            # Write the fill and extract the filled (inplace=True)
            sig = sig1.copy()
            sig[0, s] = BAD_DATA
            ex = np.zeros((s.sum()), dtype=sig.dtype)
            fillfunc([r], sig, 10, 2, True, ex)
            np.testing.assert_allclose(sig1[0, s], sig[0, s], rtol=tolerance)
            np.testing.assert_allclose(BAD_DATA, ex, rtol=tolerance)


class TestBufferWrapper(unittest.TestCase):
    def test_00(self):
        for array_shape, pattern in [
            ((2, 4), (2, 4)),
            ((2, 4), (-1, -1)),
            ((2, 4), (-2,)),
            ((2, 4), ( -1, -1, -2,)),
            ((2, 4), (-2, -1, -1)),
            ((2, 4), (-1, -2, -1)),
        ]:
            a = np.zeros(array_shape)
            so3g.test_buffer_wrapper(a, list(pattern))

        for array_shape, pattern in [
            ((2, 4), (-2, -2)),
            ((2, 4), (2,)),
            ((2, 4), (4,)),
            ((2, 4), (2, 3)),
            ((2, 4), (2, 4, -1)),
            ((2, 4), (2, -1, -1)),
            ((2, 4), (-1, -1, -1, -2)),
        ]:
            a = np.zeros(array_shape)
            with self.assertRaises(RuntimeError):
                so3g.test_buffer_wrapper(a, list(pattern))


class TestJumps(unittest.TestCase):
    """
    Test functions used for jump finding.
    """

    def test_00_clean_flag(self):
        a = np.zeros((1, 100), "int32")
        a[:, 10] = 1
        a[:, 47] = 1
        a[:, 0:8] = 1
        a[:, 50:60] = 1
        a[:, 95:] = 1
        b = a.copy(order="C")
        so3g.clean_flag(b, 5)
        flagged = np.where(b)[1]
        # fmt: off
        np.testing.assert_array_equal(
            flagged,
            np.array([0, 1, 2, 3, 4, 5, 6, 7, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 95, 96, 97, 98, 99,]),
        )
        # fmt: on

    def test_01_matched_jump(self):
        a = np.zeros((1, 100), dtype="float32", order="C")
        a[0, 50:] += 1
        b = np.zeros((1, 100), dtype="int32", order="C")
        min_size = np.ones(1, dtype="float32", order="C") * 0.5
        so3g.matched_jumps(a, b, min_size, 10)
        flagged = np.where(b)[1]
        np.testing.assert_array_equal(flagged, np.arange(46, 53))

    def test_02_quantized_jump(self):
        a = np.zeros((1, 100), dtype="float32", order="C")
        a[0, 50:] += 1
        b = np.zeros((1, 100), dtype="float32", order="C")
        so3g.find_quantized_jumps(a, b, np.zeros(1, "float32", "C"), 10, 0.5)
        flagged = np.where(b)[1]
        np.testing.assert_array_equal(flagged, np.arange(50, 60))

    def test_03_subtract_height(self):
        a = np.zeros((1, 100), dtype="float32", order="C")
        orig = a.copy()
        a[0, 25:] += 1
        a[0, 75:] -= 1
        h = np.zeros((1, 100), dtype="float32", order="C")
        h[0, 25] = 1
        h[0, 75] = -1
        j = so3g.proj.RangesMatrix.from_mask(h != 0)
        so3g.subtract_jump_heights(a, a, h, j)
        np.testing.assert_allclose(a, orig)


class TestBlockReduce(unittest.TestCase):
    """
    Test block reduction functions.
    """

    def test_00_mean(self):
        a = np.ascontiguousarray(np.arange(0, 100, dtype="float32")).reshape((1, -1))
        b = a.copy()
        so3g.block_moment(a, b, 10, 1, 0, 0)
        c = a.copy()
        for i in range(10):
            start = i * 10
            stop = start + 10
            c[0, start:stop] = np.mean(a[0, start:stop])
        np.testing.assert_allclose(b, c)

    def test_01_var(self):
        a = np.ascontiguousarray(np.arange(0, 100, dtype="float32")).reshape((1, -1))
        b = a.copy()
        so3g.block_moment(a, b, 10, 2, 1, 0)
        c = a.copy()
        for i in range(10):
            start = i * 10
            stop = start + 10
            c[0, start:stop] = np.var(a[0, start:stop])
        np.testing.assert_allclose(b, c)

    def test_02_ptp(self):
        a = np.ascontiguousarray(np.arange(0, 100, dtype="float32")).reshape((1, -1))
        b = a.copy()
        so3g.block_minmax(a, b, 10, 2, 0)
        c = a.copy()
        for i in range(10):
            start = i * 10
            stop = start + 10
            c[0, start:stop] = np.ptp(a[0, start:stop])
        np.testing.assert_allclose(b, c)


class TestGslInterpolate(unittest.TestCase):
    """
    Test interpolation using GSL.
    """

    def test_00_linear_interp_float32(self):
        t_start = 0
        t_end = 999
        t_size = 500

        t_interp_start = 0
        t_interp_end = 999
        t_interp_size = 2000

        ndet = 3
        dtype = "float32"
        order = "C"

        t = np.linspace(t_start, t_end, t_size, dtype=dtype)
        sig = np.array(
            [(i + 1) * np.sin(2 * np.pi * 0.01 * t + i) for i in range(ndet)],
            dtype=dtype,
            order=order,
        )

        t_interp = np.linspace(t_interp_start, t_interp_end, t_interp_size, dtype=dtype)

        f_template = interp1d(t, sig, fill_value="extrapolate")
        scipy_sig = f_template(t_interp)

        so3g_sig = np.zeros([ndet, t_interp_size], dtype=dtype, order=order)
        so3g.interp1d_linear(t, sig, t_interp, so3g_sig)

        tolerance = 1e-4
        np.testing.assert_allclose(scipy_sig, so3g_sig, rtol=tolerance)

    def test_01_linear_interp_float64(self):
        t_start = 0
        t_end = 999
        t_size = 500

        t_interp_start = 0
        t_interp_end = 999
        t_interp_size = 2000

        ndet = 3
        dtype = "float64"
        order = "C"

        t = np.linspace(t_start, t_end, t_size, dtype=dtype)
        sig = np.array(
            [(i + 1) * np.sin(2 * np.pi * 0.01 * t + i) for i in range(ndet)],
            dtype=dtype,
            order=order,
        )

        t_interp = np.linspace(t_interp_start, t_interp_end, t_interp_size, dtype=dtype)

        f_template = interp1d(t, sig, fill_value="extrapolate")
        scipy_sig = f_template(t_interp)

        so3g_sig = np.zeros([ndet, t_interp_size], dtype=dtype, order=order)
        so3g.interp1d_linear(t, sig, t_interp, so3g_sig)

        tolerance = 1e-10
        np.testing.assert_allclose(scipy_sig, so3g_sig, rtol=tolerance)

    def test_02_linear_extrapolation(self):
        t_start = 0.0
        t_end = 999.0
        t_size = 500

        t_interp_start = -10.0
        t_interp_end = 1009.0
        t_interp_size = 2000

        ndet = 3
        dtype = "float32"
        order = "C"

        t = np.linspace(t_start, t_end, t_size, dtype=dtype)
        sig = np.array(
            [(i + 1) * np.sin(2 * np.pi * 0.01 * t + i) for i in range(ndet)],
            dtype=dtype,
            order=order,
        )

        t_interp = np.linspace(t_interp_start, t_interp_end, t_interp_size, dtype=dtype)

        f_template = interp1d(t, sig, fill_value="extrapolate")
        scipy_sig = f_template(t_interp)

        so3g_sig = np.zeros((ndet, t_interp_size), dtype=dtype, order=order)
        so3g.interp1d_linear(t, sig, t_interp, so3g_sig)

        tolerance = 1e-4
        np.testing.assert_allclose(scipy_sig, so3g_sig, rtol=tolerance)

    def test_03_uneven_spacing(self):
        t_start = 0.0
        t_end = 999.0
        t_size = 500

        t_interp_start = 0.0
        t_interp_end = 999.0
        t_interp_size = 2000

        ndet = 3
        dtype = "float32"
        order = "C"

        # Generate uneven spaced time samples with power law
        t_pow = 1.3

        t = np.linspace(
            t_start ** (1 / t_pow), t_end ** (1 / t_pow), t_size, dtype=dtype
        )
        t = t**t_pow
        sig = np.array(
            [(i + 1) * np.sin(2 * np.pi * 0.01 * t + i) for i in range(ndet)],
            dtype=dtype,
            order=order,
        )

        t_interp = np.linspace(
            t_interp_start ** (1 / t_pow),
            t_interp_end ** (1 / t_pow),
            t_interp_size,
            dtype=dtype,
        )
        t_interp = t_interp**t_pow

        f_template = interp1d(t, sig, fill_value="extrapolate")
        scipy_sig = f_template(t_interp)

        so3g_sig = np.zeros((ndet, t_interp_size), dtype=dtype, order=order)
        so3g.interp1d_linear(t, sig, t_interp, so3g_sig)

        tolerance = 1e-4
        np.testing.assert_allclose(scipy_sig, so3g_sig, rtol=tolerance)

    def test_04_array_slicing(self):
        t_start = 0
        t_end = 999
        t_size = 500
        slice_offset = 100  # Sample index to start input arrays at

        t_interp_start = 0
        t_interp_end = 999
        t_interp_size = 2000
        interp_slice_offset = 1000  # Sample index to start interpolated arrays at

        ndet = 3
        dtype = "float32"
        order = "C"

        t = np.linspace(t_start, t_end, t_size, dtype=dtype)
        sig = np.array(
            [(i + 1) * np.sin(2 * np.pi * 0.01 * t + i) for i in range(ndet)],
            dtype=dtype,
            order=order,
        )

        t_interp = np.linspace(t_interp_start, t_interp_end, t_interp_size, dtype=dtype)

        f_template = interp1d(
            t[slice_offset:], sig[:, slice_offset:], fill_value="extrapolate"
        )
        scipy_sig = f_template(t_interp[interp_slice_offset:])

        so3g_sig = np.zeros((ndet, t_interp_size), dtype=dtype, order=order)
        so3g.interp1d_linear(
            t[slice_offset:],
            sig[:, slice_offset:],
            t_interp[interp_slice_offset:],
            so3g_sig[:, interp_slice_offset:],
        )

        tolerance = 1e-4
        np.testing.assert_allclose(
            scipy_sig, so3g_sig[:, interp_slice_offset:], rtol=tolerance
        )


class TestDetrend(unittest.TestCase):
    """
    Test detrending.
    """

    def test_00_mean_detrending(self):
        nsamps = 1000
        ndets = 3
        dtype = "float32"
        order = "C"

        x = np.linspace(0., 1., nsamps, dtype=dtype)
        signal = np.array([(i + 1) * np.sin(2*np.pi*x + i) for i in range(ndets)], dtype=dtype, order=order)

        signal_copy = signal.copy(order=order)
        signal_copy -= np.mean(signal_copy, axis=-1, dtype=dtype)[..., None]

        method = "mean"
        count = 0 # not used for mean detrending
        so3g.detrend(signal, method, count)

        rtol = 0
        atol = 1e-5
        np.testing.assert_allclose(signal_copy, signal, rtol=rtol, atol=atol)

    def test_01_median_detrending(self):
        nsamps = 1000
        ndets = 3
        dtype = "float32"
        order = "C"

        x = np.linspace(0, 1, nsamps, dtype=dtype)
        signal = np.array([(i + 1) * np.sin(2*np.pi*x + i) for i in range(ndets)], dtype=dtype, order=order)

        signal_copy = signal.copy(order=order)
        signal_copy -= np.median(signal_copy, axis=-1)[..., None]

        method = "median"
        count = 0 # not used for median detrending
        so3g.detrend(signal, method, count)

        rtol = 0.
        atol = 1e-5
        np.testing.assert_allclose(signal_copy, signal, rtol=rtol, atol=atol)

    def test_02_linear_detrending(self):
        nsamps = 1000
        ndets = 3
        dtype = "float32"
        order = "C"
        count = nsamps // 3

        x = np.linspace(0., 1., nsamps, dtype=dtype)
        signal = np.array([(i + 1) * np.sin(2*np.pi*x + i) for i in range(ndets)], dtype=dtype, order=order)

        signal_copy = signal.copy(order=order)

        # this is the sotodlib "linear" detrending algorithm copied exactly
        count_copy = max(1, min(count, signal_copy.shape[-1] // 2))
        slopes = signal_copy[..., -count_copy:].mean(axis=-1, dtype=dtype) - signal[
            ..., :count_copy
        ].mean(axis=-1, dtype=dtype)

        # ignore shape != 2 case as c++ approach only supports 1D or 2D
        for i in range(signal_copy.shape[0]):
            signal_copy[i, :] -= slopes[i] * x

        signal_copy -= np.mean(signal_copy, axis=-1)[..., None]

        method = "linear"
        so3g.detrend(signal, method, count)

        rtol = 0.
        atol = 1e-5
        np.testing.assert_allclose(signal_copy, signal, rtol=rtol, atol=atol)


class TestBinning(unittest.TestCase):
    """
    Test binning.
    """

    def test_00_binning_no_flags_float64(self):
        nsamps = 10
        ndets = 1
        dtype = "float64"
        order = "C"

        x_min = 0.0
        x_max = 1.0
        bins = 2

        x = np.linspace(x_min, x_max, nsamps, dtype=dtype)
        signal = np.array([(i + 1) * np.sin(2*np.pi*x + i) for i in range(ndets)],
                          dtype=dtype, order=order)
        weight = np.ones((nsamps), dtype=dtype, order=order)

        bin_edges = np.histogram_bin_edges(x, bins=bins, range=[x_min,x_max],)
        bin_centers = (bin_edges[1] - bin_edges[0])/2. + bin_edges[:-1] # edge to center
        nbins = len(bin_centers)

        def numpy_binning():
            binned_signal = np.full([ndets, nbins], np.nan)
            binned_signal_squared_mean = np.full([ndets, nbins], np.nan)
            binned_signal_sigma = np.full([ndets, nbins], np.nan)

            # get bin indices
            bin_indices = np.digitize(x, bin_edges) - 1
            bin_indices = np.clip(bin_indices, 0, nbins-1)

            bin_counts, _ = np.histogram(x, bins=bins, range=[x_min,x_max], weights = weight)
            mcnts = bin_counts > 0

            for i in range(ndets):
                binned_signal[i][mcnts] = np.bincount(bin_indices, weights=signal[i]*weight, minlength=nbins
                                                     )[mcnts]/bin_counts[mcnts]
                binned_signal_squared_mean[i][mcnts] = np.bincount(bin_indices, weights=(signal[i]*weight)**2, minlength=nbins
                                                     )[mcnts]/bin_counts[mcnts]

            binned_signal_sigma[:, mcnts] = np.sqrt(np.abs(binned_signal_squared_mean[:,mcnts] - binned_signal[:,mcnts]**2)
                                          ) / np.sqrt(bin_counts[mcnts])
            bin_counts_dets = np.tile(bin_counts, (ndets, 1))

            return binned_signal, binned_signal_sigma, bin_counts_dets


        binned_signal, binned_signal_sigma, bin_counts_dets = numpy_binning()

        # so3g
        binned_signal_so3g = np.zeros((ndets, nbins), dtype=dtype, order=order)
        binned_signal_sigma_so3g = np.zeros((ndets, nbins), dtype=dtype, order=order)
        bin_counts_so3g = np.zeros((ndets, nbins), dtype=np.int32, order=order)

        so3g.bin_signal(x, signal, weight, binned_signal_so3g, binned_signal_sigma_so3g,
                        bin_counts_so3g, bin_edges, x_min, x_max)

        tolerance = 1e-10
        np.testing.assert_allclose(binned_signal, binned_signal, atol=tolerance)
        np.testing.assert_allclose(binned_signal_sigma_so3g, binned_signal_sigma, atol=tolerance)
        np.testing.assert_allclose(bin_counts_so3g, bin_counts_dets, atol=tolerance)

    def test_01_binning_flags_float32(self):
        nsamps = 1000
        ndets = 3
        dtype = "float32"
        order = "C"

        x_min = 0.0
        x_max = 1.0
        bins = 10

        x = np.linspace(x_min, x_max, nsamps, dtype=dtype)
        signal = np.array([(i + 1) * np.sin(2*np.pi*x + i) for i in range(ndets)],
                          dtype=dtype, order=order)
        weight = np.ones((ndets, nsamps), dtype=dtype, order=order)
        flags = np.zeros((ndets, nsamps), dtype=np.int32, order=order)
        flags[0,::2] = 1
        flags[1,::3] = 1
        flags[2,::4] = 1

        bin_edges = np.histogram_bin_edges(x, bins=bins, range=[x_min,x_max],)
        bin_centers = (bin_edges[1] - bin_edges[0])/2. + bin_edges[:-1] # edge to center
        nbins = len(bin_centers)

        def numpy_binning():
            binned_signal = np.full([ndets, nbins], np.nan)
            binned_signal_squared_mean = np.full([ndets, nbins], np.nan)
            binned_signal_sigma = np.full([ndets, nbins], np.nan)

            # get bin indices
            bin_indices = np.digitize(x, bin_edges) - 1
            bin_indices = np.clip(bin_indices, 0, nbins-1)
            bin_counts_dets = np.full([ndets, nbins], np.nan)

            if flags.shape == (ndets, nsamps):
                flag_is_2d = True
                m_2d = ~flags.astype(bool)
            elif flags.shape == (nsamps, ):
                flag_is_2d = False
                m = ~flags.astype(bool)

            for i in range(ndets):
                if flag_is_2d:
                    m = m_2d[i]

                if weight.shape == (ndets, nsamps):
                    weight_det = weight[i]
                elif weight.shape == (nsamps, ):
                    weight_det = weight

                bin_counts_dets[i] = np.bincount(bin_indices[m], weights=weight_det[m], minlength=nbins)
                mcnts = bin_counts_dets[i] > 0
                binned_signal[i][mcnts] = np.bincount(bin_indices[m], weights=signal[i][m]*weight_det[m], minlength=nbins
                                                     )[mcnts]/bin_counts_dets[i][mcnts]
                binned_signal_squared_mean[i][mcnts] = np.bincount(bin_indices[m], weights=(signal[i][m]*weight_det[m])**2, minlength=nbins
                                                     )[mcnts]/bin_counts_dets[i][mcnts]
                binned_signal_sigma[i][mcnts] = np.sqrt(np.abs(binned_signal_squared_mean[i,mcnts] - binned_signal[i,mcnts]**2)
                                                     ) / np.sqrt(bin_counts_dets[i][mcnts])

            return binned_signal, binned_signal_sigma, bin_counts_dets


        binned_signal, binned_signal_sigma, bin_counts_dets = numpy_binning()

        # so3g
        binned_signal_so3g = np.zeros((ndets, nbins), dtype=dtype, order=order)
        binned_signal_sigma_so3g = np.zeros((ndets, nbins), dtype=dtype, order=order)
        bin_counts_so3g = np.zeros((ndets, nbins), dtype=np.int32, order=order)

        so3g.bin_flagged_signal(x, signal, weight, binned_signal_so3g,
                                binned_signal_sigma_so3g, bin_counts_so3g,
                                bin_edges, x_min, x_max, flags)

        tolerance = 1e-4
        np.testing.assert_allclose(binned_signal, binned_signal, atol=tolerance)
        np.testing.assert_allclose(binned_signal_sigma_so3g, binned_signal_sigma, atol=tolerance)
        np.testing.assert_allclose(bin_counts_so3g, bin_counts_dets, atol=tolerance)

    def test_02_binning_flags_float64(self):
        nsamps = 1000
        ndets = 3
        dtype = "float64"
        order = "C"

        x_min = 0.0
        x_max = 1.0
        bins = 10

        x = np.linspace(x_min, x_max, nsamps, dtype=dtype)
        signal = np.array([(i + 1) * np.sin(2*np.pi*x + i) for i in range(ndets)],
                          dtype=dtype, order=order)
        weight = np.ones((ndets, nsamps), dtype=dtype, order=order)
        flags = np.zeros((ndets, nsamps), dtype=np.int32, order=order)
        # flags[0,::2] = 1
        # flags[1,::3] = 1
        # flags[2,::4] = 1

        bin_edges = np.histogram_bin_edges(x, bins=bins, range=[x_min,x_max],)
        bin_centers = (bin_edges[1] - bin_edges[0])/2. + bin_edges[:-1] # edge to center
        nbins = len(bin_centers)

        def numpy_binning():
            binned_signal = np.full([ndets, nbins], np.nan)
            binned_signal_squared_mean = np.full([ndets, nbins], np.nan)
            binned_signal_sigma = np.full([ndets, nbins], np.nan)

            # get bin indices
            bin_indices = np.digitize(x, bin_edges) - 1
            bin_indices = np.clip(bin_indices, 0, nbins-1)
            bin_counts_dets = np.full([ndets, nbins], np.nan)

            if flags.shape == (ndets, nsamps):
                flag_is_2d = True
                m_2d = ~flags.astype(bool)
            elif flags.shape == (nsamps, ):
                flag_is_2d = False
                m = ~flags.astype(bool)

            for i in range(ndets):
                if flag_is_2d:
                    m = m_2d[i]

                if weight.shape == (ndets, nsamps):
                    weight_det = weight[i]
                elif weight.shape == (nsamps, ):
                    weight_det = weight

                bin_counts_dets[i] = np.bincount(bin_indices[m], weights=weight_det[m], minlength=nbins)
                mcnts = bin_counts_dets[i] > 0
                binned_signal[i][mcnts] = np.bincount(bin_indices[m], weights=signal[i][m]*weight_det[m], minlength=nbins
                                                     )[mcnts]/bin_counts_dets[i][mcnts]
                binned_signal_squared_mean[i][mcnts] = np.bincount(bin_indices[m], weights=(signal[i][m]*weight_det[m])**2, minlength=nbins
                                                     )[mcnts]/bin_counts_dets[i][mcnts]
                binned_signal_sigma[i][mcnts] = np.sqrt(np.abs(binned_signal_squared_mean[i,mcnts] - binned_signal[i,mcnts]**2)
                                                     ) / np.sqrt(bin_counts_dets[i][mcnts])

            return binned_signal, binned_signal_sigma, bin_counts_dets


        binned_signal, binned_signal_sigma, bin_counts_dets = numpy_binning()

        # so3g
        binned_signal_so3g = np.zeros((ndets, nbins), dtype=dtype, order=order)
        binned_signal_sigma_so3g = np.zeros((ndets, nbins), dtype=dtype, order=order)
        bin_counts_so3g = np.zeros((ndets, nbins), dtype=np.int32, order=order)

        so3g.bin_flagged_signal(x, signal, weight, binned_signal_so3g,
                                binned_signal_sigma_so3g, bin_counts_so3g,
                                bin_edges, x_min, x_max, flags)

        tolerance = 1e-10
        np.testing.assert_allclose(binned_signal, binned_signal, atol=tolerance)
        np.testing.assert_allclose(binned_signal_sigma_so3g, binned_signal_sigma, atol=tolerance)
        np.testing.assert_allclose(bin_counts_so3g, bin_counts_dets, atol=tolerance)


if __name__ == "__main__":
    unittest.main()
