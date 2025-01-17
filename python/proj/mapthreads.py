"""This submodule is for functions that produce "thread assignments",
i.e. disjoint RangesMatrix objects with shape (n_threads, n_dets,
n_samps) that are suitable for the threads= argument in Projection
code that projects from time to map space using OpenMP
parallelization.

"""

import so3g
import numpy as np

def get_num_threads(n_threads=None):
    """Utility function for computing n_threads.  If n_threads is not
    None, it is returned directly.  But if it is None, then the OpenMP
    thread count is returned.  Uses so3g.useful_info().

    """
    if n_threads is None:
        return so3g.useful_info()['omp_num_threads']
    return n_threads

def get_threads_domdir(sight, fplane, shape, wcs, tile_shape=None,
                       active_tiles=None, n_threads=None, fplane_rep=None,
                       plot_prefix=None):
    """Assign samples to threads according to the dominant scan
    direction.

    The dominant scan direction is first determined, which requires
    creating a 3-component map.  This projection operation can't be
    parallelized safely so it is wise to pass in a decimated detector
    set through offs_rep.  The second step is to assign pixels to
    threads; to do this a signal-sized array is projected from map
    space.  In the present implementation this uses all the detectors
    (but takes advantage of threads).  In principle this step could be
    done with a decimated detector set, though with care that the
    decimated subset covers the array well, spatially.

    Arguments:
      sight (CelestialSightLine): The boresight pointing.
      fplane (FocalPlane): The detector pointing
      shape (tuple): The map shape.
      wcs (wcs): The map WCS.
      tile_shape (tuple): The tile shape, if this should be done using
        tiles.
      active_tiles (list): The list of active tiles.  If None, this
        will be computed along the way.
      n_threads (int): The number of threads to target (defaults to
        OMP_NUM_THREADS).
      fplane_rep (FocalPlane): A representative set of
        detectors, for determining scan direction (if not present then
        fplane is used for this).
      plot_prefix (str): If present, pylab will be imported and some
        plots saved with this name as prefix.

    Returns:
      Ranges matrix with shape (n_thread, n_det, n_samp), that can be
      passes as the "threads" argument to Projectionist routines.

    """
    if plot_prefix:
        import pylab as pl
        if tile_shape is None:
            tile_iter = lambda maps: [('', maps[0])]
        else:
            tile_iter = lambda maps: [('_tile%04i' % i, _m)
                                      for i, _m in enumerate(maps)
                                      if _m is not None]

    n_threads = get_num_threads(n_threads)
    if fplane_rep is None:
        fplane_rep = fplane

    if tile_shape is None:
        # Let's pretend it is, though; this simplifies logic below and
        # doesn't cost much.
        tile_shape = shape
        active_tiles = [0]

    # The full assembly, for later.
    asm_full = so3g.proj.Assembly.attach(sight, fplane)

    # Get a Projectionist -- note it can be used with full or
    # representative assembly.
    pmat = so3g.proj.wcs.Projectionist.for_tiled(
        shape, wcs, tile_shape=tile_shape, active_tiles=active_tiles
    )
    if active_tiles is None:
        # This is OMPed, but it's still better to pass in active_tiles
        # if you've already computed it somehow.
        pmat.active_tiles = pmat.get_active_tiles(asm_full, assign=True)['active_tiles']
    tiling = pmat.tiling

    # For the scan direction map, use the "representative" subset
    # detectors, with polarization direction aligned parallel to
    # elevation.
    xi, eta, gamma = so3g.proj.quat.decompose_xieta(fplane_rep.quats)
    fplane_xl = so3g.proj.FocalPlane.from_xieta(xi, eta, gamma*0+90*so3g.proj.DEG)
    asm_rep = so3g.proj.Assembly.attach(sight, fplane_xl)
    sig = np.ones((fplane_xl.ndet, len(asm_rep.Q)), dtype='float32')
    scan_maps = pmat.to_map(sig, asm_rep, comps='TQU')

    # Compute the scan angle, based on Q and U weights.  This assumes
    # that +Q is parallel to map columns, and U increases along the
    # map diagonal.
    T, Q, U = np.sum([_m.reshape((3, -1)).sum(axis=-1)
                      for _m in scan_maps if _m is not None], axis=0)
    phi = np.arctan2(U, Q) / 2

    if plot_prefix:
        text = 'Qf=%.2f Uf=%.2f phi=%.1f deg' % (Q/T, U/T, phi / so3g.proj.DEG)
        for label, _m in tile_iter(scan_maps):
            for i in range(3):
                pl.imshow(_m[i], origin='lower')
                pl.title(label + ' '  + 'TQU'[i] + ' ' + text)
                pl.colorbar()
                pl.savefig(plot_prefix + '%02i%s%s.png' % (i, 'TQU'[i], label))
                pl.clf()

    # Now paint a map with a single parameter such that contours are
    # parallel to the scan direction.
    idx_maps = pmat.zeros((1,))
    lims = None
    for t in pmat.active_tiles:
        y0, x0 = tiling.tile_offset(t)
        ny, nx = tiling.tile_dims(t)
        y, x = y0 + np.arange(ny), x0 + np.arange(nx)
        idx_maps[t] = y[:,None]*np.sin(phi) - x[None,:]*np.cos(phi)
        _lo, _hi = idx_maps[t].min(), idx_maps[t].max()
        if lims is None:
            lims = _lo, _hi
        else:
            lims = min(_lo, lims[0]), max(_hi, lims[1])

    if plot_prefix:
        for label, _m in tile_iter(idx_maps):
            pl.imshow(_m, origin='lower')
            pl.colorbar()
            pl.savefig(plot_prefix + '10param%s.png' % label)
            pl.clf()

    # Histogram the parameter, weighted by the weights map -- this
    # tells us the number of hits for ranges of the parameter value.
    n_bins = 200
    bins = np.linspace(lims[0], lims[1], n_bins+1)
    H = np.zeros(n_bins)
    for t in pmat.active_tiles:
        H += np.histogram(idx_maps[t].ravel(), bins=bins,
                          weights=scan_maps[t][0].ravel())[0]
    del scan_maps

    # Turn histogram into a CDF.  The CDF is non-decreasing, but could
    # have some flat parts in it, so bias the whole curve to guarantee
    # it's strictly increasing.
    H = np.hstack((0, H.cumsum()))
    bias = .00001
    H = H / H.max() * (1-bias) + np.linspace(0, bias, len(H))

    # Create n_threads "super_bins" that contain roughly equal weight.
    xs = np.linspace(0, 1, n_threads + 1)
    ys = np.interp(xs, H, bins)
    superbins = np.hstack((bins[0], ys[1:-1], bins[-1]))

    # Create maps where the pixel value is a superbin index.
    for t in pmat.active_tiles:
        temp = np.zeros(idx_maps[t].shape)
        for i in range(n_threads):
            s = (superbins[i] <= idx_maps[t])*(idx_maps[t] < superbins[i+1])
            temp[s] = i
        idx_maps[t][:] = temp
        idx_maps[t].shape = (1, ) + tuple(idx_maps[t].shape)

    # Turn the superbin index maps into thread assignments.
    threads = pmat.assign_threads_from_map(asm_full, idx_maps, n_threads=n_threads)

    if plot_prefix:
        pl.plot(bins, H)
        for _x, _y in zip(superbins[1:-1], xs[1:-1]):
            pl.plot([_x, _x], [-1, _y], c='k')
            pl.plot([superbins[0], superbins[-1]], [_y, _y], c='gray')
        pl.xlim(superbins[0], superbins[-1])
        pl.ylim(-.01, 1.01)
        pl.savefig(plot_prefix + '20histo.png')
        pl.clf()

        for label, _m in tile_iter(idx_maps):
            pl.imshow(_m[0], origin='lower')
            pl.colorbar()
            pl.savefig(plot_prefix + '30thread%s.png' % label)
            pl.clf()

    return threads

