"""This submodule is for functions that produce disjoint RangesMatrix
objects suitable for use with the Projection code and OpenMP.

"""

import so3g
import numpy as np
from pixell import enmap

def get_num_threads(n_threads=None):
    if n_threads is None:
        return so3g.useful_info()['omp_num_threads']
    return n_threads


def get_threads_domdir(sight, offs, shape, wcs, n_threads=None,
                       offs_rep=None, plot_prefix=None):
    """Assign samples to threads according to the dominant scan
    direction.

    Arguments:
      sight (CelestialSightLine): The boresight pointing.
      offs (array of quaternions): The detector pointing.
      shape (tuple): The map shape.
      wcs (wcs): The map WCS.
      n_threads (int): The number of threads to target.
      offs_rep (array of quaternions): A representative set of
        detectors, for determining scan direction and relative weights
        (if not present then offs is used for this).
      plot_prefix (str): If present, pylab will be imported and some
        plots saved with this name as prefix.

    Returns:
      Ranges matrix with shape (n_thread, n_det, n_samp), that can be
      passes as the "threads" argument to Projectionist routines.

    """
    n_threads = get_num_threads(n_threads)
    if plot_prefix:
        import pylab as pl
    if offs_rep is None:
        offs_rep = offs

    xi, eta, gamma = so3g.proj.quat.decompose_xieta(offs)
    # Align the alternative detectors parallel to lines of constant el.
    offs_xl = np.array(so3g.proj.quat.rotation_xieta(xi, eta, gamma*0 + 90*so3g.proj.DEG))
    pmat   = so3g.proj.wcs.Projectionist.for_geom(shape, wcs)
    ass    = so3g.proj.Assembly.attach(sight, offs_xl)

    sig = np.ones((len(offs_xl), len(ass.Q)), dtype='float32')
    m = pmat.to_map(sig, ass, comps='TQU')
    del sig

    # Dominant angle?  This assumes that +Q is parallel to map
    # columns, and U increases along the map diagonal.
    T, Q, U = [_m.sum() for _m in m]
    phi = np.arctan2(U, Q) / 2

    if plot_prefix:
        text = 'Qf=%.2f Uf=%.2f phi=%.1f deg' % (Q/T, U/T, phi / so3g.proj.DEG)
        for i in range(3):
            pl.imshow(m[i], origin='lower')
            pl.title('TQU'[i] + ' ' + text)
            pl.colorbar()
            pl.savefig(plot_prefix + '%02i.png' % i)
            pl.clf()

    # Index the map in stripes parallel to phi
    y, x = np.arange(m.shape[-2]), np.arange(m.shape[-1])
    idx = y[:,None]*np.sin(phi) - x[None,:]*np.cos(phi)

    if plot_prefix:
        pl.imshow(idx, origin='lower')
        pl.colorbar()
        pl.savefig(plot_prefix + '10.png')
        pl.clf()

    # Bins.
    n_bins = 100
    bins = np.linspace(idx.min(), idx.max(), n_bins+1)
    H, _ = np.histogram(idx.ravel(), weights=m[0].ravel(), bins=bins)
    H = H.cumsum()
    H /= H.max()

    # Threadify.
    superbins = [bins[0]]
    for i in range(1, n_threads):
        j = (H * n_threads >= i).nonzero()[0][0]
        superbins.append(bins[j])

    superbins.append(bins[-1])

    # Convert the idx map.
    tidx_map = enmap.zeros(idx.shape, wcs)
    for i in range(n_threads):
        s = (superbins[i] <= idx)*(idx < superbins[i+1])
        tidx_map[s] = i

    # Turn that into threads!    
    threads = pmat.assign_threads_from_map(ass, tidx_map[None])

    if plot_prefix:
        pl.plot(bins[:-1], H)
        for b in superbins[1:-1]:
            pl.axvline(b)
        pl.savefig(plot_prefix + '12.png')
        pl.clf()
        pl.imshow(tidx_map, origin='lower')
        pl.colorbar()
        pl.savefig(plot_prefix + '13.png')
        pl.clf()

    return threads

