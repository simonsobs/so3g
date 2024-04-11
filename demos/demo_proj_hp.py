import numpy as np
import so3g
import healpy as hp
from matplotlib import pyplot as plt

def main():
    ######## Set up a simple TOD simulation ########
    ## Scanning params ##
    scan_rate = 1. # deg/s
    sample_rate = 50 # 1/s
    field_width = 30 # deg
    el0 = 30
    tmax = 4000

    ## Focal plane params ##
    ndets = 1000
    rmax = 0.1 # max val for xi, eta

    ## Other TOD params ##
    isignal = 0 # What coord to use as signal. 0=ra, 1=dec

    ## Mapmaking params ##
    nside = 512
    nside_tile = 16
    det_weights = None

    ## Simulate scanning and make TOD ##
    print("Simulating TOD")
    # Constant el scan
    time, az, el = sim_scan(tmax, sample_rate, scan_rate, field_width, el0)
    # Scalar signal equal to the RA. Dets are randomly placed on a circular focal plane
    signal, asm = make_tod(time, az, el, ndets, rmax, isignal)
    signal = np.asarray(signal, dtype='float32') # This is required

    ######## Do simple 1-component mapmaking  ########
    print("Making untiled map")
    ## Make an untiled map ##
    proj = so3g.proj.ProjectionistHealpix.for_healpix(nside, nside_tile=None)
    proj.assign_threads(asm, 'simple') ## Assign threads
    rhs = proj.to_map(signal, asm, det_weights=det_weights, comps='T')
    weights = proj.to_weights(asm, det_weights=det_weights, comps='T')
    # Now simple combination
    iweights = np.zeros_like(weights)
    idx = (weights != 0)    
    iweights[idx] = 1/(weights[idx])
    imap = (rhs * iweights[0])[0]

    ## Make a tiled map ##
    print("Making tiled map")
    proj = so3g.proj.ProjectionistHealpix.for_healpix(nside, nside_tile=nside_tile)
    proj.assign_threads(asm, 'tiles') ## Assign threads
    # Output of to_map and to_weights are len(nTile) lists of None (inactive tiles) or
    # arrays of the tile shape. Use untile_healpix to recover the full map
    rhs = untile_healpix(proj.to_map(signal, asm, det_weights=det_weights, comps='T'))
    weights = untile_healpix(proj.to_weights(asm, det_weights=det_weights, comps='T'))
    iweights = np.zeros_like(weights)
    idx = (weights != 0)    
    iweights[idx] = 1/(weights[idx])
    imap_tiled = (rhs * iweights[0])[0]

    ## Check they are the same ##
    assert np.array_equal(imap, imap_tiled)

    ## Plot ##
    imap[imap==0] = hp.UNSEEN
    hp.mollview(imap, nest=True)    
    plt.show()

    ## Plot gnomview ##
    ras, decs, _ = so3g.proj.quat.decompose_lonlat(asm.Q)
    mean_ra, mean_dec = np.mean(ras), np.mean(decs)
    width = max(np.max(ras)-mean_ra, mean_ra-np.min(ras))
    center = np.rad2deg([mean_ra, mean_dec])
    reso = 1.5 # arcmin
    width_pix = int(np.rad2deg(width) * 60. / reso) * 2
    hp.gnomview(imap, nest=True, rot=center, xsize=width_pix*1.4, reso=reso)
    plt.show()
    

######## Helper functions for tiled healpix maps ########

def untile_healpix(tiled):
    compressed = untile_healpix_compressed(tiled)
    active_tiles = get_active_tiles(tiled)
    return decompress_healpix(compressed, active_tiles)

def untile_healpix_compressed(tiled):
    """
    Get a compressed healpix map by removing Nones from a tiled map
    Input: list of tiles. None or (..., npix) in NEST ordering
    Out:   *np.array (nActiveTile, ..., npix) of all active tiles
    """
    return np.array([tiled[ii] for ii in range(len(tiled)) if tiled[ii] is not None])

def get_active_tiles(tiled):
    return np.array([tile is not None for tile in tiled])

def decompress_healpix(compressed, active_tiles):
    """Decompress a healpix map
    Input: See outputs of untile_healpix_compressed
    Output: np.array: Full hp map in nest
    """
    tile_shape = compressed[0].shape
    npix_per_tile = tile_shape[-1]
    super_shape = tile_shape[:-1]
    npix = len(active_tiles) * npix_per_tile # ntiles * npix_per_tile
    out = np.zeros(super_shape + (npix,))
    tile_inds = [ii for ii in range(len(active_tiles)) if active_tiles[ii]]
    for ii in range(len(compressed)):
        tile_ind = tile_inds[ii]
        out[..., npix_per_tile * tile_ind : npix_per_tile * (tile_ind+1)] = compressed[ii]
    return out

######## Helper functions to simulate a simple focal plane and TOD ########

def simulate_detectors(ndets, rmax, seed=0):
    np.random.seed(seed)
    phi = np.random.random(ndets) * (2*np.pi)
    rr =  np.sqrt(np.random.random(ndets)) * rmax # P(r) \propto r
    gamma = np.random.random(ndets) * (2*np.pi)
    xi = rr * np.cos(phi)
    eta = rr * np.sin(phi)
    return xi, eta, gamma

def make_fp(xi, eta, gamma):
    fp = so3g.proj.quat.rotation_xieta(xi, eta, gamma)
    so3g_fp = so3g.proj.FocalPlane()
    for i, q in enumerate(fp):
        so3g_fp[f'a{i}'] = q
    return so3g_fp

def extract_coord(sight, fp, isignal=0, groupsize=100):
    coord_out = []
    for ii in range(0, len(fp), groupsize):
        coords = np.array(sight.coords(fp[ii:ii+groupsize])) # (groupsize, nsamples, (ra, dec, cos(phi), sin(phi)))
        coord_out.append(coords[..., isignal])
    coord_out = np.concatenate(coord_out) # ndets, nsamples
    return coord_out

def make_signal(ndets, rmax, sight, isignal, seed=0):
    xi, eta, gamma = simulate_detectors(ndets, rmax, seed)
    fp = so3g.proj.quat.rotation_xieta(xi, eta, gamma)
    signal = extract_coord(sight, fp, isignal)
    return signal

def make_tod(time, az, el, ndets, rmax, isignal, seed=0):
    sight = so3g.proj.CelestialSightLine.naive_az_el(time, az, el)
    signal = make_signal(ndets, rmax, sight, isignal, seed)
    so3g_fp = make_fp(*simulate_detectors(ndets, rmax, seed))
    asm = so3g.proj.Assembly.attach(sight, so3g_fp)
    return signal, asm

def sim_scan(tmax, sample_rate, scan_rate, field_width, el0):
    time = np.arange(0, tmax, 1/sample_rate)
    az0 = time * scan_rate - field_width
    az = np.abs(az0 % (field_width * 2) - field_width)
    el = np.ones_like(az) * el0
    return time, az, el

##
if __name__ == '__main__':
    main()
