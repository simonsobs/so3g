import numpy as np

from spt3g import core


def get_unix_time(g3_time):
    """Convert a G3Time or G3VectorTime time object to a unix timestamp or
    numpy vector (double) of unix timestamps."""
    if isinstance(g3_time, core.G3Time):
        return g3_time.time / core.G3Units.seconds
    if isinstance(g3_time, core.G3VectorTime):
        output = np.array([t.time for t in g3_time], dtype=float)
        output /= core.G3Units.seconds
        return output


def get_g3_time(unix_time):
    """Convert a double or numpy array of floats to G3Time or
    G3VectorTime."""
    src = None
    if isinstance(unix_time, core.G3VectorDouble):
        src = (np.array(unix_time) * core.G3Units.seconds).astype('int')
    elif isinstance(unix_time, np.ndarray) and unix_time.ndim == 1:
        src = (unix_time * core.G3Units.seconds).astype('int')
    if src is not None:
        return core.G3VectorTime([core.G3Time(t) for t in src])
    return core.G3Time(int(unix_time * core.G3Units.seconds))
