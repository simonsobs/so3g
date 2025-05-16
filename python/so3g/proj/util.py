import numpy as np


def ces(el, az0, throw, v_scan, t):
    """Generate a CES scan pattern.

    Arguments:
      el: The elevation angle.
      az0: The central azimuth angle.
      throw: The azimuth half-scan amplitude.
      v_scan: The scan speed.
      t: A vector of timestamps.

    Returns:
      Vectors (az,el).  Each is the same length as t.

    The angle and time units do not matter as long as they are
    consistent.  For example, use degrees for the angles, seconds for
    time, and deg/s for v_scan.

    """
    phase = (t-t[0]) * v_scan % (4*throw)
    az = phase
    az[phase>2*throw] = 4*throw - phase[phase>2*throw]
    return az - throw + az0, el + az*0
