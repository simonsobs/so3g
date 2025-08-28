import numpy as np

from spt3g.core import quat, G3VectorQuat


"""We are using the spt3g quaternion containers,
i.e. cu3g.G3VectorQuat and cu3g.quat.  One way these are nice is that
they have accelerated quaternion arithmetic with operator overloading.
The component ordering is (1,i,j,k).  We are restricted to 0d or 1d,
and that's fine."""

DEG = np.pi/180


def euler(axis, angle):
    """
    The quaternion representing of an Euler rotation.
    
    For example, if axis=2 the computed quaternion(s) will have
    components:

        q = (cos(angle/2), 0, 0, sin(angle/2))
    
    Parameters
    ----------
    axis : {0, 1, 2}
        The index of the cartesian axis of the rotation (x, y, z).
    angle : float or 1-d float array
        Angle of rotation, in radians.
    
    Returns
    -------
    quat or G3VectorQuat, depending on ndim(angle).
    """
    # Either angle or axis or both can be vectors.
    angle = np.asarray(angle)
    shape = np.broadcast(axis, angle).shape + (4,)
    c, s = np.cos(angle/2), np.sin(angle/2)
    q = np.zeros(shape)
    q[..., 0] = c
    q[..., axis+1] = s
    if len(shape) == 1:
        return quat(*q)
    return G3VectorQuat(q)


def rotation_iso(theta, phi, psi=None):
    """Returns the quaternion that composes the Euler rotations:
   
        Qz(phi) Qy(theta) Qz(psi)

    Note arguments are in radians.
    """
    output = euler(2, phi) * euler(1, theta)
    if psi is None:
        return output
    return output * euler(2, psi)


def rotation_lonlat(lon, lat, psi=0., azel=False):
    """Returns the quaternion that composes the Euler rotations:
        
        Qz(lon) Qy(pi/2 - lat) Qz(psi)
    
    Note the three angle arguments are in radians.

    If azel is True, then the sign of lon is flipped (as though lon
    and lat were azimuth and elevation).

    """
    if azel:
        return rotation_iso(np.pi/2 - lat, -lon, psi)
    return rotation_iso(np.pi/2 - lat, lon, psi)

def rotation_xieta(xi, eta, gamma=0):
    """Returns the quaternion that rotates the center of focal plane to
    (xi, eta, gamma).  This is equivalent to composed Euler rotations:

        Qz(phi) Qy(theta) Qz (psi)

    where

        xi = - sin(theta) * sin(phi)
        eta = - sin(theta) * cos(phi)
        gamma = psi + phi

    Note arguments are in radians.

    """
    phi = np.arctan2(-xi, -eta)
    theta = np.arcsin(np.sqrt(xi**2 + eta**2))
    psi = gamma - phi
    return rotation_iso(theta, phi, psi)

def decompose_iso(q):
    """Decomposes the rotation encoded by q into the product of Euler
    rotations:
    
        q = Qz(phi) Qy(theta) Qz(psi)
    
    Parameters
    ----------
    q : quat or G3VectorQuat
        The quaternion(s) to be decomposed.
    
    Returns
    -------
    (theta, phi, psi) : tuple of floats or of 1-d arrays
        The rotation angles, in radians.
    """

    if isinstance(q, quat):
        a,b,c,d = q.a, q.b, q.c, q.d
    else:
        a,b,c,d = np.transpose(q)

    psi = np.arctan2(a*b+c*d, a*c-b*d)
    phi = np.arctan2(c*d-a*b, a*c+b*d)

    # There are many ways to get theta; this is probably the most
    # expensive but the most robust.
    theta = 2 * np.arctan2((b**2 + c**2)**.5, (a**2 + d**2)**.5)

    return (theta, phi, psi)

def decompose_lonlat(q, azel=False):
    """Like decompose_iso, but returns (lon, lat, psi) assuming that the
    input quaternion(s) were constructed as rotation_lonlat(lon, lat,
    psi).

    If azel is True, returns (-lon, lat, psi), so that first two
    coords are easily interpreted as an azimuth and elevation.

    """
    theta, phi, psi = decompose_iso(q)
    if azel:
        return (-phi, np.pi/2 - theta, psi)
    return (phi, np.pi/2 - theta, psi)

def decompose_xieta(q):
    """Like decompose_iso, but returns (xi, eta, gamma) assuming that the
    input quaternion(s) were constructed as rotation_xieta(xi, eta,
    gamma).

    """
    if isinstance(q, quat):
        a,b,c,d = q.a, q.b, q.c, q.d
    else:
        a,b,c,d = np.transpose(q)
    return (
        2 * (a * b - c * d),
        -2 * (a * c + b * d),
        np.arctan2(2 * a * d, a * a - d * d)
    )
