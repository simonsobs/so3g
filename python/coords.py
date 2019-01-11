"""
Coordinate system support, especially for using quaternions in the
context of pointing and astrometry.

Common sense alone cannot guide one in the choice of sign conventions
and the definitions of the cartesian axes in various systems.  Consult
the "coord_sys" documents in the SO repo:

   https://github.com/simonsobs/tod2maps_docs/
"""

import so3g

from spt3g import core
from spt3g import coordinateutils as cu3g

import numpy as np

def q_euler(axis, angle):
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
    c, s = np.cos(angle/2), np.sin(angle/2)
    shape = np.shape(c) + (4,)
    q = np.zeros(shape)
    q[...,0] = c
    q[...,axis+1] = s
    if len(shape) == 1:
        return cu3g.quat(*q)
    return cu3g.G3VectorQuat(q)

def decompose(q, axes=[0,1,2], cycle=None):
    """
    Decompose the rotation encoded by q into the product:

       q = q0 * q_euler(axes[0], angles[0]) *
                ...
                q_euler(axes[n-1], angles[n-1])

    Where n = len(axes).  Typically this would be used with axes=(some
    permutation of [0,1,2]).  But any list of axes can be used, with
    the caveat that it may be necessary to specify cycle='even' or
    cycle='odd' in some cases (see Notes).

    Parameters
    ----------
    q : quat or G3VectorQuat
        The quaternion(s) to be decomposed.
    axes : list drawn from {0,1,2}
        The axes of the rotations to be extracted.
    cycle : {None,'even','odd'}
        Chooses between equations that assume "even" vs. "odd"
        permutation of axes=[0,1,2].  Often optional; see notes.

    Returns
    -------
    q0 : quat or G3VectorQuat
        The residual quaternion, as described in the equation above.
    angles : list of float
        The euler angles, as described in the equation above.

    Notes
    -----
    Because rotations are not commutative, there is no general way to
    "project out" a single eulier rotation from a general rotation.
    The 'cycle' argument allows one to select between forms of the
    equations that are appropriate for even or odd permutations of
    axes=[0,1,2].  This is used in the internal recursion, and
    similarly it may be necessary if a user wants to pass in only one
    axis at a time.

    For example, the following should be equivalent (note that
    axes=[1,0,2] is an odd permutation):

    q0, (beta, alpha, gamma) = decompose(q, [1,0,2])

    q2, gamma = decompose(q,  2, 'odd')
    q1, alpha = decompose(q2, 0, 'odd')
    a0, beta  = decompose(q1, 1, 'odd')
    """
    
    if cycle is None:
        assert(np.ndim(axes) > 0)  # Choose a cycle sign if in single axis mode.
        cycle, i = 0, len(axes)-1
        while cycle == 0 and i > 0:
            cycle = (axes[i] - axes[i-1])
            i -= 1
        assert(cycle != 0)         # Must have 2 different axes to guess cycle.
        cycle = {1: 'even', 2: 'odd'}[cycle % 3]
    if np.ndim(axes) > 0:
        angles = []
        for ax in axes[::-1]:
            q, angle = decompose(q, ax, cycle=cycle)
            angles.insert(0, angle)
        return q, np.array(angles)
    qa = np.asarray(q)
    if cycle == 'even':
        i0, i1, i2, i3 = 0, 1 + axes, 1 + (axes+2)%3, 1 + (axes+1)%3
        angle = np.arctan2(qa[...,i0]*qa[...,i1] - qa[...,i2]*qa[...,i3],
                           0.5 - qa[...,i1]**2 - qa[...,i2]**2)
    elif cycle == 'odd':
        i0, i1, i2, i3 = 0, 1 + axes, 1 + (axes+1)%3, 1 + (axes+2)%3
        angle = np.arctan2(qa[...,i0]*qa[...,i1] + qa[...,i2]*qa[...,i3],
                           0.5 - qa[...,i1]**2 - qa[...,i2]**2)
    else:
        raise ValueError('cycle argument must be "even", "odd", or None.')
    return q * q_euler(axes, -angle), angle

def _is_null_rotation(q, tol=1e-6):
    """
    Returns true, element by element, if q is equal (within tol) to
    the identity rotation.
    """
    q = np.asarray(q)
    q_mag = ((q - [1,0,0,0])**2).sum(axis=-1)**.5
    return q_mag < tol

if __name__ == '__main__':
    #
    # This demo is destined for unit-testing...
    #

    # Check that our Euler operations are invertible, per axis.
    psi = .1
    for axis in [0,1,2]:
        print('Checking axis=%i, angle=%.3f' % (axis, psi))
        q = q_euler(axis, np.array([psi]))
        for cycle in ['even', 'odd']:
            q1, psi1 = decompose(q, axis, cycle)
            print(' ... %s cycle form recovers: %.3f' % (cycle, psi1))
            assert(_is_null_rotation(q1) and abs(psi1 - psi) < 1e-6)
        
    # Check that composed rotations are decomposable.
    angles = [.3, .2, .1]
    print('Checking all 3-axis permutations with angles %s' % angles)
    for sign in [1,-1]:
        for start in [0,1,2]:
            axes = [(start+sign*i)%3 for i in range(3)]
            q = cu3g.G3VectorQuat([cu3g.quat(1,0,0,0)])
            for ax,an in zip(axes, angles):
                qmod = q_euler(ax, an)
                q = q * qmod
            q0, angles = decompose(q, axes)
            print(' perm %s recovers: %s' % (axes, angles.T))
            assert(_is_null_rotation(q0))
