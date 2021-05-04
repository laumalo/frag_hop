"""
This module containts all the functions related to handle geometry operations on
the molecules trajectories.
"""
import numpy as np

def rotation_matrix_from_vectors(vec1, vec2):
    """
    Finds the rotation matrix that aligns vec1 to vec2
    Parameters
    ----------
    vec1: np.array
        A 3d "source" vector
    vec2: np.array
        A 3d "destination" vector

    Returns
    -------
    rotation_matrix: np.array
        A transform matrix (3x3) which when applied to vec1,
        aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), \
           (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + \
        kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def rotation_matrix_axis(radi, u):
    """
    Finds the general rotation matrix of a certain angle along an axis
    defined by a vector.

    Parameters
    ----------
    radi : float
        Rotation angle, in radians.
    u : np.array
        A 3d vector defining the rotation axis.

    Returns
    -------
    rotation_matrix: np.array
        A transform matrix (3x3) which when applied to a point, it rotates
        a certain angle along an axis.
    """
    ux, uy, uz = u
    rotation_matrix =  np.array([[np.cos(radi) + ux**2 * (1 - np.cos(radi)),
                      ux * uy * (1 - np.cos(radi)) - uz * np.sin(radi),
                      ux * uz * (1 - np.cos(radi)) + uy * np.sin(radi)],
                     [uy * ux * (1 - np.cos(radi)) + uz * np.sin(radi),
                         np.cos(radi) + uy**2 * (1 - np.cos(radi)),
                      uy * uz * (1 - np.cos(radi)) - ux * np.sin(radi)],
                     [uz * ux * (1 - np.cos(radi)) - uy * np.sin(radi),
                      uz * uy * (1 - np.cos(radi)) + ux * np.sin(radi),
                         np.cos(radi) + uz**2 * (1 - np.cos(radi))]],
                    dtype=np.double)
    return rotation_matrix

def distance(x, y):
    """
    It computes the distance between two points.

    Parameters
    ----------
    x : np.array
        Point
    y : np.array
        Point

    Returns
    -------
    dist : float
        Distance between the two points
    """
    import math
    return math.sqrt((x[0] - y[0]) ** 2 +
                     (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2)
