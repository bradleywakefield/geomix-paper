"""Grouped Vecchia approximation for scalable GP likelihood evaluation.

References
----------
Katzfuss, M. & Guinness, J. (2021). A general framework for Vecchia
approximations of Gaussian processes. *Statistical Science*, 36(1), 124-141.
"""

import numpy as np
from scipy.spatial import KDTree


def find_vecchia_neighbours(coords: np.ndarray, m: int) -> np.ndarray:
    """Find the *m* nearest conditioning neighbours for each location.

    Locations are ordered by their index, so location ``i`` is conditioned on
    the ``min(i, m)`` previously ordered locations that are closest in space.

    Parameters
    ----------
    coords : np.ndarray of shape (n, d)
        Spatial coordinates.
    m : int
        Maximum number of conditioning neighbours.

    Returns
    -------
    np.ndarray of shape (n, m)
        Indices of conditioning neighbours; padded with ``-1`` where fewer
        than ``m`` predecessors exist.
    """
    n = coords.shape[0]
    neighbours = np.full((n, m), -1, dtype=int)
    for i in range(1, n):
        candidate_idx = list(range(i))
        if len(candidate_idx) <= m:
            neighbours[i, : len(candidate_idx)] = candidate_idx
        else:
            subset = coords[candidate_idx]
            _, local_idx = KDTree(subset).query(coords[i], k=m)
            neighbours[i] = [candidate_idx[j] for j in local_idx]
    return neighbours
