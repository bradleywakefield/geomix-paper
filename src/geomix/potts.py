"""Potts Markov random field for spatially coherent stratigraphy."""

import numpy as np


def potts_energy(labels: np.ndarray, adjacency: list[list[int]], beta: float) -> float:
    """Compute the Potts MRF energy for a given label configuration.

    The energy is defined as

    .. math::

        U(\\mathbf{z}) = -\\beta \\sum_{(i,j) \\in \\mathcal{E}}
        \\mathbf{1}[z_i = z_j]

    Parameters
    ----------
    labels : np.ndarray of shape (n,)
        Integer class labels for each node.
    adjacency : list of lists
        ``adjacency[i]`` is a list of node indices neighbouring node ``i``.
    beta : float
        Spatial smoothing parameter (inverse temperature).

    Returns
    -------
    float
        Potts energy (lower is more spatially coherent).
    """
    energy = 0.0
    for i, neighbours in enumerate(adjacency):
        for j in neighbours:
            if j > i:  # count each edge once
                energy -= beta * float(labels[i] == labels[j])
    return energy


def build_grid_adjacency(nx: int, ny: int) -> list[list[int]]:
    """Build a 4-connectivity adjacency list for a 2-D regular grid.

    Parameters
    ----------
    nx : int
        Number of grid nodes in the x direction.
    ny : int
        Number of grid nodes in the y direction.

    Returns
    -------
    list of lists
        Adjacency list of length ``nx * ny``.
    """
    adjacency: list[list[int]] = [[] for _ in range(nx * ny)]
    for ix in range(nx):
        for iy in range(ny):
            node = ix * ny + iy
            if ix > 0:
                adjacency[node].append((ix - 1) * ny + iy)
            if ix < nx - 1:
                adjacency[node].append((ix + 1) * ny + iy)
            if iy > 0:
                adjacency[node].append(ix * ny + iy - 1)
            if iy < ny - 1:
                adjacency[node].append(ix * ny + iy + 1)
    return adjacency
