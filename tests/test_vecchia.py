"""Tests for the Vecchia approximation module."""

import numpy as np

from geomix.vecchia import find_vecchia_neighbours


def test_first_location_has_no_neighbours():
    coords = np.random.default_rng(1).standard_normal((10, 2))
    neighbours = find_vecchia_neighbours(coords, m=3)
    assert np.all(neighbours[0] == -1)


def test_neighbour_count_bounded_by_m():
    coords = np.random.default_rng(2).standard_normal((20, 3))
    m = 5
    neighbours = find_vecchia_neighbours(coords, m=m)
    assert neighbours.shape == (20, m)


def test_neighbours_are_predecessors():
    coords = np.random.default_rng(3).standard_normal((15, 2))
    neighbours = find_vecchia_neighbours(coords, m=4)
    for i, row in enumerate(neighbours):
        for idx in row:
            if idx != -1:
                assert idx < i, f"Location {i} has non-predecessor neighbour {idx}"
