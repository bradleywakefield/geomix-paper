"""Tests for the Potts MRF module."""

import numpy as np
import pytest

from geomix.potts import build_grid_adjacency, potts_energy


def test_build_grid_adjacency_shape():
    adjacency = build_grid_adjacency(3, 4)
    assert len(adjacency) == 12


def test_build_grid_adjacency_corner_node():
    # Corner node (0, 0) should have exactly 2 neighbours
    adjacency = build_grid_adjacency(3, 3)
    assert len(adjacency[0]) == 2


def test_build_grid_adjacency_interior_node():
    # Interior node (1, 1) in a 3x3 grid → index 4 → 4 neighbours
    adjacency = build_grid_adjacency(3, 3)
    assert len(adjacency[4]) == 4


def test_potts_energy_all_same():
    # All labels identical → maximum agreement → most negative energy
    adjacency = build_grid_adjacency(2, 2)
    labels_same = np.zeros(4, dtype=int)
    labels_diff = np.array([0, 1, 0, 1])
    assert potts_energy(labels_same, adjacency, beta=1.0) < potts_energy(
        labels_diff, adjacency, beta=1.0
    )


def test_potts_energy_zero_beta():
    adjacency = build_grid_adjacency(2, 2)
    labels = np.array([0, 1, 2, 0])
    assert potts_energy(labels, adjacency, beta=0.0) == pytest.approx(0.0)
