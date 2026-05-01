"""Tests for GP kernel functions."""

import numpy as np
import pytest

from geomix.gp import matern52_kernel, squared_exponential_kernel


@pytest.fixture
def simple_coords():
    rng = np.random.default_rng(0)
    return rng.standard_normal((5, 2))


def test_se_kernel_shape(simple_coords):
    K = squared_exponential_kernel(simple_coords, simple_coords)
    assert K.shape == (5, 5)


def test_se_kernel_diagonal(simple_coords):
    K = squared_exponential_kernel(simple_coords, simple_coords, variance=2.0)
    np.testing.assert_allclose(np.diag(K), 2.0)


def test_se_kernel_symmetric(simple_coords):
    K = squared_exponential_kernel(simple_coords, simple_coords)
    np.testing.assert_allclose(K, K.T)


def test_matern52_shape(simple_coords):
    K = matern52_kernel(simple_coords, simple_coords)
    assert K.shape == (5, 5)


def test_matern52_diagonal(simple_coords):
    K = matern52_kernel(simple_coords, simple_coords, variance=3.0)
    np.testing.assert_allclose(np.diag(K), 3.0)


def test_matern52_symmetric(simple_coords):
    K = matern52_kernel(simple_coords, simple_coords)
    np.testing.assert_allclose(K, K.T, atol=1e-12)
