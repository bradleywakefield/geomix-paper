"""Class-specific Gaussian process models for geotechnical properties."""

import numpy as np


def squared_exponential_kernel(
    X1: np.ndarray,
    X2: np.ndarray,
    length_scale: float = 1.0,
    variance: float = 1.0,
) -> np.ndarray:
    """Compute the squared-exponential (RBF) covariance matrix.

    Parameters
    ----------
    X1 : np.ndarray of shape (n, d)
    X2 : np.ndarray of shape (m, d)
    length_scale : float
    variance : float

    Returns
    -------
    np.ndarray of shape (n, m)
    """
    diff = X1[:, None, :] - X2[None, :, :]
    sq_dist = np.sum(diff**2, axis=-1)
    return variance * np.exp(-0.5 * sq_dist / length_scale**2)


def matern52_kernel(
    X1: np.ndarray,
    X2: np.ndarray,
    length_scale: float = 1.0,
    variance: float = 1.0,
) -> np.ndarray:
    """Compute the Matérn 5/2 covariance matrix.

    Parameters
    ----------
    X1 : np.ndarray of shape (n, d)
    X2 : np.ndarray of shape (m, d)
    length_scale : float
    variance : float

    Returns
    -------
    np.ndarray of shape (n, m)
    """
    diff = X1[:, None, :] - X2[None, :, :]
    r = np.sqrt(np.sum(diff**2, axis=-1))
    sqrt5 = np.sqrt(5.0) * r / length_scale
    return variance * (1.0 + sqrt5 + sqrt5**2 / 3.0) * np.exp(-sqrt5)
