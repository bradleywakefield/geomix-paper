"""Plotting utilities for GeoMix results."""

import matplotlib.pyplot as plt
import numpy as np


def plot_stratigraphy_map(
    coords_2d: np.ndarray,
    labels: np.ndarray,
    n_classes: int,
    title: str = "Stratigraphic Map",
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot a 2-D stratigraphic label map.

    Parameters
    ----------
    coords_2d : np.ndarray of shape (n, 2)
        Easting and northing coordinates.
    labels : np.ndarray of shape (n,)
        Integer class labels.
    n_classes : int
        Total number of stratigraphic classes.
    title : str
        Axes title.
    ax : plt.Axes, optional
        Existing axes to plot on. A new figure is created if ``None``.

    Returns
    -------
    plt.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        coords_2d[:, 0],
        coords_2d[:, 1],
        c=labels,
        cmap="tab10",
        vmin=0,
        vmax=n_classes - 1,
        s=5,
    )
    plt.colorbar(sc, ax=ax, label="Stratigraphic class")
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_title(title)
    return ax
