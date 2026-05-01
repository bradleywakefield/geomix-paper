"""Top-level GeoMix model combining all components."""


class GeoMix:
    """Bayesian hierarchical model for stratigraphic and geotechnical inference.

    Parameters
    ----------
    n_classes : int
        Number of stratigraphic classes (strata).
    potts_beta : float
        Spatial smoothing parameter for the Potts MRF.
    gp_kernel : str
        Covariance kernel for the class-specific Gaussian processes.
        One of ``"matern32"``, ``"matern52"``, ``"rbf"``.
    vecchia_m : int
        Number of nearest neighbours used in the grouped Vecchia approximation.
    """

    def __init__(
        self,
        n_classes: int = 3,
        potts_beta: float = 1.0,
        gp_kernel: str = "matern52",
        vecchia_m: int = 10,
    ) -> None:
        self.n_classes = n_classes
        self.potts_beta = potts_beta
        self.gp_kernel = gp_kernel
        self.vecchia_m = vecchia_m

    def fit(self, seismic_labels, geotechnical_data, coords):
        """Fit the model using MCMC.

        Parameters
        ----------
        seismic_labels : array-like of shape (n_locations,)
            Seismic-derived stratigraphic labels.
        geotechnical_data : array-like of shape (n_measurements, n_properties)
            Observed geotechnical measurements.
        coords : array-like of shape (n_locations, 3)
            Spatial coordinates (easting, northing, depth).
        """
        raise NotImplementedError("Model fitting is not yet implemented.")

    def predict(self, coords_new):
        """Predict stratigraphy and geotechnical properties at new locations.

        Parameters
        ----------
        coords_new : array-like of shape (n_new, 3)
            Spatial coordinates for prediction.

        Returns
        -------
        dict
            Dictionary with keys ``"stratigraphy"`` and ``"geotechnical"``,
            each containing posterior predictive samples.
        """
        raise NotImplementedError("Prediction is not yet implemented.")
