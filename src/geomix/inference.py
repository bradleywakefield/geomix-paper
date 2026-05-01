"""Inference routines for the GeoMix model."""


def run_mcmc(model, n_samples: int = 1000, n_warmup: int = 500, seed: int = 0):
    """Run MCMC inference for a fitted GeoMix model.

    Parameters
    ----------
    model : GeoMix
        An initialised (but not yet fitted) GeoMix model instance.
    n_samples : int
        Number of posterior samples to draw after warm-up.
    n_warmup : int
        Number of warm-up / burn-in samples (discarded).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict
        ArviZ ``InferenceData`` object (or equivalent) containing posterior
        samples.
    """
    raise NotImplementedError("MCMC inference is not yet implemented.")
