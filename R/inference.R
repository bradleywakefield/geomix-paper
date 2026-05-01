#' Run MCMC inference for a GeoMix model
#'
#' @param model An initialised (but not yet fitted) GeoMix model object.
#' @param n_samples Integer number of posterior samples to draw after warm-up.
#' @param n_warmup Integer number of warm-up / burn-in samples (discarded).
#' @param seed Integer random seed for reproducibility.
#' @return A posterior samples object (TODO: define format, e.g., \code{posterior::draws}).
#' @export
run_mcmc <- function(model, n_samples = 1000L, n_warmup = 500L, seed = 0L) {
  # TODO: implement MCMC inference (e.g., via Stan / rstan / cmdstanr)
  stop("run_mcmc() is not yet implemented.")
}
