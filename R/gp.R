#' Squared-exponential (RBF) covariance function
#'
#' @param X1 Numeric matrix of shape (n, d).
#' @param X2 Numeric matrix of shape (m, d).
#' @param length_scale Positive numeric length-scale parameter.
#' @param variance Positive numeric marginal variance parameter.
#' @return Numeric matrix of shape (n, m).
#' @export
squared_exponential_kernel <- function(X1, X2,
                                       length_scale = 1.0,
                                       variance = 1.0) {
  # TODO: implement squared-exponential covariance
  stop("squared_exponential_kernel() is not yet implemented.")
}

#' Matérn 5/2 covariance function
#'
#' @param X1 Numeric matrix of shape (n, d).
#' @param X2 Numeric matrix of shape (m, d).
#' @param length_scale Positive numeric length-scale parameter.
#' @param variance Positive numeric marginal variance parameter.
#' @return Numeric matrix of shape (n, m).
#' @export
matern52_kernel <- function(X1, X2,
                            length_scale = 1.0,
                            variance = 1.0) {
  # TODO: implement Matérn 5/2 covariance
  stop("matern52_kernel() is not yet implemented.")
}
