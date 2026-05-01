#' Build a 2-D grid adjacency list
#'
#' Returns a list of length \code{nrow * ncol} where each element is an
#' integer vector of neighbour indices (4-connectivity).
#'
#' @param nrow Number of rows in the grid.
#' @param ncol Number of columns in the grid.
#' @return List of integer vectors.
#' @export
build_grid_adjacency <- function(nrow, ncol) {
  # TODO: implement 4-connected grid adjacency
  stop("build_grid_adjacency() is not yet implemented.")
}

#' Compute the Potts energy for a given label configuration
#'
#' @param adjacency Adjacency list as returned by \code{build_grid_adjacency()}.
#' @param labels Integer vector of class labels (length = number of nodes).
#' @param beta Inverse temperature (non-negative numeric).
#' @return Scalar energy value.
#' @export
potts_energy <- function(adjacency, labels, beta = 1.0) {
  # TODO: implement Potts energy
  stop("potts_energy() is not yet implemented.")
}
