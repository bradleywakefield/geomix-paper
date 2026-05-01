#' Find Vecchia nearest-neighbour conditioning sets
#'
#' For each location \code{i} (in the order they appear in \code{coords}),
#' returns the indices of up to \code{m} preceding locations that are nearest
#' in Euclidean distance. The first location has no predecessors; its row is
#' filled with \code{NA}.
#'
#' @param coords Numeric matrix of coordinates, shape (n, d).
#' @param m Integer maximum number of neighbours per location.
#' @return Integer matrix of shape (n, m) with \code{NA} for missing neighbours.
#' @export
find_vecchia_neighbours <- function(coords, m = 10L) {
  # TODO: implement Vecchia nearest-neighbour search
  stop("find_vecchia_neighbours() is not yet implemented.")
}
