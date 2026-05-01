library(testthat)
library(geomix)

test_that("find_vecchia_neighbours raises not-implemented error", {
  # TODO: replace with real tests once find_vecchia_neighbours() is implemented
  coords <- matrix(rnorm(20), nrow = 10, ncol = 2)
  expect_error(
    find_vecchia_neighbours(coords, m = 3L),
    "not yet implemented"
  )
})
