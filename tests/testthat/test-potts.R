library(testthat)
library(geomix)

test_that("build_grid_adjacency raises not-implemented error", {
  # TODO: replace with real tests once build_grid_adjacency() is implemented
  expect_error(
    build_grid_adjacency(3L, 4L),
    "not yet implemented"
  )
})

test_that("potts_energy raises not-implemented error", {
  # TODO: replace with real tests once potts_energy() is implemented
  expect_error(
    potts_energy(list(), integer(0)),
    "not yet implemented"
  )
})
