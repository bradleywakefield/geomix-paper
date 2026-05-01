library(testthat)
library(geomix)

test_that("squared_exponential_kernel raises not-implemented error", {
  # TODO: replace with real tests once squared_exponential_kernel() is implemented
  expect_error(
    squared_exponential_kernel(matrix(0, 2, 2), matrix(0, 2, 2)),
    "not yet implemented"
  )
})

test_that("matern52_kernel raises not-implemented error", {
  # TODO: replace with real tests once matern52_kernel() is implemented
  expect_error(
    matern52_kernel(matrix(0, 2, 2), matrix(0, 2, 2)),
    "not yet implemented"
  )
})
