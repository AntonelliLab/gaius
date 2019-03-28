# LIBS
library(gaius)
library(testthat)

# Data ----
data("alignment")
data("alignment_list")
data("supermatrix")

# Running ----
context('Testing \'print-methods\'')
test_that('print.alignment() works', {
  expect_null(print(alignment))
})
test_that('print.alignment_list() works', {
  expect_null(print(alignment_list))
})
test_that('print.supermatrix() works', {
  expect_null(print(supermatrix))
})
test_that('print.supermatrices() works', {
  expect_null(print(supermatrices))
})
