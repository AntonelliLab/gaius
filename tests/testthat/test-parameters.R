# LIBS
library(gaius)
library(testthat)

# Running ----
context('Testing \'parameters\'')
test_that('default_pset() works', {
  expect_true(default_pset())
})
test_that('default_pget() works', {
  expect_equal(default_pget(), gaius:::default_parameters)
})
test_that('pset() and pget() work', {
  # handles both numeric and integer types
  expect_true(pset(val = 200, parameter = 'min_nbps'))
  expect_true(pset(val = 200L, parameter = 'min_nbps'))
  # errors if character
  expect_error(pset(val = "200", parameter = 'min_nbps'))
  expect_equal(pget(parameter = 'min_nbps'), 200)
})
test_that("ptype_match() works", {
  # both numeric
  expect_true(gaius:::ptype_match(1, 10))
  expect_true(gaius:::ptype_match(1, 10L))
  # different types
  expect_false(gaius:::ptype_match(1, "10"))
  expect_false(gaius:::ptype_match(1, FALSE))
})
