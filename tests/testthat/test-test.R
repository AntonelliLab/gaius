# LIBS
library(gaius)
library(testthat)

# Running ----
context('Testing \'test\'')
test_that('datadir_get() works', {
  expect_true(dir.exists(gaius:::datadir_get()))
})
