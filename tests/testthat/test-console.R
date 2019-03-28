# LIBS
library(gaius)
library(testthat)

# Running ----
context('Testing \'console\'')
test_that('char() works', {
  expect_is(object = gaius:::char('2'), class = 'character')
})
test_that('elem() works', {
  expect_is(object = gaius:::elem('2'), class = 'character')
})
test_that('stat() works', {
  expect_is(object = gaius:::stat('2'), class = 'character')
})
test_that('func() works', {
  expect_is(object = gaius:::func('2'), class = 'character')
})
test_that('obj() works', {
  expect_is(object = gaius:::obj('2'), class = 'character')
})
test_that('sq() works', {
  expect_is(object = gaius:::sq('atcg'), class = 'character')
})
test_that('cat_line() works', {
  expect_null(gaius:::cat_line('2'))
})
