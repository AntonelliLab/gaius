# LIBS
library(gaius)
library(testthat)

# Vars ----
data_dir <- gaius:::datadir_get()
alignment_dir <- file.path(data_dir, 'alignment.fasta')
alignments_dir <- file.path(data_dir, 'alignments')
alignments_files <- file.path(alignments_dir, list.files(path = alignments_dir,
                                                         pattern = '.fasta'))
# sample for speed
alignments_files <- sample(alignments_files, 5)
data('alignment')
data('alignment_list')
data('supermatrix')
data('supermatrices')

# Running ----
context('Testing \'read-write\'')
test_that('format_guess() works', {
  tester <- function(expctd, exts) {
    res <- vapply(X = paste0('seqs.', exts), FUN = gaius:::format_guess,
                  FUN.VALUE = character(1))
    all(res == expctd)
  }
  expect_true(tester('fasta', c('fa', 'fasta')))
  expect_true(tester('clustal', c('clus', 'aln')))
  expect_true(tester('msf', 'msf'))
  expect_true(tester('mase', 'mase'))
  expect_true(tester('phylip', 'phy'))
  expect_error(gaius:::format_guess('seqs.txt'))
})
test_that('.alignment_read() works', {
  res <- gaius:::.alignment_read(flpth = alignments_files[[1]],
                                 format = 'fasta')
  expect_true(inherits(res, 'alignment'))
  res <- gaius:::.alignment_read(flpth = alignments_files[[1]],
                                 format = 'guess')
  expect_true(inherits(res, 'alignment'))
})
test_that('alignment_read() works', {
  res <- alignment_read(flpths = alignments_files, format = 'fasta')
  expect_true(inherits(res, 'alignment_list'))
  res <- alignment_read(flpths = alignments_files)
  expect_true(inherits(res, 'alignment_list'))
  expect_error(alignment_read(flpths = alignments_files, format = 'clustal'))
  expect_error(alignment_read(flpths = alignments_files, format = 'notaformat'))
})
test_that('.sequences_write() works', {
  tester <- function(x) {
    flpth <- tempfile(fileext = '.fasta')
    on.exit(if(res) file.remove(flpth))
    res <- gaius:::.sequences_write(flpth = flpth, x = x, width = 80)
    res
  }
  expect_true(tester(x = alignment))
  expect_true(tester(x = supermatrix))
  expect_error(tester(x = supermatrices))
})
test_that('sequences_write() works', {
  tester <- function(x) {
    flpth <- tempfile(fileext = '.fasta')
    on.exit(if (res) file.remove(flpth))
    res <- sequences_write(flpth = flpth, x = x, width = 80)
    res
  }
  tester2 <- function(x) {
    flpth <- file.path(tempdir(), 'tester2')
    dir.create(flpth)
    on.exit(if (res) unlink(x = flpth, recursive = TRUE, force = TRUE))
    res <- sequences_write(flpth = flpth, x = x, width = 80)
    res
  }
  expect_true(tester(x = alignment))
  expect_true(tester(x = supermatrix))
  expect_true(tester2(x = alignment_list))
  expect_true(tester2(x = supermatrices))
})
