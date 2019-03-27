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

# Running ----
context('Testing \'matrix\'')
test_that('alignments_read() works', {
  alignment_list <- alignment_read(flpths = alignments_files)
  expect_true(inherits(x = alignment_list, what = 'alignment_list'))
  expect_true(inherits(x = alignment_list[[1]], what = 'alignment'))
})
test_that('gapmatrix_get() works', {
  alignment_list <- alignment_read(flpths = alignment_dir)
  alignment <- alignment_list[[1]]
  gpmtrx <- gaius:::gapmatrix_get(alignment = alignment)
  expect_equal(sum(alignment[1, 1:10] == '-'), sum(gpmtrx[1, 1:10]))
  expect_true(inherits(x = gpmtrx, what = 'gapmatrix'))
  expect_true(is.numeric(gaius:::pmissing.gapmatrix(gpmtrx)))
})
test_that('sequences_select() works', {
  alignment_list <- alignment_read(flpths = alignments_files)
  # random selection of names
  nms <- sample(x = rownames(alignment_list[[1]]), size = 10)
  selected <- gaius:::sequences_select(alignment_list = alignment_list,
                                       nms = nms)
  expect_true(inherits(x = selected, what = 'alignment_list'))
  expect_true(inherits(x = selected[[1]], what = 'alignment'))
  # ensure 10 nms
  expect_true(all(vapply(X = selected, FUN = nrow,
                         FUN.VALUE = integer(1)) == 10))
  expect_true(all(vapply(X = selected, FUN = function(x) {
    all(nms %in% rownames(x))
    }, logical(1))))
})
test_that('sequences_filter() works', {
  alignment_list <- alignment_read(flpths = alignments_files)
  nms <- sample(x = rownames(alignment_list[[1]]), size = 10)
  selected <- gaius:::sequences_select(alignment_list = alignment_list,
                                       nms = nms)
  # uber strict
  filtered <- gaius:::sequences_filter(alignment_list = selected, cutoff = 1.0,
                                       min_nbps = 10000)
  expect_true(length(filtered) == 0)
  # lenient
  filtered <- gaius:::sequences_filter(alignment_list = selected, cutoff = 0.5,
                                       min_nbps = 200)
  # cutoff is above 0.5
  ppresent <- vapply(X = filtered, FUN = function(x) {
    1 - gaius:::pmissing.gapmatrix(gaius:::gapmatrix_get(alignment = x))
  }, FUN.VALUE = numeric(1))
  expect_true(all(ppresent >= 0.5))
  # min bp is above 200
  nbps <- vapply(X = filtered, FUN = ncol, FUN.VALUE = integer(1))
  expect_true(all(nbps >= 200))
})
test_that('supermatrix_get() works', {
  alignment_list <- alignment_read(flpths = alignments_files)
  nms <- sample(x = rownames(alignment_list[[1]]), size = 10)
  selected <- gaius:::sequences_select(alignment_list = alignment_list,
                                       nms = nms)
  filtered <- gaius:::sequences_filter(alignment_list = selected, cutoff = 0.5,
                                       min_nbps = 200)
  supermatrix <- gaius:::supermatrix_get(alignment_list = filtered)
  expect_true(inherits(x = supermatrix, what = 'supermatrix'))
  expect_true(length(attr(supermatrix, 'genes')) ==
                length(attr(supermatrix, 'nbps')))
  expect_true(all(names(supermatrix) %in% nms))
})
test_that('drop_tips() works', {
  # alignment_list <- alignment_read(flpths = alignments_files)
  # nms <- sample(x = rownames(alignment_list[[1]]), size = 10)
  # selected <- gaius:::sequences_select(alignment_list = alignment_list,
  #                                      nms = nms)
  # filtered <- gaius:::sequences_filter(alignment_list = selected,
  #                                      cutoff = 0.5, min_nbps = 200)
  # supermatrix <- gaius:::supermatrix_get(alignment_list = filtered)
  # saveRDS(object = supermatrix, file = file.path(data_dir, 'supermatrix.rda'),
  #         compress = 'xz')
  supermatrix <- readRDS(file = file.path(data_dir, 'supermatrix.rda'))
  supermatrix_shrunk <- gaius:::drop_tips(supermatrix = supermatrix,
                                          cutoff = 0.99)
  expect_true(length(supermatrix) > length(supermatrix_shrunk))
  # TODO: create pgaps() and ngaps()
})
test_that('supermatrices_get() works', {
  alignment_list <- alignment_read(flpths = alignments_files)
  # create fake groups
  all_nms <- sample(unique(unlist(lapply(X = alignment_list, FUN = rownames))))
  ngroups <- 3
  size <- floor(length(all_nms)/(ngroups + 1))
  groups <- vector(mode = 'list', length = ngroups)
  for (i in seq_len(ngroups)) {
    groups[[i]] <- sample(x = all_nms, size = size)
    all_nms <- all_nms[!all_nms %in% unlist(groups)]
  }
  names(groups) <- paste0('g', seq_along(groups))
  groups[['unmatched']] <- all_nms
  # get supermatrices
  supermatrices <- supermatrices_get(groups = groups,
                                     alignment_list = alignment_list,
                                     column_cutoff = 0.1)
  expect_true(inherits(x = supermatrices, what = 'supermatrices'))
})
