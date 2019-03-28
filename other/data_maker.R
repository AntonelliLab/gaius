# make test and example data ----
library(gaius)

# alignment ----
make_alignment <- function() {
  bps <- c('-', 'a', 't', 'c', 'g')
  sq <- rep(rep(x = bps, each = 5), 10)
  nrands <- floor(length(sq) * .1)
  res <- matrix(data = NA, nrow = 20, ncol = length(sq))
  for (i in seq_len(nrow(res))) {
    randomise <- sample.int(n = length(sq), size = nrands)
    res[i, ] <- sq
    res[i, randomise] <- sample(bps, size = nrands, replace = TRUE)
  }
  rownames(res) <- paste0('t', seq_len(nrow(res)))
  class(res) <- 'alignment'
  res
}
alignment <- make_alignment()
save(object = alignment, file = file.path('data', 'alignment.rda'),
     compress = 'xz')

# alignment_list ----
make_alignment_list <- function() {
  res <- list()
  for (i in 1:5) {
    res[[i]] <- make_alignment()
  }
  names(res) <- paste0('g', seq_len(length(res)))
  class(res) <- 'alignment_list'
  res
}
alignment_list <- make_alignment_list()
save(object = alignment_list, file = file.path('data', 'alignment_list.rda'),
     compress = 'xz')

# supermatrix ----
supermatrix <- gaius:::supermatrix_get(make_alignment_list())
save(object = supermatrix, file = file.path('data', 'supermatrix.rda'),
     compress = 'xz')

# supermatrices ---
groups <- list('n1' = paste0('t', 1:4),
               'n2' = paste0('t', 5:8),
               'unmatched' = paste0('t', 9:10))
supermatrices <- supermatrices_get(groups = groups,
                                   alignment_list = alignment_list,
                                   column_cutoff = 0.1, tip_cutoff = 0.1,
                                   min_ntips = 2, min_ngenes = 2,
                                   min_nbps = 200)
save(object = supermatrices, file = file.path('data', 'supermatrices.rda'),
     compress = 'xz')
