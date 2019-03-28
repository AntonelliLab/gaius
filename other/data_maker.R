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
# TODO