#' @export
print.alignment <- function(x, ...) {
  max_ntips <- 10
  max_nclmns <- options()$width
  ntips <- ifelse(nrow(x) >= max_ntips, max_ntips, length(x))
  nclmns <- ifelse(ncol(x) >= max_nclmns, max_nclmns, ncol(x))
  part <- x[seq_len(ntips), seq_len(nclmns)]
  max_tpnm <- max(nchar(rownames(part)) + 1)
  nseq <- max_nclmns - max_tpnm - 8
  cat_line('An ', obj('alignment'),  ': ', stat(nrow(x)), ' tips x ',
           stat(ncol(x)), ' bps')
  for (i in seq_len(ntips)) {
    tip_nm <- paste0('$', rownames(part)[i])
    nspaces <- max_tpnm - nchar(tip_nm)
    spacer <- paste0(rep(' ', nspaces, replace = TRUE), collapse = '')
    s <- paste0(part[i, seq_len(nseq)], collapse = '')
    cat_line('... ', elem(tip_nm), spacer, ' ', sq(s), ' ...')
  }
  if (nrow(x) > ntips) {
    cat_line('... with ', stat(nrow(x) - ntips), ' more tip(s)')
  }
}

#' @export
print.alignment_list <- function(x, ...) {
  cat_line('List of ', stat(length(x)), ' alignments')
}
  
#' @export
print.supermatrix <- function(x, ...) {
  # stats
  nbps <- sum(attr(x, 'nbps'))
  ntips <- length(x)
  # seq print limits
  max_ntips <- 10
  max_nclmns <- options()$width
  nlines <- ifelse(ntips >= max_ntips, max_ntips, ntips)
  part <- x[seq_len(nlines)]
  max_tpnm <- max(nchar(names(part)) + 1)
  nseq <- max_nclmns - max_tpnm - 8
  # gap calc
  nmssng <- sum(vapply(X = gregexpr(pattern = '-', text = x), FUN = length,
                       FUN.VALUE = integer(1)))
  total_nbps <- nbps * length(x)
  pmssng <- signif(x = nmssng * 100/total_nbps, digits = 2)
  # first line
  cat_line('A ', obj('supermatrix'),  ': ', stat(ntips), ' tips x ',
           stat(nbps), ' bps x ', stat(length(attr(x, 'genes'))), ' genes, ',
           stat(pmssng), '% gaps')
  # seq
  for (i in seq_len(nlines)) {
    tip_nm <- paste0('$', names(part)[i])
    nspaces <- max_tpnm - nchar(tip_nm)
    spacer <- paste0(rep(' ', nspaces, replace = TRUE), collapse = '')
    s <- substr(part[[i]], start = 1, stop = nseq)
    cat_line('... ', elem(tip_nm), spacer, ' ', sq(s), ' ...')
  }
  if (ntips > nlines) {
    cat_line('... with ', stat(ntips - nlines), ' more tip(s)')
  }
}

#' @export
print.supermatrices <- function(x, ...) {
  cat('supermatrices object of [', length(x), ']', sep = '')
}