#' @export
print.alignment <- function(x, ...) {
  max_ntips <- 10
  max_nclmns <- options()$width
  ntips <- ifelse(nrow(x) >= max_ntips, max_ntips, nrow(x))
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
  cat_line(obj('alignment_list'), ' object, with ', stat(length(x)), ' ',
           obj('alignment'), ' objects\n')
  if (length(x) > 0) {
    nm <- names(x)[[1]]
    if (!is.null(nm)) {
      cat_line(elem(paste0('[["', nm, '"]] ...')))
    } else {
      cat_line(elem('[[1]] ...'))
    }
    print(x[[1]])
    cat_line()
    nm <- names(x)[[length(x)]]
    if (!is.null(nm)) {
      cat_line(elem(paste0('... [["', nm, '"]]')))
    } else {
      cat_line(elem(paste0('... [[', length(x), ']]')))
    }
  }
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
  cat_line(obj('supermatices'), ' list, with ', stat(length(x)), ' ',
           obj('supermatrix'), ' objects\n')
  if (length(x) > 0) {
    nm <- names(x)[[1]]
    if (!is.null(nm)) {
      cat_line(elem(paste0('[["', nm, '"]] ...')))
    } else {
      cat_line(elem('[[1]] ...'))
    }
    print(x[[1]])
    cat_line()
    nm <- names(x)[[length(x)]]
    if (!is.null(nm)) {
      cat_line(elem(paste0('... [["', nm, '"]]')))
    } else {
      cat_line(elem(paste0('... [[', length(x), ']]')))
    }
  }
}

#' @export
print.matched_names <- function(x) {
  format_text <- function(y) {
    n <- ifelse(length(y) > 3, 3, length(y))
    paste0(paste0(y[1:n], collapse = ', '), ' ...')
  }
  cat('[', length(x$alignment), '] names matched:\n',
      '... from `$alignment` ', format_text(x$alignment),
      '\n... to `$tree` ', format_text(x$tree),
      '\n[', length(x$unmatched), '] unmatched.\n', sep = '')
}

#' @export
print.groups <- function(x) {
  nmatched <- vapply(X = x, FUN = length, integer(1))
  nunmatched <- nmatched[['unmatched']]
  nmatched <- sum(nmatched[names(nmatched) != 'unmatched'])
  cat('Tips group object of [', length(x), '] elements:\n', '... [', nmatched,
      '] tips in mono\n... [', nunmatched, '] tips in super',
      sep = '')
}

#' @export
print.alignment_info <- function(x) {
  cat('Alignment (', stat(x$ntips), ') tips.\n', sep = '')
}

