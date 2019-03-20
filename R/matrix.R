#' @name alignments_read
#' @title Read in sequence alignments
#' @description Return alignment object from reading in sequences from
#' a FASTA file.
#' @param flpth File path to .fasta
#' @return alignment
#' @export
alignments_read <- function(flpth) {
  all_data <- readLines(flpth)
  alignments <- list()
  for (i in seq_along(all_data)) {
    bit <- all_data[[i]]
    if (grepl(pattern = '^>', x = bit)) {
      nm <- sub(pattern = '^>', '', x = bit)
      alignments[[nm]] <- NULL
    } else {
      bit <- strsplit(x = bit, split = '')[[1]]
      alignments[[nm]] <- c(alignments[[nm]], bit)
    }
  }
  lngths <- vapply(X = alignments, FUN = length, FUN.VALUE = integer(1))
  if (all(lngths != lngths[[1]])) {
    stop('Not an alignment: sequences have different lengths.')
  }
  res <- matrix(unlist(alignments), ncol = length(alignments[[1]]),
                nrow = length(alignments),
                byrow = TRUE)
  rownames(res) <- names(alignments)
  class(res) <- 'alignment'
  res
}

#' @export
print.alignment <- function(x) {
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
    sq <- paste0(part[i, seq_len(nseq)], collapse = '')
    cat_line('... ', elem(tip_nm), spacer, ' ', seq(sq), ' ...')
  }
  if (nrow(x) > ntips) {
    cat_line('... with ', stat(nrow(x) - ntips), ' more tip(s)')
  }
}

#' @name gapmatrix_get
#' @title Generate a gaps matrix
#' @description Return a matrix of TRUE and FALSE corresponding to '-' in a
#' alignment.
#' @return gapmatrix
gapmatrix_get <- function(alignments) {
  res <- alignments == '-'
  rownames(res) <- rownames(alignments)
  class(res) <- 'gapmatrix'
  res
}

#' @name sequences_select
#' @title Select sequences from alignments
#' @description Return a alignment of just sequences corresponding to
#' nms.
#' @return alignment
sequences_select <- function(alignments_list, nms) {
  alignments_get <- function(i) {
    alignments <- alignments_list[[i]]
    nbp <- ncol(alignments)
    res <- matrix(data = '-', nrow = length(nms), ncol = nbp)
    rownames(res) <- nms
    present <- rownames(alignments)[rownames(alignments) %in% nms]
    res[present, ] <- alignments[present, ]
    class(res) <- 'alignment'
    res
  }
  res <- lapply(seq_along(alignments_list), alignments_get)
  names(res) <- names(alignments_list)
  res
}

#' @name sequences_filter
#' @title Filter sequecnes from alignments
#' @description Return a filtered alignment. Columns will be dropped
#' that have fewer proportion of non-gaps then cutoff and genes will be dropped
#' that are less than min_nbps long.
#' @return alignment
sequences_filter <- function(alignments_list, cutoff = 0.9, min_nbps = 200) {
  calc <- function(alignments) {
    gapmatrix <- gapmatrix_get(alignments = alignments)
    pclmn <- 1 - (colSums(gapmatrix)/nrow(gapmatrix))
    keep_clmns <- pclmn >= cutoff
    alignments <- alignments[ ,keep_clmns]
    class(alignments) <- 'alignment'
    alignments
  }
  alignments_list <- lapply(X = alignments_list, FUN = calc)
  nbps <- vapply(X = alignments_list, FUN = ncol, FUN.VALUE = integer(1))
  alignments_list[nbps >= min_nbps]
}

#' @name supermatrix_get
#' @title Convert list of alignments into a supermatrix
#' @description Return a supermatrix object from a list of alignments. The
#' supermatrix object is smaller.
#' @return supermatrix
supermatrix_get <- function(alignments_list) {
  stick_together <- function(i) {
    res <- unlist(lapply(X = alignments_list, FUN = function(x) x[i, ]))
    paste0(res, collapse = '')
  }
  nbps <- vapply(X = alignments_list, FUN = ncol, FUN.VALUE = integer(1))
  ntips <- nrow(alignments_list[[1]])
  res <- lapply(seq_len(ntips), stick_together)
  names(res) <- rownames(alignments_list[[1]])
  attr(res, 'genes') <- names(alignments_list)
  attr(res, 'nbps') <- nbps
  class(res) <- 'supermatrix'
  res
}

#' @export
print.supermatrix <- function(x) {
  nmssng <- sum(vapply(X = gregexpr(pattern = '-', text = x), FUN = length,
                       FUN.VALUE = integer(1)))
  total_nbps <- sum(attr(x, 'nbps')) * length(x)
  pmssng <- signif(x = nmssng * 100/total_nbps, digits = 2)
  cat('supermatrix:\n', '... [', length(attr(x, 'genes')), '] genes\n',
      '... [', length(x), '] tips\n',
      '... [', sum(attr(x, 'nbps')), '] bps long\n',
      '... [', pmssng, '%] gaps\n', sep = '')
}

#' @name drop_tips
#' @title Drop tips from a supermatrix
#' @description Return a supermatrix object with tips dropped that have more
#' non-gaps than cutoff.
#' @return supermatrix
drop_tips <- function(supermatrix, cutoff = 0.5) {
  nmissing <- vapply(X = gregexpr(pattern = '-', text = supermatrix),
                     FUN = length, FUN.VALUE = integer(1))
  nbps <- sum(attr(supermatrix, 'nbps'))
  pmissing <- 1 - (nmissing/nbps)
  keep_tips <- pmissing >= cutoff
  res <- supermatrix[keep_tips]
  attr(res, 'genes') <- attr(supermatrix, 'genes')
  attr(res, 'nbps') <- attr(supermatrix, 'nbps')
  class(res) <- 'supermatrix'
  res
}

#' @name supermatrices_get
#' @title Create supermatrices from list of alignments and groups object
#' @description Return a supermatrices object from a list of alignments
#' and a groups object.
#' @param groups groups object, determines names of monophyletic groups
#' @param alignments list of alignments
#' @param column_cutoff Min. prop. of non-gaps in an alignment column
#' @param tip_cutoff Min. prop. of non-gaps in an alignment row
#' @param min_ngenes Min. no. of genes in a supermatrix
#' @param min_nbps Min. length (in base pairs) for a gene in a supermatrix
#' @return supermatrix
#' @export
supermatrices_get <- function(groups, alignments, column_cutoff = .5,
                              tip_cutoff = column_cutoff, min_ntips = 5,
                              min_ngenes = 5, min_nbps = 200) {
  res <- list()
  all_tips <- character(0)
  for (grp_id in names(groups)) {
    nms <- groups[[grp_id]]
    # select sequences from alignments
    alignments_list <- sequences_select(nms = nms,
                                         alignments_list = alignments)
    # filter selected sequences
    alignments_list <- sequences_filter(alignments_list = alignments_list,
                                         cutoff = column_cutoff,
                                         min_nbps = min_nbps)
    if (length(alignments_list) == 0) {
      next
    }
    # merge into supermatrix
    supermatrix <- supermatrix_get(alignments_list = alignments_list)
    # drop tips
    supermatrix <- drop_tips(supermatrix = supermatrix, cutoff = tip_cutoff)
    if (length(supermatrix) >= min_ntips &
        length(attr(supermatrix, 'genes')) >= min_ngenes) {
      res[[grp_id]] <- supermatrix
      # record tips
      all_tips <- c(all_tips, names(supermatrix))
    }
  }
  attr(res, 'tips') <- all_tips
  class(res) <- 'supermatrices'
  res
}

#' @export
print.supermatrices <- function(x) {
  cat('supermatrices object of [', length(x), ']', sep = '')
}
