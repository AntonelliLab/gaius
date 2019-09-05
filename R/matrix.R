#' @name gapmatrix_get
#' @title Generate a gaps matrix
#' @description Return a matrix of TRUE and FALSE corresponding to '-' in a
#' alignment.
#' @param alignment alignment object
#' @return gapmatrix
gapmatrix_get <- function(alignment) {
  res <- alignment == '-'
  rownames(res) <- rownames(alignment)
  class(res) <- 'gapmatrix'
  res
}

pmissing <- function(x) {
  UseMethod('pmissing', x)
}

#' @rdname gapmatrix_get
#' @param x gapmatrix object
pmissing.gapmatrix <- function(x) {
  sum(x)/(nrow(x) * ncol(x))
}

#' @name sequences_select
#' @title Select sequences from alignments
#' @description Return a alignment of just sequences corresponding to
#' nms.
#' @param alignment_list alignment list object
#' @param nms names of tips to be selected
#' @return alignment_list
sequences_select <- function(alignment_list, nms) {
  alignment_get <- function(i) {
    alignment <- alignment_list[[i]]
    nbp <- ncol(alignment)
    res <- matrix(data = '-', nrow = length(nms), ncol = nbp)
    rownames(res) <- nms
    present <- rownames(alignment)[rownames(alignment) %in% nms]
    res[present, ] <- alignment[present, ]
    class(res) <- 'alignment'
    res
  }
  res <- lapply(seq_along(alignment_list), alignment_get)
  names(res) <- names(alignment_list)
  class(res) <- 'alignment_list'
  res
}

#' @name sequences_filter
#' @title Filter sequecnes from alignments
#' @description Return a filtered alignment. Columns will be dropped
#' that have fewer proportion of non-gaps then cutoff and genes will be dropped
#' that are less than min_nbps long.
#' @param alignment_list alignment list object
#' @param cutoff Proportion of non-gaps
#' @param min_nbps Minimum sequence length
#' @return alignment_list
sequences_filter <- function(alignment_list, cutoff = 0.9, min_nbps = 200) {
  calc <- function(alignment) {
    gapmatrix <- gapmatrix_get(alignment = alignment)
    pclmn <- 1 - (colSums(gapmatrix)/nrow(gapmatrix))
    keep_clmns <- pclmn >= cutoff
    alignment <- alignment[ ,keep_clmns]
    class(alignment) <- 'alignment'
    alignment
  }
  alignment_list <- lapply(X = alignment_list, FUN = calc)
  nbps <- vapply(X = alignment_list, FUN = ncol, FUN.VALUE = integer(1))
  res <- alignment_list[nbps >= min_nbps]
  class(res) <- 'alignment_list'
  res
}

#' @name supermatrix_get
#' @title Convert list of alignments into a supermatrix
#' @description Return a supermatrix object from a list of alignments. The
#' supermatrix object is smaller.
#' @param alignment_list alignment list object
#' @return supermatrix
supermatrix_get <- function(alignment_list) {
  stick_together <- function(i) {
    res <- unlist(lapply(X = alignment_list, FUN = function(x) x[i, ]))
    paste0(res, collapse = '')
  }
  nbps <- vapply(X = alignment_list, FUN = ncol, FUN.VALUE = integer(1))
  ntips <- nrow(alignment_list[[1]])
  res <- lapply(seq_len(ntips), stick_together)
  names(res) <- rownames(alignment_list[[1]])
  gene_nms <- names(alignment_list)
  if (is.null(gene_nms)) {
    gene_nms <- as.character(seq_along(alignment_list))
  }
  attr(res, 'genes') <- gene_nms
  attr(res, 'nbps') <- nbps
  class(res) <- 'supermatrix'
  res
}

#' @name drop_tips
#' @title Drop tips from a supermatrix
#' @description Return a supermatrix object with tips dropped that have more
#' non-gaps than cutoff.
#' @param supermatrix supermatrix object
#' @param cutoff Proportion of non-gap cutoff
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
#' @param alignment_list list of alignments
#' @param column_cutoff Min. prop. of non-gaps in an alignment column
#' @param tip_cutoff Min. prop. of non-gaps in an alignment row
#' @param min_ntips Min. no. of tips in a supermatrix
#' @param min_ngenes Min. no. of genes in a supermatrix
#' @param min_nbps Min. length (in base pairs) for a gene in a supermatrix
#' @return supermatrices
#' @export
supermatrices_get <- function(alignment_list, groups = NULL, column_cutoff = .5,
                              tip_cutoff = column_cutoff, min_ntips = 5,
                              min_ngenes = 5, min_nbps = 200) {
  res <- list()
  all_tips <- character(0)
  if (is.null(groups)) {
    nms <- unique(unname(unlist(lapply(X = alignment_list, FUN = rownames))))
    groups <- list('all' = nms)
    class(groups) <- 'groups'
  }
  for (grp_id in names(groups)) {
    nms <- groups[[grp_id]]
    # select sequences from alignments
    selected <- sequences_select(nms = nms, alignment_list = alignment_list)
    # filter selected sequences
    filtered <- sequences_filter(alignment_list = selected,
                                 cutoff = column_cutoff, min_nbps = min_nbps)
    if (length(filtered) == 0) {
      next
    }
    # merge into supermatrix
    supermatrix <- supermatrix_get(alignment_list = filtered)
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

# TODO: update docs
#' @rdname supermatrices_get
#' @export
sift <- function(x, ...) {
  UseMethod('sift', x)
}

#' @rdname supermatrices_get
#' @param x supermatrices
#' @param keep Names of matrices to keep, character vector.
#' @param drop Names of matrices to drop, character vector
#' @return supermatrices
#' @export
sift.supermatrices <- function(x, keep = NULL, drop = NULL) {
  check <- function(nms) {
    pull <- nms %in% names(x)
    if (!all(pull)) {
      nms[!pull]
      msg <- paste0('Names not found:\n',
                    paste0(char(nms[!pull]), collapse = ', '))
      stop(msg)
    }
    names(x) %in% nms
  }
  if (!is.null(keep)) {
    x <- x[check(nms = keep)]
  }
  if (!is.null(drop)) {
    x <- x[!check(nms = drop)]
  }
  class(x) <- 'supermatrices'
  x
}