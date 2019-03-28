#' @name nms_from_fasta
#' @title Read tip labels from sequence file
#' @description Return .
#' @param flpth File path to .fasta file
#' @return character vector
nms_from_fasta <- function(flpth) {
  lines <- readLines(con = flpth)
  nms <- lines[grepl(pattern = '^>', x = lines)]
  sub(pattern = '^>', replacement = '', nms)
}


alignment_info_get <- function(flpth) {
  nms <- nms_from_fasta(flpth = flpth)
  res <- list(nms = nms, ntips = length(nms), flpth = flpth)
  class(res) <- 'alignment_info'
  res
}

nms_from_alignments <- function(flpths) {
  calc <- function(flpth) {
    alignment <- alignment_info_get(flpth = flpth)
    alignment[['nms']]
  }
  nms <- lapply(X = flpths, FUN = calc)
  unique(unlist(nms))
}

name_match <- function(alignment_names, tree_names, max_dist = .1,
                       alignment_patterns = alignment_names,
                       tree_patterns = tree_names) {
  # Checks
  if (any(duplicated(alignment_names)) | any(duplicated(tree_names))) {
    stop('Duplicated names.')
  }
  if (length(alignment_names) != length(alignment_patterns) |
      length(tree_names) != length(tree_patterns)) {
    stop('`*_names` and `*_patterns` must be same lengths.')
  }
  calc <- function(i) {
    # Calc prop. Levenshtein distance, select tree name with lowest distance
    alignment_pattern <- alignment_patterns[[i]]
    dists <- utils::adist(x = alignment_pattern, y = tree_patterns,
                          partial = TRUE)[1, ]
    pdists <- dists/nchar(alignment_pattern)
    pssbls <- which(pdists < max_dist)
    npssbls <- length(pssbls)
    if (npssbls > 1) {
      res <- tree_names[pssbls[which.min(pdists[pssbls])]]
    } else if (npssbls == 1) {
      res <- tree_names[pssbls]
    } else {
      res <- ''
    }
    res
  }
  matched_treenms <- vapply(X = seq_along(alignment_names), FUN = calc,
                            FUN.VALUE = character(1))
  unmatched <- alignment_names[matched_treenms == '']
  matched <- alignment_names[matched_treenms != '']
  matched_treenms <- matched_treenms[matched_treenms != '']
  res <- list('alignment' = matched, 'tree' = matched_treenms,
              'unmatched' = unmatched)
  class(res) <- 'matched_names'
  res
}

groups_get <- function(matched_names, tree, max_size, min_size) {
  group_get <- function(tree, groups) {
    ptids <- treeman::getNdSlt(tree = tree, slt_nm = 'ptid', id = tree@root)
    ptids <- ptids[!ptids %in% tree@tips]
    for (ptid in ptids) {
      subtree <- treeman::getSubtree(tree = tree, id = ptid)
      pssbls <- unname(nms[nms %in% subtree@tips])
      size <- length(pssbls)
      if (size <= max_size & size >= min_size) {
        groups[[ptid]] <- pssbls
      } else {
        if (subtree@ntips > min_size) {
          groups <- group_get(tree = subtree, groups = groups)
        }
      }
    }
    groups
  }
  nms <- matched_names$tree
  alignment_names <- matched_names$alignment
  groups <- group_get(tree = tree, groups = list())
  res <- lapply(X = groups, function(x) alignment_names[nms %in% x])
  unmatched <- c(matched_names[['unmatched']],
                 alignment_names[!alignment_names %in% unlist(res)])
  res[['unmatched']] <- unmatched
  class(res) <- "groups"
  res
}
