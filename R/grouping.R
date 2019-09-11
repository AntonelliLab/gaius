# Private ----
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

group_search <- function(matching_tree_names, tree_file) {
  group_get <- function(tree, groups) {
    ptids <- treeman::getNdSlt(tree = tree, slt_nm = 'ptid', id = tree@root)
    ptids <- ptids[!ptids %in% tree@tips]
    for (ptid in ptids) {
      subtree <- treeman::getSubtree(tree = tree, id = ptid)
      pssbls <- unname(matching_tree_names[matching_tree_names %in%
                                             subtree@tips])
      size <- length(pssbls)
      if (size <= pget('max_ntips') & size >= pget('min_ntips')) {
        groups[[ptid]] <- pssbls
      } else {
        if (subtree@ntips > pget('min_ntips')) {
          groups <- group_get(tree = subtree, groups = groups)
        }
      }
    }
    groups
  }
  tree <- treeman::readTree(file = tree_file)
  group_get(tree = tree, groups = list())
}

# Public ----
#' @export
names_from_alignments <- function(flpths) {
  calc <- function(flpth) {
    alignment <- alignment_info_get(flpth = flpth)
    alignment[['nms']]
  }
  nms <- lapply(X = flpths, FUN = calc)
  unique(unlist(nms))
}

#' @export
names_from_tree <- function(flpth) {
  # TODO: avoid reading in tree?
  tree <- treeman::readTree(file = flpth)
  gsub(pattern = '\\s', replacement = "_", x = tree@tips)
}

#' @export
name_match <- function(alignment_names, tree_names,
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
    pssbls <- which(pdists < pget('max_name_dist'))
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

#' @export
groups_get <- function(tree_file, alignment_files = NULL,
                       matched_names = NULL) {
  if (is.null(matched_names) & !is.null(alignment_files)) {
    # keep only the first part of the word for the patterns
    tree_names <- names_from_tree(flpth = tree_file)
    tree_patterns <- sub(pattern = '_.*$', replacement = '', x = tree_names)
    alignment_names <- names_from_alignments(flpths = alignment_files)
    alignment_patterns <- sub(pattern = '\\s.*$', replacement = '',
                              x = alignment_names)
    matched_names <- name_match(alignment_names = alignment_names,
                                tree_names = tree_names,
                                alignment_patterns = alignment_patterns,
                                tree_patterns = tree_patterns)
  }
  nms <- matched_names[['tree']]
  alignment_names <- matched_names[['alignment']]
  if (length(nms) == 0 & length(alignment_names) == 0) {
    warning('No names matched.')
  }
  groups <- group_search(matching_tree_names = nms, tree_file = tree_file)
  res <- lapply(X = groups, function(x) alignment_names[nms %in% x])
  unmatched <- c(matched_names[['unmatched']],
                 alignment_names[!alignment_names %in% unlist(res)])
  res[['unmatched']] <- unmatched
  class(res) <- "groups"
  res
}
