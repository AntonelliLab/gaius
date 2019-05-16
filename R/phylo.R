#' @name phylo
#' @title Generate phylogenies from supermatrices
#' @description Create phylogenies from supermatrices using method specified.
#' @param supermatrices supermatrices object
#' @param method phylogeny generation method
#' @return filepaths to phylogenies, character vector
#' @export
phylo <- function(supermatrices, method = c('raxml', 'astral')) {
  method <- match.arg(method)
  phylo_program <- switch(method, raxml = raxml, astral = astral)
  phylo_program(supermatrices = supermatrices)
}

raxml <- function(supermatrices) {
  outsider_install(repo = 'dombennett/om..raxml', service = 'github')
  foo <- outsider::module_import('raxml', repo = 'dombennett/om..raxml')
  alignment_files <- paste0(alignment_files, '_mafft_alignment.fasta')
  for (i in seq_along(sequence_files)) {
    tmpfl <- tempfile(pattern = 'supermatrix', fileext = '.fasta')
    sequences_write(x = supermatrices, flpth = tmpfl)
    foo('--auto', sequence_files[[i]], '>', alignment_files[[i]])
  }
  alignment_files
}

astral <- function() {
  outsider_install(repo = 'dombennett/om..mafft', service = 'github')
  foo <- outsider::module_import('mafft', repo = 'dombennett/om..mafft')
  alignment_files <- paste0(alignment_files, '_mafft_alignment.fasta')
  for (i in seq_along(sequence_files)) {
    foo('--auto', sequence_files[[i]], '>', alignment_files[[i]])
  }
  alignment_files
}
