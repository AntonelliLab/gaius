#' @name align
#' @title Generate alignments from sequences
#' @description Create alignments from sequence files using method specified.
#' @param sequence_files file paths to sequence files
#' @param method alignment method
#' @return filepaths to alignments, character vector
#' @export
align <- function(sequence_files, method = c('mafft', 'pasta')) {
  # TODO: Multi-processing? Should be at outsider level
  method <- match.arg(method)
  align_program <- switch(method, mafft = mafft, pasta = pasta)
  alignment_files <- tools::file_path_sans_ext(sequence_files)
  align_program(sequence_files = sequence_files,
                alignment_files = alignment_files)
}

mafft <- function(sequence_files, alignment_files) {
  outsider_install(repo = 'dombennett/om..mafft', service = 'github')
  foo <- outsider::module_import('mafft', repo = 'dombennett/om..mafft')
  alignment_files <- paste0(alignment_files, '_mafft_alignment.fasta')
  for (i in seq_along(sequence_files)) {
    foo(c('--auto', sequence_files[[i]], '>', alignment_files[[i]]))
  }
  alignment_files
}

pasta <- function(sequence_files, alignment_files) {
  # TODO
  # outsider_install(repo = 'dombennett/om..pasta', service = 'github')
  NULL
}
