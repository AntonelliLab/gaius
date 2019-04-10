# TODO: add additional arguments for msaR

#' @name viz
#' @title Visualise an alignment or supermatrix
#' @description Plot a Shiny multiple alignment object.
#' @param x alignment object
#' @return NULL
NULL

#' @export
viz <- function(x) {
  UseMethod('viz', x)
}

#' @rdname viz
#' @export
viz.alignment <- function(x) {
  nms <- rownames(x)
  x <- lapply(X = seq_len(nrow(x)), FUN = function(y) x[y, ])
  msar_wrapper(sequences = x, nms = nms)
}

#' @rdname viz
#' @export
viz.supermatrix <- function(x) {
  nms <- names(x)
  x <- lapply(X = x, FUN = function(y) strsplit(x = y, split = "")[[1]])
  msar_wrapper(sequences = x, nms = nms)
}

#' @name msar_wrapper
#' @title Visualise sequences
#' @description Plot a Shiny multiple alignment object.
#' @param sequences list of biological sequences as separate characters
#' @param nms names of sequences in list
#' @return NULL
msar_wrapper <- function(sequences, nms) {
  if (!requireNamespace("msaR", quietly = TRUE)) {
    msg <- paste0("Package ", obj("msaR"), " required. Run ",
                  elem("install.packages(\"msaR\")"), ' or similar.')
    stop(msg, call. = FALSE)
  }
  # write out, then read in
  tmpfl <- tempfile(fileext = '.fasta')
  seqinr::write.fasta(sequences = sequences, names = nms, file.out = tmpfl)
  print(msaR::msaR(msa = tmpfl))
  on.exit(file.remove(tmpfl))
  invisible(NULL)
}


