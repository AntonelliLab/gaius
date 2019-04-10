# Private ----
#' @name format_guess
#' @title Guess the sequence file format
#' @description Return sequence file format from a file path based on the
#' file's extension. If no format found, errors.
#' @param flpth File path to sequence file. Character.
#' @return character
format_guess <- function(flpth) {
  fmt <- NULL
  ext <- tools::file_ext(x = flpth)
  if (grepl(pattern = '(clus|aln)', x = ext, ignore.case = TRUE)) {
    fmt <- "clustal"
  } else if (grepl(pattern = 'msf', x = ext, ignore.case = TRUE)) {
    fmt <- "msf"
  } else if (grepl(pattern = 'mase', x = ext, ignore.case = TRUE)) {
    fmt <- "mase"
  } else if (grepl(pattern = 'phy', x = ext, ignore.case = TRUE)) {
    fmt <- "phylip"
  } else if (grepl(pattern = 'fa', x = ext, ignore.case = TRUE)) {
    fmt <- "fasta"
  } else {
    msg <- paste0('Hmmmm.... I don\'t know this file extension for ',
                  char(flpth))
    stop(msg, call. = FALSE)
  }
  fmt
}

#' @name .alignment_read
#' @title Read alignments
#' @description Return "alignment" from a filepath.
#' @details If \code{format} is 'guess', file format will be determined from
#' file extension.
#' @param flpth File path to output file. Character.
#' @param format Sequence file format. Character.
#' @seealso \code{\link[seqinr]{read.alignment}}
#' @return alignment
.alignment_read <- function(flpth, format) {
  if (format == 'guess') {
    fmt <- format_guess(flpth)
  } else {
    fmt <- format
  }
  temp <- seqinr::read.alignment(file = flpth, format = fmt)
  res <- t(vapply(X = temp$seq, FUN = function(x) {
    strsplit(x = x, split = '')[[1]]
  }, FUN.VALUE = character(nchar(temp$seq[[1]]))))
  rownames(res) <- temp[['nam']]
  class(res) <- 'alignment'
  res
}

#' @name .sequences_write
#' @title Write sequences
#' @description Return "alignment" or "supermatrix"  object to .FASTA file.
#' @param flpth File path to output file. Character.
#' @param x Sequence object, either "alignment" or "supermatrix".
#' @param width Number of sequence columns. Integer.
#' @return Logical
.sequences_write <- function(flpth, x, width) {
  sqs_from_alignment <- function(x) {
    sqs <- split(x, rep(seq_len(nrow(x)), each = ncol(x)))
    nms <- rownames(x)
    list('sqs' = sqs, 'nms' = nms)
  }
  sqs_from_supermatrix <- function(x) {
    sqs <- lapply(X = x, FUN = function(x) strsplit(x = x, split = '')[[1]])
    nms <- names(x)
    list('sqs' = sqs, 'nms' = nms)
  }
  # TODO: make this switch based on inherits()
  tmp <- switch(methods::is(x)[[1]], alignment = sqs_from_alignment(x),
                supermatrix = sqs_from_supermatrix(x),
                stop(paste0('Hmm.... ', obj('x'),
                            ' doesn\t seem to be a sequence object.')))
  seqinr::write.fasta(sequences = tmp[['sqs']], names = tmp[['nms']],
                      file.out = flpth, open = 'w', nbchar = width)
  invisible(file.exists(flpth))
}

# Public ----
#' @name alignment_read
#' @title Read in sequence alignment(s)
#' @description Return alignment list object from reading in sequences from
#' alignment file(s). Use \code{format} to specify alignment file type. If no
#' \code{format} given, the function will guess based on file extension.
#' @details Available file formats: "clustal", "msf", "mase", "phylip" and 
#' "fasta".
#' @param flpths File path(s) to alignment file(s)
#' @param format Sequnece file format. Default 'guess'. Character.
#' @return alignment list
#' @seealso \code{\link[seqinr]{read.alignment}}
#' @export
alignment_read <- function(flpths, format = c('guess', "clustal", "msf", "mase",
                                              "phylip", "fasta")) {
  format <- match.arg(format)
  res <- lapply(X = flpths, FUN = .alignment_read, format = format)
  class(res) <- 'alignment_list'
  names(res) <- tools::file_path_sans_ext(basename(path = flpths))
  res
}

#' @name sequences_write
#' @title Write sequences to file
#' @description Takes a sequence-based object and writes the sequences to file
#' in .FASTA format. For objects with multiple sequence groups, "supermatrices"
#' and "alignment_list", the \code{flpth} should point to a directory not a
#' file. Files are then created in this directory based on the group name for
#' each sequence group.
#' @param x sequences as "supermatrix", "supermatrices", "alignment" or
#' "alignment_list" objects.
#' @param flpth File path to output file or folder. Character.
#' @param width Max number of base pairs per line, default 80. Integer.
#' @return Logical
#' @export
sequences_write <- function(x, flpth, width = 80) {
  if (inherits(x = x, what = c('alignment', 'supermatrix'))) {
    res <- .sequences_write(flpth = flpth, x = x, width = width)
  } else if (inherits(x = x, what = c('alignment_list', 'supermatrices'))) {
    if (!dir.exists(flpth)) {
      stop(paste0(char(flpth), ' doesn\'t exist.'))
    }
    flpths <- file.path(flpth, paste0(names(x), '.fasta'))
    res <- vapply(X = seq_along(x), FUN = function(i) {
      .sequences_write(flpth = flpths[[i]], x = x[[i]], width = width)
    }, FUN.VALUE = logical(1))
  } else {
    msg <- paste0(obj(x), ' is not a sequence object')
    stop(msg, call. = FALSE)
  }
  all(res)
}
