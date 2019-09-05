#' @name partition_file
#' @title Generate partition file for RAxML
#' @description Write out partition text file for RAxML from supermatrix
#' information.
#' @param nbp Number of base pairs per cluster, character vector
#' @param type Type of cluster, e.g. DNA, gene, character vector
#' @param flpth File path for resulting file
#' @export
partition_file <- function(nbp, type = 'DNA, gene', flpth) {
  gene <- strt <- 1
  prttn_txt <- ''
  if (length(type) == 1) {
    type <- rep(type, length(nbp))
  } else {
    if (length(type) != length(nbp)) {
      stop('type and nbp must be same length or type of length 1')
    }
  }
  for (i in seq_along(nbp)) {
    end <- nbp[[i]] + strt - 1
    prttn_txt <- paste0(prttn_txt, type[[i]],
                        gene, ' = ', strt, '-',
                        end, '\n')
    strt <- end + 1
    gene <- gene + 1
  }
  cat(prttn_txt, file = flpth)
}