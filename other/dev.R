# TODO:
# - alignment function
# - phylogeny function
# - supertree function


readSqs <- function(fl) {
  all_data <- readLines(fl)
  sqs <- list()
  for(i in seq_along(all_data)) {
    bit <- all_data[[i]]
    if(grepl(pattern='^>', x=bit)) {
      nm <- sub(pattern='^>', '', x=bit)
      sqs[[nm]] <- ''
    } else {
      sqs[[nm]] <- paste0(sqs[[nm]], bit)
    }
  }
  sqs
}

writeSqs <- function(sqs, fl) {
  fasta <- ''
  for(i in seq_along(sqs)) {
    sq <- sqs[[i]]
    dfln <- names(sqs)[[i]]
    fasta <- paste0(fasta, '>', dfln, '\n',
                    sq, '\n\n')
  }
  cat(fasta, file=fl)
}

partition <- function(lngths, fl) {
  gene <- strt <- 1
  prttn_txt <- ''
  for(lngth in lngths) {
    end <- lngth + strt - 1
    prttn_txt <- paste0(prttn_txt, 'DNA, gene',
                   gene, ' = ', strt, '-',
                   end, '\n')
    strt <- end + 1
    gene <- gene + 1
  }
  cat(prttn_txt, file=fl)
}

wd <- '~/Desktop/supermatrix_testing/'
fl <- file.path(wd, 'alignment1.fasta')

sqs1 <- readSqs(file.path(wd, 'alignment1.fasta'))
sqs2 <- readSqs(file.path(wd, 'alignment2.fasta'))
lngths <- c(nchar(sqs1[[1]]), nchar(sqs2[[1]]))
partition(lngths, fl=file.path(wd, 'partition.txt'))
fllr1 <- paste0(rep('-', lngths[[1]]), collapse='')
fllr2 <- paste0(rep('-', lngths[[2]]), collapse='')
all_nms <- unique(c(names(sqs1), names(sqs2)))
all_nms <- sort(all_nms)
sprmtrx <- vector('list', length=length(all_nms))
names(sprmtrx) <- all_nms
for(nm in all_nms) {
  sq1 <- sqs1[[nm]]
  sq1 <- ifelse(is.null(sq1), fllr1, sq1)
  sq2 <- sqs2[[nm]]
  sq2 <- ifelse(is.null(sq2), fllr2, sq2)
  sprmtrx[[nm]] <- paste0(sq1, sq2)
}
all(sapply(sprmtrx, nchar) == nchar(sprmtrx[[1]]))
sprmtrx <- sprmtrx[!grepl('^Arecaceae', names(sprmtrx))]
writeSqs(sprmtrx, fl=file.path(wd, 'supermatrix.fasta'))
