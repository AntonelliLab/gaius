# 
# flpth <- gaius:::datadir_get('alignment.fasta')
# 
# fmt <- "fasta"
# 
# format_guess <- function(flpth) {
#   ext <- tools::file_ext(x = flpth)
#   if (grepl(pattern = '(clus|aln)', x = ext, ignore.case = TRUE)) {
#     fmt <- "clustal"
#   } else if (grepl(pattern = 'msf', x = ext, ignore.case = TRUE)) {
#     fmt <- "msf"
#   } else if (grepl(pattern = 'mase', x = ext, ignore.case = TRUE)) {
#     fmt <- "mase"
#   } else if (grepl(pattern = 'phy', x = ext, ignore.case = TRUE)) {
#     fmt <- "phylip"
#   } else {
#     
#   }
# }
# 
# if (!is.null(format)) {
#   fmt <- format
# }
# temp <- seqinr::read.alignment(alignment, format = fmt)
# temp$nam <- do.call("rbind", lapply(strsplit(temp$nam, " "), 
#                                     "[[", 1))
# output <- data.frame(strsplit(gsub("[\r\n]", "", unlist(temp$seq)), 
#                               split = ""))
# colnames(output) <- temp$nam
# output