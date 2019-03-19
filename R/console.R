char <- function(x) {
  crayon::green(encodeString(x, quote = "'"))
}

stat <- function(...) {
  crayon::blue(...)
}

func <- function(x) {
  crayon::red(paste0(x, '()'))
}

obj <- function(x) {
  crayon::red(encodeString(x, quote = "'"))
}

seq <- function(x) {
  x <- strsplit(x = x, split = '')[[1]]
  x[x == 'a'] <- crayon::bgBlue('a')
  x[x == 't'] <- crayon::bgGreen('t')
  x[x == 'c'] <- crayon::bgRed('c')
  x[x == 'g'] <- crayon::bgYellow('g')
  paste0(x, collapse = '')
}

cat_line <- function(...) {
  cat(paste0(..., "\n"), sep = "")
}