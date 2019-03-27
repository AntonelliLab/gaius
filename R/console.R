char <- function(x) {
  crayon::green(encodeString(x, quote = "'"))
}

elem <- function(x) {
  crayon::silver(crayon::italic(x))
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

sq <- function(x) {
  x <- strsplit(x = x, split = '')[[1]]
  pull <- !x %in% c('a', 't', 'c', 'g')
  x[pull] <- crayon::bgBlack(crayon::white(x[pull]))
  x[x == 'a'] <- crayon::bgBlue(crayon::white('a'))
  x[x == 't'] <- crayon::bgGreen(crayon::white('t'))
  x[x == 'c'] <- crayon::bgRed(crayon::white('c'))
  x[x == 'g'] <- crayon::bgYellow(crayon::white('g'))
  x <- paste0(x, collapse = '')
  x
}

cat_line <- function(...) {
  cat(paste0(..., "\n"), sep = "")
}