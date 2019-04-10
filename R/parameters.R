default_parameters <- list('max_name_dist' = 0.15, 'column_cutoff' = 0.5,
                           'tip_cutoff' = 0.1, 'min_ntips' = 5L,
                           'max_ntips' = 100L, 'min_ngenes' = 5L,
                           'min_nbps' = 250L)

#' @name parameters
#' @title Set and get parameters
#' @description Set new and default parameters. Look up currently set
#' parameters.
#' @param parameter Name of parameter
NULL

#' @rdname parameters
#' @export
default_pset <- function() {
  options(gaius_parameters = default_parameters)
  invisible(TRUE)
}

#' @rdname parameters
#' @export
default_pget <- function() {
  default_parameters
}

#' @rdname parameters
#' @param val Value to be assigned to parameter
#' @export
pset <- function(val, parameter = names(default_pget())) {
  if (length(parameter) != length(val)) {
    stop(paste0(char(val), " and ", char(parameter), " not the same length."))
  }
  nm <- match.arg(arg = parameter, several.ok = TRUE)
  prmtrs <- options()[['gaius_parameters']]
  for (i in seq_along(nm)) {
    expected_class <- class(default_parameters[[nm[[i]]]])
    if (ptype_match(val1 = val[[i]], val2 = default_parameters[[nm[[i]]]])) {
      prmtrs[[nm[[i]]]] <- val[[i]]
    } else {
      msg <- paste0(char(nm[[i]]), ' must be a ', obj(expected_class),
                    ' object')
      stop(msg)
    }
  }
  options(gaius_parameters = prmtrs)
  invisible(TRUE)
}

#' @rdname parameters
#' @export
pget <- function(parameter = names(default_pget())) {
  nm <- match.arg(arg = parameter, several.ok = TRUE)
  prmtrs <- options()[['gaius_parameters']]
  if (length(nm) == 1) {
    res <- prmtrs[[nm]]
  } else {
    res <- prmtrs[nm]
  }
  res
}

ptype_match <- function(val1, val2) {
  is.numeric(val1) == is.numeric(val2) &
    is.character(val1) == is.character(val2) &
    is.logical(val1) == is.logical(val2)
}

.onAttach <- function(...) {
  default_pset()
}