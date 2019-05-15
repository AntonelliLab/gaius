outsider_install <- function(repo, service) {
  if (!outsider::is_module_installed(repo = repo)) {
    message(paste0('Outsider module ', char(repo), ' not available.',
                   '\nAttempting to install automatically.',
                   ' Note: user permission required.'))
    res <- outsider::module_install(repo = repo, service = service,
                                    force = FALSE)
    if (!res) {
      warning(paste0('Failed to install ', char(repo)))
    }
  }
  invisible(TRUE)
}