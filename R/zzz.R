.onLoad <- function(libname, pkgname) {
  kw <- tryCatch(ergm::ergm_keyword(), error = function(e) NULL)
  if (is.data.frame(kw) && !("signed" %in% kw$name)) {
    try(ergm::ergm_keyword(name = "signed", short = "signed",
                           description = "for signed networks built by as_signed_network()",
                           popular = FALSE, package = "Bstergm"), silent = TRUE)
  }
}
