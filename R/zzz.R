# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("INLAjoint", pkgname, libname)
# }

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    # Disable INLA's version checking to avoid network calls
    try({
      INLA::inla.setOption(inla.timeout = 0)
      suppressMessages(attachNamespace("INLA"))
    }, silent = TRUE)
  }
}

# Avoid check note
if(getRversion() >= "2.15.1") utils::globalVariables(c(".", "V2"))

# Internal INLA functions
# These are needed for predictions functionality
.inla_tempdir_safe <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA package is required for predictions.")
  }
  INLA::inla.tempdir()
}

.inla_run_many_safe <- function(n, wd, num.threads = 1, cleanup = TRUE, verbose = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA package is required for predictions.")
  }
  INLA::inla.run.many(n, wd, num.threads = num.threads, cleanup = cleanup, verbose = verbose)
}

INLAjointStartupMessage <- function()
{

  msg <- c(paste0("INLAjoint version ",
utils::packageVersion("INLAjoint")),
"\nType 'citation(\"INLAjoint\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- INLAjointStartupMessage()
  msg[1] <- paste("Package 'INLAjoint' version", utils::packageVersion("INLAjoint"))
  packageStartupMessage(msg)

  # Check INLA
  if (!requireNamespace("INLA", quietly = TRUE)) {
    packageStartupMessage("Note: INLA package is required but not installed.")
    packageStartupMessage("Install with: install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/testing'), dep=TRUE)")
  }

  invisible()
}


