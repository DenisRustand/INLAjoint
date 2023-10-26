# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("INLAjoint", pkgname, libname)
# }

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
  invisible()
}


