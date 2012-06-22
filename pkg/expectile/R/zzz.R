# Allows conflicts. For more information, see library() and
# conflicts() in R base.
.conflicts.OK <- TRUE


#.onAttach <- function(libname, pkgname) {
.First.lib <- function(libname, pkgname) {
  library.dynam(pkgname, pkgname, libname, now=FALSE);

  pkgD <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pkgD$Version, " (", pkgD$Date, ")",
      " successfully loaded. See ?", pkgname, " for help.\n", sep="");
}


###########################################################################
# HISTORY:
# 2012-03-23
# o Now .onAttach() uses packageStartupMessage() instead of cat().
# 2008-04-12
# o Now we're finally are using a dynamically loaded library; no more
#   calling external binaries.
# 2007-06-12
# o Replaced .First.lib() with .onAttach().
# 2007-05-20
# o WORKAROUND: Put quotation marks around default 'cfit' command.
###########################################################################
