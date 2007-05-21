# Allows conflicts. For more information, see library() and
# conflicts() in R base.
.conflicts.OK <- TRUE


.First.lib <- function(libname, pkgname) {
  pkg <- Package(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));

  # Assert that the binaries install successfully
  if (.Platform$OS.type == "windows") {
    binfile <- "cfit.exe";
  } else {
    binfile <- "cfit";
  }
  pathname <- system.file("bin", binfile, package=pkgname);
  if (is.null(pathname)) {
    stop("Could not locate 'bin/", binfile, "' in package '", getName(pkg), "' (v", getVersion(pkg), "). It seems like the installation of the package failed.  Please report this to ", getMaintainer(pkg), ".");
  }

  # cfit system command
  cfitPath <- system.file("bin", package=pkgname);
  pathname <- file.path(cfitPath, "cfit");

  # Set the default 'cfit' command (within quotation marks)
  cmd <- sprintf("\"%s\"", pathname);
  options("cfit"=cmd);

  cat(getName(pkg), " v", getVersion(pkg), " (", getDate(pkg), ")",
      " successfully loaded. See ?", pkgname, " for help.\n", sep="");
}

###########################################################################
# HISTORY:
# 2007-05-20
# o WORKAROUND: Put quotation marks around default 'cfit' command.
###########################################################################
