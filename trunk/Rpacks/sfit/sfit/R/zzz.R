# Allows conflicts. For more information, see library() and
# conflicts() in R base.
.conflicts.OK <- TRUE


.First.lib <- function(libname, pkgname) {
  # Add this package's cdf/ directory to the CDF path
  cfitPath <- system.file("bin", package=pkgname);
  options("cfit"=file.path(cfitPath, "cfit"));

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

  cat(getName(pkg), " v", getVersion(pkg), " (", getDate(pkg), ")",
      " successfully loaded. See ?", pkgname, " for help.\n", sep="");
}
