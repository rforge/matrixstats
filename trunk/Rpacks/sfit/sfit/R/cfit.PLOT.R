setMethodS3("getEdges", "cfit", function(object, ...) {
  u <- list();
  for (ii in 1:(nrow(object)-1)) {
    for (jj in (ii+1):nrow(object)) {
      idx <- c(ii,jj);;
      u <- c(u, list(object[idx,]));
    }
  }

  u;
}, private=TRUE) # getEdges.cfit()


setMethodS3("points", "cfit", function(x, dim=c(1,2), ...) {
  # To please R CMD check
  object <- x;

  xyz <- object[,dim];
  points(xyz, ...);
}) # points.cfit()


setMethodS3("lines", "cfit", function(x, dim=c(1,2), ...) {
  # To please R CMD check
  object <- x;

  u <- getEdges(object);

  for (ii in seq(u)) {
    xy <- u[[ii]][,dim];
    lines(xy, ...);
  }
}) # lines.cfit()


setMethodS3("points3d", "cfit", function(object, dim=c(1,2,3), ...) {
  xyz <- object[,dim];
  points3d(xyz, ...);
}) # points3d.cfit()


setMethodS3("lines3d", "cfit", function(object, dim=c(1,2,3), ...) {
  u <- getEdges(object);

  for (ii in seq(u)) {
    xyz <- u[[ii]][,dim];
    lines3d(xyz, ...);
  }
}) # lines3d.cfit()


###########################################################################
# HISTORY:
# 2007-06-10
# o BUG FIX: lines3d() for 'cfit' queried non-existing objects.
# 2006-05-07
# o Created.  For now, these functions are only for internal use and
#   for the examples.
###########################################################################
