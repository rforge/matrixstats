setMethodS3("getEdges", "cfit", function(object, ...) {
  u <- list();
  for (i in 1:(nrow(object)-1)) {
    for (j in (i+1):nrow(object)) {
      idx <- c(i,j);;
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

  for (i in seq(u)) {
    xy <- u[[i]][,dim];
    lines(xy, ...);
  }
}) # lines.cfit()


setMethodS3("points3d", "cfit", function(object, dim=c(1,2,3), ...) {
  xyz <- object[,dim];
  points3d(xyz, ...);
}) # points3d.cfit()


setMethodS3("lines3d", "cfit", function(object, dim=c(1,2,3), ...) {
  u <- getEdges(object);

  for (i in seq(u)) {
    xyz <- u[[i]][,dim];
    lines3d(xyz, col=cols[iter]);
  }
}) # lines3d.cfit()


###########################################################################
# HISTORY:
# 2006-05-07
# o Created.  For now, these functions are only for internal use and
#   for the examples.
###########################################################################
