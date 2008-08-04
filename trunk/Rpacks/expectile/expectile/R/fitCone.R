setMethodS3("fitCone", "matrix", function(y, ...) {
  sfit2(y, ..., fitcone=TRUE);
})

setMethodS3("fitCone", "data.frame", function(y, ...) {
  fitCone(as.matrix(y), ...);
})


###########################################################################
# HISTORY:
# 2008-04-12
# o Added fitCone() and fitSimplex(), which are simply wrappers to sfit().
###########################################################################
