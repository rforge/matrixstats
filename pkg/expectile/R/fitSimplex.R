setMethodS3("fitSimplex", "matrix", function(y, ...) {
  sfit2(y, ..., fitcone=FALSE);
})



###########################################################################
# HISTORY:
# 2008-04-12
# o Added fitCone() and fitSimplex(), which are simply wrappers to sfit().
###########################################################################
