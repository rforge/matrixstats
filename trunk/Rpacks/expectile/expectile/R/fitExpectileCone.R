setMethodS3("fitCone", "matrix", function(y, ...) {
  fitExpectileCone(y, ...);
})

setMethodS3("fitExpectileCone", "matrix", function(y, ...) {
  fit <- sfit2(y, ..., fitcone=TRUE);
  class(fit) <- c("ExpectileCone", class(fit));
  fit;
})


setMethodS3("points", "ExpectileCone", function(x, ...) {
  # To please R CMD check
  fit <- x;

  X <- fit$X;
  X <- t(X);
  points(X, ...);
})


setMethodS3("lines", "ExpectileCone", function(x, ...) {
  # To please R CMD check
  fit <- x;

  X <- fit$X;
  xy <- (X[,c(3,1,2)]-X[,1])+X[,1];
  xy <- t(xy);
  lines(xy, ...);
})


setMethodS3("drawApex", "ExpectileCone", function(fit, ...) {
  X <- fit$X;
  X <- X[,1];
  X <- t(X);
  points(X, ...);
})



setMethodS3("radials", "ExpectileCone", function(fit, ...) {
  X <- fit$X;
  usr <- matrix(par("usr"), nrow=2);
  stretch <- max(apply(usr, MARGIN=2, FUN=diff));
  xy <- stretch*(X[,c(3,1,2)]-X[,1]) + X[,1];
  xy <- t(xy);
  lines(xy, ...);
})



###########################################################################
# HISTORY:
# 2008-08-31
# o fitCone() calls fitExpectileCone() as is.
# o Added fitExpectileCone(), which returns an ExpectileCone object.
# o Added points(), lines() and radials() to ExpectileCone.
# 2008-04-12
# o Added fitCone() and fitSimplex(), which are simply wrappers to sfit().
###########################################################################
