#########################################################################/**
# @set "class=matrix"
# @RdocMethod sfit2
# @alias sfit2
# @alias fitExpectileCone
# @aliasmethod fitExpectileCone
# @alias fitSimplex
# @aliasmethod fitSimplex
# @alias fitCone
# @aliasmethod fitCone
#
# @title "Fit a simplex or polyhedral cone to multivariate data"
#
# \description{
#   @get "title" by decomposing data PxN @matrix
#   \eqn{Y = X B + E}, where 
#   \eqn{X} is a PxM @matrix, 
#   \eqn{B} is a "mostly non-negative" MxN @matrix, and
#   \eqn{E} is a PxN @matrix of noise, all with \eqn{M-1 \leq P}.
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A PxN @matrix (or @data.frame) containing P variables and 
#     N observations in \eqn{R^N}.}
#   \item{M}{Number of vertices, M-1 <= P.}.
#   \item{w}{An optional @vector in [0,1] of length N specifying weight
#     for each observation.}
#   \item{lambda}{A scalar vertex assigment parameter in [1,Inf).}
#   \item{alpha}{A @double in [0,1] specifying the desired expectile.}
#   \item{family}{A @character string specifying the ....}
#   \item{robustConst}{A @double constant multiplier of MAR scale estimate.}
#   \item{tol}{A positive @double tolerance for expectile estimation.}
#   \item{maxIter}{The maximum number of iterations in estimation step.}
#   \item{Rtol}{A postive @double tolerance in linear solve, 
#      before a vertex is ignored.} 
#   \item{priorX, priorW}{(Optional) Prior simplex PxM @matrix and 
#      M vertex weights.  An @Inf weight corresponds to a fixed vertex.
#      If @NULL, no priors are used.
#   }
#   \item{initX}{(Optional) An initial simplex PxM @matrix ('X').
#      If @NULL, the initial simplex is estimated automatically.}
#   \item{fitCone}{If @TRUE, the first vertex is treated as an apex and
#     the opposite face has its own residual scale estimator.}
#   \item{verbose}{if @TRUE, iteration progress is printed to standard error.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list structure elements:
#   \item{X}{the fitted simplex, as a PxM @matrix.}
#   \item{B}{Affine coefficients, as an MxN @matrix.}
# }
#
# \details{
#   Given multidimensional data matrix Y with P rows (variables)
#   and N columns (observations), decompose Y into two matrices,
#   X (P-by-M) and B (M-by-N) as
#     \eqn{Y = X B + E},
#   where P may be larger than M-1.
#
#   In simplex fitting mode, \eqn{B_j} for each observation
#   sums to one, and mostly non-negative. The columns of X are the 
#   estimated vertices of the simplex enclosing most points.
#   
#   In cone fitting mode, the first column of X is apex of the cone, while
#   the others are directions of the rays emanating from the apex, with
#   the vector norms standardized to one. The first row of B is
#   always equal to one, and the remaining rows are mostly non-negative.
#   They don't necessarily sum to one.
# } 
#
# \examples{@include "..\incl\fitCone.matrix.Rex"}
#
# \author{
#   Algorithm and native code by Pratyaksha (Asa) Wirapati.
#   R interface by Henrik Bengtsson.
# }
#
# \references{
#  [1] P. Wirapati, & T. Speed, \emph{Fitting polyhedrial cones and
#     simplices to multivariate data points}, Walter and Eliza Hall Institute 
#     of Medical Research, December 30, 2001.\cr
#  [2] P. Wirapati and T. Speed, \emph{An algorithm to fit a simplex 
#     to a set of multidimensional points}, Walter and Eliza Hall Institute 
#     of Medical Research, January 15, 2002.\cr
# }
#
#*/#########################################################################
setMethodS3("sfit2", "matrix", function(y, M=dim(y)[1]+1, w=rep(1,dim(y)[2]),
            lambda=2, alpha=0.05, 
            family=c("biweight", "huber", "normal"), robustConst=4.685,
            tol=0.001, maxIter=60, Rtol=1e-7, 
            priorX=NULL, priorW=NULL,
            initX=NULL,
            fitCone=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dimY <- dim(y);
  P <- dimY[1];
  N <- dimY[2];

  # Argument 'y':
  y <- as.double(y);
  dim(y) <- dimY;

  # Argument 'M':
  if(P < M-1) {
    throw("Too many vertices (P=", P, ") for the data dimension (M=", M, "): P < M-1");
  }
  dimX <- c(P,M);

  # Argument 'lambda':
  lambda <- as.double(lambda);
  if (length(lambda) != 1)
    throw("Argument 'lambda' must a scalar.");
  if (!is.finite(lambda))
    throw("Argument 'lambda' is out of range: ", lambda);

  # Argument 'alpha':
  alpha <- as.double(alpha);
  if (length(alpha) != 1)
    throw("Argument 'alpha' must a scalar.");
  if (!is.finite(alpha) || alpha < 0 || alpha > 1)
    throw("Argument 'alpha' is out of range [0,1]: ", alpha);

  # Argument 'w':
  if (is.null(w)) {
    w <- rep(1, times=N);
  } else if (is.numeric(w)) {
    if (length(w) == 1) {
      w <- rep(w, times=N);
    } else if (length(w) != N) {
      throw("The length of argument 'w' does not equal the number of variables: ", length(w), " != ", N);
    }
  } else {
    throw("Argument 'w' must be numeric: ", class(w)[1]);
  }

  # Argument 'family':
  family <- match.arg(family);

  # Argument 'initX':
  if (is.null(initX)) {
    X <- matrix(0.0, nrow=P, ncol=M);
    autoInit <- 1;
  } else {
    if (!all.equal(dim(initX), dimX)) {
      throw("Argument 'initX' has the incorrect dimension: ", 
            nrow(initX), "x", ncol(initX), 
            " != ", P, "x", "M, where M=", M);
    }
    X <- as.double(initX); 
    autoInit <- 0;
  }

  # Argument 'priorX':
  if (!is.null(priorX)) {
    if (!all.equal(dim(priorX), dimX)) {
      throw("Argument 'priorX' has the incorrect dimension: ", 
            nrow(priorX), "x", ncol(priorX), 
            " != ", P, "x", "M, where M=", M);
    }
    X0 <- as.double(priorX);
  } else {
    X0 <- NULL;
  }

  # Argument 'priorW':
  if (!is.null(priorW)) {
    if (length(priorW) != M) {
      throw("The length of argument 'priorW' does not match the number of vertices: ", length(priorW), " != ", M); 
    }
    wX0 <- as.double(priorW);
    if (any(is.na(wX0)))
      throw("Argument 'priorW' contains missing values.");
    if (any(wX0 < 0))
      throw("Argument 'priorW' contains negative weights.");
    # Workaround for +Inf;
    wX0[is.infinite(wX0)] <- .Machine$double.xmax;
  } else {
    wX0 <- NULL;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Fit simplex
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Encode parameters
  familyCode <- switch(family,
    normal   = 0,
    huber    = 1,
    biweight = 2
  );

  fit <- .C("Rwrapper_sfit0", as.integer(N), as.integer(P), as.double(y),
    as.double(w), as.integer(M), as.integer(autoInit), as.double(lambda),
    as.double(alpha), as.integer(familyCode), as.double(robustConst),
    as.integer(fitCone), as.integer(verbose),
    as.double(tol), as.integer(maxIter), as.double(Rtol),
    X=X, Beta=double(M*N),
    wX0=wX0, X0=X0,
    PACKAGE="expectile");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dim(fit$X) <- dimX;
  dim(fit$Beta) <- c(M,N);
  if (!is.null(fit$X0))
    dim(fit$X0) <- dimX;
  if (!is.null(fit$wX0))
    dim(fit$wX0) <- M;

  # Update the names
  names <- names(fit);
  names[1] <- "P";
  names[2] <- "N";
  names[3] <- "y";
  names[4] <- "w";
  names[5] <- "M";
  names[6] <- "autoInit";
  names[7] <- "lambda";
  names[8] <- "alpha";
  names[9] <- "familyCode";
  names[10] <- "robustConst";
  names[11] <- "fitCone";
  names[12] <- "verbose";
  names[13] <- "tol";
  names[14] <- "maxIter";
  names[15] <- "Rtol";
  names(fit) <- names;

  class(fit) <- "sfit2";

  fit;
}, protected=TRUE)


###########################################################################
# HISTORY:
# 2008-09-08
# o Added R support for fitting with prior simplex.  Added example code.
# o WP updated native sfit to accept prior simplex with weights.
# 2008-09-03
# o Updated the validation of 'initMatrix'.
# o Added validation for 'lambda'.
# 2008-04-12
# o Added fitCone() and fitSimplex(), which are simply wrappers to sfit().
###########################################################################
