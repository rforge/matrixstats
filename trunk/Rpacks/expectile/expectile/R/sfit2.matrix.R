#########################################################################/**
# @set "class=matrix"
# @RdocMethod sfit2
# @aliasmethod fitSimplex
# @aliasmethod fitCone
#
# @title "Fit a simplex or polyhedral cone to multivariate data"
#
# \description{
#   @get "title" by decomposing data @matrix
#   \eqn{Y = X B + E}, where \eqn{B} is "mostly non-negative".
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A PxN @matrix (or @data.frame) containing P variables and 
#     N observations in \eqn{R^N}.}
#   \item{M}{Number of vertices, M-1 <= P}.
#   \item{w}{An optional @vector in [0,1] of length N specifying weight
#     for each observation.}
#   \item{lambda}{Vertex assigment parameters.}
#   \item{alpha}{A @double specifying the desired expectile.}
#   \item{family}{A @character string specifying the ....}
#   \item{robustConst}{A @double constant multiplier of MAR scale estimate.}
#   \item{tol}{A @double tolerance for expectile estimation.}
#   \item{maxIter}{The maximum number of iterations in estimation step.}
#   \item{Rtol}{A @double tolerance in linear solve, 
#      before a vertex is ignored.} 
#   \item{initSimplex}{A user-supplied initial simplex, otherwise automatic.}
#   \item{fitCone}{If @TRUE, the first vertex is treated as an apex and
#     the opposite face has its own residual scale estimator.}
#   \item{verbose}{if @TRUE, iteration progress is printed to standard error.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list structure elements:
#   \itemize{
#     \item{X}{the fitted simplex, as a PxM @matrix.}
#     \item{B}{Affine coefficients, as an MxN @matrix.}
#   }
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
# \examples{@include "..\incl\sfit2.matrix.Rex"}
#
# \author{
#   Pratyaksha (Asa) Wirapati, \email{wirapati@wehi.edu.au}.
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
setMethodS3("sfit2", "matrix", function(y, M, w=rep(1,dim(y)[2]),
            lambda=2, alpha=0.05, 
            family=c("biweight", "huber", "normal"), robustConst=4.685,
            tol=0.001, maxIter=60, Rtol=1e-7, 
            initSimplex=NULL,
            fitCone=FALSE, verbose=FALSE, ...) {
  P <- dim(y)[1];
  N <- dim(y)[2];

  # Argument 'M':  
  if(P < M-1) {
    stop("too many vertices for the data dimension");
  }

  # Argument 'lambda':

  # Argument 'alpha':

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

  # Argument 'initSimplex':
  if(is.null(initSimplex)) {
    X <- matrix(0.0, nrow=P, ncol=M);
    autoInit <- 1;
  } else {
    if (!identical(dim(initSimplex), c(P,M))) {
      throw("Argument 'initSimplex' has the incorrect dimension: ", 
            nrow(initSimplex), "x", ncol(initSimplex), " != ", P, "x", "M");
    }
    X <- initSimplex; 
    autoInit <- 0;
  }


  familyCode <- switch(family,
    normal   = 0,
    huber    = 1,
    biweight = 2
  );

  fit <- .C("Rwrapper_sfit", as.integer(N), as.integer(P), as.double(y),
    as.double(w), as.integer(M), as.integer(autoInit), as.double(lambda),
    as.double(alpha), as.integer(familyCode), as.double(robustConst),
    as.integer(fitCone), as.integer(verbose),
    as.double(tol), as.integer(maxIter), as.double(Rtol),
    X=as.double(X), Beta=double(M*N),
    PACKAGE="expectile");

  dim(fit$X) <- c(P,M);
  dim(fit$Beta) <- c(M,N);

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

  fit;
}, protected=TRUE)


###########################################################################
# HISTORY:
# 2008-04-12
# o Added fitCone() and fitSimplex(), which are simply wrappers to sfit().
###########################################################################
