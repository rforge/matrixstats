%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  sfit2.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{sfit2.matrix}
\alias{sfit2.matrix}
\alias{sfit2.matrix}

 \alias{sfit2}
 \alias{fitExpectileCone}
 \alias{matrix.fitExpectileCone}
\alias{fitExpectileCone.matrix}
\alias{fitExpectileCone,matrix-method}

 \alias{fitSimplex}
 \alias{matrix.fitSimplex}
\alias{fitSimplex.matrix}
\alias{fitSimplex,matrix-method}

 \alias{fitCone}
 \alias{matrix.fitCone}
\alias{fitCone.matrix}
\alias{fitCone,matrix-method}


 \title{Fit a simplex or polyhedral cone to multivariate data}

 \description{
   Fit a simplex or polyhedral cone to multivariate data by decomposing data PxN \code{\link[base]{matrix}}
   \eqn{Y = X B + E}, where
   \eqn{X} is a PxM \code{\link[base]{matrix}},
   \eqn{B} is a "mostly non-negative" MxN \code{\link[base]{matrix}}, and
   \eqn{E} is a PxN \code{\link[base]{matrix}} of noise, all with \eqn{M-1 \leq P}.
 }

 \usage{\method{sfit2}{matrix}(y, M=dim(y)[1] + 1, w=rep(1, dim(y)[2]), lambda=2, alpha=0.05, family=c("biweight", "huber", "normal"), robustConst=4.685, tol=0.001, maxIter=60, Rtol=1e-07, priorX=NULL, priorW=NULL, initX=NULL, fitCone=FALSE, verbose=FALSE, ...)}

 \arguments{
   \item{y}{A PxN \code{\link[base]{matrix}} (or \code{\link[base]{data.frame}}) containing P variables and
     N observations in \eqn{R^N}.}
   \item{M}{Number of vertices, M-1 <= P.}.
   \item{w}{An optional \code{\link[base]{vector}} in [0,1] of length N specifying weight
     for each observation.}
   \item{lambda}{A scalar vertex assigment parameter in [1,Inf).}
   \item{alpha}{A \code{\link[base]{double}} in [0,1] specifying the desired expectile.}
   \item{family}{A \code{\link[base]{character}} string specifying the ....}
   \item{robustConst}{A \code{\link[base]{double}} constant multiplier of MAR scale estimate.}
   \item{tol}{A positive \code{\link[base]{double}} tolerance for expectile estimation.}
   \item{maxIter}{The maximum number of iterations in estimation step.}
   \item{Rtol}{A postive \code{\link[base]{double}} tolerance in linear solve,
      before a vertex is ignored.}
   \item{priorX, priorW}{(Optional) Prior simplex PxM \code{\link[base]{matrix}} and
      M vertex weights.  An \code{\link[base:is.finite]{Inf}} weight corresponds to a fixed vertex.
      If \code{\link[base]{NULL}}, no priors are used.
   }
   \item{initX}{(Optional) An initial simplex PxM \code{\link[base]{matrix}} ('X').
      If \code{\link[base]{NULL}}, the initial simplex is estimated automatically.}
   \item{fitCone}{If \code{\link[base:logical]{TRUE}}, the first vertex is treated as an apex and
     the opposite face has its own residual scale estimator.}
   \item{verbose}{if \code{\link[base:logical]{TRUE}}, iteration progress is printed to standard error.}
   \item{...}{Not used.}
 }

 \value{
   Returns a named \code{\link[base]{list}} structure elements:
   \item{X}{the fitted simplex, as a PxM \code{\link[base]{matrix}}.}
   \item{B}{Affine coefficients, as an MxN \code{\link[base]{matrix}}.}
 }

 \details{
   Given multidimensional data matrix Y with P rows (variables)
   and N columns (observations), decompose Y into two matrices,
   X (P-by-M) and B (M-by-N) as
     \eqn{Y = X B + E},
   where P may be larger than M-1.

   In simplex fitting mode, \eqn{B_j} for each observation
   sums to one, and mostly non-negative. The columns of X are the
   estimated vertices of the simplex enclosing most points.
 
   In cone fitting mode, the first column of X is apex of the cone, while
   the others are directions of the rays emanating from the apex, with
   the vector norms standardized to one. The first row of B is
   always equal to one, and the remaining rows are mostly non-negative.
   They don't necessarily sum to one.
 }

 \examples{# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example with simulated data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of observations
n <- 20000

# Offset and cross talk
a0 <- c(50,300)
A <- matrix(c(1,0.2,0.5,1), nrow=2, ncol=2)  # cross-talk

# the true signal is joint gamma
z <- matrix(rgamma(2*n, shape=0.25, scale=100), nrow=2, ncol=n)

# Observed signal plus Gaussian error
eps <- matrix(rnorm(2*n, mean=0, sd=10), nrow=2, ncol=n)
y <- A \%*\% z + a0 + eps

 
layout(matrix(1:4, nrow=2, byrow=TRUE))
par(mar=c(5,4,2,2)+0.1)
lim <- c(0,1000)
xlab <- expression(y[1])
ylab <- expression(y[2])

for (withPrior in c(FALSE, TRUE)) {
  if (withPrior) {
    priorX <- matrix(c(a0, 0,0, 0,0), nrow=2, ncol=3)
    priorW <- c(Inf,0,0)
    # Fit cone
    fit <- fitCone(y, priorX=priorX, priorW=priorW)
    stopifnot(identical(fit$X[,1], a0))
  } else {
    # Fit cone
    fit <- fitCone(y)
    fit0 <- fit
  }

  cat("Estimated cone:\n")
  print(fit$X)
 
  plot(t(y), pch=".", xlim=lim, ylim=lim, xlab=xlab, ylab=ylab)
  points(fit, pch=19, cex=1.5, col="#aaaaaa")
  radials(fit, col="#aaaaaa", lwd=2)
  drawApex(fit, pch=19, cex=1, col="tomato")
  lines(fit, col="tomato", lwd=2)
 
 
  # The rectified data points
  xlab <- expression(hat(x)[1])
  ylab <- expression(hat(x)[2])
  plot(t(fit$Beta[2:3,]), pch=".", xlab=xlab, ylab=ylab)
  points(0,0, pch=19, cex=1.5, col="tomato") # the apex
  lines(c(0,0,lim[2]), c(lim[2],0,0), lwd=2, col="tomato")
}
}

 \author{
   Algorithm and native code by Pratyaksha (Asa) Wirapati.
   R interface by Henrik Bengtsson.
 }

 \references{
  [1] P. Wirapati, & T. Speed, \emph{Fitting polyhedrial cones and
     simplices to multivariate data points}, Walter and Eliza Hall Institute
     of Medical Research, December 30, 2001.\cr
  [2] P. Wirapati and T. Speed, \emph{An algorithm to fit a simplex
     to a set of multidimensional points}, Walter and Eliza Hall Institute
     of Medical Research, January 15, 2002.\cr
 }
\keyword{methods}
