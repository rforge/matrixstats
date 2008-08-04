#########################################################################/**
# @RdocFunction loexp
#
# @title "Expectile regression for time series curve"
#
# \description{
#   @get "title" by decomposing data @matrix
#   \eqn{Y = X B + E}, where \eqn{B} is "mostly non-negative".
# }
#
# @synopsis
#
# \arguments{
#  \item{y}{A @numeric time series @vector of length T.}
#  \item{w}{A @numeric @vector of T case weights.
#     If unspecificed, all weights are set to 1.}
#  \item{sigma}{The standard deviation of Gaussian kernel.}
#  \item{polyo}{The order of polynomial, currently can be 0, 1, or 2.}
#  \item{alpha}{The desired expectile.}
#  \item{biweight}{Parameter used in Tukey's biweight function.}
#  \item{tol}{The tolerance for expectile estimation.}
#  \item{maxIter}{The maximum number of iterations in estimation step.}
# }
#
# \value{
#  A @list with components:
#   \item{intparams}{An @integer @vector parameters: return status
#     (0=,1=,2=,3=), length of input vector, \code{polyo} and \code{maxit}.}
#   \item{dblparams}{A @double @vector of parameters: 
#     sigma, alpha, biweight and tolerance.}
#   \item{y}{The input @vector of data.}
#   \item{w}{The input @vector of case weights.}
#   \item{outy}{The output @vector of expectiles for the input data \code{y}.}
#   \item{outw}{The output @vector of weights used in fitting.}
# }
#
# \author{Asa Wirapati, Mark Robinson}
#
# \references{
#   [1] ...
# }
#
# @examples "../incl/loexp.Rex"
#
#*/#########################################################################
loexp <- function(y, w=rep(1,length(y)), sigma=40, polyo=2, alpha=0.5, biweight=4.685, tol=0.0001, maxIter=50) {
  T <- length(y);

  # Argument 'w':
  if(length(w) != T) {
    stop("The number of elements in argument 'w' and argument 'y' does not match: ", length(w), " != ", T);
  }

  outy <- rep(0.0, times=T);
  outw <- rep(0.0, times=T);
  intparams <- c(-1, T, polyo, maxIter);
  dblparams <- c(sigma, alpha, biweight, tol);
  #cat("intparams in R:",intparams,"\n");
  #cat("dblparams in R:",dblparams,"\n");
  .C("call_loexp", 
     intparams=as.integer(intparams), dblparams=as.double(dblparams),
     y=as.double(y), w=as.double(w),
     outy=as.double(outy), outw=as.double(outw),
     PACKAGE="expectile");
} # loexp()


###########################################################################
# HISTORY:
# 2008-03-24
# o Created by Mark Robinson (mrobinson@wehi.edu.au).
###########################################################################



