%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  loexp.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{loexp}
\alias{loexp}


\title{Expectile regression for time series curve}

\description{
  Expectile regression for time series curve by decomposing data \code{\link[base]{matrix}}
  \eqn{Y = X B + E}, where \eqn{B} is "mostly non-negative".
}

\usage{loexp(y, w=rep(1, length(y)), sigma=40, polyo=2, alpha=0.5, biweight=4.685, tol=1e-04, maxIter=50)}

\arguments{
 \item{y}{A \code{\link[base]{numeric}} time series \code{\link[base]{vector}} of length T.}
 \item{w}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of T case weights.
    If unspecificed, all weights are set to 1.}
 \item{sigma}{The standard deviation of Gaussian kernel.}
 \item{polyo}{The order of polynomial, currently can be 0, 1, or 2.}
 \item{alpha}{The desired expectile.}
 \item{biweight}{Parameter used in Tukey's biweight function.}
 \item{tol}{The tolerance for expectile estimation.}
 \item{maxIter}{The maximum number of iterations in estimation step.}
}

\value{
 A \code{\link[base]{list}} with components:
  \item{intparams}{An \code{\link[base]{integer}} \code{\link[base]{vector}} parameters: return status
    (0=,1=,2=,3=), length of input vector, \code{polyo} and \code{maxit}.}
  \item{dblparams}{A \code{\link[base]{double}} \code{\link[base]{vector}} of parameters:
    sigma, alpha, biweight and tolerance.}
  \item{y}{The input \code{\link[base]{vector}} of data.}
  \item{w}{The input \code{\link[base]{vector}} of case weights.}
  \item{outy}{The output \code{\link[base]{vector}} of expectiles for the input data \code{y}.}
  \item{outw}{The output \code{\link[base]{vector}} of weights used in fitting.}
}

\author{Asa Wirapati, Mark Robinson}

\references{
  [1] ...
}

\examples{
data(loexp.example)

plot(ex1, pch=19, type="b")

# 50th expectile (median)
lines(loexp(ex1, alpha=0.5)$outy, lwd=4, col="blue")

# 0.5th expectile (baseline)
lines(loexp(ex1, alpha=0.005)$outy, lwd=4, col="red")

plot(ex2.y, pch=19, type="b")
# give weight=0 for 0s
lines(loexp(ex2.y, w=ex2.w, alpha=0.005)$outy, lwd=4, col="blue")
}

