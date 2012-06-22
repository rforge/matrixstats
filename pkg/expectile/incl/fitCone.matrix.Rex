# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
y <- A %*% z + a0 + eps

  
layout(matrix(1:4, nrow=2, byrow=TRUE))
par(mar=c(5,4,2,2)+0.1)
lim <- c(0,1000)
xlab <- expression(y[1])
ylab <- expression(y[2])

for (withPrior in c(FALSE, TRUE)) {
  if (withPrior) {
    priorX <- matrix(c(a0, 0,0, 0,0), nrow=2, ncol=3)
    priorW <- c(Inf,0,0)
    priorW <- c(+100,0,0)
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
