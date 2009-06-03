# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Example with simulated data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Number of observations
n <- 20000

# Offset and cross talk
a0 <- c(100,250)
A <- matrix(c(1,0.2,0.5,1), nrow=2, ncol=2)  # cross-talk

# the true signal is joint gamma
z <- matrix(rgamma(2*n, shape=0.25, scale=100), nrow=2, ncol=n)

# Observed signal plus Gaussian error
eps <- matrix(rnorm(2*n, mean=0, sd=10), nrow=2, ncol=n)
y <- A %*% z + a0 + eps

# Fit cone
fit <- fitCone(y)

lim <- c(0,1000)
xlab <- expression(y[1])
ylab <- expression(y[2])
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
