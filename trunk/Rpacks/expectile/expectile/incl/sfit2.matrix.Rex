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
fit <- fitCone(y, M=3)

layout(matrix(1:4, ncol=2, byrow=TRUE))
lim <- c(0,1000)
plot(t(y), xlim=lim, ylim=lim, pch=".")
X <- fit$X
points(t(X[,1]), pch=19, cex=1.5, col="tomato") # the apex
lines(t(X[,1] + lim[2]*cbind(X[,2], 0, X[,3])), lwd=2, col="tomato")

# The rectified data points
plot(t(fit$Beta[2:3,]), pch=".")
points(0,0, pch=19, cex=1.5, col="tomato") # the apex
lines(c(0,0,lim[2]), c(lim[2],0,0), lwd=3, col="tomato")
