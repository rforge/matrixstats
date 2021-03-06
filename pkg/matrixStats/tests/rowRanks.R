library("matrixStats")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Consistency checks
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(1)

cat("Consistency checks:\n")
K <- if (Sys.getenv("_R_CHECK_FULL_") == "") 4 else 20
for (kk in 1:K) {
  cat("Random test #", kk, "\n", sep="")

  # Simulate data in a matrix of any shape
  dim <- sample(500:2000, size=2)
  n <- prod(dim)
  x <- rnorm(n, sd=10)
  dim(x) <- dim

  # Add NAs?
  if ((kk %% 4) %in% c(3,0)) {
    cat("Adding NAs\n")
    nna <- sample(n, size=1)
    x[sample(length(x), size=nna)] <- as.double(NA)
  }

  # Integer or double?
  if ((kk %% 4) %in% c(2,0)) {
    cat("Coercing to integers\n")
    storage.mode(x) <- "integer"
  }

  # rowRanks():
  y1 <- matrixStats::rowRanks(x)
  y2 <- t(apply(x, MARGIN=1, FUN=rank, na.last="keep", ties.method="max"))
  stopifnot(identical(y1,y2))

  y3 <- matrixStats::colRanks(t(x))
  stopifnot(identical(y1,y3))

  # colRanks():
  y1 <- matrixStats::colRanks(x)
  y2 <- t(apply(x, MARGIN=2, FUN=rank, na.last="keep", ties.method="max"))
  stopifnot(identical(y1,y2))

  y3 <- matrixStats::rowRanks(t(x))
  stopifnot(identical(y1,y3))
} # for (kk ...)
