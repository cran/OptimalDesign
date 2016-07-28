od.plot <- function (w, X=NULL, main="", del=1e-3) 
{
  ## Check the input.

  if (!is.vector(w) || !is.numeric(w) || !all(is.finite(w)) || (min(w) < 0)) stop("w must be a vector of non-negative real numbers")
  if (!((is.vector(X) || is.matrix(X)) && is.numeric(X)) && !is.null(X)) stop("X must be a numeric vector, numeric matrix, or NULL")
  if (!is.vector(del) || !is.numeric(del) || (prod(is.finite(del)) == 0) || (length(del) > 1) || !(min(del) > 0)) stop("del must be a positive real number")
  n <- length(w)
  if (is.null(X)) {
    X <- as.matrix(1:n)
  } else if (is.vector(X)) {
   X <- as.matrix(X)
  }
  if (dim(X)[1] != n) stop("the sizes of w and X are not compatible")

  ## Plot the graph of the design.

  if (dim(X)[2] == 1) {
    plot(as.vector(X), w, type = "h", lwd = 2, main = main)
  } else {
    wd <- w/sum(w)
    if (max(wd) < del) { 
      print("The graph is empty, because all design weights are smaller than del.", quote=FALSE)
    } else {
      pos <- wd > del
      n.pos <- sum(pos)
      sizes <- rep(7 * sqrt(del), n)
      sizes[pos] <- 7 * sqrt(wd[pos])
      cols <- rep("#CCCCCC", n)
      cols[pos] <- rainbow(n.pos)
      plot(as.data.frame(X), col = cols, cex = sizes, pch = 19, main = main, asp = 1)
    }
  }

}
