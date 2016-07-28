od.print <- function (w, X=NULL, del=1e-6) 
{
  ## Check the input.

  if (!is.vector(w) || !is.numeric(w) || (prod(is.finite(w)) == 0) || (min(w) < 0)) stop("w must be a vector of non-negative real numbers")
  if (!((is.vector(X) || is.matrix(X)) && is.numeric(X)) && !is.null(X)) stop("X must be a numeric vector, numeric matrix, or NULL")
  if (!is.vector(del) || !is.numeric(del) || !all(is.finite(del)) || (length(del) > 1) || !(min(del) >= 0)) stop("del must be a non-negative real number")
  n <- length(w)
  if (is.null(X)) {
    X <- matrix(nrow=n, ncol=0)
  } else if (is.vector(X)) {
   X <- as.matrix(X)
  }
  if (dim(X)[1] != n) stop("the sizes of w and X are not compatible")

  ## Print a compact form of the design.

  ind.supp <- (1:n)[w >= del]
  n.supp <- length(ind.supp)
  if(n.supp == 0){
    print("The print is empty, because all design weights are smaller than del.", quote=FALSE)
  } else {
    k <- dim(X)[2]
    X.supp <- matrix(nrow=n.supp, ncol=k)
    if(k > 0) {
      X.supp[1:n.supp, 1:k] <- X[ind.supp,]
      colnames(X.supp) <- colnames(X)
      if (is.null(colnames(X.supp)))
        colnames(X.supp) <- rep("", k)
    }
    mat.print <- cbind(X.supp, w[ind.supp])
    rownames(mat.print) <- ind.supp
    colnames(mat.print)[k + 1] <- "weight"
    print(mat.print)
  }

}
