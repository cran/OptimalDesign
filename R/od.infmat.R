od.infmat <-
function(F,w){

## Verify F
if (!is.matrix(F) || !is.numeric(F) || !all(is.finite(F))) stop("F must be a matrix of real numbers.")
n <- dim(F)[1]; m <- dim(F)[2]

## Verify w
if (!is.vector(w) || !is.numeric(w) || !all(is.finite(w)) || (min(w) < 0)) stop("w must be a vector of non-negative real numbers.")
if (length(w) != n) stop("The length of w must be equal to the number dim(F)[1] of design points.")

one <- t(rep(1,m))
t((w %*% one) * F) %*% F
}
