od.crit <- function(F, w, crit="D", R=NULL, tol=1e-12){

## Verify F
if (!is.matrix(F) || !is.numeric(F) || !all(is.finite(F))) stop("F must be a matrix of real numbers.")
n <- dim(F)[1]; m <- dim(F)[2]

## Verify w
if (!is.vector(w) || !is.numeric(w) || !all(is.finite(w)) || (min(w) < 0)) stop("w must be a vector of non-negative real numbers.")
if (length(w) != n) stop("The length of w must be equal to the number dim(F)[1] of design points.")

## Verify crit
if (!is.vector(crit) || !is.character(crit) || (length(crit) != 1) || !is.element(crit, c("D", "A", "IV"))) stop("crit must be 'D', 'A', or 'IV'.")

## Verify R
if (is.null(R)) R <- 1:n
if (!is.vector(R) || !is.numeric(R) || (length(unique(R)) != length(R)) || !all(is.element(R, 1:n))) stop("R must be NULL, or a vector, the elements of which form a subset of 1,2,...,dim(F)[1].")

if (crit=="IV") {
  M <- od.infmat(F.IVtoA(F, R), w)
} else {
  M <- od.infmat(F,w)
}

if (min(abs(eigen(M)$values)) <= m*tol){
  return(0)
} else {
  if (crit=="A" || crit=="IV") return(m*sum(diag(solve(M)))^(-1))
  if (crit=="D") return((abs(det(M)))^(1/m))
}
}


