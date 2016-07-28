od.RCs <-
function(F, N, w0=NULL, crit="D", R=NULL, w1=NULL, kappa=0, tab=NULL, graph=NULL, t.max=120)
{

## Verify F
if (!is.matrix(F) || !is.numeric(F) || !all(is.finite(F))) stop("F must be a matrix of real numbers.")
n <- dim(F)[1]; m <- dim(F)[2]
if (n < m) stop("The number dim(F)[1] of design points must be greater or equal to the number dim(F)[2] of parameters of the model.")
if (m < 2) stop("The number dim(F)[2] of parameters must be at least 2. Use the procedure od.m1 for models with a one-dimensional parameter.")  

## Verify N
if (!is.vector(N) || !is.numeric(N) || !all(is.finite(N)) || (length(N) != 1) || (N <= 0) || (N != round(N))) stop("The required size N of the design must be a natural number.")
if (N < m) stop("The required size N of the design should not be smaller than the number of parameters.")

### condition number test
print(paste("Condition number of F'F is", rcond(t(F)%*%F)), quote=FALSE)

## Verify w0 
if (is.null(w0)) w0 <- rep(0, n)
if (!is.vector(w0) || !is.numeric(w0) || !all(is.finite(w0)) || (min(w0) < 0)) stop("w0 must be NULL, or a vector of non-negative real numbers.")
if (length(w0) != n) stop("The length of w0 must be equal to the number dim(F)[1] of design points.")
if (sum(w0) > N) stop("w0 must satisfy the size constraint sum(w0) <= N.")   # Pre ine ako RCs
if (sum(w0) > (N - 1)) stop("w0 must satisfy the size constraint sum(w0) <= N-1.")  # Pre RCs

## Verify w1 
if (is.null(w1)) w1 <- w0             
if (!is.vector(w1) || !is.numeric(w1) || !all(is.finite(w1)) || (min(w1) < 0)) stop("w1 must be NULL, or a vector of non-negative real numbers.")
if (length(w1) != n) stop("The length of w1 must be equal to the number dim(F)[1] of design points.")
if (sum(w1) > N) stop("w1 must satisfy the size constraint sum(w1) <= N.")
if (any(w0 > w1)) stop("w1 must satisfy the constraints w1 >= w0.")

## Verify crit
if (!is.vector(crit) || !is.character(crit) || (length(crit) != 1) || !is.element(crit, c("D", "A", "IV"))) stop("crit must be 'D', 'A', or 'IV'.")

## Verify R
if (is.null(R)) R <- 1:n
if (!is.vector(R) || !is.numeric(R) || (length(unique(R)) != length(R)) || !all(is.element(R, 1:n))) stop("R must be NULL, or a vector, the elements of which form a subset of 1,2,...,dim(F)[1].")

## Verify kappa
if (!is.vector(kappa) || !is.numeric(kappa) || !all(is.finite(kappa)) || (length(kappa) != 1) || (kappa < 0)) stop("kappa must be a non-negative real number.")

## Verify graph
if (!is.null(graph)) {
  if (is.vector(graph) && is.numeric(graph)) {
    if ((length(unique(graph)) != length(graph)) || !all(is.element(graph, 1:m))) stop("If graph is a numeric vector, its elements must form a subset of 1,2,...,dim(F)[2].")
  } else  if (is.vector(graph) && is.character(graph)) {
    if ((length(unique(graph)) != length(graph)) || !all(is.element(graph, colnames(F)))) stop("If graph is a character vector, its elements must form a subset of colnames(F).")
  } else {
    stop("graph must be a numeric or a character vector.")
  } 
}
     
## Verify tab
if (!is.null(tab)) {
  if (is.vector(tab) && is.numeric(tab)) {
    if ((length(unique(tab)) != length(tab)) || !all(is.element(tab, 1:m))) stop("If tab is a numeric vector, its elements must form a subset of 1,2,...,dim(F)[2].")
  } else if (is.vector(tab) && is.character(tab)) {
    if ((length(unique(tab)) != length(tab)) || !all(is.element(tab, colnames(F)))) stop("If tab is a character vector, its elements must form a subset of colnames(F).")
  } else {
    stop("tab must be a numeric or a character vector.")
  } 
}

## Verify t.max
if (!is.vector(t.max) || !is.numeric(t.max) || !all(is.finite(t.max)) || (length(t.max) != 1) || (t.max <= 0)) stop("t.max must be a positive real number.")


if (crit=="D") r <- od.D.RCs(F=F, N=N, w0=w0, w1=w1, t.max=t.max, kappa=kappa)  
if (crit=="A") r <- od.A.RCs(F=F, N=N, w0=w0, w1=w1, t.max=t.max, kappa=kappa)  
if (crit=="IV") r<-od.A.RCs(F=F.IVtoA(F, R), N=N, w0=w0, w1=w1, t.max=t.max, kappa=kappa)


if (!is.null(graph) && (r$Phi.best>1e-6)){
  title <- paste("od.RCs: n=", n, ", m=", m, ", N=", N ,", crit=", crit ,", Phi=", sep="")
  title <- paste(title, round(r$Phi.best, 6), ", t=", r$t.act, " secs", sep="")
  od.plot(r$w.best, F[, graph], main=title)
}

if (!is.null(tab)) {
        od.print(r$w.best, F[, tab])
    } else {
        od.print(r$w.best)
    }



res <- list(method= "Resource constraints heuristic", w.best=r$w.best, Phi.best=r$Phi.best, t.act=r$t.act)
return(res)
}

