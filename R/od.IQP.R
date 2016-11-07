od.IQP <-
function(F, b, A=NULL, w0=NULL, crit="D", R=NULL, w1=NULL, kappa=1e-9, tab=NULL, graph=NULL, t.max=120){

if(requireNamespace('gurobi',quietly=TRUE)){

## Verify F
if (!is.matrix(F) || !is.numeric(F) || !all(is.finite(F))) stop("F must be a matrix of real numbers.")
n <- dim(F)[1]; m <- dim(F)[2]
if (n < m) stop("The number dim(F)[1] of design points must be greater or equal to the number dim(F)[2] of parameters of the model.")
if (m < 2) stop("The number dim(F)[2] of parameters must be at least 2. Use the procedure od.m1 for models with a one-dimensional parameter.")

### condition number test
print(paste("Reciprocal condition number of F'F is", rcond(t(F)%*%F)), quote=FALSE)

## Verify b
if (!is.vector(b) || !is.numeric(b) || !all(is.finite(b))) stop("b must be a vector of real numbers.")

## Verify A
if (is.null(A)) {
  A <- t(as.matrix(rep(1, n)))
  if (length(b) != 1) stop("If A is NULL then b must have length 1.")
}
if (!is.matrix(A) || !is.numeric(A) || !all(is.finite(A))) stop("A must be NULL, or a matrix of real numbers.")
if (dim(A)[1] != length(b)) stop("The dimension dim(A)[1] must be equal to the number length(b) of constraints.")
if (dim(A)[2] != n) stop("The dimension dim(A)[2] must be equal to the number dim(F)[1] of design points.")

## Verify w0 
if (is.null(w0)) w0 <- rep(0, n)
if (!is.vector(w0) || !is.numeric(w0) || !all(is.finite(w0)) || (min(w0) < 0)) stop("w0 must be NULL, or a vector of non-negative real numbers.")
if (length(w0) != n) stop("The length of w0 must be equal to the number dim(F)[1] of design points.")


## Verify w1 
if (is.null(w1)) w1<-rep(0,n)             
if (!is.vector(w1) || !is.numeric(w1) || !all(is.finite(w1)) || (min(w1) < 0)) stop("w1 must be NULL, or a vector of non-negative real numbers.")
if (length(w1) != n) stop("The length of w1 must be equal to the number dim(F)[1] of design points.")


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


## Convert w0 to constraints  
augment <- function(w0, A, b){ 
	if (!all(w0==0)) {  	
		A <- rbind(A, -diag(n)) 	
		b <- c(b, -w0) 	
	} 
	return(list(A = A,b = b))
}

aug <- augment(w0, A, b)
A <- aug$A
b <- aug$b

if (crit=="D") r<-od.D.IQP(F=F, b=b, A=A, w1=w1, t.max=t.max, kappa=kappa)
if (crit=="A") r<-od.A.IQP(F=F, b=b, A=A, w1=w1, t.max=t.max, kappa=kappa)
if (crit=="IV") r<-od.A.IQP(F=F.IVtoA(F, R), b=b, A=A, w1=w1, t.max=t.max, kappa=kappa)


r<-list(method="Integer quadratic programming", w.best=r$w.best, Phi.best=r$Phi.best, status=r$status, t.act=r$t.act)

if ((!is.null(graph)) && (r$Phi.best>1e-6)){
  title <- paste("od.IQP: n=", n, ", m=", m, ", crit=", crit ,", Phi=", sep="")
  title <- paste(title, round(r$Phi.best, 6), ", t=", r$t.act, " secs", sep="")
  od.plot(r$w.best, F[, graph], main=title)
}

if (!is.null(tab)) {
        od.print(r$w.best, F[, tab])
    } else {
        od.print(r$w.best)
    }



return(r)
}else{
 print("Package gurobi is required to run this algorithm.")
 r <- list(w.best=NULL, Phi.best=0, status="Gurobi not installed.",       t.act=0)
 return(r)
}
}
