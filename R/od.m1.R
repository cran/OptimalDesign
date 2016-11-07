od.m1 <-
function (F, b, A=NULL, w0=NULL, type = "exact", kappa=1e-9, tab=NULL, graph=NULL, t.max = 120) 
{
    requireNamespace('Matrix',quietly=TRUE)
 
  if(requireNamespace('gurobi',quietly=TRUE)){
 
    ## Verify F
if (!is.matrix(F) || !is.numeric(F) || !all(is.finite(F))) stop("F must be a matrix of real numbers.")
n <- dim(F)[1]; m <- dim(F)[2]
if (m != 1) stop("This function is intended to solve problems with m=1.") 

print(paste("Reciprocal condition number of F'F is", rcond(t(F) %*% F)), quote = FALSE)

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
if (any(A %*% w0 > b)) stop("w0 must satisfy the constraints A %*% w0 <= b.")

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

## convert w0 to constraints  
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

Fp <- F + matrix(runif(n, min=-kappa, max=kappa), nrow=n)

    start <- as.numeric(proc.time()[3])
    info <- paste("Running m1 for cca", t.max, "seconds")
    info <- paste(info, " starting at ", Sys.time(), ".", sep = "")
    print(info, quote = FALSE)
    model <- list()
    model$obj <- Fp^2
    model$A <- Matrix::Matrix(A, sparse = TRUE)
    model$rhs <- b
    model$sense <- "<="
    if (type == "approximate") {
        model$vtypes <- "C"
    }
    else if (type == "exact") {
        model$vtypes <- "I"
    }
    else {
        stop("Unknown value of parameter type, use 'exact' or 'approximate'.")
    }
    model$modelsense <- "max"
    params <- list(TimeLimit = t.max, OptimalityTol = 1e-09, 
        FeasibilityTol = 1e-09)
    res <- gurobi::gurobi(model)

   status <- res$status
   w.best <- NULL; Phi.best <- 0
   w.temp <- res$x
   if(is.numeric(w.temp) && is.vector(w.temp)){
     w.temp <- pmax(0, w.temp)
     if(length(w.temp) == dim(A)[2]){
     if(sum(A %*% w.temp - b <= 0) == dim(A)[1]){
       w.best <- w.temp[1:n]
       Phi.best <- res$objval
     }
    }
  }


    w.best <- res$x[1:n]
    Phi.best <- res$objval
    t.act <- round(as.numeric(proc.time()[3]) - start, 2)
    info <- paste("m1 finished after", t.act, "seconds at", Sys.time())
    print(info, quote = FALSE)
    
    r <- list(method="Mixed integer linear programming", w.best = w.best, Phi.best = Phi.best, t.act = t.act)

if ((!is.null(graph)) && (r$Phi.best>1e-6)){
  title <- paste("od.m1: n=", n, ", m=", 1, ", Phi=", sep="")
  title <- paste(title, round(r$Phi.best, 6), ", t=", r$t.act, " secs", sep="")
  od.plot(r$w.best, F[, graph], main=title)
}

if (!is.null(tab)) {
        od.print(r$w.best, F[, tab])
    } else {
        od.print(r$w.best)
    }




 }else{
 print("Package gurobi is required to run this algorithm.")
 r <- list(method="Mixed integer linear programming", w.best=NULL, Phi.best=0, status="Gurobi not installed.", t.act=0)
 return(r)

}
}
