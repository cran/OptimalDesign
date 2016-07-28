od.D.IQP <-
function(F, b, A, w1, kappa, t.max){

# Authors: Lenka Filova, Radoslav Harman
# R Library: OptimalDesign
#
# Description:
# An implementation of the "DQ"-method for computing D-efficient exact designs
#   under general linear constraints from the paper:
#
# Harman R, Filova L (2014) "Computing Efficient Exact Designs
# of Experiments Using Integer Quadratic Programming", 
# Computational Statistics & Data Analysis, Vol. 71, pp. 1159-1167
#
# Uses libraries gurobi and Matrix.
#
# Args:
# F ... n times m matrix of regressors corresponding to
#       m>=2 model parameters, and n>=m design points.
# A,b ... the kxn matrix, and the kx1 vector of constraints Aw<=b.
# w1 ... the approximate design used for the criterion of DQ-optimality.
#   If null, w1 is computed using od.D.AA or od.D.SOCP
# graph ... a vector of regressor components to be plotted with the design.
# t.max ... approximate upper limit on the time of computation.
#
# Returns:
# w.best ... best design found within the time limit.
# Phi.best ... the value of the D-criterion of w.best.
# t.act ... the actual time taken by the computation.
#
# Works reasonably well for most problems with n<=..., m<=..., N<=...

start <- as.numeric(proc.time()[3])
info <- paste("Running D.IQP for cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

requireNamespace('Matrix',quietly=TRUE)
  
if(requireNamespace('gurobi',quietly=TRUE)){
n<-dim(F)[1]
m<-dim(F)[2]

Fp <- F + matrix(runif(n * m, min=-kappa, max=kappa), nrow=n)

if (is.vector(A)) A<-t(as.matrix(A))

A<-as.matrix(A)
tone.m <- t(rep(1, m))

if(all(w1==0)){
  if(length(b) == 1 && all(A==rep(1,n))){ 
    w1 <- od.D.AA(F=Fp, N=b, w1=rep(0,n), alg="oom", eff=0.99999, t.max=t.max / 10)$w.best
  } else {
  w1 <- od.D.SOCP(Fp, b, A, type = "approximate", kappa=0, t.max=t.max / 3)$w.best
  }
}

if(is.null(w1)) stop("Approximate design not available.")

M.app <- od.infmat(Fp, w1)
invM.app<-solve(M.app)

# Generate the arguments for the quadratic integer program
f <- rep(0, n)
H <- matrix(0, nrow=n, ncol=n)
for (i in 1:n)
  f[i] <- -2 * t(Fp[i, ]) %*% invM.app %*% Fp[i, ]
for (i in 1:n){
  for (j in i:n){
    H[j, i] <- H[i, j] <- 0.5 * (t(Fp[i, ]) %*% invM.app %*% Fp[j, ])^2
  }
}

model <- list()
model$obj <- f
model$Q <- H
model$A <- Matrix::Matrix(A, sparse=TRUE)
model$rhs <- c(b)

k <- length(b)
model$sense <- rep("<=", k)
model$lb <- rep(0, n)
model$vtypes <- "I"
 
# Use gurobi to compute the DQ-optimal design 
params <- list(TimeLimit=t.max, MIPGap=0)
result <- gurobi::gurobi(model, params)

status <- result$status
w.best <- NULL; Phi.best <- 0
w.temp <- result$x
if (is.numeric(w.temp) && is.vector(w.temp)) {
  if (length(w.temp) >= dim(A)[2]) {
     w.temp <- w.temp[1:dim(A)[2]]
     w.temp <- pmax(0, w.temp)
     if (sum(A %*% w.temp - b <= 1e-4) == dim(A)[1]) {
        w.best <- w.temp
        Phi.best <- od.crit(F, w.best, crit="D")  
     }
  }
}


t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("D.IQP finished after", t.act, "seconds at", Sys.time())
print(info, quote=FALSE)


list(w.best=w.best, Phi.best=Phi.best, status=status, t.act=t.act)
  }else{
 return("Package gurobi is required to run this algorithm")
}
}
