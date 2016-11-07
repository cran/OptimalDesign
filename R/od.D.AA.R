od.D.AA <-
function(F, N, w1, alg, eff, t.max){ 
#
# Description:
# Selected algorithms for D-optimal approximate design of experiments
# under the standard "size" constraint.
#
# Authors: Radoslav Harman, Lenka Filova
# R Library: OptimalDesign
#  
# Args:
# F ... nxm model matrix, where 
#       n>=2 is the size of the design space,
#       m>=2 is the number of model parameters.
# N ... required size of the design, that is, the number of trials.
# w1 ... the initial design. Before using, w1 is normalized.
# alg ... which method(s) will be used (in which order within one iteration).
#     d ...   Standard vertex direction method
#     o ...   Ordered vertex exchange method
#     m ...   Standard multiplicative method
#             Set w1=rep(1/n,n) pls if you use multiplicative steps only.
#     The user can "mix its own cocktail" such as alg="ooom", or alg="doom".
# graph ... a vector of regressor components to be plotted with the design.
# eff ... the efficiency of a design considered to be perfectly optimal.
# t.max ... approximate upper limit on the time of computation.
#
# Returns:
# w.best ... best design found within the time limit.
# Phi.best ... the value of the D-criterion of w.best.
# eff.best ... the lower bound on the D-efficiency of w.best.
# t.act ... the actual time taken by the computation.
#
# The time elapsed and the efficiency attained is checked after each iteration.
#
# Possible improvements in the future:
#   - Add the Newton step to optimize the support weights. 
#   - Improve initial design heuristic.
#   - Add the deletion of non-supporting points.

start <- as.numeric(proc.time()[3])
info <- paste("Running D.AA", alg, "for maximum cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

n <- dim(F)[1]        # The size of the design space
m <- dim(F)[2]        # The number of parameters of the model
next.sec <- 0         # Counter of seconds elapsed
n.iter <- 0           # Counter of iterations
del <- 1e-30          # Numerical zero
index <- 1:n
eff.inv <- 1 / eff 
one <- t(rep(1,m))
F <- as.matrix(F)

n.steps <- nchar(alg) # The number of steps in the mix
step.type <- rep("", n.steps)
for(i in 1:n.steps)
  step.type[i] <- substr(alg, i, i)

# Pre-compute elementary information matrices.
A <- array(0, dim=c(n, m, m))
for (i in 1:n)
  A[i, , ] <- F[i, ]%*%t(F[i, ])

# If w1 is not set, use a greedy initial design.
if (all(w1==0)){
  w <- rep(0, n)
  supp <- sample(1:n,m) 
  w[supp] <- 1
  F.supp <- F[supp, ]
  M <- t((w[supp] %*% one) * F.supp) %*% F.supp
  E <- 1e-6 * diag(m)
  for (k in 1:(3*m)){
    G <- F %*% solve(M + E)
    d.fun <- apply(G * F, 1, sum)
    i.best <- which.max(d.fun)
    w[i.best] <- w[i.best] + 1
    M <- M + A[i.best, , ]
  }
  w <- w / sum(w)
} else {
  w <- w1 / sum(w1)
}

# Compute the basic characteristics of w
supp <- index[w > del]
w.supp <- w[supp]
n.supp <- length(w.supp)
F.supp <- F[supp, ]
M <- t((w.supp %*% one) * F.supp) %*% F.supp
G <- F %*% solve(M)
d.fun <- apply(G * F, 1, sum) / m 
max.d <- max(d.fun)  
ind.max.d <- min(which(d.fun >= max.d - del))

finish <- (max.d < eff.inv)
while(!finish){
  n.iter <- n.iter + 1

  tm <- as.numeric(proc.time()[3]) - start
  if (tm > next.sec){
     info <- paste("Time:", round(tm, 1), "secs, Actual value:", N * det(M)^(1/m), ", Eff-bound:",  1 / max.d)
     print(info, quote=FALSE)
     next.sec <- ceiling(tm)
  }

  for(i in 1:n.steps){

    if(step.type[i] == "d"){
      # Perform the vertex direction step.
      alpha <- (max.d - 1) / (m * max.d - 1)  
      w[supp] <- (1 - alpha) * w.supp 
      w[ind.max.d] <- w[ind.max.d] + alpha
      if(!is.element(ind.max.d, supp)){
        supp <- sort(c(supp, ind.max.d))
        n.supp <- n.supp + 1
        F.supp <- F[supp, ]
      }
      w.supp <- w[supp]
      M <- (1 - alpha) * M + alpha * A[ind.max.d, , ]       
    }
 
    if(step.type[i] == "m"){
      # Perform the multiplicative step.
      w[supp] <- w.supp <- w.supp * d.fun[supp]
      M <- t((w.supp %*% one) * F.supp) %*% F.supp
    }

    if(step.type[i] == "o"){
      # Perform the weight-ordered support points exchanges with the overall variance-best point
      k <- supp[order(w[supp])]
      lx <- ind.max.d
      Alx <- A[lx, , ]
      for(j in 1:n.supp){
        kx <- k[j]
        v <- c(kx, lx)
        cv <- F[v, ] %*% solve(M, t(F[v, ])) 
        alpha <- 0.5 * (cv[1, 1] - cv[2, 2]) / (cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + del) 
        alpha <- max(min(alpha, w[lx]), -w[kx]) 
        w[kx] <- w[kx] + alpha
        w[lx] <- w[lx] - alpha  
        M <- M + alpha * (A[kx, , ] - Alx) 
      }

      supp <- index[w > del]
      w.supp <- w[supp]  
      n.supp <- length(supp)
      F.supp <- F[supp, ]
    }

    # Compute the variance function and related characteristics
    G <- F %*% solve(M)
    d.fun <- apply(G * F, 1, sum) / m 
    max.d <- max(d.fun)  
    ind.max.d <- min(which(d.fun >= max.d - del))
  }

  # Stop if the design is essentially perfectly optimal or if the time elapsed.
  if((max.d < eff.inv) || (as.numeric(proc.time()[3]) > start + t.max))
    finish <- TRUE
}
   
t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("D.AA", alg, " finished after", t.act, "seconds at", Sys.time())
info <- paste(info, "with", n.iter, "iterations.")
print(info, quote=FALSE)

print(paste("D-criterion value:", N * det(M)^(1/m)), quote=FALSE)
print(paste("Efficiency at least:", 1/max.d), quote=FALSE)



list(w.best=N * w, Phi.best=N * det(M)^(1/m), Eff.best=1/max.d, t.act=t.act)
}
