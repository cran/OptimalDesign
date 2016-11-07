od.A.AA <-
function(F, N, w1, alg, lambda, eff, t.max){ 
#
# Description:
# Selected algorithms for A-optimal approximate design of experiments
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
# lambda ... a parameter for the multiplicative algorithm.
# graph ... a vector of regressor components to be plotted with the design.
# eff ... the efficiency of a design considered to be perfectly optimal.
# t.max ... approximate upper limit on the time of computation.
#
# Returns:
# w.best ... best design found within the time limit.
# Phi.best ... the value of the A-criterion of w.best.
# eff.best ... the lower bound on the A-efficiency of w.best.
# t.act ... the actual time taken by the computation.
#
# The time elapsed and the efficiency attained is checked after each iteration.
#
# Possible improvements in the future:
#   - Add the Newton step to optimize the support weights. 
#   - Improve initial design heuristic.
#   - Add the deletion of non-supporting points.

stepsize <- function (M, w, v) {
# Computes the A-optimal weight delta for the weight exchange in the sense:
# w[v[1]] <- w[v[1]] + delta
# w[v[2]] <- w[v[2]] - delta
  M.inv <- solve(M)
  dv <- F[v, ] %*% M.inv %*% t(F[v, ]) 
  av <- F[v, ] %*% M.inv %*% M.inv %*% t(F[v, ]) 
  A <- av[1,1] - av[2,2] 
  B <- 2 * dv[1,2] * av[1,2] - dv[2,2] * av[1,1] - dv[1,1] * av[2,2]  
  C <- dv[1,1] - dv[2,2] 
  D <- dv[1,2]^2 - dv[1,1] * dv[2,2]
  G <- B * C - A * D
  k <- v[1]
  l <- v[2]

  if ((abs(G) <= 1e-10) && (abs(B) >= 1e-10)) {
    r0 <- -A / (2*B)  
    if ((-w[k] <= r0) && (r0 <= w[l])) return(r0)  
  }

  if (abs(G) >= 1e-10){  
    H <- B^2 - A * G
    if (H >= 0){
      K <- sqrt(H)
      L <- K / G  
      M <- -B / G
      r1 <- M + L
      r2 <- M - L 
      if ((-w[k] <= r1) && (r1 <= w[l])) return(r1)  
      if ((-w[k] <= r2) && (r2 <= w[l])) return(r2)
    }
  } 

  if (A >= 0) return(w[l]) else return(-w[k])
}

start <- as.numeric(proc.time()[3])
info <- paste("Running A.AA", alg, "for maximum cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

n <- dim(F)[1]        # The size of the design space
m <- dim(F)[2]        # The number of parameters of the model
next.sec <- 0         # Counter of seconds elapsed
n.iter <- 0           # Counter of iterations
del <- 1e-30          # Numerical zero
index <- 1:n
one <- t(rep(1,m))
F <- as.matrix(F)

n.steps <- nchar(alg) # The number of steps in the mix
step.type <- rep("", n.steps)
for(i in 1:n.steps)
  step.type[i] <- substr(alg, i, i)

# Pre-compute elementary information matrices.
A <- array(0, dim=c(n, m, m))
for (i in 1:n)
  A[i, , ] <- F[i, ] %*% t(F[i, ])

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
    a.fun <- apply(G * G, 1, sum) 
    i.best <- which.max(a.fun)
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
M.inv <- solve(M)
G <- F %*% M.inv 
a.fun <- apply(G * G, 1, sum)  
max.a <- max(a.fun)  
ind.max.a <- min(which(a.fun >= max.a - del))
sum.vars <- sum(diag(M.inv))

finish <- (sum.vars / max.a > eff)
while(!finish){
  n.iter <- n.iter + 1

  tm <- as.numeric(proc.time()[3]) - start
  if (tm > next.sec){
     info <- paste("Time:", round(tm, 1), "secs, Actual value:", N * m / sum(diag(M.inv)), ", Eff-bound:",  sum.vars / max.a)
     print(info, quote=FALSE)
     next.sec <- ceiling(tm)
  }

  for(i in 1:n.steps){

    if(step.type[i] == "d"){
      # Perform the vertex direction step.
      a <- -sum(diag(M.inv))
      b <- max.a
      cc<- as.numeric(t(F[ind.max.a,]) %*% M.inv %*% F[ind.max.a,])
      beta <- (sqrt((b - b * cc) / (a * cc + b)) - 1) / cc
      alpha <- beta / (1 + beta)
      w[supp] <- (1 - alpha) * w.supp 
      w[ind.max.a] <- w[ind.max.a] + alpha
      if(!is.element(ind.max.a, supp)){
        supp <- sort(c(supp, ind.max.a))
        n.supp <- n.supp + 1
        F.supp <- F[supp, ]
      }
      w.supp <- w[supp]
      M <- (1 - alpha) * M + alpha * A[ind.max.a, , ]       
    }
 
    if(step.type[i] == "m"){
      # Perform the multiplicative step.
      w.supp <- w.supp * (a.fun[supp]) ^ lambda
      w[supp] <- w.supp <- w.supp / sum(w.supp)
      M <- t((w.supp %*% one) * F.supp) %*% F.supp
    }

    if(step.type[i] == "o"){
      # Perform the weight-ordered support points exchanges with the overall variance-best point
      k <- supp[order(w[supp])]
      lx <- ind.max.a
      Alx <- A[lx, , ]
      for(j in 1:n.supp){
        kx <- k[j]
        v <- c(kx, lx)
        alpha <- stepsize(M, w, v)
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
    M.inv <- solve(M)
    G <- F %*% M.inv 
    a.fun <- apply(G * G, 1, sum)  
    max.a <- max(a.fun)  
    ind.max.a <- min(which(a.fun >= max.a - del))
    sum.vars <- sum(diag(M.inv)) 
  }

  # Stop if the design is essentially perfectly optimal or if the time elapsed.
  if((sum.vars /max.a > eff) || (as.numeric(proc.time()[3]) > start + t.max))  
    finish <- TRUE
}
   
t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("A.AA finished after", t.act, "seconds at", Sys.time())
info <- paste(info, "with", n.iter, "iterations.")
print(info, quote=FALSE)

crit.final <- N * m / sum(diag(M.inv))
eff.final <- sum.vars / max.a

print(paste("A-criterion value:", crit.final), quote=FALSE)
print(paste("Efficiency at least:", eff.final), quote=FALSE)


list(w.best=N * w, Phi.best=crit.final, Eff.best=eff.final, t.act=t.act)
}
