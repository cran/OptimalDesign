od.D.KL <-
function(F, N, w1, K, L, variant, kappa, t.max){

# Authors: Radoslav Harman, Lenka Filova
# R Library: OptimalDesign
#
# Description:
# An implementation of the KL exchange algorithm for computing D-efficient
#   exact designs under the standard "size" constraint.
# K is the (upper bound on the) number of "least promising" support points
#   of the current design, for which exchanges are attempted.
# L is the (upper bound on the) number of "most promising" candidate
#   design points for which exchanges are attempted.
# Variant "a" means that the exchanges are performed in an efficient order and
#   any improving exchange is immediately executed.
# Variant "b" means that all permissible exchanges are evaluated and
#   the best exchange (or a random one of the best exchanges) is executed.
# If the algorithm stops in a local optimum before the alloted time ellapsed,
#   the computation is restarted.
# 
# Args:
# F ... nxm model matrix, i.e., n>=2 is the size of the design space
#   and m>=1 is the number of model parameters.
# N ... required size of the design, i.e., the number of trials.
# w1 ... the initial design for the greedy phase. Before the greedy phase,
#   it is replaced by floor(w1), i.e., w1 can be an approximate design.
#   Note: w1 is used only in the first restart. 
# K, L, variant ... parameters of the procedure as explained above.
# t.max ... approximate upper limit on the time of computation.
# graph ... a vector of regressor components to be plotted with the design.
#
# Returns:
# w.best ... best design found within the time limit.
# Phi.best ... the value of the D-criterion of w.best.
# t.act ... the actual time taken by the computation.
#
# Works reasonably well for most problems with n<=1000 and m<=10.
# Consider: w1 can be an integer part of a D-efficient approximate design.
# Possible improvements:
#   - In "a" make the order of search of the KL neighbourhood more efficient.
#   - Improve the heuristics for the choice of K and L.

start <- as.numeric(proc.time()[3])
info <- paste("Running D.KL", variant, "for cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

n <- dim(F)[1]     # The size of the design space
m <- dim(F)[2]     # The number of parameters of the model
p <- ceiling(m/2)  # The size of the random part of the initial design
next.sec <- 0      # Counter of seconds elapsed
n.ex <- 0          # Counter of exchanges
n.rest <- 0        # Counter of restarts
one <- t(rep(1,m))
E <- 10e-9 * diag(m)
Fp <- F + matrix(runif(n * m, min=-kappa, max=kappa), nrow=n)

# If the parameters K,L are not set, use a heuristic assignment
if (is.null(K))
  K <- max(c(10,min(ceiling(sqrt(c(N, n))))))
if (is.null(L))
  L <- max(c(10,min(ceiling(sqrt(n)))))
print(paste("Setting K=", K, ", L=", L, sep=""), quote=FALSE)

if (variant != "a") 
  D0 <- matrix(-1, ncol=L, nrow=K)

# Pre-compute elementary information matrices.
# Note: This approach can be modified in the case of memory problems.
A <- array(0, dim=c(n, m, m))
for (i in 1:n) A[i, , ] <- Fp[i, ] %*% t(Fp[i, ])

finish.all <- FALSE
detM.best <- 0

while (!finish.all){
  n.rest <- n.rest + 1 

  # If not supplied as w1, generate a random starting design with p trials
  if (is.null(w1)) {
    w <- as.vector(table(c(sample(1:n, p, replace=TRUE), 1:n)) - 1)
  } else {
    w <- floor(w1)
    w1 <- NULL  # Start the second run from a random design
  }
  M <- t((w %*% one) * Fp) %*% Fp

  # The forward greedy phase
  sz <- sum(w)
  if (N > sz+1e-9){
    for (k in 1:(N-sz)){
      G <- F %*% solve(M + E)
      d.fun <- apply(G * Fp, 1, sum)
      i.best <- which.max(d.fun)
      w[i.best] <- w[i.best] + 1
      M <- M + A[i.best, , ]
    }
  }

  # Recalculate M afresh to avoid numeric errors
  M <- t((w %*% one) * Fp) %*% Fp
  detM <- det(M)

  if (detM.best < detM){
    w.best <- w
    M.best <- M
    detM.best <- detM
  }

  # The phase of exchanges in w
  finish.all <- finish <- as.numeric(proc.time()[3]) > start + t.max

  while (!finish){
    tm <- as.numeric(proc.time()[3]) - start
    if (tm > next.sec){
       info <- paste("Time:", round(tm, 1), "secs, Best value:", detM.best^(1/m))
       print(info, quote=FALSE)
       next.sec <- ceiling(tm)
    }
    Minv <- solve(M + E)
    G <- Fp %*% Minv
    d.fun <- apply(G * Fp, 1, sum)

    non.supp <- w < 1e-09  # Indices of non-supporting design points
    Kd.fun <- d.fun
    Kd.fun[non.supp] <- Inf
    Kord <- order(Kd.fun)
    Kact <- min(c(n-sum(non.supp), K))
    Kind <- Kord[1:Kact]
    Lact <- min(c(n, L))
    Lind <- order(d.fun)[(n-Lact+1):n]

    if (variant == "a"){

      # Attempt exchanges until improvement is found
      imp <- FALSE
      for (iL in Lact:1){
        for (iK in 1:Kact){
          M.temp <- M + A[Lind[iL], , ] - A[Kind[iK], , ]
          detM.temp <- det(M.temp)
          if(detM.temp > detM + 1e-12){
            w[Lind[iL]] <- w[Lind[iL]] + 1
            w[Kind[iK]] <- w[Kind[iK]] - 1
            M <- M.temp
            detM <- detM.temp
            n.ex <- n.ex + 1
            imp <- TRUE
            break
          }
        }
        if (imp)
          break
      }
    } else {

      # Evaluate the complete KL-neighbourhood     
      Kd.fun <- d.fun[Kind]
      Ld.fun <- d.fun[Lind]
      Cov <- G %*% t(Fp)
      KLd.fun <- Cov[Kind, Lind]

      D <- D0
      D[1:Kact, 1:Lact] <- (1 - Kd.fun) %*% t(1 + Ld.fun) + KLd.fun * KLd.fun

      max.ind <- which(D >= max(D) - 1e-9, arr.ind = TRUE) 
      l.max <- dim(max.ind)[1]
      i.row <- sample(1:l.max, 1)
      iK <- max.ind[i.row, 1]
      iL <- max.ind[i.row, 2]

      imp <- FALSE
      if (D[iK, iL] > 1 + 1e-12){ 
        w[Kind[iK]] <- w[Kind[iK]] - 1
        w[Lind[iL]] <- w[Lind[iL]] + 1
        M <- M - A[Kind[iK], , ] + A[Lind[iL], , ]
        detM <- detM * D[iK, iL]
        n.ex <- n.ex + 1
        imp <- TRUE
      } 
    }

    if (as.numeric(proc.time()[3]) > start + t.max)
      finish.all <- TRUE
    if (finish.all || !imp)
      finish <- TRUE
  }

  # Recalculate M afresh to avoid numeric errors
  M <- t((w %*% one) * Fp) %*% Fp
  detM <- det(M)
  if (detM > detM.best){
    w.best <- w
    M.best <- M
    detM.best <- detM
  }
}

t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("D.KL", variant, " finished after", t.act, "seconds at", Sys.time())
print(info, quote=FALSE)
info <- paste("with", n.rest, "restarts and", n.ex, "exchanges.")
print(info, quote=FALSE)



list(w.best=w.best, Phi.best=od.crit(F,w.best, crit="D"), t.act=t.act)
}
