od.A.KL <-
function(F, N, w1, K, L, variant, kappa, t.max){

# Authors: Radoslav Harman, Lenka Filova
# R Library: OptimalDesign
#
# Description:
# An implementation of the KL exchange algorithm for computing A-efficient
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
# Phi.best ... the value of the A-criterion of w.best.
# t.act ... the actual time taken by the computation.
#
# Works reasonably well for most problems with n<=1000 and m<=10.
# Consider: w1 can be an integer part of an A-efficient approximate design.
# Possible improvements:
#   - In "a" make the order of search of the KL neighbourhood more efficient.
#   - Improve the heuristic for the choice of K and L.

start <- as.numeric(proc.time()[3])
info <- paste("Running A.KL", variant, "for cca", t.max, "seconds")
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
  T0 <- matrix(-1, ncol=L, nrow=K)

# Pre-compute elementary information matrices.
# Note: This approach can be changed in the case of memory problems.
A <- array(0, dim=c(n, m, m))
for (i in 1:n) A[i, , ] <- Fp[i, ] %*% t(Fp[i, ])

finish.all <- FALSE
trMinv.best <- Inf

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
      G <- Fp %*% solve(M + E)
      h.fun <- apply(G * G, 1, sum) / (1 + apply(G * Fp, 1, sum)) 
      i.best <- which.max(h.fun)
      w[i.best] <- w[i.best] + 1
      M <- M + A[i.best, , ]
    }
  }

  # Recalculate M afresh to avoid numeric errors
  M <- t((w %*% one) * Fp) %*% Fp
  trMinv <- sum(diag(solve(M + E)))

  if (trMinv.best > trMinv){
    w.best <- w
    M.best <- M
    trMinv.best <- trMinv
  }

  # The phase of exchanges in w
  finish.all <- finish <- as.numeric(proc.time()[3]) > start + t.max

  while (!finish){
    tm <- as.numeric(proc.time()[3]) - start

    if(det(M.best) < 1e-12){
      Phi.actual <- 0
      } else {
      Phi.actual <- m / sum(diag(solve(M.best)))
    }

    if (tm > next.sec){
       info <- paste("Time:", round(tm, 1), "secs, Best value:", Phi.actual)
       print(info, quote=FALSE)
       next.sec <- ceiling(tm)
    }
    Minv <- solve(M + E)
    G <- Fp %*% Minv
    d.fun <- apply(G * Fp, 1, sum)
    d2.fun <- apply(G * G, 1, sum) 
    h.fun <-  d2.fun / (1 + d.fun) 

    non.supp <- w < 1e-09  # Indices of non-supporting design points
    Kh.fun <- h.fun
    Kh.fun[non.supp] <- Inf
    Kord <- order(Kh.fun)
    Kact <- min(c(n-sum(non.supp), K))
    Kind <- Kord[1:Kact]
    Lact <- min(c(n, L))
    Lind <- order(h.fun)[(n-Lact+1):n]

    if (variant == "a"){

      # Attempt exchanges until improvement is found
      imp <- FALSE
      for (iL in Lact:1){
        for (iK in 1:Kact){
          M.temp <- M + A[Lind[iL], , ] - A[Kind[iK], , ]
          trMinv.temp <- sum(diag(solve(M.temp + E)))
          if(trMinv.temp < trMinv - 1e-12){
            w[Lind[iL]] <- w[Lind[iL]] + 1
            w[Kind[iK]] <- w[Kind[iK]] - 1
            M <- M.temp
            trMinv <- trMinv.temp
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
      Kd2.fun <- d2.fun[Kind]
      Ld2.fun <- d2.fun[Lind]
      Cov <- G %*% t(Fp)
      Cov2 <- G %*% t(G)
      KLd.fun <- Cov[Kind, Lind]
      KLd2.fun <- Cov2[Kind, Lind]

      T <- T0
      # The following formulas can be obtained from the Woodbury identity.
      T[1:Kact, 1:Lact] <- (1 - Kd.fun) %*% t(Ld2.fun) - Kd2.fun %*% t(1 + Ld.fun) + 2 * KLd.fun * KLd2.fun
      T[1:Kact, 1:Lact] <- T[1:Kact, 1:Lact] / ((1 - Kd.fun) %*% t(1 + Ld.fun) + KLd.fun * KLd.fun)

      max.ind <- which(T >= max(T) - 1e-9, arr.ind = TRUE) 
      l.max <- dim(max.ind)[1]
      i.row <- sample(1:l.max, 1)
      iK <- max.ind[i.row, 1]
      iL <- max.ind[i.row, 2]

      imp <- FALSE
      if (T[iK, iL] > 1e-12){ 
        w[Kind[iK]] <- w[Kind[iK]] - 1
        w[Lind[iL]] <- w[Lind[iL]] + 1
        M <- M - A[Kind[iK], , ] + A[Lind[iL], , ]
        trMinv <- trMinv - T[iK, iL]
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
  trMinv <- sum(diag(solve(M + E)))
  if (trMinv < trMinv.best){
    w.best <- w
    M.best <- M
    trMinv.best <- trMinv
  }
}

t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("A.KL finished after", t.act, "seconds at", Sys.time())
print(info, quote=FALSE)
info <- paste("with", n.rest, "restarts and", n.ex, "exchanges.")
print(info, quote=FALSE)

Phi.final <- od.crit(F, w.best, crit="A")


list(w.best=w.best, Phi.best=Phi.final, t.act=t.act)
}
