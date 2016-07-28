od.A.RC <-
function(F, b, A, w0, w1, kappa, t.max){
#
# Authors: Radoslav Harman, Lenka Filova
# R Library: OptimalDesign
#
# Description:
# An implementation the algorithm for computing A-efficient exact designs from
#
# Harman R, Bachrata A, Filova L (2016) "Construction of efficient
# experimental designs under multiple resource constraints",
# Applied Stochastic Models in Business and Industry, Vol. 32, pp. 3-17.
#
# Args:
# F ... n times m matrix of regressors corresponding to
#       m>=2 model parameters, and n>=m design points.
# A,b ... the kxn matrix, and the kx1 vector of constraints Aw<=b.
#   A must have non-negative elements, b must have positive elements.
#   Each column of A must have at least one strictly positive element.
# w0 ... the nx1 vector representing the design to be augmented.
#   Design w0 must satisfy A*w0<=b.
# w1 ... the initial design, i.e., the vector of size n of the
#   numbers of trials in individual design points.
#   If w1=NULL, the algorithm uses a random initial design.
#   Note that the floor of an optimal approximate design is often a good start.
# graph ... a vector of regressor components to be plotted with the design.
# t.max ... approximate upper limit on the time of computation.
#
# Returns:
# w.best ... the best design found within the time limit.
# Phi.best ... the value of the A-criterion of w.best.
# t.act ... the actual time taken by the computation.
#
# Works reasonably well for most problems with n<=200, m<=10, N<=100.
# For difficult problems set t.max to be more than n.
#
# Possible improvements:
#   - Speed-up the attribute computation method (e.g. a linear one).
#   - Choose a faster greed-evaluation heuristic (e.g. a direct one).
#   - Allow excursions outside the feasible set.
#   - Create a more efficient restart heuristic. 
#   - Make a parallel/C++ implementation.
#   - Improve the efficiency of the random step.

start <- as.numeric(proc.time()[3])
info <- paste("Running A.RC for cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

n.round <- 12    # The number of significant digits used for attributes
eps <- 1e-9      # Used to correct possible round-off errors
V.max <- 500000  # Maximal size of the tabu vector
back.max<-20     # Maximal number of backward steps in one excursion
V <- NULL        # Initialization of the tabu vector

if (is.vector(A)) A<-t(as.matrix(A))


n <- dim(F)[1]   # The number of model parameters
m <- dim(F)[2]   # The size of the design space
next.sec <- 0    # Counter of the seconds elapsed from the start
n.evals <- 0     # Counter of the number of criterion evaluations

Fp <- F + matrix(runif(n * m, min=-kappa, max=kappa), nrow=n)

tone.m <- t(rep(1, m))
one.n <- rep(1, n)
tone.n <- t(one.n)
E <- eps * diag(m)
b <- b + eps

Phi <- function(w){ 
  n.evals <<- n.evals + 1
  m / sum(diag(solve(t((w %*% tone.m) * Fp) %*% Fp + E)))
}


Phi.exact <- function(w){ 
  n.evals <<- n.evals + 1
  M <- t((w %*% tone.m) * F) %*% F
  if(rcond(M) > 1e-15){
    val <- max(c(0, m / sum(diag(solve(M)))))
  } else {
    val <- 0
  }
  val
}

########################################################################
# - Vector w represents the current design. The design has its
#   characteristics res, d and attrb that are always updated with w.
# - Vector res=b-A*w determines the residual amounts of the k resources.
# - Vector d represents numbers of possible additional trials in pts of w.
# - attrb is the "characteristic attribute" of the design w to build V.
# - val is the "local heuristic evaluation" of the design w.

compute.d <- function(res.act){  
  kr <- (res.act %*% tone.n) / A
  kr[is.na(kr)] <- Inf
  floor(apply(kr, 2, min) + eps)
}

compute.gamma<-function(numer,denom){
# Computes the minimum of the vector numer/denom from all components where
# denom is not zero. If all denoms are zero, the result is zero.
  if (!any(denom > eps)){
    gamma <- 0
  } else {
    ratio <- numer / denom
    gamma <- min(ratio[denom > eps])
  }
  gamma
}

explore.up <- function(){
# Returns the index of the "best" nonV upper neighbor of w
# If there is no nonV upper neigbor, returns 0.
  wei <- w 
  val <- -one.n                  
  for (i in 1:n){
    if (d[i] > 0){  
      wei[i] <- w[i] + 1 
      attrei <- signif(Phi(wei), n.round) 
      if (!is.element(attrei, V)){  
        dei <- compute.d(res - A[, i])
        gammaei <- compute.gamma(res - A[, i], A %*% dei)
        val[i] <- Phi(wei + gammaei * dei)            
      }
      wei[i] <- w[i]
    }
  }
  up.index <- which.max(val)
  max.val <- val[up.index] 
  if (max.val < -eps){            
    up.index<-0
  }
  up.index
}

explore.down <- function(){         
# Returns the index of the "best" nonV lower neighbor of w
# If there is no nonV lower neigbor, returns 0.
  wei <- w
  val <- -one.n
  for (i in 1:n){
    if (w[i] > w0[i]){  
      wei[i] <- w[i] - 1 
      attrei <- signif(Phi(wei), n.round)
      if (!is.element(attrei, V)){
        dei <- compute.d(res + A[, i])
        gammaei <- compute.gamma(res + A[, i], A %*% dei)   
        val[i] <- Phi(wei + gammaei * dei)
      }
      wei[i] <- w[i]
    }
  }
  down.index <- which.max(val)
  max.val <- val[down.index]
  if (max.val < (-eps)){
    down.index <- 0
  }
  down.index
}

random.step<-function(){ # Efficiency of this can be improved.
  success <- FALSE
  while(!success){
    dir <- sample(c(FALSE, TRUE), 1, prob=c(0.2, 0.8))
    if(dir){
      up.index <- sample(n, 1) 
      if (d[up.index] > 0){
        w[up.index] <<- w[up.index] + 1
        res <<- res - A[, up.index]
        d <<- compute.d(res)
        attrb <<- signif(Phi(w), n.round)
        success <- TRUE
      }
    } else {
      down.index <- sample(n, 1)
      if (w[down.index] > w0[down.index]){
        w[down.index] <<- w[down.index] - 1
        res <<- res + A[,down.index]
        d <<- compute.d(res)
        attrb <<- signif(Phi(w), n.round)
        success <- TRUE
        back.no <<- back.no + 1
      }
    }
  }
}

#  MAIN BODY

finish <- sat <- FALSE
back.no <- 0
start <- as.numeric(proc.time()[3])

if (is.null(w0))
  w0 <- rep(0, n)

# If w1 is null, use a random initialization.
if (is.null(w1)){
  w <- w0
  res <- b - A %*% w
  d <- compute.d(res)    
  maximal <- FALSE
  while(!maximal){
    val <- (d > 0)
    if(sum(val) == 0){
      maximal <- TRUE           
    } else {
      i <- sample(1:n, 1, prob=val)  
      w[i] <- w[i]+1
      res <- b - A %*% w
      d <- compute.d(res)
    }
  }
} else {
  w <- w1
  res <- b - A %*% w
  d <- compute.d(res)
}
attrb <- signif(Phi(w), n.round)

w.best <- w
res.best <- res
d.best <- d
attr.best <- attrb
Phi.best <- Phi.exact(w.best)

while(!finish){ 

  tm <- as.numeric(proc.time()[3]) - start
  if (tm > next.sec){
     info <- paste("Time:", round(tm, 1), "secs, Best value:", m / sum(diag(solve(t((w.best %*% tone.m) * F) %*% F + E))) )
     print(info, quote=FALSE)
     next.sec <- ceiling(tm)
  }

  if(!is.element(attrb, V)){
    V <- union(attrb, V)
    up.index <- explore.up()
    if(up.index){           
      w[up.index] <- w[up.index] + 1
      res <- res - A[,up.index]
      d <- compute.d(res)
      attrb <- signif(Phi(w), n.round)
    } else {
      if (!any(as.logical(d))){   # Zefektivnit  
        Phi.w.cmp <- Phi.exact(w)
        if (Phi.w.cmp > Phi.best){
          w.best <- w
          res.best <- res
          d.best <- d
          attr.best <- attrb
          Phi.best <- Phi.w.cmp
          back.no <- 0
        }
      }
      down.index <- explore.down()
      if(down.index){      
        w[down.index] <- w[down.index] - 1
        res <- res + A[, down.index]
        d <- compute.d(res)
        attrb <- signif(Phi(w),n.round)
        back.no <- back.no + 1
      } else {
        random.step()
      }
    }   
  } else {
    down.index <- explore.down()
    if(down.index){
      w[down.index] <- w[down.index] - 1
      res <- res + A[, down.index]
      d <- compute.d(res)
      attrb <- signif(Phi(w), n.round)
      back.no <- back.no + 1
    } else {
      up.index <- explore.up()
      if(up.index){
        w[up.index] <- w[up.index] + 1
        res <- res - A[, up.index]
        d <- compute.d(res)
        attrb <- signif(Phi(w), n.round)
      } else {
        random.step()
      }
    }
  }
   
  if (back.no>back.max){
    w <- w.best
    res <- res.best
    d <- d.best
    attrb <- attr.best
    back.no <- 0
  }
       
  ptm <- as.numeric(proc.time()[3]) - start       
  if((ptm > t.max) || (length(V) >= V.max))
    finish <- TRUE
}

t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("A.RC finished after", t.act, "seconds at", Sys.time())
print(info, quote=FALSE)
info <- paste("with", n.evals, "criterion evaluations.")
print(info, quote=FALSE)


list(w.best=w.best, Phi.best=Phi.best, t.act=t.act)
}
