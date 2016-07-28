od.A.RCs <-
function(F, N, w0, w1, kappa, t.max){
#
# Authors: Radoslav Harman, Lenka Filova
# R Library: OptimalDesign
#
# Description:
# An implementation of a special case of the algorithm from
#
# Harman R, Bachrata A, Filova L (2016) "Construction of efficient
# experimental designs under multiple resource constraints",
# Applied Stochastic Models in Business and Industry, Vol. 32, pp. 3-17.
#
# This special case computes A-efficient exact designs
# under the standard "size" constraint.
# 
# Args:
# F ... n times m matrix of regressors corresponding to
#       m>=2 model parameters, and n>=m design points.
# N ... required size of the design (the number of trials of the experiment).
# w0 ... the design to be augmented, i.e., the vector of size n of the
#   minimum numbers of trials in individual design points, w0>=0, sum(w0)<N.
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
# Works reasonably well for most problems with n<=500, m<=10, N<=100.
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
info <- paste("Running A.RCs for cca", t.max, "seconds")
info <- paste(info, " starting at ", Sys.time(), ".", sep="") 
print(info, quote=FALSE)

n.round <- 12    # The number of significant digits used for attributes
eps <- 1e-9     # Used to correct possible round-off errors
V.max <- 500000  # Maximal size of the tabu vector
e.max <- 200     # Maximum length of the excursion
V <- NULL        # Initialization of the tabu vector

n <- dim(F)[1]   # The number of model parameters
m <- dim(F)[2]   # The size of the design space
next.sec <- 0    # Counter of the seconds elapsed from the start
n.evals <- 0     # Counter of the number of criterion evaluations

Fp <- F + matrix(runif(n * m, min=-kappa, max=kappa), nrow=n)

tone.m <- t(rep(1, m))
one.n <- rep(1, n)
E <- eps * diag(m)

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
#   characteristics res and Phi.w that are always updated with w.
# - Number res=N-sum(w) determines the residual number of trials.
# - Number Phi.w is the attribute of w, i.e., rounded criterion value.
    
explore.up <- function(){
# Returns the index of the "best" nonV upper neighbor of w
# If there is no nonV upper neigbor, returns 0.
  wei <- w 
  greed <- -one.n
  shift <- (res - 1) / n * one.n
  if (res > eps){
    for (i in 1:n){  
      wei[i] <- w[i] + 1 
      Phi.wei <- signif(Phi(wei), n.round)
      if (!is.element(Phi.wei, V)){
        greed[i] <- Phi(wei + shift)
      }
      wei[i] <- w[i]
    }
  }
  up.index <- which.max(greed)
  max.greed <- greed[up.index]
  if (max.greed < -eps){
    up.index <- 0
  }
  up.index
}

explore.down <- function(){
# Returns the index of the "best" nonV lower neighbor of w
# If there is no nonV lower neigbor, returns 0.
  wei <- w
  greed <- -one.n
  shift <- (res - 1) / n * one.n
  for (i in 1:n){
    if (w[i] > w0[i]){  
      wei[i] <- w[i] - 1 
      Phi.wei <- signif(Phi(wei), n.round)
      if (!is.element(Phi.wei, V)){
        greed[i] <- Phi(wei + shift)
      }
      wei[i] <- w[i]
    }
  }
  down.index <- which.max(greed)
  max.greed <- greed[down.index]
  if (max.greed < 0){
    down.index <- 0
  }
  down.index
}

random.step <- function(){  # Efficiency of this can be improved.
  success <- FALSE
  while(!success){
    dir <- sample(c(FALSE, TRUE), 1, prob=c(0.2, 0.8))
    if(dir && (res > eps)){
      up.index <- sample(n, 1) 
      w[up.index] <<- w[up.index] + 1
      res <<- res - 1
      Phi.w <<- signif(Phi(w), n.round)
      success <- TRUE
    } else {
      down.index <- sample(n, 1)
      if (w[down.index] > w0[down.index]){
        w[down.index] <<- w[down.index] - 1
        res <<- res + 1
        Phi.w <<- signif(Phi(w), n.round)
        success <- TRUE
      }
    }
  }
}

#  MAIN BODY

finish <- sat <- FALSE
e.length <- 1 
start <- as.numeric(proc.time()[3])

if (is.null(w0))
  w0 <- rep(0, n)

# If w1 is not supplied, use a completely random initialization.
if (is.null(w1)){
  w <- w0
  res <- N - sum(w)
  if(res > eps){
    for(i in 1:res){
      i <- sample(1:n, size=1)  
      w[i] <- w[i] + 1
    }
  }
  res <- 0
} else {
  w <- w1
  res <- N - sum(w) 
}
Phi.w <- signif(Phi(w), n.round)

w.best <- w 
res.best <- res  
Phi.best <- Phi.exact(w.best)

while(!finish){

  tm <- as.numeric(proc.time()[3]) - start
  if (tm > next.sec){
     info <- paste("Time:", round(tm, 1), "secs, Best value:", Phi.best)
     print(info, quote=FALSE)
     next.sec <- ceiling(tm)
  }

  if(!is.element(Phi.w, V)){
    V <- union(Phi.w, V)
    up.index <- explore.up()
    if(up.index){     
      w[up.index] <- w[up.index] + 1
      res <- res - 1
      Phi.w <- signif(Phi(w), n.round)
    } else {
      if (res < eps){
        sat <- TRUE
        Phi.w.cmp <- Phi.exact(w)
        if (Phi.w.cmp > Phi.best){
          w.best <- w
          Phi.best <- Phi.w.cmp
          e.length <- 0
        }
      }
      down.index <- explore.down()
      if(down.index){       
        w[down.index] <- w[down.index] - 1
        res <- res + 1
        Phi.w <- signif(Phi(w), n.round)
      } else {
        random.step()
      }
    } 
  } else {
    down.index <- explore.down()
    if(down.index){
      w[down.index] <- w[down.index] - 1
      res <- res + 1
      Phi.w <- signif(Phi(w), n.round)
    } else {
      up.index <- explore.up()
      if(up.index){
        w[up.index] <- w[up.index] + 1
        res <- res - 1
        Phi.w <- signif(Phi(w), n.round)
      } else {
        random.step()
      }    
    }
  }

  e.length <- e.length + 1
  if (sat && (e.length > e.max)){
    w <- w.best
    res <- 0
    Phi.w <- signif(Phi(w), n.round)
    e.length <- 0
  }

  ptm <- as.numeric(proc.time()[3]) - start       
  if((ptm > t.max) || (length(V) >= V.max))
    finish <- TRUE
}

t.act <- round(as.numeric(proc.time()[3]) - start, 2)
info <- paste("A.RCs finished after", t.act, "seconds at", Sys.time())
print(info, quote=FALSE)
info <- paste("with", n.evals, "criterion evaluations.")
print(info, quote=FALSE)


list(w.best=w.best, Phi.best=Phi.best, t.act=t.act)
}
