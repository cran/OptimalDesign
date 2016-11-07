\name{od.RCs}
\alias{od.RCs}

\title{Efficient exact size-constrained design using the RC heuristic}

\description{
  Computes an efficient exact design under standard (size) constraints using the RC heuristic.
} 

\usage{
  od.RCs(F, N, w0=NULL, crit="D", R=NULL, w1=NULL, kappa=0, 
         tab=NULL, graph=NULL, t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. The rows of \code{F} represent the \code{m}-dimensional
         regressors corresponding to \code{n} design points. It is assumed that \code{n>=m>=2}. Cf. \link{od.m1} for models
        with \code{1}-dimensional regressors.}
  \item{N}{The required size of the design, i.e., the sum of the components of the resulting design. It is assumed that \code{N>=m}.}
  \item{w0}{The \code{n} times \code{1} vector representing the design to be augmented. Design \code{w0} must satisfy
        \code{sum(w0)<=N-1}. The argument \code{w0} can also be \code{NULL}. In that case \code{w0} is set to the vector 
        of zeros.}
  \item{crit}{The optimality criterion. Possible values are \code{"D", "A", "IV"}.}
  \item{R}{The region of summation for the IV-optimality criterion. The argument \code{R} must be a subvector of \code{1:n},
        or \code{NULL}. If \code{R=NULL}, the procedure uses \code{R=1:n}. Argument \code{R} is ignored if \code{crit="D"},
        or if \code{crit="A"}.}
  \item{w1}{The \code{n} times \code{1} nonnegative vector representing the initial design. Design \code{w1} must satisfy
        \code{w0<=w1} and \code{sum(w1)<=N}. Argument \code{w1} can also be \code{NULL}. In that case the procedure
        sets \code{w1} to be \code{w0}.}
  \item{kappa}{A small non-negative perturbation parameter.}
  \item{tab}{A vector determining the regressor components to be printed with the resulting design.
        This argument should be a subvector of \code{1:n}, or a subvector of \code{colnames(F)}, or it can be \code{NULL}.
        If \code{tab=NULL}, the design is not printed.}
  \item{graph}{A vector determining the regressor components to be plotted with the resulting design.
        This argument should be a subvector of \code{1:n}, or a subvector of \code{colnames(F)}, or it can be \code{NULL}.
        If \code{graph=NULL}, the resulting design is not visualized.}
  \item{t.max}{The time limit for the computation.}
}

\details{
  This is a numerically more efficient implementation of the general algorithm \link{od.RC}
  for the case of the standard (size) constraint, which is a special case of the resource constraints.
  See the reference for details. 
  
  The information matrix of \code{w1} should preferably have the reciprocal condition number of at least \code{1e-5}. Note that the floor of
  an optimal approximate design of size \code{N} (computed using \link{od.AA}) is often a good initial design. Alternatively,
  the initial design can be the result of a different optimal design procedure, such as \link{od.KL}. In any case, the
  model should be non-singular in the sense that there \emph{exists} an exact design \code{w} with a well conditioned information
  matrix. If this requirement is not satisfied, the computation may fail, or it may produce a deficient design. 

  If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
  (see the help page of the function \link{od.crit}) is non-singular, preferably with a reciprocal condition number of at least \code{1e-5}.
  If this requirement is not satisfied, the computation may fail, or it may produce a deficient design.
  
  The perturbation parameter \code{kappa} can be used to add \code{n*m} iid random numbers from the uniform distribution
  in \code{[-kappa,kappa]} to the elements of \code{F} before the optimization is executed. This can be helpful for
  generating a random design from the potentially large set of optimal or nearly-optimal designs. However,
  the RC heuristic uses a tabu principle based on the criterion values of designs, therefore in some problems a nonzero
  \code{kappa} can be detrimental to the optimization process.
  
  The procedure always returns a permissible design, but in some cases, especially if \code{t.max} is too small,
  the resulting design can be inefficient. The performance depends on the problem and on the hardware used, but in most
  cases the function can compute a nearly-optimal exact design for a problem with a hundred design points
  within minutes of computing time. Because this is a heuristic method, we advise the user to verify the quality of the
  resulting design by comparing it to the result of an alternative method (such as \link{od.KL}) and/or by computing its
  efficiency relative to the corresponding optimal approximate design (computed using \link{od.AA}).
}

\value{
  A list with the following components:
   \item{method}{The method used for computing the design \code{w.best}.} 
   \item{w.best}{The best design found.}
   \item{Phi.best}{The value of the criterion of \code{w.best}.}
   \item{t.act}{The actual time taken by the computation.}
}

\references{
  Harman R, Bachrata A, Filova L (2016): Heuristic construction of exact experimental designs under multiple resource constraints, Applied Stochastic Models in Business and Industry, Volume 32, pp. 3-17
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
  \code{\link{od.RC}, \link{od.KL}, \link{od.AA}}
}

\examples{
# We will compute the D-optimal design for the block model with additive 
# effects of blocks and treatments. Suppose that we have 15 blocks of 
# size two and 10 treatments. The problem is equivalent to the problem 
# of D-optimality in the standard model with 45 design points, 15 trials, 
# and the matrix of regressors F.block computed as below 
# (see the reference for details).

k <- 0
T1 <- T2 <- rep(0,45)
F.aux <- matrix(0, nrow=45, ncol=10)
for(i in 1:9){
  for(j in (i + 1):10){
    k <- k + 1
    T1[k] <- i
    T2[k] <- j
    F.aux[k, i] <- 1
    F.aux[k, j] <- (-1)
  }
}
F.block <-  F.aux \%*\% eigen(matrix(1, ncol=10, nrow=10))$vectors[, 2:10]

# Compute the D-optimal exact design of size 15 using the RCs procedure.
res <- od.RCs(F.block, 15, crit = "D", t.max=2.5)
res$Phi.best

# Solve the same problem using the KL procedure to check the obtained 
# criterion value.
od.KL(F.block, 15, crit = "D", t.max=1.5)$Phi.best

# Display the treatments in the 15 blocks of the obtained design.
data.frame(T1=T1[as.logical(res$w.best)], T2=T2[as.logical(res$w.best)])

# Note that the concurrence graph of the design is the Petersen graph.
}
