\name{od.m1}
\alias{od.m1}

\title{Optimal design for models with a one-dimensional parameter}

\description{
  Computes an optimal approximate or exact design under general linear constraints for a model with a one-dimensional parameter. 
}

\usage{
  od.m1(F, b, A=NULL, w0=NULL, type="exact", kappa=1e-9, 
        tab=NULL, graph=NULL, t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{1} matrix of regressors corresponding to \code{n} design points
          and \code{1} model parameter.}
  \item{b, A}{The real vector of length \code{k} and the \code{k} times \code{n} matrix of reals numbers.
         The linear constraints \code{A\%*\%w<=b, w0<=w} define the set of permissible designs \code{w} (where \code{w0} is a
         described below.) The argument \code{A} can also be \code{NULL}; in that case \code{b} must be a positive number
         and \code{A} is set to the \code{1} times \code{n} matrix of ones.}
  \item{w0}{The non-negative vector of length \code{n} representing the design to be augmented. This argument can
         also be \code{NULL}; in that case, \code{w0} is set to the vector of zeros.}
  \item{type}{Specifies whether exact or approximate design is to be computed.}
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
   For the case of \code{m=1} parameter, the problem of optimal design is much simpler than for \code{m>=2}. First, for \code{m=1}
   all standardized information criteria coincide, therefore the D-, A-, and IV-optimal designs are the same,
   and we can call them just "optimal". Second, the information matrix is a real number, therefore we can call it just "information".
   Under the most common size constraint, the optimal design is any design supported on the set of maxima of \code{F[1,1]^2,...,F[n,1]^2},
   which is straightforward to find. Under a non-standard linear constraint, however, the problem becomes a less trivial knapsack problem,
   which is here solved by the integer linear programming solver of gurobi.
  
   The model should be non-singular in the sense that there exists an exact design \code{w} satisfying the constraints \code{0<=w0<=w} and
   \code{A\%*\%w<=b}, with a nonzero information. If this requirement is not satisfied, the computation may fail, or produce a deficient design. 
   
   If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
   (in this case a real number; see the help page of the function \link{od.crit}) is non-zero. If this requirement is not satisfied, the
   computation may fail, or it may produce a deficient design.
   
   The perturbation parameter \code{kappa} can be used to add \code{n*m} iid random numbers from the uniform distribution
   in \code{[-kappa,kappa]} to the elements of \code{F} before the optimization is executed. This can be helpful for
   increasing the numerical stability of the computation or for generating a random design from the potentially large set of optimal or
   nearly-optimal designs.
   
   The performance depends on the problem and on the hardware used, but in most cases the function can compute an optimal exact design for
   a problem with a thousand design points within seconds of computing time.
}

\value{
  A list with the following components: 
  \item{method}{The method used for computing the design \code{w.best}.}
  \item{w.best}{the best permissible design found, or \code{NULL}. The value of \code{w.best} will be \code{NULL} if the computation fails. 
      This can happen, if no permissible solution is found within the time limit, no permissible solution exists, or the problem is unbounded; see the
      \code{status} variable for more details. Note that even if \code{w.best} is a permissible design, then it still can have a singular information
      matrix; cf. the \code{Phi.best} variable.}
  \item{Phi.best}{The value of the criterion of optimality of the design \code{w.best}. If \code{w.best} has a singular information matrix or if the
      computation fails, the value of \code{Phi.best} will be \code{0}.}
  \item{status}{The status variable of the gurobi optimization procedure; see the gurobi solver documentation for details.}
  \item{t.act}{The actual time taken by the computation.}
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\examples{
if(require("gurobi")){
# We will demonstrate the procedure on a simple randomly generated 
# knapsack problem. Here, the squares of elements of F correspond 
# to the values of n available items, the elements of the 1 times n 
# matrix A correspond to the weights of the n items, and the real 
# number b is the upper limit on the total weight of the items that 
# can  be put into the knapsack. 
# The resulting binary "optimal design" determines which of the items 
# should we put into the knapsack to steal the highest possible value.

n <- 200 # There are this many items to choose from.
F.square <- matrix(sample(1:10, n, replace=TRUE), ncol=1) 
# Generate random prices of items.
A <- matrix(sample(1:10, n, replace=TRUE), nrow=1) 
# Generate random weights of items in kgs.
A <- rbind(A, diag(n))  
# We assume that there is just one copy of each item.
b <- c(n / 4, rep(1,n)) 
# The capacity of the knapsack is n/4 kgs.

# Compute the optimal design, which in this case determines how many 
# (0 or 1) of each of the n items should we put into the knapsack.
od.m1(sqrt(F.square), b, A)

# Note: one can compare the result with a specialized function 
# as follows:
# library(adagio); knapsack(A[1,], F.square[,1], n / 4)

# However, od.m1 is more general than the standard knapsack functions.
# Suppose, for instance, that the uncle asks that we must be sure to 
# take the items number 1, 13 and 66. We will compute the most valuable 
# selection of items that fit into our knapsack and contain the 
# required items.
w0 <- rep(0, n); w0[c(1, 13, 66)] <- 1
od.m1(sqrt(F.square), b, A, w0, t.max=2)
}
}

