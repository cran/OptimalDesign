\name{od.IQP}
\alias{od.IQP}

\title{Efficient exact design using integer quadratic programming}

\description{
  Computes an efficient exact design under general linear constraints using the approach of integer quadratic programming.
}

\usage{
  od.IQP(F, b, A=NULL, w0=NULL, crit="D", R=NULL, w1=NULL, 
         kappa=1e-9, tab=NULL, graph=NULL, t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. The rows of \code{F} represent the \code{m}-dimensional
        regressors corresponding to \code{n} design points. It is assumed that \code{n>=m>=2}. Use \link{od.m1} for models
        with 1-dimensional regressors.}
  \item{b, A}{The real vector of length \code{k} and the \code{k} times \code{n} matrix of reals numbers.
         The linear constraints \code{A\%*\%w<=b, w0<=w} define the set of permissible designs \code{w} (where \code{w0} is a
         described below.) The argument \code{A} can also be \code{NULL}; in that case \code{b} must be a positive number
         and \code{A} is set to the \code{1} times \code{n} matrix of ones.}
  \item{w0}{The non-negative vector of length \code{n} representing the design to be augmented. This argument can
         also be \code{NULL}; in that case, \code{w0} is set to the vector of zeros.}
  \item{crit}{The optimality criterion. Possible values are \code{"D", "A", "IV"}.}
  \item{R}{The region of summation for the IV-optimality criterion. The argument \code{R} must be a subvector of \code{1:n},
        or \code{NULL}. If \code{R=NULL}, the procedure uses \code{R=1:n}. Argument \code{R} is ignored if \code{crit="D"},
        or if \code{crit="A"}.}
  \item{w1}{The non-negative vector of length \code{n} representing the optimal or nearly-optimal approximate design
        for the design problem (including the design constraints). The argument \code{w1} can also be \code{NULL}.
        In that case the procedure computes \code{w1} using an appropriate procedure.}
  \item{kappa}{A small non-negative perturbation parameter.}
  \item{tab}{A vector determining the regressor components to be printed with the resulting design.
        This argument should be a subvector of \code{1:n}, or a subvector of \code{colnames(F)}, or it can be \code{NULL}.
        If \code{tab=NULL}, the design is not printed.}
  \item{graph}{A vector determining the regressor components to be plotted with the resulting design.
        This argument should be a subvector of \code{1:n}, or a subvector of \code{colnames(F)}, or it can be \code{NULL}.
        If \code{graph=NULL}, the resulting design is not visualized.}
  \item{t.max}{The time limit for the computation.}
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

\details{
   The procedure computes an efficient exact design by means of integer quadratic programming.
   The idea is to use a quadratic criterion, which approximates the target criterion in the
   neighborhood of the information matrix of the optimal approximate design. See the reference
   for details.
   
   The model should be non-singular in the sense that there exists an exact design \code{w} satisfying the constraints \code{0<=w0<=w} and
   \code{A\%*\%w<=b}, with a non-singular information matrix, preferably with the reciprocal condition number of at least \code{1e-5}. If this requirement
   is not satisfied, the computation may fail, or it may produce a deficient design. 
   
   If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
   (see the help page of the function \code{od.crit}) is non-singular, preferably with a reciprocal condition number of at least \code{1e-5}.
   If this requirement is not satisfied, the computation may fail, or it may produce a deficient design.
   
   The perturbation parameter \code{kappa} can be used to add \code{n*m} iid random numbers from the uniform distribution
   in \code{[-kappa,kappa]} to the elements of \code{F} before the optimization is executed. This can be helpful for
   increasing the numerical stability of the computation or for generating a random design from the potentially large set of optimal or
   nearly-optimal designs.
   
   The performance depends on the problem and on the hardware used, but in most cases the function can compute an optimal or nearly-optimal
   exact design for a problem with a thousand design points within minutes of computing time. Although mostly reliable and fast, there are some
   pathological problems, even small in size, where the resulting design is deficient. We advise the user to verify the
   quality of the resulting design by comparing it to the result of an alternative method (such as \code{od.RC} or \code{od.MISOCP})
   and/or by computing its efficiency relative to the corresponding optimal approximate design (e.g., by using \code{od.SOCP}).
}

\references{
  Harman R., Filova L. (2014): Computing efficient exact designs of experiments using integer quadratic programming, Computational Statistics & Data Analysis, Volume 71, pp. 1159-1167
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
  \code{\link{od.RC}, \link{od.MISOCP}, \link{od.SOCP}}
}

\examples{
if(require("gurobi")){
# Consider the full quadratic model on the domain {-1,-0.9,...,0.9,1}^2, 
# that is 2 factors, each with 11 levels -1,-0.9,...,0.9,1. Suppose that 
# the cost of an observation for the combination (i1,i2) of factors is 
# equal to (i1+1.1)+(i2+1.1) price units. Our time and budget constraints 
# dictate that we cannot perform more than 18 observations and spend more
# than 28 price units. Moreover, we do not wish to perform more than one 
# observation under the same combination of factors. Let us compute 
# an IV-efficient design under these constraints.

# Create the matrix of regressors of the model.
F.quad <- F.cube(~x1 + x2 + I(x1^2) + I(x2^2) + I(x1 * x2), 
          c(-1, -1), c(1, 1), c(11, 11))

# Create the matrix A and the vector b such that our constraints on the 
# feasible design w are equivalent to A * w <= b. 
A.quad <- matrix(0, nrow=123, ncol=121)
for (i in 1:121){
  A.quad[1, i] <- 1
  A.quad[2, i] <- sum(F.quad[i, 2:3]) + 2.2
  A.quad[i + 2, i] <- 1
}
b.quad <- c(18, 28, rep(1, 121))

# Compute an IV-efficient design under the constraints defined above.
resIV <- od.IQP(F.quad, b.quad, A.quad, crit="IV", graph = c("x1","x2"), 
                t.max = 60)

# Verify the quality of the resulting design by computing its efficiency 
# with respect to the IV-optimal approximate design.
resIV.approx <- od.IQP(F.quad, b.quad, A.quad, crit="IV")
resIV$Phi.best / resIV.approx$Phi.best
}
}
