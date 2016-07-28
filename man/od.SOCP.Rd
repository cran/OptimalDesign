\name{od.SOCP}
\alias{od.SOCP}

\title{Optimal approximate design using second-order cone programming}

\description{
  Computes an efficient approximate experimental design under general linear constraints using the approach of second-order cone programming.
}

\usage{
  od.SOCP(F, b, A=NULL, w0=NULL, crit="D", R=NULL, kappa=1e-9, 
          tab=NULL, graph=NULL, t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. The rows of \code{F} represent the \code{m}-dimensional
         regressors corresponding to \code{n} design points. It is assumed that \code{n>=m>=2}. Use \link{od.m1} for models
         with 1-dimensional regressors. For D-optimality, the current implementation supports the models with \code{m<=10}. }
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
   The procedure computes an efficient approximate design by converting the optimal design problem to 
   a specific problem of second-order cone programming; see the reference for details. The advantage of this
   approach is the possibility to construct approximate designs under a general system of linear constraints. In particular,
   this function provides means of computing informative lower bounds on the efficiency of the linearly constrained
   exact designs computed using methods such as \link{od.RC}, \link{od.IQP}, and \link{od.MISOCP}.
  
   The model should be non-singular in the sense that there exists an approximate design \code{w} satisfying the constraints \code{0<=w0<=w} and
   \code{A\%*\%w<=b}, with a non-singular information matrix, preferably with the condition number of at least \code{1e-5}. If this requirement
   is not satisfied, the computation may fail, or it may produce a deficient design. 
   
   If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
   (see the help page of the function \link{od.crit}) is non-singular, preferably with a condition number of at least \code{1e-5}.
   If this requirement is not satisfied, the computation may fail, or it may produce a deficient design.
   
   The perturbation parameter \code{kappa} can be used to add \code{n*m} iid random numbers from the uniform distribution
   in \code{[-kappa,kappa]} to the elements of \code{F} before the optimization is executed. This can be helpful for
   increasing the numerical stability of the computation or for generating a random design from the potentially large set of optimal or
   nearly-optimal designs.
   
   The performance depends on the problem and on the hardware used, but in most cases the function can compute an optimal or nearly-optimal
   exact design for a problem with a thousand design points within minutes of computing time. If the only constraint on the design is the
   standard constraint on the size, the function \link{od.AA} should be a preferred choice.
}

\references{
  Sagnol G, Harman R (2015): Computing exact D-optimal designs by mixed integer second-order cone programming. The Annals of Statistics, Volume 43, Number 5, pp. 2198-2224.
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
  \code{\link{od.AA}, \link{od.MISOCP}, \link{od.IQP}, \link{od.RC}}
}

\examples{
if(require("gurobi")){
# Suppose that we model the mean values of the observations of 
# a circadian rhythm by a third-degree trigonometric model on 
# a discretization of the interval [0, 2*pi]. We would like to 
# construct a D-efficient design. 
# However, the distance of successive times of observations should not be 
# smaller than the (1/72)-th of the interval (20 minutes, if the interval 
# represents one day). Also, we cannot perform more than 8 observations
# in total, and more than 4 observations in the interval [2*pi/3, 2*pi] 
# (i.e., during the 16 non-working hours).

# Create the matrix of regressors. 
F.trig <- F.cube(~I(cos(x1)) + I(sin(x1)) + 
          I(cos(2 * x1)) + I(sin(2 * x1)) +
          I(cos(3 * x1)) + I(sin(3 * x1)), 
          1 / 144, 2 * pi, 288) 

# Create the constraints.
b.trig <- c(rep(1, 285), 12, 4)
A.trig <- matrix(0, nrow=287, ncol=288)
for(i in 1:285) A.trig[i, i:(i+3)] <- 1  
A.trig[286,] <- 1; A.trig[287, 97:288] <- 1

# Compute the D-optimal approximate design under the constraints.
res.trig <- od.SOCP(F.trig, b.trig, A.trig, crit="D")

# Inspect the resulting approximate design.
od.plot(res.trig$w.best)
od.print(round(res.trig$w.best,2))

# It is clear that a very efficient exact design of size 8 satisfying 
# the constraints performs the observations in design points 
# 1, 34, 63, 96, 134, 173, 212, 251, i.e.
w.exact <- rep(0, 288)
w.exact[c(1, 34, 63, 96, 134, 173, 212, 251)] <- 1

# Indeed, the efficiency of this exact design relative to the optimal 
# approximate design is:
od.crit(F.trig, w.exact) / od.crit(F.trig, res.trig$w.best)

# Of course, it is also possible to directly use an exact-design function 
# such as od.MISOCP for this problem, or it is possible to use 
# the optimal approximate design to decrease the support size of 
# the candidate design points set, and then use an exact-design 
# procedure.

# See also the examples in the help files of functions od.RC, od.IQP 
# and od.MISOCP, where od.SOCP is used to compute informative lower 
# bounds on the efficiencies of exact designs.
}
}
