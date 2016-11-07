\name{od.MISOCP}
\alias{od.MISOCP}

\title{Efficient exact design using mixed integer second-order cone programming}

\description{
  Computes an efficient exact experimental design under general linear constraints using the approach of mixed integer second-order cone programming.
}

\usage{
  od.MISOCP(F, b, A=NULL, w0=NULL, crit="D", R=NULL, kappa=1e-9, 
            tab=NULL, graph=NULL, t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. The rows of \code{F} represent the \code{m}-dimensional
         regressors corresponding to \code{n} design points. It is assumed that \code{n>=m>=2}. Use \link{od.m1} for models
         with 1-dimensional regressors. For D-optimality, the current implementation supports the models with \code{m<=10}. }
  \item{b, A}{The real vector of length \code{k} and the \code{k} times \code{n} matrix of reals numbers.
         The linear constraints \code{A\%*\%w<=b, w0<=w} define the set of permissible designs \code{w} (where \code{w0} is a
         described below.) The argument \code{A} can also be \code{NULL}; in that case \code{b} must be a non-negative number
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
   The procedure computes an efficient exact design by converting the optimal design problem to 
   a specific problem of mixed integer second-order cone programming; see the reference for details. The advantage of this
   approach is the possibility to construct exact designs under a general system of linear constraints.
  
   The model should be non-singular in the sense that there exists an exact design \code{w} satisfying the constraints \code{0<=w0<=w} and
   \code{A\%*\%w<=b}, with a non-singular information matrix, preferably with the reciprocal condition number of at least \code{1e-5}. If this requirement
   is not satisfied, the computation may fail, or it may produce a deficient design. 
   
   If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
   (see the help page of the function \link{od.crit}) is non-singular, preferably with a reciprocal condition number of at least \code{1e-5}.
   If this requirement is not satisfied, the computation may fail, or it may produce a deficient design.
   
   The perturbation parameter \code{kappa} can be used to add \code{n*m} iid random numbers from the uniform distribution
   in \code{[-kappa,kappa]} to the elements of \code{F} before the optimization is executed. This can be helpful for
   increasing the numerical stability of the computation or for generating a random design from the potentially large set of optimal or
   nearly-optimal designs.
   
   The performance strongly depends on the problem and on the hardware used, but in most cases the function can compute an optimal or nearly-optimal
   exact design for a problem with a hundred design points within minutes of computing time. We advise the user to verify the quality of the
   resulting design by comparing it to the result of an alternative method (such as \link{od.IQP} and \link{od.RC}) and/or by computing its
   efficiency relative to the corresponding optimal approximate design (computed using \link{od.SOCP}).
   In the special case of the single constraint on the size, it is generally more efficient to use the function \link{od.KL}, or the function
  \link{od.RCs}.
}

\references{
  Sagnol G, Harman R (2015): Computing exact D-optimal designs by mixed integer second-order cone programming. The Annals of Statistics, Volume 43, Number 5, pp. 2198-2224.
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
  \code{\link{od.IQP}, \link{od.RC}, \link{od.SOCP}, \link{od.KL}, \link{od.RCs}}
}

\examples{
if(require("gurobi")){
# Consider a dose-response study where both efficacy and toxicity are 
# observed as 0/1 outcomes for each patient, where the probability of 
# the outcome 1 under the dose x is modeled by the logistic function: 
# exp(ae+be*x)/(1+exp(ae+be*x)) for the efficacy, and 
# exp(at+bt*x)/(1+exp(at+bt*x)) for the toxicity. We can choose the 
# doses x in the range 1 mg to 150 mgs. 
# The aim is to estimate the parameters ae,be using the D-optimal
# design, or using the A-optimal design. 

# Because this is a non-linear model, the optimal designs will depend 
# on the unknown values of the parameters. We will use the approach of 
# local optimality with the following nominal values of the parameters: 
tle <- c(-10, 0.2)
tlt <- c(-20, 0.2)

# It is simple to show that the localized information matrix for 
# (ae,be) is the information matrix of the standard model with 
# the following regressors:
F.logistic <- matrix(0, nrow=150, ncol=2)
for (i in 1:150) 
  F.logistic[i, ] <- 
        c(sqrt(exp(tle[1]+tle[2]*i))/(1+exp(tle[1]+tle[2]*i)), 
        i*sqrt(exp(tle[1]+tle[2]*i))/(1+exp(tle[1]+tle[2]*i)))
                       
# The constraints on the experiment are twofold: We can have at most 
# N=100 subjects and we also require that the expected number of 
# "failed" trials is at most 10. A trial is considered to be a failure 
# if it leads to either a toxic response, or if it is not efficacious. 
# These constraints can be expressed as A*w<=b:
efficacy.prob <- function(x) 
    exp(tle[1]+tle[2]*x)/(1+exp(tle[1]+tle[2]*x))
toxicity.prob <- function(x) 
    exp(tlt[1]+tlt[2]*x)/(1+exp(tlt[1]+tlt[2]*x))
failure.prob <- function(x) 
    1 - (1 - toxicity.prob(x)) * efficacy.prob(x)
b <- c(100, 10); A <- rbind(rep(1,150), failure.prob(1:150))

# Now we can compute the designs:
res.D <- od.MISOCP(F.logistic, b, A, crit="D")
res.A <- od.MISOCP(F.logistic, b, A, crit="A")

# Let us verify the quality of the designs by computing their efficiency 
# relative to the approximate optimal designs:
res.D$Phi.best / od.SOCP(F.logistic, b, A, crit="D")$Phi.best
res.A$Phi.best / od.SOCP(F.logistic, b, A, crit="A")$Phi.best

# We can plot the failure probability curve (red), the toxicity 
# probability curve (black), the efficacy probability curve (green), 
# the D-optimal design (orange) and the A-optimal design (blue):
plot(failure.prob(1:150), type="l", 
     ylab="probability / proportion of subjects", 
     lwd=3, col="red")
lines(toxicity.prob(1:150), col="black")
lines(efficacy.prob(1:150), col="green")
lines(res.D$w.best/100, type="h", col="orange")
lines(res.D$w.best/100, type="h", col="orange", lwd=2)
lines(res.A$w.best/100, type="h", col="blue", lwd=2)

# Note that both designs perform the observations generally at two 
# different levels, one of which is a dose which leads to a 50-percent 
# failure of efficacy. Based on these designs, none of the patients 
# are put on dangerously high levels of doses.
}
}

