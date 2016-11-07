\name{od.AA}
\alias{od.AA}

\title{Optimal approximate size-constrained design}

\description{
  Computes an optimal approximate design under the standard (size) constraint using a mixture of methods.
}

\usage{
  od.AA(F, N=1, crit="D", R=NULL, w1=NULL, alg="oom", lambda=0.9, 
        tab=NULL, graph=NULL, eff=1-1e-09, t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. The rows of \code{F} represent the \code{m}-dimensional
         regressors corresponding to \code{n} design points. It is assumed that \code{n>=m>=2}. Cf. \link{od.m1} for models
        with \code{1}-dimensional regressors.}
  \item{N}{The required size of the design, i.e., the sum of the components of the resulting design. If \code{N=1} the
        result will be in the standard form of a normalized (probability) measure representing proportions of trials
        in individual design points.}
  \item{crit}{The optimality criterion. Possible values are \code{"D", "A", "IV"}.}
  \item{R}{The region of summation for the IV-optimality criterion. The argument \code{R} must be a subvector of \code{1:n},
        or \code{NULL}. If \code{R=NULL}, the procedure uses \code{R=1:n}. Argument \code{R} is ignored if \code{crit="D"},
        or if \code{crit="A"}.}
  \item{w1}{The real vector of length \code{n} with non-negative components, representing the initial design.
        Before using, \code{w1} is normalized. The information matrix of \code{w1} must be nonsingular.
        The argument \code{w1} can also be \code{NULL}; in that case the function generates a random initial design.
        Set all components of \code{w1} to strictly positive numbers if you use multiplicative steps only.}
  \item{alg}{The methods (\code{"d", "o", "m"}) in the order used within one iteration (see Details).}
  \item{lambda}{A real number in the interval \code{(0,1)} representing the parameter of the multiplicative algorithm
        used for computing A- and IV-optimal designs. Parameter \code{lambda} is ignored if \code{crit="D"}.}
  \item{tab}{A vector determining the regressor components to be printed with the resulting design.
        This argument should be a subvector of \code{1:n}, or a subvector of \code{colnames(F)}, or it can be \code{NULL}.
        If \code{tab=NULL}, the design is not printed.}
  \item{graph}{A vector determining the regressor components to be plotted with the resulting design.
        This argument should be a subvector of \code{1:n}, or a subvector of \code{colnames(F)}, or it can be \code{NULL}.
        If \code{graph=NULL}, the resulting design is not visualized.}
  \item{eff}{The efficiency for the stopping rule in the interval \code{(0,1)}. The algorithm will be stopped if a
        lower bound on the efficiency of the actual design is equal or greater than \code{eff}.}
  \item{t.max}{The time limit for the computation.}
}

\details{
  The function implements three algorithms for the computation of optimal approximate designs:
  the standard vertex-direction method (\code{"d"}), the standard multiplicative method (\code{"m"}), and a novel ordered
  vertex exchange method (\code{"o"}). The methods can be freely combined. For instance, setting \code{method="doom"}
  means that within each iteration the program will execute one vertex direction step, then two ordered
  vertex exchange steps, and then one multiplicative step. The justification of this approach is similar to that of 
  Yu (2011); see the references.
  
  The information matrix of \code{w1} should have the reciprocal condition number of at least \code{1e-5}. Even if no initial design is
  provided, the model should be non-singular in the sense that there \emph{exists} an approximate design \code{w} with a well
  conditioned information matrix. If this requirement is not satisfied, the computation may fail, or it may produce a deficient design. 

  If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
  (see the help page of the function \link{od.crit}) is non-singular, preferably with a reciprocal condition number of at least \code{1e-5}.
  If this requirement is not satisfied, the computation may fail, or it may produce a deficient design.
  
  The performance of the algorithms depends on the problem, on the choice of the parameters, and on the hardware used. However,
  if the method \code{"o"} is used at least once within each iteration, the function tends to reliably compute at least
  \code{99.9\%}-efficient approximate design for a problem with a ten thousand design points within seconds of computing time.
}

\value{
  A list with the following components: 
   \item{method}{The method used for computing the design \code{w.best}.}
   \item{w.best}{The best design found.}
   \item{Phi.best}{The value of the criterion of \code{w.best}.}
   \item{eff.best}{The lower bound on the efficiency of \code{w.best}.}
   \item{t.act}{The actual time taken by the computation.}
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
  \code{\link{od.SOCP}, \link{od.KL}, \link{od.RCs}}
}

\references{
Yu Y (2011): "D-optimal designs via a cocktail algorithm", Statistics and Computing, Vol. 21, No. 4, pp. 475-481
}

\examples{
# Create the regressors matrix for the third-degree trigonometric model 
# on a "partial circle".
F.trig <- F.cube(~I(cos(x1)) + I(sin(x1)) + 
          I(cos(2 * x1)) + I(sin(2 * x1)) + 
          I(cos(3 * x1)) + I(sin(3 * x1)), 
          -pi/2, pi/2, 200)

# Compute the D-, A-, and IV-optimal standardized approximate designs.
res.trig.D <- od.AA(F.trig, 1, alg = "doom", crit = "D", 
              graph=1:7, t.max=4)
res.trig.A <- od.AA(F.trig, 1, alg = "doom", crit = "A", 
              graph=1:7, t.max=4)
res.trig.IV <- od.AA(F.trig, 1, alg = "doom", crit = "IV", 
               graph=1:7, t.max=4)

# Create a model with 2000 regressors generated from the 5 dimensional 
# multivariate normal distribution. (A random linear regression model 
# with 2000 design points and 5 parameters.)
F.norm <- matrix(rnorm(10000), ncol=5)

# Compute the D-, and A-optimal standardized approximate designs. 
res.norm.D <- od.AA(F.norm, 1, alg = "om", crit = "D", 
                    graph=1:5, t.max=4)
res.norm.A <- od.AA(F.norm, 1, alg = "om", crit = "A", 
                    graph=1:5, t.max=4)
}