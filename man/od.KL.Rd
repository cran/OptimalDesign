\name{od.KL}
\alias{od.KL}

\title{Efficient exact size-constrained design using the KL exchange heuristic}

\description{
  Computes an efficient exact design under the standard (size) constraint using the KL exchange heuristic. 
}

\usage{
  od.KL(F, N, crit="D", R=NULL, w1=NULL, K=NULL, L=NULL, 
        variant="a", kappa=0, tab=NULL, graph=NULL, 
        t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. The rows of \code{F} represent the \code{m}-dimensional
         regressors corresponding to \code{n} design points. It is assumed that \code{n>=m>=2}. Cf. \link{od.m1} for models
        with \code{1}-dimensional regressors.}
  \item{N}{The required size of the design, i.e., the sum of the components of the resulting design. It is assumed that 
         \code{N>=m}.}
  \item{crit}{The optimality criterion. Possible values are \code{"D", "A", "IV"}.}
  \item{R}{The region of summation for the IV-optimality criterion. The argument \code{R} must be a subvector of \code{1:n},
        or \code{NULL}. If \code{R=NULL}, the procedure uses \code{R=1:n}. Argument \code{R} is ignored if \code{crit="D"},
        or if \code{crit="A"}.}
  \item{w1}{The real vector of length \code{n} with non-negative integer components used as an initial design for the
        forward phase of the algorithm. Required: \code{sum(w1)<=N}. The argument \code{w1} can also be \code{NULL};
        in that case the procedure generates a random initial design.}
  \item{K, L}{The natural numbers representing tuning parameters of the procedure; see Details. One or both arguments \code{K,L}
       can also be \code{NULL}; then the function selects its own value of the parameter(s).}
  \item{variant}{The variant of the exchange heuristic; implemented values are \code{"a", "b"}.}
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
  This implementation of the KL algorithm is generally based on the ideas described in Atkinson et al. (2007);
  see the references. 
  
  The information matrix of \code{w1} should preferably have the reciprocal condition number of at least \code{1e-5}. Even if no initial design is
  provided, the model should be non-singular in the sense that there \emph{exists} an exact design \code{w} of size \code{N} with a well
  conditioned information matrix. If this requirement is not satisfied, the computation may fail, or it may produce a deficient design. 

  If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
  (see the help page of the function \link{od.crit}) is non-singular, preferably with a reciprocal condition number of at least \code{1e-5}.
  If this requirement is not satisfied, the computation may fail, or it may produce a deficient design.
  
  The tuning parameter \code{K} is the (upper bound on the) number of "least promising" support points
  of the current design, for which exchanges are attempted. The tuning parameter \code{L} is the (upper
  bound on the) number of "most promising" candidate design points for which exchanges are attempted.

  Variant \code{"a"} means that the exchanges are performed in an efficient order and any improving exchange
  is immediately executed. Variant \code{"b"} means that all permissible exchanges are evaluated and
  the best one of the exchanges (or a random one of the best exchanges) is executed.

  If the algorithm stops in a local optimum before the allotted time elapsed, the computation is restarted.
  The result is the best design found within all restarts.
  
  The perturbation parameter \code{kappa} can be used to add \code{n*m} iid random numbers from the uniform distribution
  in \code{[-kappa,kappa]} to the elements of \code{F} before the optimization is executed. This can be helpful for
  generating a random design from the potentially large set of optimal or nearly-optimal designs.
  
  The performance of the function depends on the problem, on the chosen parameters, and on the hardware used, but in most
  cases the function can compute a nearly-optimal exact design for a problem with a thousand design points within seconds
  of computing time. Because this is only a heuristic, we advise the user to verify the quality of the resulting design by
  comparing it to the result of an alternative method (such as \link{od.RCs}) and/or by computing its efficiency relative
  to the corresponding optimal approximate design (computed using \link{od.AA}).
}

\value{
  A list with the following components: 
   \item{method}{The method used for computing the design \code{w.best}.}
   \item{w.best}{The best design found.}
   \item{Phi.best}{The value of the criterion of \code{w.best}.}
   \item{t.act}{The actual time taken by the computation.}
}

\references{
  Atkinson AC, Donev AN, Tobias RD (2007): Optimum experimental designs, with SAS. Vol. 34. Oxford: Oxford University Press.
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
  \code{\link{od.RCs}, \link{od.AA}}
}

\examples{
# Consider the quadratic Scheffe mixture model with 3 mixture components, 
# each one with permissible levels 0,0.02,...,0.7, i.e., without high 
# proportions of components.
# We will calculate an A-efficient exact design of size 18.

# Compute the regressors for the mixture model without constraints 
# on the region.
# (Note: Here, the constraints are on the design region; it is generally 
# much simpler to  find an optimal design under non-standard constraints 
# on the design region than find an optimal design under non-standard 
# constraints on the design itself.)
F.scheffe <- F.simplex(~x1 + x2 + x3 + I(x1 * x2) + I(x1 * x3) + 
             I(x2 * x3) - 1, 3, 51)

# Remove the trials with high values of the mixture components.
# The resulting design space will have 966 design points.
F.scheffe <- F.scheffe[apply(F.scheffe[, 1:3], 1, max) <= 0.7,]

# Compute an A-efficient exact designs with 18 observations.
res.exact <- od.KL(F.scheffe, 18, crit = "A", tab=1:3, 
                   graph=1:3, t.max=4)

# Verify the quality of the resulting design by computing its efficiency 
# relative to the A-optimal approximate size-constrained design.
res.approx <- od.AA(F.scheffe, 18, crit = "A", eff=1-1e-9)
res.exact$Phi.best / res.approx$Phi.best
}

