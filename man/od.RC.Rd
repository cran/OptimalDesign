\name{od.RC}
\alias{od.RC}

\title{Efficient exact design using the RC heuristic}

\description{
  Computes an efficient exact design under linear resource constraints using the RC heuristic.
}

\usage{
  od.RC(F, b, A=NULL, w0=NULL, crit="D", R=NULL, w1=NULL, 
        kappa=0, tab=NULL, graph=NULL, t.max=120)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. The rows of \code{F} represent the \code{m}-dimensional
         regressors corresponding to \code{n} design points. It is assumed that \code{n>=m>=2}. Use \link{od.m1} for models
        with 1-dimensional regressors.}
  \item{b, A}{The vector of length \code{k} with positive real components and the \code{k} times \code{n} matrix of non-negative  
         reals numbers. Each column of \code{A} must have at least one strictly positive element. The linear constraints
         \code{A\%*\%w<=b, w0<=w} define the set of permissible designs \code{w} (where \code{w0} is a described below.)
         The argument \code{A} can also be \code{NULL}; in that case \code{b} must be a positive number and \code{A}
        is set to the \code{1} times \code{n} matrix of ones.}
  \item{w0}{The \code{n} times \code{1} vector of nonnegative numbers representing the design to be augmented. Design
         \code{w0} must satisfy \code{A*w0<=b} and it is also required that \code{w0} is not the only permissible exact design.
         The argument \code{w0} can also be \code{NULL}; in that case, \code{w0} is set to the vector of zeros.}
  \item{crit}{The optimality criterion. Possible values are \code{"D", "A", "IV"}.}
  \item{R}{The region of summation for the IV-optimality criterion. The argument \code{R} must be a subvector of \code{1:n},
        or \code{NULL}. If \code{R=NULL}, the procedure uses \code{R=1:n}. Argument \code{R} is ignored if \code{crit="D"},
        or if \code{crit="A"}.}
  \item{w1}{The \code{n} times \code{1} nonnegative vector which represents the initial design. The design \code{w1} must satisfy
        \code{w0<=w1} and \code{A*w1<=b}. The argument \code{w1} can also be \code{NULL}. In that case the procedure sets
        \code{w1} to be \code{w0}.}
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
  This is an implementation of the algorithm proposed by Harman et al. employing the tabu search principle, and related to 
  the DETMAX procedure; see References. The inequalities \code{A\%*\%w<=b},
  \code{w0<=w} with the specific properties mentioned above, form the so-called resource constraints which encompass many
  practical restrictions on the design, always permit a feasible solution, and lead to a bounded set of feasible solutions.

  The information matrix of \code{w1} should preferably have the condition number of at least \code{1e-5}. Note that the floor
  of an optimal approximate design (computed using \link{od.SOCP}) is sometimes a good initial design. Alternatively,
  the initial design can be the result of a different optimal design procedure, such as \link{od.IQP}. Even if no initial design
  is provided, the model should be non-singular in the sense that there \emph{exists} an exact design \code{w} with a well
  conditioned information matrix, satisfying all constraints. If this requirement is not satisfied, the computation may fail,
  or it may produce a deficient design. 

  If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L}
  (see the help page of the function \link{od.crit}) is non-singular, preferably with a condition number of at least \code{1e-5}.
  If this requirement is not satisfied, the computation may fail, or produce a deficient design.
  
  The perturbation parameter \code{kappa} can be used to add \code{n*m} iid random numbers from the uniform distribution
  in \code{[-kappa,kappa]} to the elements of \code{F} before the optimization is executed. This can be helpful for
  generating a random design from the potentially large set of optimal or nearly-optimal designs. However,
  the RC heuristic uses a tabu principle based on the criterion values of designs, therefore in some problems a nonzero
  \code{kappa} can be detrimental to the optimization process.
  
  The procedure always returns a permissible design, but in some cases, especially if \code{t.max} is too small,
  the resulting design can be inefficient. The performance depends on the problem and on the hardware used, but in most
  cases the function can compute a nearly-optimal exact design for a problem with a hundred design points
  within minutes of computing time. Because this is a heuristic method, we advise the user to verify the quality of the
  resulting design by comparing it to the result of an alternative method (such as \link{od.IQP} and \link{od.MISOCP})
  and/or by computing its efficiency relative to the corresponding optimal approximate design (computed using \link{od.SOCP}).
  In the special case of the single constraint on the size, it is generally more efficient to use the functions \link{od.KL}, or
  \link{od.RCs}.
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
  \code{\link{od.IQP}, \link{od.MISOCP}, \link{od.SOCP}, \link{od.KL}, \link{od.RCs}}
}

\examples{
\donttest{
# Consider the spring balance weighing model with 6 items of unknown 
# weight. Suppose that we have already performed one weighing of each 
# item separately. We will compute an A-efficient augmentation designs 
# with 20 additional weighings. Then we will compute A-efficient designs 
# under the additional restriction that no item can be used more 
# than 8 times.

# Create the matrix of regressors of the model and the design 
# to be augmented.
F.sbw <- F.cube(~x1 + x2 + x3 + x4 + x5 + x6 - 1, rep(0, 6), 
                rep(1, 6), rep(2, 6))
w0 <- rep(0, 64); w0[apply(F.sbw, 1, sum)==1] <- 1

# Compute an exact A-efficient augmentation design with 26 
# total weighings.
b.sbw <- 26; A.sbw <- matrix(1, nrow=1, ncol=64)
res <- od.RC(F.sbw, b.sbw, A.sbw, w0=w0, crit="A", tab=1:6, t.max = 30)

# There are many A-optimal designs for this problem, which is possible 
# to see by running the line above several times with a very small 
# non-zero kappa. Note that each of the A-optimal experiments uses each 
# item exactly 11 times. Suppose that we can use each item at 
# most 8 times.

# Create the constraints A.eight * w <= b.eight representing 
# the restriction that no item can be used more than eight times 
# in the entire experiment.
b.eight <- rep(8, 6); A.eight <- t(F.sbw)

# Compute an exact A-efficient design with 26 total weighings under 
# all constraints defined above.
b.sbw <- c(b.eight, 26); A.sbw <- rbind(A.eight, rep(1,64))
res <- od.RC(F.sbw, b.sbw, A.sbw, w0=w0, crit="A", tab=1:6, t.max = 60)

# To find a lower bound on the true efficiency of the resulting design, 
# let us compare it to the A-optimal approximate design.
res.approx <- od.SOCP(F.sbw, b.sbw, A.sbw, w0=w0, crit="A", 
                      tab=1:6, t.max = 20)
res$Phi.best / res.approx$Phi.best
}
}




