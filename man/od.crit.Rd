\name{od.crit}
\alias{od.crit}

\title{Optimality criterion}

\description{
  Computes the value of the criterion of D-, A-, or IV-optimality of a given design.
}

\usage{
  od.crit(F, w, crit="D", R=NULL, tol=1e-12)
}

\arguments{
  \item{F}{The \code{n} times \code{m} matrix of real numbers. Rows of \code{F} represent the \code{m}-dimensional
         regressors corresponding to \code{n} design points.}
  \item{w}{The non-negative vector of length \code{n} representing the design.}
  \item{crit}{The optimality criterion. Possible values are \code{"D", "A","IV"}. }
  \item{R}{The region of summation for the IV-optimality criterion. The argument \code{R} must be a subvector of \code{1:n},
        or \code{NULL}. If \code{R=NULL}, the procedure uses \code{R=1:n}. Argument \code{R} is ignored if \code{crit="D"},
        or if \code{crit="A"}.}
  \item{tol}{A small positive number to determine singularity of the information matrix.}
}

\value{
  The value of the concave, positive homogeneous version of the selected real-valued
  criterion applied to the information matrix of the design \code{w} in the linear regression
  model with \code{m}-dimensional regressors \code{F[1,],...,F[n,]} corresponding to \code{n} design points.
}

\details{
  Let \code{w} be a design with information matrix \code{M}, let \code{n} be the number of design points and let \code{m} be the number of
  parameters of the model.

  For \code{w}, the value of the criterion of D-optimality is computed as \code{(det(M))^(1/m)} and the value of the criterion
  of A-optimality is computed as \code{m/trace(M.inv)}, where \code{M.inv} is the inverse of \code{M}.
  
  The IV-optimal design, sometimes called I-optimal or V-optimal, minimizes the integral of the variances
  of the BLUEs of the response surface over a region \code{R}, or the sum of the variances over \code{R}, if \code{R} is finite;
  see Section 10.6 in Atkinson et al. Let the matrix \code{L} be the integral (or the sum) of \code{F[x,]\%*\%t(F[x,])} over \code{x} in \code{R}.
  If the criterion of IV-optimality is selected, the region \code{R} should be chosen such that the associated matrix \code{L} is non-singular.
  Then, let \code{L=t(C)\%*\%C} be the Cholesky decomposition of \code{L}. The design \code{w} is IV-optimal in the model given by \code{F},
  if and only if \code{w} is A-optimal for the model with the regressors matrix \code{F\%*\%C.inv}, where \code{C.inv} is the inverse of \code{C}. 
  
  For the purpose of this package, the value of the IV-criterion for \code{w} is \code{m/trace(N.inv)},
  where \code{N.inv} is the inverse of the information matrix of \code{w} in the model given by regressors matrix \code{F\%*\%C.inv}, and 
  every computational problem of IV-optimality is converted to the corresponding problem of A-optimality. The argument
  \code{R} is assumed to be a subset of \code{1:n}. If the application requires that \code{R} is not a subset of the set of design points,
  the user should compute the matrix \code{C}, transform the model as described above, and use the procedures for A-optimality.

  If the information matrix is singular, the value of all three criteria is zero. An information matrix is considered
  to be singular, if its minimal eigenvalue is smaller than \code{m*tol}. 
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\references{
  Atkinson AC, Donev AN, Tobias RD (2007): Optimum Experimental Designs, with SAS, Oxford University Press, Oxford
}

\seealso{
\code{\link{od.infmat}, \link{od.print}, \link{od.plot}}
}

\examples{
# The matrix of regressors for the spring balance weighing model with 
# 6 weighed items.
F.sbw <- F.cube(~x1 + x2 + x3 + x4 + x5 + x6 - 1, rep(0, 6), 
                rep(1, 6), rep(2, 6))

# The value of all 3 optimality criteria for the design of size 15 
# that weighs each pair of items exactly once.
w2 <- rep(0, 64); w2[apply(F.sbw, 1, sum)==2] <- 1
od.crit(F.sbw, w2, "D")
od.crit(F.sbw, w2, "A")
od.crit(F.sbw, w2, "IV")

# The value of all 3 optimality criteria for the design of size 15 that 
# weighs each quadruple of items exactly once.
w4 <- rep(0, 64); w4[apply(F.sbw, 1, sum)==4] <- 1
od.crit(F.sbw, w4, "D")
od.crit(F.sbw, w4, "A")
od.crit(F.sbw, w4, "IV")
}

