\name{F.cube}
\alias{F.cube}

\title{Model matrix on a cube}

\description{
  Creates the matrix of regressors for a model on a \code{d}-dimensional rectangular grid of points,
  i.e., for a factorial model with \code{d} quantitative factors.
}

\usage{
  F.cube(formula, lower, upper, n.levels)
}

\arguments{
  \item{formula}{The formula of the model.}
  \item{lower}{The real vector of length \code{d} determining the smallest values of factors.}
  \item{upper}{The real vector of length \code{d} determining the largest values of factors.}
  \item{n.levels}{The vector of length \code{d} consisting of the numbers greater than or equal to 2, determining the numbers of levels of each factor.}
}

\value{
  The \code{n} times \code{m} matrix of regressors with the \code{m}-dimensional rows corresponding to 
  the regressors of the \code{n} design points (\code{n} is equal to \code{prod(n.levels)}).
}

\details{
   The rules for creating the model formula are the same as in the \code{lm} function but:
        1) the formula must not contain the dependent variable;
        2) the \code{d} factors (independent variables) must be labeled \code{x1, x2}, ...
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
  \code{\link{F.simplex}}
}

\examples{
# The matrix of regressors for the trigonometric model of the second 
# degree on the "full circle" discretized to 24 equidistant points.
F.cube(~I(cos(x1)) + I(sin(x1)) + I(cos(2*x1)) + I(sin(2*x1)), 
       0, (46 / 24) * pi, 24)

# The matrix of regressors of the cubic model on a discretized 
# interval [-1,1].
F.cube(~x1 + I(x1^2) + I(x1^3), -1, 1, 21)

# The matrix of regressors of the 3-factor model with the intercept, 
# the linear and the second-order interaction terms; 
# the first 2 factors have levels -1, 0, 1, and the third factor has 
# levels -1, 0.5, 0, 0.5, 1.
F.cube(~x1*x2 + x1*x3 + x2*x3, rep(-1, 3), rep(1, 3), c(3, 3, 5))
}

