\name{F.simplex}
\alias{F.simplex}

\title{Model matrix on a simplex}

\description{
  Creates the matrix of regressors for a model on a discretized \code{(d-1)}-dimensional  
probability simplex, i.e., for a mixture model with \code{d} mixture components.
}

\usage{
  F.simplex(formula, n.factors, n.levels)
}

\arguments{
  \item{formula}{The formula of the model.}
  \item{n.factors}{The integer number greater than or equal to 2, determining the number of factors.}
  \item{n.levels}{The integer number greater than or equal to 2, determining the numbers of levels of each factor (all factors have the same number of levels).}
}

\value{
  The \code{n} times \code{m} matrix of regressors corresponding to \code{m} model parameters and \code{n} design points,
  where \code{n} is equal to \code{(n.factors + n.levels - 2)} choose \code{(n.factors - 1)}.
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
  \code{\link{F.cube}}
}

\examples{
  # The matrix of regressors for the Scheffe quadratic mixture model 
  # with 3 mixture components, each with levels {0, 0.25, 0.5, 0.75, 1}.
  F.simplex(~x1 + x2 + x3 + I(x1*x2) + I(x1*x3) + I(x2*x3) - 1, 3, 5)
  
  # The matrix of regressors for a "special cubic" model with 3 mixture 
  # components, each with levels {0, 0.2, 0.4, 0.6, 0.8, 1}.
  F.simplex(~x1 + x2 + x3 + I(x1*x2) + I(x1*x3) + I(x2*x3) + 
            I(x1*x2*x3) - 1, 3, 6)
  
  # Note that one must be careful when choosing a model for a mixture 
  # experiment: Let us compute the matrix of regressors of the simple 
  # linear mixture model with 4 mixture components, each with levels 
  # {0, 0.5, 1}.
  F1 <- F.simplex(~x1 + x2 + x3 + x4, 4, 3)
  
  # The model has only 4 parameters and as many as 10 design points, 
  # but there is no design that guarantees estimability of the 
  # parameters. This can be shown by evaluating:
  det(od.infmat(F1, rep(1, 10)))
}

