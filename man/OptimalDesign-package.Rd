\name{OptimalDesign-package}
\alias{OptimalDesign-package}
\alias{OptimalDesign}
\docType{package}

\title{
  OptimalDesign
}

\description{
  \code{OptimalDesign} is a toolbox for constructing optimal (or efficient) experimental designs on finite sets of candidate design points. It provides algorithms for computing \code{D}-, \code{A}-, \code{I}-, and \code{c}-optimal approximate and exact designs of experiments for regression models with uncorrelated observations.

  The package can be used for a wide range of applications, including optimal screening, response-surface, mixture, and weighing designs, as well as designs for standard nonlinear, generalized linear, and survival models via the approach of local optimality.

  Although the design space (the set of candidate design points) is required to be finite, this is typically not a limitation in practice. In many applications the design space is intrinsically finite. Even when a continuous region is assumed theoretically, implementable factor levels are always constrained to a finite set, by instrument resolution and operating procedures.

  In this package, a design on a finite candidate set is represented by a vector \code{w} of replications of the candidate design points. The components of \code{w} are non-negative numbers, one for each candidate point. For exact designs these replications are integers, i.e., they can be directly used. For approximate designs the replication numbers may be non-integer, that is, they must be properly rounded before the design is implemented. Optimal approximate designs also have an auxiliary utility for the construction of exact designs, for example by providing a bound on the efficiency of any exact design.

  The model is represented by the matrix \code{Fx} of candidate regressors: each row of \code{Fx} corresponds to one candidate design point, and each column to one model parameter.

  Several functions allow the user to impose linear constraints directly on the design vector \code{w} (for example, on the total cost if the trials at different design points have different costs, on upper bounds on replication numbers of individual design points or entire strata of the design set, and so on).

  Note that constraints restricting the set of candidate points themselves can be handled simply by removing the impermissible rows from \code{Fx} before computing the design.

  A typical workflow is as follows.

  1) Specify the candidate design space and construct the matrix \code{Fx}, either manually or using helper functions such as \code{Fx_cube}, \code{Fx_simplex}, \code{Fx_blocks}, \code{Fx_glm}, \code{Fx_dose}, or \code{Fx_survival}.

  2) Choose an optimality criterion (\code{"D"}, \code{"A"}, \code{"I"}, \code{"C"} or \code{"c"}) and, if the \code{c}- or \code{C}-criterion is used, specify a linear combination of parameters via the vector \code{h}. Note that the \code{C}-criterion is a regularized version of the standard \code{c}-criterion.

  3) Use \code{od_REX} for standard approximate designs, \code{od_KL} for standard exact designs constrained only by the size (i.e., the total number of trials), or \code{od_RC}, \code{od_AQUA}, or \code{od_MISOCP} to compute an efficient design under more general linear constraints.

  4) Possibly inspect and visualize the resulting design using \code{optcrit}, \code{infmat}, \code{od_print}, \code{od_plot}, \code{dirder}, \code{effbound}, or \code{varfun}.

  Several functions offer additional options that are useful in applications. For example, setting \code{bin = TRUE} in \code{od_KL}, \code{od_AQUA}, \code{od_MISOCP}, or \code{od_RC} enforces a binary design, in which each candidate design point can be used at most once (non-replicable design points). The argument \code{w0} allows the user to fix a lower bound on the replication numbers (also called protected runs), i.e., to compute an optimal augmentation of a given design.
}

\author{
  Radoslav Harman, Lenka Filova
}

\keyword{ Optimal Design }
\keyword{ A-optimality }
\keyword{ D-optimality }
\keyword{ I-optimality }
\keyword{ c-optimality }
