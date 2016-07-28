\name{od.plot}
\alias{od.plot}

\title{Plot of a design}

\description{
  Plots a graph of a given design.
}

\usage{
  od.plot(w, X=NULL, main="", del=0.001) 
}

\arguments{
  \item{w}{The non-negative vector of length \code{n} representing the design.}
  \item{X}{An \code{n} times \code{d} real matrix giving the \code{d} coordinates of each of the \code{n} design points.
           The argument \code{X} can also be \code{NULL}; in that case, the design weights will be plotted against the
           index vector \code{1:n}.}
  \item{main}{A character string used for the main title of the plot.}
  \item{del}{A small positive number. Each design weight smaller than \code{del} will be treated as zero.}
}

\details{
  If supplied, the rows of \code{X} represent coordinates of the design points in the \code{d}-dimensional design region.
  The procedure generates a plot of all basic two-dimensional projections of the design \code{w}, such that
  the circles corresponding to the same design points have the same color, and the areas of the circles
  are approximately proportional to the weights of the design points.
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
\code{\link{od.crit}, \link{od.infmat}, \link{od.print}}
}

\examples{
# Plot of a uniform design on a discretized line segment.
od.plot(rep(1/21, 21), seq(-1, 1, by=0.1), "Discrete uniform design")

# Plot of a non-uniform simplex lattice design.
F <-  F.simplex(~x1 + x2 + x3 + I(x1*x2) + I(x1*x3) + I(x2*x3) - 1, 3, 4)
od.plot(c(4, 1, 1, 4, 1, 2, 1, 1, 1, 4), F[,c(1, 2, 3)], 
        "Non-uniform simplex lattice design")

# Plot of a Hadamard design on the cube.
F <- F.cube(~x1 + x2 + x3 - 1, rep(-1, 3), rep(1, 3), rep(11, 3))
w <- rep(0, 11^3); w[c(1, 121, 11^3-10, 11^3-110)] <- 1
od.plot(w, F, "Hadamard design (standard view)")
U <- qr.Q(qr(matrix(rnorm(9), ncol=3)))
od.plot(w, F \%*\% U, "Hadamard design (randomly rotated view)")
}
