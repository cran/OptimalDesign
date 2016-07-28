\name{od.print}
\alias{od.print}

\title{Print of a design}

\description{
  Prints a table of a given design.
}

\usage{
  od.print(w, X=NULL, del=1e-6) 
}

\arguments{
  \item{w}{The non-negative vector of length \code{n} representing the design.}
  \item{X}{An \code{n} times \code{d} real matrix giving the \code{d} coordinates of each of the \code{n} design points. 
           The argument \code{X} can also be \code{NULL}; in that case, the non-zero design weights will be printed with their index in \code{1:n}.}
  \item{del}{A small positive number. Each design weight smaller than \code{del} will be treated as zero, therefore it will not be displayed in the table.}
}

\details{
  If supplied, the rows of \code{X} represent coordinates of the design points in the \code{d}-dimensional design region.
  The procedure prints a compact table of the design \code{w}, excluding the coordinates smaller than \code{del}.
  This is often useful, because the support size of the experimental design (the number of non-zero design weights)
  is usually much smaller than the size \code{n} of the design space.
}

\author{ 
  Radoslav Harman, Lenka Filova
}

\seealso{
\code{\link{od.crit}, \link{od.infmat}, \link{od.plot}}
}

\examples{
# Print the uniform approximate design on 24 design points in [0, 2*pi].
w1 <- rep(1/24, 24); od.print(w1, seq(0, (46 / 24)* pi, length=24))

# Print an approximate design on 24 design points in [0, 2*pi] 
# with 8 support points.
w2 <- rep(c(1/8, 0, 0), 8); 
      od.print(w2, seq(0, (46 / 24)* pi, length=24))

# Print an approximate design on 24 design points in [0, 2*pi] 
# with 3 support points.
w3 <- rep(c(1/3, rep(0,7)), 3); 
      od.print(w3, seq(0, (46 / 24)* pi, length=24))

# For the trigonometric model of the first degree on the "full circle" 
# discretized to 24 equispaced points, all three designs have the same 
# information matrix:
F.trig <- F.cube(~I(cos(x1)) + I(sin(x1)), 0, (46 / 24)* pi, 24)
round(od.infmat(F.trig, w1), 6)
round(od.infmat(F.trig, w2), 6)
round(od.infmat(F.trig, w3), 6)
}
