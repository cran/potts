\name{packPotts}
\alias{packPotts}
\alias{inspectPotts}
\alias{unpackPotts}
\title{Transform Potts Model Data}
\description{
    transform Potts model data from integer matrix to raw vector and
    vice versa.
}
\usage{
packPotts(x, ncolor)
inspectPotts(raw)
unpackPotts(raw)
}
\arguments{
  \item{x}{integer matrix containing Potts model data.  Colors are coded
      from one to \code{ncolor}.}
  \item{ncolor}{integer scalar, number of colors.}
  \item{raw}{vector of type \code{"raw"}.}
}
\value{
  for \code{packPotts} a vector of type \code{"raw"}.
  for \code{inspectPotts} a list containing components
  \code{ncolor}, \code{nrow}, and \code{ncol}.
  for \code{unpackPotts} an integer matrix.
}
\examples{
x <- matrix(sample(4, 2 * 3, replace = TRUE), nrow = 2)
x
foo <- packPotts(x, ncolor = 4)
foo
inspectPotts(foo)
unpackPotts(foo)
}
\keyword{misc}
