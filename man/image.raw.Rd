\name{image.raw}
\alias{image.raw}
\title{Plot Potts Model Data}
\description{
    plot Potts model data.
}
\usage{
\method{image}{raw}(x, col = c("white", "red", "blue", "green",
    "black", "cyan", "yellow", "magenta"), \dots)
}
\arguments{
  \item{x}{an R vector of class \code{"raw"} that encodes a realization
      of a Potts model, typically the output of \code{\link{packPotts}}
      or of \code{\link{potts}}.}
  \item{col}{a vector of colors.  Must be as many as number of colors of
      Potts model.}
  \item{\dots}{other arguments passed to \code{image.default}.}
}
\section{Bugs}{Too slow for large images.  Needs to be rewritten for
efficient plotting.}
\seealso{\code{\link{potts}}}
\keyword{misc}
