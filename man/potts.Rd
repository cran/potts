\name{potts}
\alias{potts}
\alias{potts.raw}
\alias{potts.potts}
\title{Potts Models}
\description{
    Simulate Potts model using Swendsen-Wang algorithm.
}
\usage{
potts(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"), debug = FALSE,
    outfun = NULL, ...)
}
\arguments{
  \item{obj}{an R vector of class \code{"raw"} that encodes a realization
      of a Potts model, typically the output of \code{\link{packPotts}}.
      Alternatively, an object of class \code{"potts"} from a
      previous invocation of this function can be supplied,
      in which case any missing arguments
      are taken from this object.}
  \item{param}{numeric, canonical parameter vector.  Last component must
      nonnegative (see Details below).}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{boundary}{type of boundary conditions.
      The value of this argument can be abbreviated.}
  \item{debug}{return additional debugging information.}
  \item{outfun}{controls the output.  If a function, then the batch means
      of \code{outfun(tt, \dots)} are returned.  The argument \code{tt}
      is the canonical statistic of the Potts model having the same length
      as the argument \code{param} of this function.  If \code{NULL}, the
      batch means of the canonical statistic are returned.}
  \item{...}{additional arguments for \code{outfun}.}
}
\details{
Runs a Swendsen-Wang algorithm producing a Markov chain with equilibrium
distribution having the specified Potts model.  The state of a Potts model
is a collection of random variables taking values in a finite set.  Here
the finite set is 1, \ldots, \code{ncolor} and the elements are called
\dQuote{colors}.  The random variables are associated with the nodes of
a rectangular lattice, represented by \code{\link{unpackPotts}} as a matrix.
In keeping with calling the values \dQuote{colors}, the random variables
themselves are often called \dQuote{pixels}.  The probability model is an
exponential family with canonical statistic vector of length \code{ncolor + 1}.
The first \code{ncolor} components are the counts of the number of pixels
of each color.  The last component is the number of pairs of neighboring
pixels colored the same.  The corresponding canonical parameter, last
component of the canonical parameter vector (argument \code{param})
must be nonnegative for the Swendsen-Wang algorithm to work (Potts models
are defined for negative dependence parameter, but can't be simulated
using this algorithm).

In the default boundary specification (\code{"torus"}), also called toroidal
or periodic boundary conditions, the vertical edges of the pixel matrix 
are considered glued together, as are the horizontal edges.
Thus corresponding pixels in the first and last rows are considered neighbors,
as are corresponding pixels in the first and last columns.  In the other
boundary specifications there is no such gluing: pixels in the the relative
interiors of the first and last rows and first and last columns have only
three neighbors, and the four corner pixels have only two neighbors.

In the \code{"torus"} and \code{"free"} boundary specifications, all pixels
are counted in determining the color count canonical statistics, which thus
range from zero to \code{nrow * ncol}, where \code{nrow} and \code{ncol} are
the number of rows and columns of the pixel matrix.
In the \code{"condition"} boundary specification, all pixels in the first
and last rows and first and last columns are fixed (conditioned on), and only
the random pixels are counted in determining the color count canonical
statistics, which thus range from zero to \code{(nrow - 2) * (ncol - 2)}.

In the \code{"torus"} boundary specification, all pixels have four neighbors,
so the neighbor pair canonical statistic ranges from zero
to \code{2 * nrow * ncol}.
In the \code{"free"} boundary specification, pixels in the interior have four
neighbors,  those in the relative interior of edges have three, and those
in the corners have two, so the neighbor pair canonical statistic ranges from
zero to \code{nrow * (ncol - 1) + (nrow - 1) * ncol}.
In the \code{"condition"} boundary specification, only neighbor pairs in which
at least one pixel is random are counted, so the neighbor pair canonical
statistic ranges from zero
to \code{(nrow - 2) * (ncol - 1) + (nrow - 1) * (ncol - 2)}.
}
\value{
  an object of class \code{"potts"},
  which is a list containing at least the following components:
  \item{initial}{initial state of Markov chain in the format output
  by \code{\link{packPotts}}.}
  \item{final}{final state of Markov chain in the same format.}
  \item{initial.seed}{value of \code{.Random.seed} before the run.}
  \item{final.seed}{value of \code{.Random.seed} after the run.}
  \item{time}{running time of Markov chain from \code{\link{system.time}}.}
  \item{param}{canonical parameter vector.}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{boundary}{the argument \code{boundary}.}
  \item{batch}{an \code{nbatch} by \code{nout} matrix, where \code{nout} is
      the length of the result returned by \code{outfun} or \code{length(param)}
      if \code{outfun == NULL}; each row is the batch means for the result
      of \code{outfun} or the canonical statistic vector for one batch
      of Markov chain iterations.}
}
\examples{
ncolor <- as.integer(4)
beta <- log(1 + sqrt(ncolor))
theta <- c(rep(0, ncolor), beta)

nrow <- 100
ncol <- 100
x <- matrix(1, nrow = nrow, ncol = ncol)
foo <- packPotts(x, ncolor)

out <- potts(foo, theta, nbatch = 10)
out$batch
\dontrun{image(out$final)}
}
\keyword{misc}
