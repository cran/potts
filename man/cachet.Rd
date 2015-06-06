\name{Cache T}
\alias{generate_t_cache}
\alias{gengridcache}
\alias{gensingleton}
\alias{singleton}
\alias{gentwopixel}
\alias{twopixel}
\alias{twopixel.nonoverlap}
\alias{genfourpixel}
\alias{fourpixel}
\alias{fourpixel.nonoverlap}
\alias{genthreebythree}
\alias{ninepixel.nonoverlap}
\alias{genfourbyfour}
\alias{sixteenpixel.nonoverlap}
\title{Cache calculated Canonical Statistics for Potts Models.}
\description{
  Variety of functions to support caching of calculated canonical
  statistics for Potts Models.  There is some attempt at being 'smart'
  with when to regenerate the statistics.
}
\usage{
generate_t_cache(x, ncolor, t_stat, sizeA, npixel, f,
                 fapply=lapply, gridcache=NULL)

gengridcache(ncolor, sizeCA, ncol)

gensingleton(ncolor)

singleton(x, ncolor, a, idx, gridcache=NULL)

gentwopixel(ncolor)

twopixel(x, ncolor, a, idx, gridcache=NULL)

twopixel.nonoverlap(x, ncolor, a, idx, gridcache=NULL)

genfourpixel(ncolor)

fourpixel(x, ncolor, a, idx, gridcache=NULL)

fourpixel.nonoverlap(x, ncolor, a, idx, gridcache=NULL)

genthreebythree(ncolor)

ninepixel.nonoverlap(x, ncolor, a, idx, gridcache=NULL)

genfourbyfour(ncolor)

sixteenpixel.nonoverlap(x, ncolor, a, idx, gridcache=NULL)
}
\arguments{

  \item{t_stat}{numerical vector of length \code{ncolor}.  Contains the
                canonical statistic for the whole image.}

  \item{sizeA}{numerical.  The number of elements in \eqn{\mathcal{A}}{sA}.}

  \item{sizeCA}{numerical.  The number of elements in \eqn{C^A}{C^A}.}
  
  \item{npixel}{numerical.  The number of pixels in one element of
  \eqn{\mathcal{A}}{sA}.}

  \item{f}{function.  Takes arguments \code{x}, \code{ncolor}, \code{a},
           \code{idx} and \code{ncolor}.  Returns value of
           \code{t_calc_innergrid} with window \eqn{A_a}{A(a)} replaced
           by the \eqn{\texttt{idx}^{th}}{\code{idx}-th}
	   element of \eqn{C^A}{C^A}.}

  \item{fapply}{function.  It should behave exactly as lapply does.
                You can use this argument to enable parallel computing.}
  
  \item{gridcache}{list.  Optional.  If non-null, it is a list of the
                   elements of \eqn{C^A}{C^A}.}

  \item{x}{numeric, 2 dimensional matrix, elements in 1, \dots, \code{ncolor},
           representing a Potts model.}

  \item{ncolor}{numeric.  Number of colors in this Potts Model.}

  \item{ncol}{numeric.  Gives the number of columns in a rectangular
  window.}

  \item{a}{numeric.  Indicates which member of \eqn{\mathcal{A}}{sA} is being
    referenced.}

  \item{idx}{numeric.  Indicates which element of \eqn{C^A}{C^A} is
  being referenced.}
  
}
\details{

  For a description of notation and terminology, see
  \code{\link{composite.ll}}.

  This set of functions is used to generate cached calculations of the
  canonical statistic of a Potts model suitable for passing into
  \code{composite.ll} or \code{gr.composite.ll}.

  All of the calculations using \code{\link{composite.ll}} and these
  caching functions need one of the color components to be dropped for
  the model to be identifiable.  For simplicity, the first color is
  dropped by \code{generate_t_cache}.  In computing the composite log
  likelihood for a Potts model with \code{ncolor} colors, we are
  interested in many calculations across \eqn{C^A}{C^A}, the set of all
  permutations of colors across a window.  These functions facilitate
  those calculations.  \code{gridcache} is a list of \eqn{C^A}{C^A}.

  \code{generate_t_cache} is the main function, and the others are
  intended to be used in conjunction with it.  \code{generate_t_cache}
  creates a list of arrays.  Each array represents one window in the
  image, and each row of the array contains the value of
  \eqn{t(x)}{t(x)} (with one component dropped) found by replacing the
  pixels in that window with one of the elements of \eqn{C^A}{C^A}.

  \code{gengridcache} can generate the \code{gridcache} for any rectangular
  window, give the number of colors, size of \eqn{C^A}{C^A}, and number
  of columns in the window.  \code{gensingleton}, \code{gentwopixel},
  \code{genfourpixel}, \code{genthreebythree} and \code{genfourbyfour}
  are all just simple wrappers for \code{gengridcache}.

  \code{singleton}, \code{twopixel}, \code{twopixel.nonoverlap},
  \code{fourpixel}, \code{fourpixel.nonoverlap},
  \code{ninepixel.nonoverlap} and \code{sixteenpixel.nonoverlap} are
  intended to be passed to \code{generate_t_cache} in the argument
  \code{f}.  They are used to calculate \eqn{t(ca_{idx} \cup X \textbackslash
  A_a)}{\code{t( ca(idx) union X \\ A(a) )}} for the
  \eqn{\texttt{idx}^{th}}{idx-th} element of \eqn{C^{A_a}}{C^A(a)}.

  Functions that have \code{overlap} and \code{nonoverlap} versions
  generate a overlapping and nonoverlapping set of windows respectively.
  
  \code{singleton} is for a single pixel window (Besag or MPLE).

  \code{twopixel} does a two horizontal pixel window.

  \code{fourpixel} does a two by two pixel window.

  \code{ninepixel} does a three by three pixel window.

  \code{sixteenpixel} does a four by four pixel window.
  
}
\value{
  
  Functions that start with \code{gen} return a \code{list} of the
  elements of \eqn{C^A}{C^A}.  

  The other functions (e.g. \code{twopixel}, \code{fourpixel}, \dots)
  return the result of replacing the \code{a}-th window of \code{x} with
  the \code{idx}-th element of \eqn{C^A}{C^A} and calculating
  \code{calc_t_innergrid} for that window.
  
}
\seealso{
  \code{\link{composite.ll}}, \code{\link{calc_t}}.
}

\examples{

ncolor <- 4
beta   <- log(1+sqrt(ncolor))
theta  <- c(rep(0,ncolor), beta)

nrow <- 32
ncol <- 32

x <- matrix(sample(ncolor, nrow*ncol, replace=TRUE), nrow=nrow, ncol=ncol)
foo <- packPotts(x, ncolor)
out <- potts(foo, theta, nbatch=10)
x <- unpackPotts(out$final)

t_stat <- calc_t(x, ncolor)
t_cache_mple <- generate_t_cache(x, ncolor, t_stat, nrow*ncol, 1,
                                 singleton)

\dontrun{
# use multicore to speed things up.
library(multicore)
t_cache_mple <- generate_t_cache(x, ncolor, t_stat, nrow*ncol, 1,
                                 singleton, fapply=mclapply)
}

}
\keyword{misc}
