\name{composite.ll}
\alias{composite.ll}
\alias{gr.composite.ll}
\title{Composite Log Likelihood for Potts Models}
\description{
  Calculate Composite Log Likelihood (CLL) and the gradient of the CLL
  for Potts models.
}
\usage{
composite.ll(theta, t_stat, t_cache=NULL, fapply=lapply)
gr.composite.ll(theta, t_stat, t_cache=NULL, fapply=lapply) 
}
\arguments{

  \item{theta}{numeric canonical parameter vector.  The CLL will be
               evaluated at this point.  It is assumed that the
               component corresponding to the first color has been
               dropped.}

  \item{t_stat}{numeric, canonical statistic vector.  The value of the
                canonical statistic for the full image.}

  \item{t_cache}{list of arrays.  \code{t_cache[[i]][j,]} = the value of
    \eqn{t}{\code{t}} with window \eqn{A_i}{A(i)} replaced by the
    \eqn{j^{th}}{j-th} element of \eqn{C^A}{C^A}.}
	       
  \item{fapply}{function.  Expected to function as \code{lapply} does.  Usefull
    for enabling parallel processing.  E.g. use the \code{mclapply} function
    from the \pkg{multicore} package.}

}
\details{

  For the given value of \code{theta} \code{composite.ll} and
  \code{gr.composite.ll} calculate the CLL and the gradient of the CLL
  respectively for a realized Potts model represented by \code{t_stat}
  and \code{t_cache}.

  \eqn{\mathcal{A}}{sA} is the set of all \emph{windows} to be used in
  calculating the Composite Log Likelihood (CLL) for a Potts model.  A
  \dfn{window} is a collection of adjacent pixels on the lattice of the
  Potts model.  \eqn{A}{A} is used to represent a generic window in
  \eqn{\mathcal{A}}{sA -- meaning `script A'}, the code in this package
  expects that all the windows in \eqn{\mathcal{A}}{sA} have the same
  size and shape.  \eqn{|A|}{|A|} is used to denote the size, or number
  of pixels in a window.  Each pixel in a Potts takes on a value in
  \eqn{C}{C}, the set of possible colors.  For simplicity, this
  implementation takes \eqn{C = \{1,\dots,\texttt{ncolor}\}}{C =
  \{1,\dots,\code{ncolor}\}}.  Elements of \eqn{C}{C} will be referenced
  using \eqn{c_j}{c(j)} with \eqn{j \in \{1,\dots,\texttt{ncolor}\}}{j
  in \{1,\dots,\code{ncolor}\}}.  \eqn{C^A}{C^A} is used to denote all
  the permutations of \eqn{C}{C} across the window \eqn{A}{A}, and
  \eqn{|C|^{|A|}}{|C|^|A|} is used to denote the size of \eqn{C^A}{C^A}.
  In an abuse of notation, we use \eqn{A_a}{A(a)} to refer to the
  \eqn{a^{th}}{a-th} element of \eqn{\mathcal{A}}{sA}.  No ordinal or
  numerical properties of \eqn{\mathcal{A}}{sA}, \eqn{C}{C} or
  \eqn{C^A}{C^A} are used, only that each element in the sets are
  referenced by one and only one indexing value.


  
}
\value{
  \code{composite.ll} returns CLL evaluated at \code{theta}.

  \code{gr.composite.ll} returns a numeric vector of length
  \code{length(theta)} containing the gradient of the CLL at \code{theta}.
}
\seealso{
  \code{\link{generate_t_cache}}, \code{\link{calc_t}}.
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
t_cache_two <- generate_t_cache(x, ncolor, t_stat, nrow*ncol/2, 2,
                                twopixel.nonoverlap)

composite.ll(theta[-1], t_stat, t_cache_mple)
gr.composite.ll(theta[-1], t_stat, t_cache_mple)

\dontrun{
optim.mple <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                    t_stat, t_cache_mple, method="BFGS",
                    control=list(fnscale=-1))
optim.mple$par

optim.two <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                   t_stat, t_cache_two, method="BFGS",
                   control=list(fnscale=-1))
optim.two$par
}

\dontrun{
# or use mclapply to speed things up.
library(multicore)
optim.two <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                   t_stat, t_cache_two, mclapply, method="BFGS",
                   control=list(fnscale=-1))
optim.two$par

}

}
\keyword{misc}
