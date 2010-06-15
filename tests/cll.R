
library(potts)

set.seed(42)

ncolor <- as.integer(4)
beta <- log(1 + sqrt(ncolor))
theta <- c(rep(0, ncolor), beta)
nrow <- 32
ncol <- 32

# create potts image
x <- matrix(sample(ncolor, nrow*ncol, replace=TRUE), 
            nrow = nrow, ncol = ncol)
out <- potts(packPotts(x, ncolor), theta, nbatch=1000, blen=1)
x <- unpackPotts(out$final)

# create cache
t_stat <- calc_t(x, ncolor)
t_stat
t_cache_mple <- generate_t_cache(x, ncolor, t_stat, nrow*ncol, 1, 
                                 singleton)

theta.initial <- rep(1, ncolor)
optim.mple <- optim(theta.initial, composite.ll, gr=gr.composite.ll, 
                    t_stat, t_cache_mple, method="BFGS", 
                    control=list(fnscale=-1))
