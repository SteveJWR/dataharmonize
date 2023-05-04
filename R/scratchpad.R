
X.ref = sim.data
y.ref = sim.data$Y
z.ref = sim.data$Z
n.impute = n.impute
Y.train = Y.train
Z.train = Z.train
cond.y = cond.y
cond.z = cond.z
mu.y = mu.y
mu.z = mu.z
ref.cols = ref.cols
ker.set = ker.set
R.bins = 1000
threshold = 5 *
  10^(-5)
max.iter = 50
verbose = F
init.latents = NA
latent.set.covariates = NA
init.latent.set.y = NA
init.latent.set.z = NA

ImputeOutcomes <- function (X.ref, y.ref, z.ref, n.impute, Y.train, Z.train, cond.y,
                            cond.z, mu.y, mu.z, ref.cols, ker.set, R.bins = 1000, threshold = 5 *
                              10^(-5), max.iter = 50, verbose = F, init.latents = NA,
                            latent.set.covariates = NA, init.latent.set.y = NA, init.latent.set.z = NA)
{
  if (missing(ref.cols) | missing(ker.set)) {
    if (missing(n.impute)) {
      stop("must indicate number of imputations")
    }
    if (missing(cond.y) | missing(cond.y)) {
      stop("must specify conditional distributions")
    }
    if (missing(mu.y) | missing(mu.z)) {
      stop("must specify regularization parameters")
    }
    if (missing(Y.train) | missing(Z.train)) {
      stop("must specify training data")
    }
    Ny <- length(cond.y(0.5)) - 1
    Nz <- length(cond.z(0.5)) - 1
    A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
    A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
    outcome.cols.y = stringr::str_detect(colnames(Y.train),
                                         "y\\d+") | stringr::str_detect(colnames(Y.train),
                                                                        "y")
    outcome.cols.z = stringr::str_detect(colnames(Z.train),
                                         "z\\d+") | stringr::str_detect(colnames(Z.train),
                                                                        "z")
    y.tr <- Y.train[, outcome.cols.y]
    z.tr <- Z.train[, outcome.cols.z]
    cond.block <- array(NA, c(nrow(X.ref), Ny + 1, Nz + 1))
    cond.block.Boot <- array(NA, c(nrow(X.ref), Ny + 1, Nz +
                                     1))
    p.hat.y <- compute_edf(y.tr, Ny)
    p.hat.z <- compute_edf(z.tr, Nz)
    mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y, mu.y,
                                  threshold = threshold, max.iter = max.iter)
    mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z, mu.z,
                                  threshold = threshold, max.iter = max.iter)
    mixture.y <- mix.y$latent
    mixture.z <- mix.z$latent
    p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,
                                                       cond.z, mixture.y, mixture.z)
    for (id in 1:nrow(X.ref)) {
      cond.block[id, , ] <- p.z.cond.y.slice
    }
    cond.mat <- matrix(NA, nrow(X.ref), Nz + 1)
    for (j in seq(length(y.ref))) {
      if (!is.na(y.ref[j])) {
        cond.mat[j, ] <- cond.block[j, y.ref[j] + 1,
        ]
      }
      else {
        cond.mat[j, ] <- rep(1/(Nz + 1), Nz + 1)
      }
    }
    Z.imp <- matrix(NA, length(y.ref), n.impute)
    for (k in 1:length(y.ref)) {
      if (!is.na(y.ref[k])) {
        p.vec <- as.numeric(cond.mat[k, ])
        Z.imp.row <- sample(seq(0, Nz), n.impute, prob = p.vec,
                            replace = T)
        Z.imp[k, ] <- Z.imp.row
      }
      else {
        Z.imp[k, ] <- rep(NA, n.impute)
      }
    }
    na.idx <- which(is.na(Z.imp[, 1]))
    for (i in na.idx) {
      Z.imp[i, ] <- z.ref[i]
    }
  }
  else {
    if (length(ref.cols) == 1) {
      X.un <- data.frame(tmp = unique((X.ref[, ref.cols])))
    }
    else {
      X.un <- unique_X(X.ref[, ref.cols])
    }
    colnames(X.un) = ref.cols
    if (length(ref.cols) != length(ker.set)) {
      stop("ref.cols and ker.set must be the same length")
    }
    if (missing(n.impute)) {
      stop("must indicate number of imputations")
    }
    if (missing(cond.y) | missing(cond.y)) {
      stop("must specify conditional distributions")
    }
    if (missing(mu.y) | missing(mu.z)) {
      stop("must specify regularization parameters")
    }
    if (missing(Y.train) | missing(Z.train)) {
      stop("must specify training data")
    }
    Ny <- length(cond.y(0.5)) - 1
    Nz <- length(cond.z(0.5)) - 1
    A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
    A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
    outcome.cols.y = stringr::str_detect(colnames(Y.train),
                                         "y\\d+") | stringr::str_detect(colnames(Y.train),
                                                                        "y")
    outcome.cols.z = stringr::str_detect(colnames(Z.train),
                                         "z\\d+") | stringr::str_detect(colnames(Z.train),
                                                                        "z")
    y.tr <- as.matrix(Y.train[, outcome.cols.y])
    z.tr <- as.matrix(Z.train[, outcome.cols.z])
    X.train.Y <- as.data.frame(Y.train[, ref.cols])
    X.train.Z <- as.data.frame(Z.train[, ref.cols])
    colnames(X.train.Y) = ref.cols
    colnames(X.train.Z) = ref.cols
    cond.block <- array(NA, c(nrow(X.ref), Ny + 1, Nz + 1))
    for (i in seq(nrow(X.un))) {
      if (verbose) {
        cat(paste("Unique mixture estimate:", i, "/",
                  nrow(X.un)), end = "\r")
      }
      x <- as.numeric(X.un[i, ])
      weights.y <- weight_vec(x, X.train.Y, ker.set)
      weights.z <- weight_vec(x, X.train.Z, ker.set)
      p.hat.y <- compute_edf(y.tr[, 1], Ny, weights.y)
      p.hat.z <- compute_edf(z.tr[, 1], Nz, weights.z)
      if (init.latents) {
        j = which(!colSums(t(latent.set.covariates) !=
                             x))
        init.latent.y = init.latent.set.y[j, ]
        init.latent.z = init.latent.set.z[j, ]
        mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y,
                                      mu.y, threshold = threshold, max.iter = max.iter,
                                      init.latent = init.latent.y)
        mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z,
                                      mu.z, threshold = threshold, max.iter = max.iter,
                                      init.latent = init.latent.z)
      }
      else {
        mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y,
                                      mu.y, threshold = threshold, max.iter = max.iter)
        mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z,
                                      mu.z, threshold = threshold, max.iter = max.iter)
      }
      mixture.y <- mix.y$latent
      mixture.z <- mix.z$latent
      p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,
                                                         cond.z, mixture.y, mixture.z, R.bins)
      if (length(ref.cols) > 1) {
        match.idx <- which(apply(X.ref[, ref.cols], 1,
                                 function(z) return(all(z == x))))
      }
      else {
        match.idx <- which(X.ref[, ref.cols] == x)
      }
      for (id in match.idx) {
        cond.block[id, , ] <- p.z.cond.y.slice
      }
    }
    cond.mat <- matrix(NA, nrow(X.ref), Nz + 1)
    for (j in seq(length(y.ref))) {
      if (!is.na(y.ref[j])) {
        cond.mat[j, ] <- cond.block[j, y.ref[j] + 1,
        ]
      }
      else {
        cond.mat[j, ] <- rep(1/(Nz + 1), Nz + 1)
      }
    }
    Z.imp <- matrix(NA, length(y.ref), n.impute)
    for (k in 1:length(y.ref)) {
      if (!is.na(y.ref[k])) {
        p.vec <- as.numeric(cond.mat[k, ])
        Z.imp.row <- sample(seq(0, Nz), n.impute, prob = p.vec,
                            replace = T)
        Z.imp[k, ] <- Z.imp.row
      }
      else {
        Z.imp[k, ] <- rep(NA, n.impute)
      }
    }
    na.idx <- which(is.na(Z.imp[, 1]))
    for (i in na.idx) {
      Z.imp[i, ] <- z.ref[i]
    }
  }
  return(Z.imp)
}







# saveRDS(mixture.y, "data/mixture_y.rds")
# saveRDS(mixture.z, "data/mixture_z.rds")
mixture.y <- readRDS("data/mixture_y.rds")
mixture.z <- readRDS("data/mixture_z.rds")



library(reshape2)
library(ggplot2)
library(Matrix)


p.gz <- as.matrix(p.gamma.zeta)
longData<-melt(as.matrix(p.gz.coarse - p.gz.fine))
#longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))


plot(rowMeans(p.gz.coarse - p.gz.fine))
