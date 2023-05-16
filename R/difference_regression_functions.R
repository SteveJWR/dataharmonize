

categorize <- function(data){
  tmp <- data %>% mutate(group = ifelse(((sex == 1) & (educ < 16)), 1,
                                        ifelse(((sex == 2) & (educ < 16)), 2,
                                               ifelse(((sex == 1) & (educ >= 16)), 3, 4))))
  #tmp <- tmp[, !names(tmp) %in% c('sex','educ')]
  tmp$group <- factor(tmp$group, levels = c(1,2,3,4))
  return(tmp)
}


lag_pair_visit_label <- function(X,time.lag, time.window = as.Date("2001/01/01") - as.Date("2000/01/01")){

  X.out <- X %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(date, by_group = T) %>%
    dplyr::mutate(date_diff = date - min(date))

  # only take first day or time lagged day
  X.out <- X.out %>%
    dplyr::filter((date_diff >= time.lag & date_diff <= time.lag + time.window) | date_diff == 0)

  X.out <- X.out %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(date, by_group = T) %>%
    dplyr::mutate(visit = row_number())

  X.out <- X.out %>%
    dplyr::filter(visit <= 2)

  X.out <- X.out %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(num_visits = max(visit))
  X.out <- X.out %>% dplyr::filter(num_visits == 2)
  return(X.out)
}



ZScoreMatchDifferences <- function(X.ref,y.train,z.train, Ny = 30,Nz = 30) {
  mu.y <- mean(y.train$y, na.rm = T)
  sd.y <- sd(y.train$y, na.rm = T)

  mu.z <- mean(z.train$z, na.rm = T)
  sd.z <- sd(z.train$z, na.rm = T)


  #quantile.vec <- rep(NA, Ny + 1)
  conversion.vec <- rep(NA, Ny + 1)

  for(y in seq(0,Ny)){
    q <- round((sd.z/sd.y)*(y - mu.y) + mu.z)
    conversion.vec[y + 1] = max(min(q,Nz),0)
  }

  X.ref.new <- X.ref
  X.ref.new <- X.ref.new %>% arrange(id,visit)
  idx.set.1 <- 2*(1:(nrow(X.ref.new)/2) ) - 1
  idx.set.2 <- 2*(1:(nrow(X.ref.new)/2) )

  y.out <- X.ref.new$y
  z.out <- X.ref.new$z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]

  z1.vec <- z.out[idx.set.1]
  z2.vec <- z.out[idx.set.2]
  z.diff <- z2.vec - z1.vec

  X.out <- X.ref.new[idx.set.1,]
  return(list("X" = X.out, "Z" = z.diff))
}

QuantileMatchDifferences <-  function(X.ref,y.train,z.train, Ny = 30,Nz = 30) {

  #quantile.vec <- rep(NA, Ny + 1)
  conversion.vec <- rep(NA, Ny + 1)

  for(y in seq(0,Ny)){
    p <- mean(y.train$y <= y)
    q <- round(quantile(z.train$z, p))
    q <- min(q,Nz)
    conversion.vec[y + 1] = q
  }

  X.ref.new <- X.ref
  X.ref.new <- X.ref.new %>% arrange(id,visit)
  idx.set.1 <- 2*(1:(nrow(X.ref.new)/2) ) - 1
  idx.set.2 <- 2*(1:(nrow(X.ref.new)/2) )

  y.out <- X.ref.new$y
  z.out <- X.ref.new$z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]

  z1.vec <- z.out[idx.set.1]
  z2.vec <- z.out[idx.set.2]
  z.diff <- z2.vec - z1.vec

  X.out <- X.ref.new[idx.set.1,]
  return(list("X" = X.out, "Z" = z.diff))
}


ImputeOutcomeDifferences <- function(X.ref, y.ref, z.ref, n.impute, Y.train, Z.train, cond.y,
                                     cond.z, mu.y, mu.z, ref.cols, ker.set, R.bins = 1000, threshold = 5 *
                                       10^(-5), max.iter = 50, verbose = F, init.latents = F,
                                     latent.set.covariates = NA, init.latent.set.y = NA, init.latent.set.z = NA){

  X.ref.sort <- X.ref %>% arrange(id, visit)
  y.ref.sort <- X.ref.sort$y
  z.ref.sort <- X.ref.sort$z

  Z.imp <- ImputeOutcomes(X.ref.sort, y.ref.sort, z.ref.sort, n.impute, Y.train, Z.train, cond.y,
                          cond.z, mu.y, mu.z, ref.cols, ker.set, R.bins = R.bins, threshold = threshold,
                          max.iter = max.iter, verbose = verbose, init.latents = init.latents,
                          latent.set.covariates = latent.set.covariates,
                          init.latent.set.y = init.latent.set.y, init.latent.set.z = init.latent.set.z)

  id.set <- unique(X.ref$id)

  num.ids = nrow(X.ref)/2
  idx.1 <- seq(1,2*num.ids, 2)
  idx.2 <- seq(2,2*num.ids, 2)
  complete.vec <- c()
  # idx.1 = which(X.ref$id %in% id.set & X.ref$visit == 1)
  # idx.2 = which(X.ref$id %in% id.set & X.ref$visit == 2)
  #
  z.obs1 = !is.na(X.ref.sort$z[idx.1])
  z.obs2 = !is.na(X.ref.sort$z[idx.2])
  complete.vec <- 1*z.obs1*z.obs2



  #TODO: figure out how to make this faster??
  X.out <- X.ref.sort[idx.1,]
  colnames(X.out) <- colnames(X.ref)
  X.out$complete = complete.vec
  Z.diff <- Z.imp[idx.2,] - Z.imp[idx.1,]
  return(list("X" = X.out, "Z" = Z.diff))
}








ImputationRegressionDifferencesGLMBootstrap <- function(formula, X.ref, y.ref, z.ref,
                                                        n.impute, Y.train, Z.train,
                                                        cond.y, cond.z, mu.y, mu.z,
                                                        ref.cols, ker.set, R.bins = 1000,
                                                        threshold = 5 * 10^(-5), max.iter = 3,
                                                        B.boot = 200, verbose = F)
{
  if (length(ref.cols) == 1) {
    X.un <- data.frame(tmp = unique((X.ref[, ref.cols])))
  }
  else {
    X.un <- unique_X(X.ref[, ref.cols])
  }
  colnames(X.un) = ref.cols
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
  A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
  outcome.cols.y = stringr::str_detect(colnames(Y.train), "y\\d+") |
    stringr::str_detect(colnames(Y.train), "y")
  outcome.cols.z = stringr::str_detect(colnames(Z.train), "z\\d+") |
    stringr::str_detect(colnames(Z.train), "z")
  y.tr <- as.matrix(Y.train[, outcome.cols.y])
  z.tr <- as.matrix(Z.train[, outcome.cols.z])
  X.train.Y <- as.data.frame(Y.train[, ref.cols])
  X.train.Z <- as.data.frame(Z.train[, ref.cols])
  colnames(X.train.Y) = ref.cols
  colnames(X.train.Z) = ref.cols
  latent.set.covariates = X.un
  init.latent.set.y = matrix(NA, nrow = nrow(X.un), ncol = R.bins)
  init.latent.set.z = matrix(NA, nrow = nrow(X.un), ncol = R.bins)
  for (i in seq(nrow(X.un))) {
    if (verbose) {
      cat(paste("Unique mixture estimate:", i, "/", nrow(X.un)),
          end = "\r")
    }
    x <- as.numeric(X.un[i, ])
    weights.y <- weight_vec(x, X.train.Y, ker.set)
    weights.z <- weight_vec(x, X.train.Z, ker.set)
    p.hat.y <- compute_edf(y.tr[, 1], Ny, weights.y)
    p.hat.z <- compute_edf(z.tr[, 1], Nz, weights.z)
    mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y, mu.y,
                                  threshold = threshold)
    mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z, mu.z,
                                  threshold = threshold)
    mixture.y <- mix.y$latent
    mixture.z <- mix.z$latent
    init.latent.set.y[i, ] = mixture.y
    init.latent.set.z[i, ] = mixture.z
  }
  idx.y = seq(nrow(Y.train))
  idx.z = seq(nrow(Z.train))


  coef.boot = NULL
  var.boot = NULL
  for (b in seq(B.boot)) {
    boot.idx.y = sample(idx.y, replace = T)
    boot.idx.z = sample(idx.z, replace = T)
    Y.train.boot = Y.train[boot.idx.y, ]
    Z.train.boot = Z.train[boot.idx.z, ]
    imp.model <- ImputeOutcomeDifferences(X.ref, y.ref, z.ref, n.impute,
                                          Y.train.boot, Z.train.boot, cond.y, cond.z, mu.y,
                                          mu.z, ref.cols, ker.set, R.bins, verbose = F, max.iter = max.iter,
                                          init.latents = TRUE, latent.set.covariates = latent.set.covariates,
                                          init.latent.set.y = init.latent.set.y, init.latent.set.z = init.latent.set.z)

    Z.impute = imp.model$Z
    X.frame = as.data.frame(imp.model$X)
    res.tmp <- ImputationRegressionGLM(formula, X.frame,
                                       Z.impute, fit.cc = F)
    if (is.null(coef.boot)) {
      coef.boot <- matrix(NA, ncol = length(res.tmp$coefficients),
                          nrow = B.boot)
      coef.boot[b, ] <- res.tmp$coefficients
    }
    else {
      coef.boot[b, ] <- res.tmp$coefficients
    }
    if (is.null(var.boot)) {
      var.boot <- matrix(NA, ncol = length(res.tmp$variance),
                         nrow = B.boot)
      var.boot[b, ] <- res.tmp$variance
    }
    else {
      var.boot[b, ] <- res.tmp$variance
    }
    if (verbose) {
      m1 = (round(20 * b/B.boot))
      m2 = 20 - m1
      progress.bar = paste0("|", strrep("=", m1), strrep("-",
                                                         m2), "|")
      cat(paste0("Bootstrap:", b, "/", B.boot, "  ", progress.bar),
          end = "\r")
    }
  }
  beta.est <- colMeans(coef.boot, na.rm = T)
  total.var.block <- matrix(data = colMeans(var.boot, na.rm = T),
                            nrow = length(beta.est), ncol = length(beta.est)) + var(coef.boot,
                                                                                    na.rm = T)
  z.scores <- beta.est/sqrt(diag(total.var.block))
  p.vals.imp <- 2 * pnorm(-abs(z.scores))
  out.list <- list(coefficients = beta.est, variance = total.var.block,
                   `z-scores` = z.scores, `p-values` = p.vals.imp)
  return(out.list)
}














