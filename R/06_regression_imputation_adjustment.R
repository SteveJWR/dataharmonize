rm(list = ls())

library(dplyr)
library(dnoiseR)
library(ggplot2)
library(sandwich)
source("R/paper_functions.R")
#source("R/01_functions.R")
# model 1
# model 2



slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}


include.bootstrap = T
B.boot = 50
n.sims = 5 # TODO: Change this to 500

#n.set <- c(1000,2000)
#n.set <- c(1000, 2000, 5000)


n.true <- 10**(6)

Ny <- 30 # TODO: return to 15
Nz <- 30

R.bins = 1000

cond.y <- generate_cond_binomial(Ny)
cond.z <- generate_cond_binomial(Nz)
# h1 <- 2
# h2 <- 0.5
# cond.y <- conditional_mkm(Ny, gaussian_kernel, h1)
# cond.z <- conditional_mkm(Nz, gaussian_kernel, h2)

#rho.grid = c(0.3,0.3,0.7,0.7)



# uniform age sample
#X <- round(runif(n.true, min = 54.5, max = 80.5))
x.grid = c(55,63,72,80)
x.grid = seq(55,80)
X <- rep(x.grid, n.true/4)
X <- rep(x.grid, round(n.true/length(x.grid)))


n.set <- c(100,200,500,1000, 2000, 5000)
n.set <- round(c(100,200,500,1000, 2000, 5000)/13)*13 # smoothed it over

rho.grid = 0.8*(x.grid < 61) +  0.2*(x.grid >= 61)

#X <- rep(c(55,80), n.true/2)
lat.mu <- (-1/10*(X - 67.5)) +  0.5
lat.mu <- -(1/10*(X - 67.5))^2 +  0.5



run.true.pop = F
if(run.true.pop){
  #U <- runif(n.true)
  #gamma.true.seq <- dnoiseR::logistic(lat.mu + (U <= 1/2)*(2 + rnorm(length(X),sd = 1)) -(U > 1/2)*(2 + rnorm(length(X),sd = 1)))
  gamma.true.seq <- logistic(lat.mu +rnorm(length(X),sd = 1))
  #gamma.true.seq <- logistic(lat.mu +rnorm(length(X),sd = 1))
  # plot(hist(gamma.true.seq))
  # plot(hist(sqrt(gamma.true.seq)))
  #z <- simulate_test_cond(rep(c(2,2),n.true/2),cond.z,gamma.true.seq)
  z <- rbinom(length(gamma.true.seq),Nz, gamma.true.seq)
  #z <- simulate_test_cond_uni(cond.z ,gamma.true.seq)
  dat.pop <- data.frame("Z" = z, "age" = X)

  mod.full.pop <- glm(Z ~ age, data = dat.pop)
  beta.true <- mod.full.pop$coefficients[2]

  #plot(hist(z,breaks = 36))
} else {
  beta.true <- 0
}


#unif.1 <- scale_kernel(uniform_kernel,1)
gaussian_kernel.2 <- scale_kernel(gaussian_kernel,2)
ref.cols <- "age"
ker.set <- list(gaussian_kernel.2)

compute.true.latent = T
if(compute.true.latent){
  x.grid = c(55,63,72,80)
  x.grid = seq(55,80)
  X <- rep(x.grid, n.true/4) # true.x
  X <- rep(x.grid, round(n.true/length(x.grid)))

  lat.mu <- -(1/10*(X - 67.5))^2 +  0.5
  #lat.mu2 <- (-1/5*(X - 67.5) +  2.5)*sin(X)

  #U <- runif(n.true)
  #gamma2 <- logistic(lat.mu2 +  (U <= 1/2)*(2 + rnorm(length(X),sd = 1)) -(U > 1/2)*(2 + rnorm(length(X),sd = 1)))

  gamma2 <- logistic(lat.mu +rnorm(length(X),sd = 1))
  #gamma.true.seq <- logistic(lat.mu +rnorm(length(X),sd = 1))#gamma2 <- logistic(lat.mu2 +rnorm(length(X),sd = 1))
  # one to one mapping of latent traits
  gamma1 <- sqrt(gamma2) #gamma2**(2)
  #gamma1 <- gamma2*(1/5) + 4/5
  # true.mix.y <- hist(gamma1, breaks = seq(0,1,length.out = 1001))$density
  # true.mix.z <- hist(gamma2, breaks = seq(0,1,length.out = 1001))$density
  # true.mix.y <- true.mix.y/sum(true.mix.y)
  # true.mix.z <- true.mix.z/sum(true.mix.z)
  # plot(hist(gamma1, breaks = 50))
  # plot(hist(gamma2, breaks = 50))
  # y1 <- rbinom(length(gamma1),Ny, gamma1)
  # y2 <- rbinom(length(gamma1),Ny, gamma1)
  # z1 <- rbinom(length(gamma2),Nz, gamma2)
  # z2 <- rbinom(length(gamma2),Nz, gamma2)

  # Y <- matrix(c(y1,y2), ncol= 2)
  # Z <- matrix(c(z1,z2), ncol= 2)
  round.X <- round(X)
  X.ref = data.frame("age" = round.X)
  X.un <- data.frame(tmp = sort(unique((X.ref[, ref.cols]))))
  colnames(X.un) = "age"
  mixture.y.set <- matrix(NA,nrow = nrow(X.un), ncol = R.bins)
  mixture.z.set <- matrix(NA,nrow = nrow(X.un), ncol = R.bins)
  for (j in seq(nrow(X.un))) {
    x <- as.numeric(X.un[j, ])
    match.id = which(x == round.X)
    gamma1.counts = round((R.bins - 1)*gamma1[match.id]) + 1
    mixture.y = compute_edf(gamma1.counts, R.bins)
    mixture.y = mixture.y[2:(R.bins + 1)]

    gamma2.counts = round((R.bins - 1)*gamma2[match.id]) + 1
    mixture.z = compute_edf(gamma2.counts, R.bins)
    mixture.z = mixture.z[2:(R.bins + 1)]

    mixture.y.set[j,] = mixture.y
    mixture.z.set[j,] = mixture.z
  }
}
# plot(mixture.y.set[1,])
# plot(mixture.y.set[2,])
# plot(mixture.y.set[3,])
# plot(mixture.y.set[4,])
#
#
# plot(mixture.z.set[1,])
# plot(mixture.z.set[2,])
# plot(mixture.z.set[3,])
# plot(mixture.z.set[4,])

#unif.1 <- scale_kernel(uniform_kernel,1)
gaussian_kernel.2 <- scale_kernel(gaussian_kernel,2)
ref.cols <- "age"
ker.set <- list(gaussian_kernel.2)



# deviation results

dev.cc <- matrix(NA, nrow = n.sims, ncol = length(n.set))
dev <- matrix(NA, nrow = n.sims, ncol = length(n.set))
dev.cov.adj <- matrix(NA, nrow = n.sims, ncol = length(n.set))
dev.true.latent <- matrix(NA, nrow = n.sims, ncol = length(n.set))
dev.z.score <- matrix(NA, nrow = n.sims, ncol = length(n.set))
dev.quantile <- matrix(NA, nrow = n.sims, ncol = length(n.set))
dev.bootstrap <- matrix(NA, nrow = n.sims, ncol = length(n.set))
# coverage results

cover.cc <- matrix(NA, nrow = n.sims, ncol = length(n.set))
cover <- matrix(NA, nrow = n.sims, ncol = length(n.set))
cover.cov.adj <- matrix(NA, nrow = n.sims, ncol = length(n.set))
cover.true.latent <- matrix(NA, nrow = n.sims, ncol = length(n.set))
cover.z.score <- matrix(NA, nrow = n.sims, ncol = length(n.set))
cover.quantile <- matrix(NA, nrow = n.sims, ncol = length(n.set))
cover.bootstrap <- matrix(NA, nrow = n.sims, ncol = length(n.set))

load.prior.results = F
if(load.prior.results){
  dev.cc <- readRDS("data/06_dev_cc.rds")
  dev <- readRDS("data/06_dev.rds")
  dev.cov.adj <- readRDS("data/06_dev_cov_adj.rds")
  dev.true.latent <- readRDS("data/06_dev_true_latent.rds")
  dev.z.score <- readRDS("data/06_dev_zscore.rds")
  dev.quantile <- readRDS("data/06_dev_quantile.rds")
  dev.bootstrap <- readRDS("data/06_dev_bootstrap.rds")

  cover.cc <- readRDS( "data/06_cover_cc.rds")
  cover <- readRDS("data/06_cover.rds")
  cover.cov.adj <- readRDS("data/06_cover_cov_adj.rds")
  cover.true.latent <- readRDS("data/06_cover_true_latent.rds")
  cover.z.score <- readRDS("data/06_cover_zscore.rds")
  cover.quantile <- readRDS("data/06_cover_quantile.rds")
  dev.bootstrap <- readRDS("data/06_dev_bootstrap.rds")
}

#unif.1 <- scale_kernel(uniform_kernel,1)



ref.cols <- "age"
gaussian_kernel.1 <- scale_kernel(gaussian_kernel,5/(n.set[1]**(1/5)))
gaussian_kernel.2 <- scale_kernel(gaussian_kernel,5/(n.set[2]**(1/5)))
gaussian_kernel.3 <- scale_kernel(gaussian_kernel,5/(n.set[3]**(1/5)))
gaussian_kernel.4 <- scale_kernel(gaussian_kernel,5/(n.set[4]**(1/5)))
gaussian_kernel.5 <- scale_kernel(gaussian_kernel,5/(n.set[5]**(1/5)))
gaussian_kernel.6 <- scale_kernel(gaussian_kernel,5/(n.set[6]**(1/5)))

ker.set.1 <- list(gaussian_kernel.1)
ker.set.2 <- list(gaussian_kernel.2)
ker.set.3 <- list(gaussian_kernel.3)
ker.set.4 <- list(gaussian_kernel.4)
ker.set.5 <- list(gaussian_kernel.5)
ker.set.6 <- list(gaussian_kernel.6)

ker.list <- list(ker.set.1,
                 ker.set.2,
                 ker.set.3,
                 ker.set.4,
                 ker.set.5,
                 ker.set.6)
# gaussian_kernel.10 <- scale_kernel(gaussian_kernel,10)
# ref.cols <- "age"
# ker.set <- list(gaussian_kernel.10)


mu.y = 0.0 #0.1 #0.3
mu.z = 0.0 #0.1 #0.3
n.impute = 50
fmla <- formula(outcome ~ age )
# we want to plot coverage and bias of conversion

# TODO: Try again in case gamma2 is bimodal

set.seed(id) # finding the correct id's for the simulations.
sim.start = 1

for(sim in sim.start:n.sims){
  cat(paste0("Number of sims: ", sim, "/",n.sims), end = "\n")
  if(sim %% 20 == 0){

    saveRDS(dev.cc, paste0("data/06_dev_cc",id, ".rds"))
    saveRDS(dev, paste0("data/06_dev",id, ".rds"))
    saveRDS(dev.cov.adj, paste0("data/06_dev_cov_adj",id, ".rds"))
    saveRDS(dev.true.latent, paste0("data/06_dev_true_latent",id, ".rds"))
    saveRDS(dev.z.score, paste0("data/06_dev_zscore",id, ".rds"))
    saveRDS(dev.quantile, paste0("data/06_dev_quantile",id, ".rds"))

    saveRDS(cover.cc, paste0("data/06_cover_cc",id, ".rds"))
    saveRDS(cover, paste0("data/06_cover",id, ".rds"))
    saveRDS(cover.cov.adj, paste0("data/06_cover_cov_adj",id, ".rds"))
    saveRDS(cover.true.latent, paste0("data/06_cover_true_latent",id, ".rds"))
    saveRDS(cover.z.score, paste0("data/06_cover_zscore",id, ".rds"))
    saveRDS(cover.quantile, paste0("data/06_cover_quantile",id, ".rds"))
    if(include.bootstrap){
      saveRDS(dev.bootstrap, paste0("data/06_dev_bootstrap",id, ".rds"))
      saveRDS(cover.bootstrap, paste0("data/06_cover_bootstrap",id, ".rds"))
    }
  }

  for(i in seq(length(n.set))){
    #print(i)
    n = n.set[i]
    ker.set = ker.list[[i]]
    #X <- sample(x.grid, n, replace = T)
    X <- rep(x.grid, each = n/4) # grid for different observed X values
    X <- rep(x.grid, round(n/length(x.grid)))

    #TODO: ensure that this is set at the beginning with an x.grid parameter
    #X <- rep(c(55,80), n/2)

    lat.mu2 <- -(1/10*(X - 67.5))^2 +  0.5
    #lat.mu2 <- (-1/5*(X - 67.5) +  2.5)*sin(X)

    # missingness pattern
    idx1 <- c()
    idx2 <- c()
    for(k in seq(length(x.grid))){
      x = x.grid[k]
      rho = rho.grid[k]
      idx = which(X == x)
      idx1 = c(idx1,sample(idx, size = round(rho*length(idx)), replace = F))
      idx2 = c(idx2,setdiff(idx, idx1))
    }


    U <- runif(n)
    gamma2 <- logistic(lat.mu2 +rnorm(length(X),sd = 1))
    #gamma.true.seq <- logistic(lat.mu +rnorm(length(X),sd = 1))#gamma2 <- logistic(lat.mu2 +rnorm(length(X),sd = 1))
    # one to one mapping of latent traits
    gamma1 <- sqrt(gamma2)
    #gamma1 <- gamma2*(1/5) + 4/5
    # true.mix.y <- hist(gamma1, breaks = seq(0,1,length.out = 1001))$density
    # true.mix.z <- hist(gamma2, breaks = seq(0,1,length.out = 1001))$density
    # true.mix.y <- true.mix.y/sum(true.mix.y)
    # true.mix.z <- true.mix.z/sum(true.mix.z)
    # plot(hist(gamma1, breaks = 50))
    # plot(hist(gamma2, breaks = 50))
    # y1 <- rbinom(length(gamma1),Ny, gamma1)
    # y2 <- rbinom(length(gamma1),Ny, gamma1)
    # z1 <- rbinom(length(gamma2),Nz, gamma2)
    # z2 <- rbinom(length(gamma2),Nz, gamma2)

    # Y <- matrix(c(y1,y2), ncol= 2)
    # Z <- matrix(c(z1,z2), ncol= 2)
    # Y <- simulate_test_cond(obs.set = rep(2,length(gamma1)),cond.y ,gamma1)
    # Z <- simulate_test_cond(obs.set = rep(2,length(gamma2)),cond.z ,gamma2)
    Y <- simulate_test_cond(obs.set = rep(1,length(gamma1)),cond.y ,gamma1)
    Z <- simulate_test_cond(obs.set = rep(1,length(gamma2)),cond.z ,gamma2)
    y1 <- Y[,1]
    z1 <- Z[,1]

    sim.data <- data.frame("Y" = Y[,1],"Z" = Z[,1], "age" = X)
    sim.data$Y[idx2] = NA
    sim.data$Z[idx1] = NA
    Y <- as.matrix(Y[idx1,1])
    X.y <- X[idx1]
    Z <- as.matrix(Z[idx2,1])
    X.z <- X[idx2]
    Y.train <- cbind(Y, X.y)
    Z.train <- cbind(Z, X.z)
    # colnames(Y.train) = c("y1", "y2","age") # do we need any of these?
    # colnames(Z.train) = c("z1", "z2","age")
    colnames(Y.train) = c("y", "age") # do we need any of these?
    colnames(Z.train) = c("z", "age")
    cc.sim.data <- data.frame("Z" = Z[,1], "age" = X.z)

    mod.full.pop <- glm(Z ~ age, data = cc.sim.data)
    beta.hat <- mod.full.pop$coefficients[2]


    # ensure in the R package we have the prefixes as y and z
    colnames(Y) <- c("y")
    colnames(Z) <- c("z")

    Y = as.data.frame(Y)
    Z = as.data.frame(Z)
    list.zscore <- ZScoreConversion(sim.data,Y,Z, Ny, Nz)
    list.quantile <- QuantileConversion(sim.data,Y,Z, Ny, Nz)


    z.score.sim.data <- cbind(list.zscore$Z, list.zscore$X)
    colnames(z.score.sim.data)[1] = "outcome"

    quantile.sim.data <- cbind(list.quantile$Z, list.quantile$X)
    colnames(quantile.sim.data)[1] = "outcome"

    #naive z matching
    fit.z.score <- glm(fmla, data = z.score.sim.data)
    z.score.coefs <- fit.z.score$coefficients
    #z.score.coefs



    #quantile matching
    fit.quantile <- glm(fmla, data = quantile.sim.data)
    quantile.match.coefs <- fit.quantile$coefficients
    #quantile.match.coefs


    impute.sim <- ImputeOutcomes(sim.data,sim.data$Y,sim.data$Z,n.impute,
                                 Y,Z,cond.y,cond.z,
                                 mu.y,mu.z,R.bins = 1000)

    #TODO: Ensure the conditional probabilities are matching correctly.
    impute.sim.cov.adj <- ImputeOutcomes(sim.data,sim.data$Y,sim.data$Z,n.impute,
                                         Y.train,Z.train,cond.y,cond.z,
                                         mu.y,mu.z,ref.cols, ker.set, R.bins = 1000)

    impute.sim.true.latents <- ImputeOutcomesTrueLatent(sim.data,sim.data$Y,sim.data$Z,n.impute,
                                                        cond.y,cond.z,
                                                        ref.cols, mixture.y.set, mixture.z.set)


    X.frame <- data.frame("age" = X)
    X.frame$complete <- 1*!is.na(sim.data$Z)

    imp.reg <- ImputationRegressionGLM(fmla, X.frame, impute.sim, fit.cc = T)
    imp.reg.cov.adj <- ImputationRegressionGLM(fmla, X.frame, impute.sim.cov.adj, fit.cc = T)
    imp.reg.true.latent <- ImputationRegressionGLM(fmla, X.frame, impute.sim.true.latents, fit.cc = T)




    beta.impute <-imp.reg$coefficients[2]
    beta.cov.adj <- imp.reg.cov.adj$coefficients[2]
    beta.true.latent <- imp.reg.true.latent$coefficients[2]

    beta.hat.cc <-imp.reg$`cc-coefficients`[2]
    beta.z.score <-z.score.coefs[2]
    beta.quantile <-quantile.match.coefs[2]


    imp.reg$`cc-coefficients`[2] - 2.0*sqrt(imp.reg$`cc-variance`[2,2])  <= beta.true & imp.reg$`cc-coefficients`[2] + 2.0*sqrt(imp.reg$`cc-variance`[2,2]) >= beta.true
    imp.reg$coefficients[2] - 2.0*sqrt(imp.reg$variance[2,2])  <= beta.true & imp.reg$coefficients[2] + 2.0*sqrt(imp.reg$variance[2,2]) >= beta.true
    imp.reg.cov.adj$coefficients[2] - 2.0*sqrt(imp.reg.cov.adj$variance[2,2])  <= beta.true & imp.reg.cov.adj$coefficients[2] + 2.0*sqrt(imp.reg.cov.adj$variance[2,2]) >= beta.true


    dev.cc[sim,i] = beta.hat.cc - beta.true
    dev[sim,i] = beta.impute - beta.true
    dev.cov.adj[sim,i] = beta.cov.adj - beta.true
    dev.true.latent[sim,i] = beta.true.latent - beta.true
    dev.z.score[sim,i] = beta.z.score - beta.true
    dev.quantile[sim,i] = beta.quantile - beta.true


    cover.cc[sim,i] = imp.reg$`cc-coefficients`[2] - 2.0*sqrt(imp.reg$`cc-variance`[2,2])  <= beta.true & imp.reg$`cc-coefficients`[2] + 2.0*sqrt(imp.reg$`cc-variance`[2,2]) >= beta.true
    cover[sim,i] = imp.reg$coefficients[2] - 2.0*sqrt(imp.reg$variance[2,2])  <= beta.true & imp.reg$coefficients[2] + 2.0*sqrt(imp.reg$variance[2,2]) >= beta.true
    cover.cov.adj[sim,i] = imp.reg.cov.adj$coefficients[2] - 2.0*sqrt(imp.reg.cov.adj$variance[2,2])  <= beta.true & imp.reg.cov.adj$coefficients[2] + 2.0*sqrt(imp.reg.cov.adj$variance[2,2]) >= beta.true
    cover.true.latent[sim,i] = imp.reg.true.latent$coefficients[2] - 2.0*sqrt(imp.reg.true.latent$variance[2,2])  <= beta.true & imp.reg.true.latent$coefficients[2] + 2.0*sqrt(imp.reg.true.latent$variance[2,2]) >= beta.true
    cover.z.score[sim,i] = beta.z.score - 2.0*sqrt(sandwich(fit.z.score)[2,2]) <= beta.true & beta.z.score + 2.0*sqrt(sandwich(fit.z.score)[2,2]) >=  beta.true
    cover.quantile[sim,i] = beta.quantile - 2.0*sqrt(sandwich(fit.quantile)[2,2]) <= beta.true & beta.quantile + 2.0*sqrt(sandwich(fit.quantile)[2,2]) >=  beta.true


    if(include.bootstrap){
      bootstrap.results <- ImputationRegressionGLMBootstrap(fmla, sim.data,sim.data$Y,sim.data$Z,n.impute,
                                                            Y.train,Z.train,cond.y,cond.z,
                                                            mu.y,mu.z,ref.cols, ker.set, R.bins = 1000, B.boot = B.boot, verbose = T)

      dev.bootstrap[sim,i] = bootstrap.results$coefficients[2] - beta.true
      cover.bootstrap[sim,i] = bootstrap.results$coefficients[2] - 2.0*sqrt(bootstrap.results$variance[2,2])  <= beta.true & bootstrap.results$coefficients[2] + 2.0*sqrt(bootstrap.results$variance[2,2]) >= beta.true

    }


  }

}




########

saveRDS(dev.cc, paste0("data/06_dev_cc",id, ".rds"))
saveRDS(dev, paste0("data/06_dev",id, ".rds"))
saveRDS(dev.cov.adj, paste0("data/06_dev_cov_adj",id, ".rds"))
saveRDS(dev.true.latent, paste0("data/06_dev_true_latent",id, ".rds"))
saveRDS(dev.z.score, paste0("data/06_dev_zscore",id, ".rds"))
saveRDS(dev.quantile, paste0("data/06_dev_quantile",id, ".rds"))

saveRDS(cover.cc, paste0("data/06_cover_cc",id, ".rds"))
saveRDS(cover, paste0("data/06_cover",id, ".rds"))
saveRDS(cover.cov.adj, paste0("data/06_cover_cov_adj",id, ".rds"))
saveRDS(cover.true.latent, paste0("data/06_cover_true_latent",id, ".rds"))
saveRDS(cover.z.score, paste0("data/06_cover_zscore",id, ".rds"))
saveRDS(cover.quantile, paste0("data/06_cover_quantile",id, ".rds"))

if(include.bootstrap){
  saveRDS(dev.bootstrap, paste0("data/06_dev_bootstrap",id, ".rds"))
  saveRDS(cover.bootstrap, paste0("data/06_cover_bootstrap",id, ".rds"))
}



make.plots = F

if(make.plots){
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200

  # update this section to concatenate the results
  dev.cc <- readRDS("data/06c_dev_cc1.rds")
  dev <- readRDS("data/06c_dev1.rds")
  dev.cov.adj <- readRDS("data/06c_dev_cov_adj1.rds")
  dev.true.latent <-  readRDS("data/06c_dev_true_latent1.rds")
  dev.z.score <- readRDS("data/06c_dev_zscore1.rds")
  dev.quantile <- readRDS("data/06c_dev_quantile1.rds")
  dev.bootstrap <- readRDS("data/06c_dev_bootstrap1.rds")

  cover.cc <- readRDS("data/06c_cover_cc1.rds")
  cover <- readRDS("data/06c_cover1.rds")
  cover.cov.adj <- readRDS("data/06c_cover_cov_adj1.rds")
  cover.true.latent <-  readRDS("data/06c_cover_true_latent1.rds")
  cover.z.score <- readRDS("data/06c_cover_zscore1.rds")
  cover.quantile <- readRDS("data/06c_cover_quantile1.rds")
  cover.bootstrap <- readRDS("data/06c_cover_bootstrap1.rds")
  for(j in seq(2,200)){

    dev.cc.tmp <- readRDS(paste0("data/06c_dev_cc",j,".rds"))
    dev.tmp <- readRDS(paste0("data/06c_dev",j,".rds"))
    dev.cov.adj.tmp <- readRDS(paste0("data/06c_dev_cov_adj",j,".rds"))
    dev.true.latent.tmp <- readRDS(paste0("data/06c_dev_true_latent",j,".rds"))
    dev.z.score.tmp <- readRDS(paste0("data/06c_dev_zscore",j,".rds"))
    dev.quantile.tmp <- readRDS(paste0("data/06c_dev_quantile",j,".rds"))
    dev.bootstrap.tmp <- readRDS(paste0("data/06c_dev_bootstrap",j,".rds"))

    cover.cc.tmp <- readRDS(paste0("data/06c_cover_cc",j,".rds"))
    cover.tmp <- readRDS(paste0("data/06c_cover",j,".rds"))
    cover.cov.adj.tmp <- readRDS(paste0("data/06c_cover_cov_adj",j,".rds"))
    cover.true.latent.tmp <- readRDS(paste0("data/06c_cover_true_latent",j,".rds"))
    cover.z.score.tmp <- readRDS(paste0("data/06c_cover_zscore",j,".rds"))
    cover.quantile.tmp <- readRDS(paste0("data/06c_cover_quantile",j,".rds"))
    cover.bootstrap.tmp <- readRDS(paste0("data/06c_cover_bootstrap",j,".rds"))

    dev.cc <- rbind(dev.cc, dev.cc.tmp)
    dev <- rbind(dev, dev.tmp)
    dev.cov.adj <- rbind(dev.cov.adj, dev.cov.adj.tmp)
    dev.true.latent <- rbind(dev.true.latent,dev.true.latent.tmp)
    dev.z.score <- rbind(dev.z.score, dev.z.score.tmp)
    dev.quantile <- rbind(dev.quantile, dev.quantile.tmp)
    dev.bootstrap <- rbind(dev.bootstrap, dev.bootstrap.tmp)

    cover.cc <- rbind(cover.cc, cover.cc.tmp)
    cover <- rbind(cover, cover.tmp)
    cover.cov.adj <- rbind(cover.cov.adj, cover.cov.adj.tmp)
    cover.true.latent <- rbind(cover.true.latent,cover.true.latent.tmp)
    cover.z.score <- rbind(cover.z.score, cover.z.score.tmp)
    cover.quantile <- rbind(cover.quantile, cover.quantile.tmp)
    cover.bootstrap <- rbind(cover.bootstrap, cover.bootstrap.tmp)
  }

  n.sims = nrow(cover.cc)
  n.set <- round(c(100,200,500,1000, 2000, 5000)/13)*13 # smoothed it over


  cc.mean.bias <- colMeans(dev.cc, na.rm = T)
  mean.bias <- colMeans(dev, na.rm = T)
  mean.bias.cov.adj <- colMeans(dev.cov.adj, na.rm = T)
  mean.bias.true.latent <- colMeans(dev.true.latent, na.rm = T)
  mean.bias.bootstrap <- colMeans(dev.bootstrap, na.rm = T)
  z.score.mean.bias <- colMeans(dev.z.score, na.rm = T)
  quantile.bias <- colMeans(dev.quantile, na.rm = T)

  cc.rmse <- sqrt(colMeans(abs(dev.cc)^2, na.rm = T))
  rmse <-sqrt( colMeans(abs(dev)^2, na.rm = T))
  rmse.cov.adj <-sqrt( colMeans(abs(dev.cov.adj)^2, na.rm = T))
  rmse.true.latent <-sqrt( colMeans(abs(dev.true.latent)^2, na.rm = T))
  rmse.bootstrap <-sqrt( colMeans(abs(dev.bootstrap)^2, na.rm = T))
  z.score.rmse <- sqrt(colMeans(abs(dev.z.score)^2, na.rm = T))
  quantile.rmse <-sqrt( colMeans(abs(dev.quantile)^2, na.rm = T))

  cc.rmse.sd <- colSDs(dev.cc, na.rm = T)/sqrt(n.sims)
  rmse.sd <- colSDs(dev, na.rm = T)/sqrt(n.sims)
  rmse.cov.adj.sd <- colSDs(dev.cov.adj, na.rm = T)/sqrt(n.sims)
  rmse.true.latent.sd <- colSDs(dev.true.latent, na.rm = T)/sqrt(n.sims)
  rmse.bootstrap.sd <- colSDs(dev.bootstrap, na.rm = T)/sqrt(n.sims)
  z.score.rmse.sd <- colSDs(dev.z.score, na.rm = T)/sqrt(n.sims)
  quantile.rmse.sd <- colSDs(dev.quantile, na.rm = T)/sqrt(n.sims)


  res.data <- data.frame("method" = c(rep("Complete Case", length(n.set)),
                                      rep("DNOISE", length(n.set)),
                                      rep("DNOISE (Cov. Adj.)", length(n.set)),
                                      rep("DNOISE (Oracle)", length(n.set)),
                                      rep("DNOISE (Bootstrap)", length(n.set)),
                                      rep("Z Score", length(n.set)),
                                      rep("Quantile", length(n.set))),
                         "n" = c(n.set,n.set,n.set,
                                 n.set,n.set,n.set, n.set),
                         "bias" = c(cc.mean.bias,mean.bias, mean.bias.cov.adj,
                                    mean.bias.true.latent, mean.bias.bootstrap, z.score.mean.bias,quantile.bias),
                         "rmse" = c(cc.rmse,rmse,rmse.cov.adj,
                                    rmse.true.latent, rmse.bootstrap,
                                    z.score.rmse,quantile.rmse),
                         "rmse_sd" = c(cc.rmse.sd,rmse.sd,rmse.cov.adj.sd,
                                       rmse.true.latent.sd,
                                       rmse.true.latent.sd,
                                       z.score.rmse.sd,quantile.rmse.sd))


  plt.bias <- ggplot(res.data, aes(x = log(n), y = bias, color = method)) +
    geom_line() +
    ggtitle("Bias of Regression Estimate") #+
  #geom_line(aes(x = n, y = rmse, color = method)) #+
  #geom_errorbar(aes(ymin = bias - 2*rmse, ymax = bias + 2*rmse))

  plt.bias

  png(filename = "plots/sim_binomial_bias.png",
      width = png.width, height = png.height, res = png.res)

  plt.bias
  # Close the pdf file
  dev.off()

  plt.rmse <- ggplot(res.data, aes(x = log(n), y = log(rmse), color = method)) +
    geom_line() +
    geom_errorbar(aes(ymin = log(rmse - 2*rmse_sd), ymax = log(rmse + 2*rmse_sd)))+
    ggtitle("RMSE of Regression Estimate")

  plt.rmse
  png(filename = "plots/sim_binomial_rmse.png",
      width = png.width, height = png.height, res = png.res)

  plt.rmse
  # Close the pdf file
  dev.off()



  cc.coverage <- colMeans(cover.cc, na.rm = T)
  coverage <- colMeans(cover, na.rm = T)
  coverage.cov.adj <- colMeans(cover.cov.adj, na.rm = T)
  coverage.true.latent <- colMeans(cover.true.latent, na.rm = T)
  coverage.bootstrap <- colMeans(cover.bootstrap, na.rm = T)
  z.score.coverage <- colMeans(cover.z.score, na.rm = T)
  quantile.coverage  <- colMeans(cover.quantile, na.rm = T)

  cov.vec <- c(cc.coverage,coverage,coverage.cov.adj,
               coverage.true.latent, coverage.bootstrap, z.score.coverage,quantile.coverage)
  cov.error <- sqrt(cov.vec*(1 - cov.vec)/sum(!is.na(cover[,1])))

  cov.data <- data.frame("method" = c(rep("Complete Case", length(n.set)),
                                      rep("DNOISE", length(n.set)),
                                      rep("DNOISE (Cov. Adj.)", length(n.set)),
                                      rep("DNOISE (Oracle)", length(n.set)),
                                      rep("DNOISE (Bootstrap)", length(n.set)),
                                      rep("Z Score", length(n.set)),
                                      rep("Quantile", length(n.set))),
                         "n" = c(n.set,n.set,n.set,n.set, n.set, n.set, n.set),
                         "coverage" = cov.vec,
                         "error" = cov.error)


  #cov.data <- cov.data %>% filter(n != 500)
  plt.coverage <- ggplot(cov.data, aes(x = log(n), y = coverage, color = method)) +
    geom_line(position=position_dodge(width=0.2)) +
    geom_errorbar(aes(ymin = coverage - 2*error, ymax = coverage + 2*error), width=0.5,
                  linewidth=0.5, position=position_dodge(width=0.2)) +
    geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
    ggtitle("Coverage of Regression Estimate")

  plt.coverage

  png(filename = "plots/sim_binomial_coverage.png",
      width = png.width, height = png.height, res = png.res)

  plt.coverage
  # Close the pdf file
  dev.off()

}









#### AFTER HERE, REMOVE



#
#
#
#
#
#
#
# ### Current setup is pretty good,
# #TODO: Add this idea of using quantile matching
# # or a generalized optimal transport method.
#
# # colSDs <- function(X, ...){
# #   X.means <- colMeans(X, ...)
# #   v.out <- sqrt(colMeans((X.means - X)^2, ...))
# #   return(v.out)
# # }
#
#
#
#
#
# # expectile curve,
#
#
#
#
# Ny = 200
# Nz = 200
# n = 500000
# X <- round(runif(n, min = 54.5, max = 80.5))
#
# lat.mu2 <- -1/20*(X - 67.5) +  2.5
# #lat.mu2 <- (-1/5*(X - 67.5) +  2.5)*sin(X)
#
#
#
# U <- runif(n)
# gamma2 <- logistic(lat.mu2 +  (U <= 1/2)*(2 + rnorm(length(X),sd = 1)) -(U > 1/2)*(2 + rnorm(length(X),sd = 1)))
# gamma2 <- logistic(lat.mu2 +  (U <= 1/2)*( rnorm(length(X),sd = 1)) -(U > 1/2)*( rnorm(length(X),sd = 1)))
# #gamma2 <- logistic(lat.mu2 +rnorm(length(X),sd = 1))
# # one to one mapping of latent traits
# gamma1 <- sqrt(gamma2)
# #gamma1 <- gamma2*(1/5) + 4/5
# # true.mix.y <- hist(gamma1, breaks = seq(0,1,length.out = 1001))$density
# # true.mix.z <- hist(gamma2, breaks = seq(0,1,length.out = 1001))$density
# # true.mix.y <- true.mix.y/sum(true.mix.y)
# # true.mix.z <- true.mix.z/sum(true.mix.z)
# # plot(hist(gamma1, breaks = 50))
# # plot(hist(gamma2, breaks = 50))
# y1 <- rbinom(length(gamma1),Ny, gamma1)
# z1 <- rbinom(length(gamma2),Nz, gamma2)
#
#
#
#
# quantile.vec <- rep(NA,Ny + 1)
# expectile.vec <- rep(NA,Ny + 1)
#
#
#
#
# for(y.idx in seq(0,Ny)){
#   print(y.idx)
#   which.idx <- which(y1 == y.idx)
#
#   cond.mean <- mean(z1[which.idx], na.rm = T)
#   expectile.vec[y.idx + 1] = cond.mean
#
#   p <- mean(y1 <= y.idx)
#   #q <- ceiling(quantile(z1, p, type = 3))
#   q <- quantile(z1, p, na.rm = T)
#
#   quantile.vec[y.idx + 1] = q
# }
#
#
# plot(seq(0,Ny), quantile.vec, type = "l", col = "blue")
# lines(seq(0,Ny), expectile.vec, type = "l", col = "green")
# legend(60, 60, legend=c("Quantile (OT)", "Conditional Expectation"),
#        col=c("blue", "green"), lty=c(1,1), cex=0.8)
#
#
#
#
#
#
#
# library(mvtnorm)
# theta.seq <- c(-1,-0.5,0,.3,.5,.8,.9,1)
#
#
# n = 1000000
# Ny = 60
# Nz = 60
# conditional.expect <- matrix(NA, nrow = length(theta.seq), ncol = Ny + 1)
#
# for(i in seq(length(theta.seq))){
#   theta = theta.seq[i]
#   sigma <- matrix(c(1,theta,theta,1), 2,2)
#   X <- rmvnorm(n, sigma = sigma)
#   p1 <- pnorm(X[,1])
#   p2 <- pnorm(X[,2])
#   y <- qbinom(p1,Ny, .5)
#   z <- qbinom(p2,Nz, .8)
#   expectile.vec <- rep(NA,Ny + 1)
#   for(y.idx in seq(0,Ny)){
#     print(y.idx)
#     which.idx <- which(y == y.idx)
#
#     cond.mean <- mean(z[which.idx], na.rm = T)
#     expectile.vec[y.idx + 1] = cond.mean
#   }
#   conditional.expect[i,] = expectile.vec
# }
#
#
#
#
# quantile.vec <- rep(NA,Ny + 1)
# for(y.idx in seq(0,Ny)){
#   print(y.idx)
#   which.idx <- which(y == y.idx)
#   p <- mean(y <= y.idx)
#   q <- quantile(z, p, na.rm = T)
#   quantile.vec[y.idx + 1] = q
# }
#
#
#
# lty=2
#
# plot(seq(0,Ny), quantile.vec, type = "l", col = "blue",xlim = c(0,Ny*4/2), ylim = c(30,Nz), xlab = "y", ylab = "z")
# lines(seq(0,Ny), conditional.expect[1,], type = "l", col = "red", lty = 1)
# lines(seq(0,Ny), conditional.expect[2,], type = "l", col = "red", lty = 2)
# lines(seq(0,Ny), conditional.expect[3,], type = "l", col = "purple", lty = 1)
# lines(seq(0,Ny), conditional.expect[4,], type = "l", col = "green", lty = 1)
# lines(seq(0,Ny), conditional.expect[5,], type = "l", col = "green", lty = 2)
# lines(seq(0,Ny), conditional.expect[6,], type = "l", col = "green", lty = 3)
# lines(seq(0,Ny), conditional.expect[7,], type = "l", col = "green", lty = 4)
# lines(seq(0,Ny), conditional.expect[8,], type = "l", col = "darkgreen", lty = 5)
# legend(Ny, Nz, legend=c("Quantile (OT)", "Conditional (theta = -1)",
#                         "Conditional (theta = -.5)",
#                         "Conditional (theta = 0)",
#                         "Conditional (theta = 0.3)",
#                         "Conditional (theta = 0.5)",
#                         "Conditional (theta = 0.8)",
#                         "Conditional (theta = 0.9)",
#                         "Conditional (theta = 1)"),
#        col=c("blue","red", "red", "purple",
#              rep("green", 4), "darkgreen"), lty=c(1,1,2,1,1:6), cex=0.8)
#
#



