rm(list = ls())

# TODO: Change this format for the imputation

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



n.sims = 5 # TODO: Change this to 500


Ny <- 35
Nz <- 35

R.bins = 1000

cond.y <- generate_cond_binomial(Ny)
cond.z <- generate_cond_binomial(Nz)

# conditional distributions with misspecified models
#generate_mkm_list()

h.small = 1
h.med = 3
h.large = 9

ker <- gaussian_kernel

cond.y.small.h <- conditional_mkm(Ny,ker,h.small)
cond.z.small.h <- conditional_mkm(Nz,ker,h.small)

cond.y.med.h <- conditional_mkm(Ny,ker,h.med)
cond.z.med.h <- conditional_mkm(Nz,ker,h.med)

cond.y.large.h <- conditional_mkm(Ny,ker,h.large)
cond.z.large.h <- conditional_mkm(Nz,ker,h.large)


# h1 <- 2
# h2 <- 0.5
# cond.y <- conditional_mkm(Ny, gaussian_kernel, h1)
# cond.z <- conditional_mkm(Nz, gaussian_kernel, h2)

#rho.grid = c(0.3,0.3,0.7,0.7)


#n.true = 10000
# uniform age sample
#X <- round(runif(n.true, min = 54.5, max = 80.5))
#x.grid = c(55,63,72,80)

#X <- rep(x.grid, n.true/4)
#X <- rep(x.grid, round(n.true/length(x.grid)))

#
x.grid = seq(55,80)

n.set <- c(100,200,500,1000, 2000, 5000)
n.set <- round(c(100,200,500,1000, 2000, 5000)/13)*13 # smoothed it over

# missing at random piece
rho.grid = 0.8*(x.grid < 61) +  0.2*(x.grid >= 61)
rho.grid = rep(0.5, length(x.grid))

#X <- rep(c(55,80), n.true/2)
#lat.mu <- (-1/10*(X - 67.5)) +  0.5
#lat.mu <- -(1/10*(X - 67.5))^2 +  0.5


#unif.1 <- scale_kernel(uniform_kernel,1)
gaussian_kernel.2 <- scale_kernel(gaussian_kernel,2)
ref.cols <- "age"
ker.set <- list(gaussian_kernel.2)



# results

res <- matrix(NA, nrow = n.sims, ncol = length(n.set))
res.no.reg <- matrix(NA, nrow = n.sims, ncol = length(n.set))

res.cov.adj <- matrix(NA, nrow = n.sims, ncol = length(n.set))
res.cov.adj.small.h <- matrix(NA, nrow = n.sims, ncol = length(n.set))
res.cov.adj.med.h <- matrix(NA, nrow = n.sims, ncol = length(n.set))
res.cov.adj.large.h <- matrix(NA, nrow = n.sims, ncol = length(n.set))

res.cov.adj.no.reg <- matrix(NA, nrow = n.sims, ncol = length(n.set))
res.cov.adj.small.h.no.reg <- matrix(NA, nrow = n.sims, ncol = length(n.set))
res.cov.adj.med.h.no.reg <- matrix(NA, nrow = n.sims, ncol = length(n.set))
res.cov.adj.large.h.no.reg <- matrix(NA, nrow = n.sims, ncol = length(n.set))

res.z.score <- matrix(NA, nrow = n.sims, ncol = length(n.set))


load.prior.results = F
if(load.prior.results){
  res  <- readRDS("data/05_pred.rds")
  res.no.reg  <- readRDS("data/05_pred_no_reg.rds")

  res.cov.adj <- readRDS("data/05_pred_cov_adj.rds")
  res.cov.adj.small.h <- readRDS("data/05_pred_cov_adj_small_h.rds")
  res.cov.adj.med.h <- readRDS("data/05_pred_cov_adj_med_h.rds")
  res.cov.adj.large.h <- readRDS("data/05_pred_cov_adj_large_h.rds")

  res.cov.adj.no.reg <- readRDS("data/05_pred_cov_adj_no_reg.rds")
  res.cov.adj.small.h.no.reg <- readRDS("data/05_pred_cov_adj_small_h_no_reg.rds")
  res.cov.adj.med.h.no.reg <- readRDS("data/05_pred_cov_adj_med_h_no_reg.rds")
  res.cov.adj.large.h.no.reg <- readRDS("data/05_pred_cov_adj_large_h_no_reg.rds")

  res.z.score <- readRDS("data/05_pred_zscore.rds")
}


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


# set mu to reduce as a function of the sample
mu.y = 100 # reduce /n
mu.z = 100 # reduce /n


set.seed(id) # finding the correct id's for the simulations.
sim.start = 1

for(sim in sim.start:n.sims){
  cat(paste0("Number of sims: ", sim, "/",n.sims), end = "\n")
  if(sim %% 20 == 0){


    saveRDS(res, paste0("data/05_pred",id, ".rds"))
    saveRDS(res.no.reg, paste0("data/05_pred_no_reg",id, ".rds"))

    saveRDS(res.cov.adj, paste0("data/05_pred_cov_adj",id, ".rds"))
    saveRDS(res.cov.adj.small.h, paste0("data/05_pred_cov_adj_small_h",id, ".rds"))
    saveRDS(res.cov.adj.med.h, paste0("data/05_pred_cov_adj_med_h",id, ".rds"))
    saveRDS(res.cov.adj.large.h, paste0("data/05_pred_cov_adj_large_h",id, ".rds"))

    saveRDS(res.cov.adj.no.reg, paste0("data/05_pred_cov_adj_no_reg",id, ".rds"))
    saveRDS(res.cov.adj.small.h.no.reg, paste0("data/05_pred_cov_adj_small_h_no_reg",id, ".rds"))
    saveRDS(res.cov.adj.med.h.no.reg, paste0("data/05_pred_cov_adj_med_h_no_reg",id, ".rds"))
    saveRDS(res.cov.adj.large.h.no.reg, paste0("data/05_pred_cov_adj_large_h_no_reg",id, ".rds"))

    saveRDS(res.z.score, paste0("data/05_pred_zscore",id, ".rds"))
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


    #
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
    z.test <- sim.data$Z[idx1]
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

    #TODO: add all the other kernel and regularized versions

    mean.y <- mean(Y.train[,1])
    mean.z <- mean(Z.train[,1])

    sd.y <- sd(Y.train[,1])
    sd.z <- sd(Z.train[,1])



    naive.pred <- normalScoreConversionProb(sim.data$Y[idx1],mean.y,sd.y,mean.z,sd.z, Nz)

    cond.pred <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                        Y.train,Z.train,cond.y,cond.z,
                                        mu.y/n,mu.z/n,R.bins = 1000)

    cond.pred.no.reg <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                               Y.train,Z.train,cond.y,cond.z,
                                               0,0, R.bins = 1000)




    cond.pred.cov.adj <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                Y.train,Z.train,cond.y,cond.z,
                                                mu.y/n,mu.z/n,ref.cols, ker.set, R.bins = 1000)

    cond.pred.cov.adj.no.reg <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                       Y.train,Z.train,cond.y,cond.z,
                                                       0,0,ref.cols, ker.set, R.bins = 1000)





    cond.pred.small.h <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                Y.train,Z.train,cond.y.small.h,cond.z.small.h,
                                                mu.y/n,mu.z/n,ref.cols, ker.set, R.bins = 1000)

    cond.pred.small.h.no.reg <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                       Y.train,Z.train,cond.y.small.h,cond.z.small.h,
                                                       0,0,ref.cols, ker.set, R.bins = 1000)





    cond.pred.med.h <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                Y.train,Z.train,cond.y.med.h,cond.z.med.h,
                                                mu.y/n,mu.z/n,ref.cols, ker.set, R.bins = 1000)

    cond.pred.med.h.no.reg <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                       Y.train,Z.train,cond.y.med.h,cond.z.med.h,
                                                       0,0,ref.cols, ker.set, R.bins = 1000)




    cond.pred.large.h <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                Y.train,Z.train,cond.y.large.h,cond.z.large.h,
                                                mu.y/n,mu.z/n,ref.cols, ker.set, R.bins = 1000)

    cond.pred.large.h.no.reg <- predictedDistributions(sim.data[idx1,],sim.data$Y[idx1],
                                                       Y.train,Z.train,cond.y.large.h,cond.z.large.h,
                                                       0,0,ref.cols, ker.set, R.bins = 1000)




    res.z.score[sim,i] <- empirical_cross_entropy(naive.pred,z.test)

    res[sim,i]  <- empirical_cross_entropy(cond.pred,z.test)
    res.no.reg[sim,i]  <- empirical_cross_entropy(cond.pred.no.reg,z.test)

    res.cov.adj[sim,i] <- empirical_cross_entropy(cond.pred.cov.adj,z.test)
    res.cov.adj.no.reg[sim,i] <- empirical_cross_entropy(cond.pred.cov.adj.no.reg,z.test)

    res.cov.adj.small.h[sim,i] <- empirical_cross_entropy(cond.pred.small.h,z.test)
    res.cov.adj.small.h.no.reg[sim,i] <- empirical_cross_entropy(cond.pred.small.h.no.reg,z.test)

    res.cov.adj.med.h[sim,i] <- empirical_cross_entropy(cond.pred.med.h,z.test)
    res.cov.adj.med.h.no.reg[sim,i] <- empirical_cross_entropy(cond.pred.med.h.no.reg,z.test)

    res.cov.adj.large.h[sim,i] <- empirical_cross_entropy(cond.pred.large.h,z.test)
    res.cov.adj.large.h.no.reg[sim,i] <- empirical_cross_entropy(cond.pred.large.h.no.reg,z.test)


  }
}




########

saveRDS(res, paste0("data/05_pred",id, ".rds"))
saveRDS(res.no.reg, paste0("data/05_pred_no_reg",id, ".rds"))

saveRDS(res.cov.adj, paste0("data/05_pred_cov_adj",id, ".rds"))
saveRDS(res.cov.adj.small.h, paste0("data/05_pred_cov_adj_small_h",id, ".rds"))
saveRDS(res.cov.adj.med.h, paste0("data/05_pred_cov_adj_med_h",id, ".rds"))
saveRDS(res.cov.adj.large.h, paste0("data/05_pred_cov_adj_large_h",id, ".rds"))

saveRDS(res.cov.adj.no.reg, paste0("data/05_pred_cov_adj_no_reg",id, ".rds"))
saveRDS(res.cov.adj.small.h.no.reg, paste0("data/05_pred_cov_adj_small_h_no_reg",id, ".rds"))
saveRDS(res.cov.adj.med.h.no.reg, paste0("data/05_pred_cov_adj_med_h_no_reg",id, ".rds"))
saveRDS(res.cov.adj.large.h.no.reg, paste0("data/05_pred_cov_adj_large_h_no_reg",id, ".rds"))

saveRDS(res.z.score, paste0("data/05_pred_zscore",id, ".rds"))



make.plots = F

#TODO: Finish this set of plots
if(make.plots){
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200

  # update this section to concatenate the results
  dev.cc <- readRDS("data/06_dev_cc1.rds")
  dev <- readRDS("data/06_dev1.rds")
  dev.cov.adj <- readRDS("data/06_dev_cov_adj1.rds")
  dev.true.latent <-  readRDS("data/06_dev_true_latent1.rds")
  dev.z.score <- readRDS("data/06_dev_zscore1.rds")
  dev.quantile <- readRDS("data/06_dev_quantile1.rds")
  dev.bootstrap <- readRDS("data/06_dev_bootstrap1.rds")

  cover.cc <- readRDS("data/06_cover_cc1.rds")
  cover <- readRDS("data/06_cover1.rds")
  cover.cov.adj <- readRDS("data/06_cover_cov_adj1.rds")
  cover.true.latent <-  readRDS("data/06_cover_true_latent1.rds")
  cover.z.score <- readRDS("data/06_cover_zscore1.rds")
  cover.quantile <- readRDS("data/06_cover_quantile1.rds")
  cover.bootstrap <- readRDS("data/06_cover_bootstrap1.rds")

  #TODO: Finish the set of plots
  saveRDS(res, paste0("data/05_pred",id, ".rds"))
  saveRDS(res.no.reg, paste0("data/05_pred_no_reg",id, ".rds"))

  saveRDS(res.cov.adj, paste0("data/05_pred_cov_adj",id, ".rds"))
  saveRDS(res.cov.adj.small.h, paste0("data/05_pred_cov_adj_small_h",id, ".rds"))
  saveRDS(res.cov.adj.med.h, paste0("data/05_pred_cov_adj_med_h",id, ".rds"))
  saveRDS(res.cov.adj.large.h, paste0("data/05_pred_cov_adj_large_h",id, ".rds"))

  saveRDS(res.cov.adj.no.reg, paste0("data/05_pred_cov_adj_no_reg",id, ".rds"))
  saveRDS(res.cov.adj.small.h.no.reg, paste0("data/05_pred_cov_adj_small_h_no_reg",id, ".rds"))
  saveRDS(res.cov.adj.med.h.no.reg, paste0("data/05_pred_cov_adj_med_h_no_reg",id, ".rds"))
  saveRDS(res.cov.adj.large.h.no.reg, paste0("data/05_pred_cov_adj_large_h_no_reg",id, ".rds"))

  saveRDS(res.z.score, paste0("data/05_pred_zscore",id, ".rds"))
  for(j in seq(2,200)){

    dev.cc.tmp <- readRDS(paste0("data/06_dev_cc",j,".rds"))
    dev.tmp <- readRDS(paste0("data/06_dev",j,".rds"))
    dev.cov.adj.tmp <- readRDS(paste0("data/06_dev_cov_adj",j,".rds"))
    dev.true.latent.tmp <- readRDS(paste0("data/06_dev_true_latent",j,".rds"))
    dev.z.score.tmp <- readRDS(paste0("data/06_dev_zscore",j,".rds"))
    dev.quantile.tmp <- readRDS(paste0("data/06_dev_quantile",j,".rds"))
    dev.bootstrap.tmp <- readRDS(paste0("data/06_dev_bootstrap",j,".rds"))

    cover.cc.tmp <- readRDS(paste0("data/06_cover_cc",j,".rds"))
    cover.tmp <- readRDS(paste0("data/06_cover",j,".rds"))
    cover.cov.adj.tmp <- readRDS(paste0("data/06_cover_cov_adj",j,".rds"))
    cover.true.latent.tmp <- readRDS(paste0("data/06_cover_true_latent",j,".rds"))
    cover.z.score.tmp <- readRDS(paste0("data/06_cover_zscore",j,".rds"))
    cover.quantile.tmp <- readRDS(paste0("data/06_cover_quantile",j,".rds"))
    cover.bootstrap.tmp <- readRDS(paste0("data/06_cover_bootstrap",j,".rds"))

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
                                      rep("DNOISE (cov.adj.)", length(n.set)),
                                      rep("DNOISE (T.L.)", length(n.set)),
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
    geom_line()  #+
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
    geom_errorbar(aes(ymin = log(rmse - 2*rmse_sd), ymax = log(rmse + 2*rmse_sd)))

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
                                      rep("DNOISE (cov. adj.)", length(n.set)),
                                      rep("DNOISE (T.L.)", length(n.set)),
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
    geom_hline(yintercept=0.95, linetype='dotted', col = 'black')

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



