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


# Latent piece wise Linear
latent.pwl = T

n.sims = 5 # TODO: Change this to 500


Ny <- 30 #
Nz <- 30 #

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
# We include the covariate shift in the prediction

rho.grid = 0.8*(x.grid < 61) +  0.2*(x.grid >= 61)

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
# is it something with the amount of smoothing?
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
    if(latent.pwl){
      lat.mu2 <- 0.5 + pmin( -(1/100*(X - 71)),  -(1/5*(X - 71)))
    }
    # lat.mu2 <- (-1/5*(X - 67.5) +  2.5)*sin(X)

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

    gamma1 <- sqrt(gamma2)

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


    colnames(Y.train) = c("y", "age")
    colnames(Z.train) = c("z", "age")



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




if(latent.pwl){
  saveRDS(res, paste0("data/05_pwl_pred",id, ".rds"))
  saveRDS(res.no.reg, paste0("data/05_pwl_pred_no_reg",id, ".rds"))

  saveRDS(res.cov.adj, paste0("data/05_pwl_pred_cov_adj",id, ".rds"))
  saveRDS(res.cov.adj.small.h, paste0("data/05_pwl_pred_cov_adj_small_h",id, ".rds"))
  saveRDS(res.cov.adj.med.h, paste0("data/05_pwl_pred_cov_adj_med_h",id, ".rds"))
  saveRDS(res.cov.adj.large.h, paste0("data/05_pwl_pred_cov_adj_large_h",id, ".rds"))

  saveRDS(res.cov.adj.no.reg, paste0("data/05_pwl_pred_cov_adj_no_reg",id, ".rds"))
  saveRDS(res.cov.adj.small.h.no.reg, paste0("data/05_pwl_pred_cov_adj_small_h_no_reg",id, ".rds"))
  saveRDS(res.cov.adj.med.h.no.reg, paste0("data/05_pwl_pred_cov_adj_med_h_no_reg",id, ".rds"))
  saveRDS(res.cov.adj.large.h.no.reg, paste0("data/05_pwl_pred_cov_adj_large_h_no_reg",id, ".rds"))

  saveRDS(res.z.score, paste0("data/05_pwl_pred_zscore",id, ".rds"))

} else {
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


make.plots = F

#TODO: Finish this set of plots
if(make.plots){
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200


  #TODO: Finish the set of plots
  if(latent.pwl){
    res <- readRDS(paste0("data/05_pwl_pred",1, ".rds"))
    res.no.reg <- readRDS(paste0("data/05_pwl_pred_no_reg",1, ".rds"))

    res.cov.adj <- readRDS(paste0("data/05_pwl_pred_cov_adj",1, ".rds"))
    res.cov.adj.small.h <- readRDS(paste0("data/05_pwl_pred_cov_adj_small_h",1, ".rds"))
    res.cov.adj.med.h <- readRDS(paste0("data/05_pwl_pred_cov_adj_med_h",1, ".rds"))
    res.cov.adj.large.h <- readRDS(paste0("data/05_pwl_pred_cov_adj_large_h",1, ".rds"))

    res.cov.adj.no.reg <- readRDS(paste0("data/05_pwl_pred_cov_adj_no_reg",1, ".rds"))
    res.cov.adj.small.h.no.reg <- readRDS(paste0("data/05_pwl_pred_cov_adj_small_h_no_reg",1, ".rds"))
    res.cov.adj.med.h.no.reg <- readRDS(paste0("data/05_pwl_pred_cov_adj_med_h_no_reg",1, ".rds"))
    res.cov.adj.large.h.no.reg <- readRDS(paste0("data/05_pwl_pred_cov_adj_large_h_no_reg",1, ".rds"))

    res.z.score <- readRDS(paste0("data/05_pwl_pred_zscore",1, ".rds"))
  } else {
    res <- readRDS(paste0("data/05_pred",1, ".rds"))
    res.no.reg <- readRDS(paste0("data/05_pred_no_reg",1, ".rds"))

    res.cov.adj <- readRDS(paste0("data/05_pred_cov_adj",1, ".rds"))
    res.cov.adj.small.h <- readRDS(paste0("data/05_pred_cov_adj_small_h",1, ".rds"))
    res.cov.adj.med.h <- readRDS(paste0("data/05_pred_cov_adj_med_h",1, ".rds"))
    res.cov.adj.large.h <- readRDS(paste0("data/05_pred_cov_adj_large_h",1, ".rds"))

    res.cov.adj.no.reg <- readRDS(paste0("data/05_pred_cov_adj_no_reg",1, ".rds"))
    res.cov.adj.small.h.no.reg <- readRDS(paste0("data/05_pred_cov_adj_small_h_no_reg",1, ".rds"))
    res.cov.adj.med.h.no.reg <- readRDS(paste0("data/05_pred_cov_adj_med_h_no_reg",1, ".rds"))
    res.cov.adj.large.h.no.reg <- readRDS(paste0("data/05_pred_cov_adj_large_h_no_reg",1, ".rds"))

    res.z.score <- readRDS(paste0("data/05_pred_zscore",1, ".rds"))
  }



  for(j in seq(2,200)){

    if(latent.pwl){
      res.tmp <- readRDS(paste0("data/05_pwl_pred",j, ".rds"))
      res.no.reg.tmp <- readRDS(paste0("data/05_pwl_pred_no_reg",j, ".rds"))

      res.cov.adj.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj",j, ".rds"))
      res.cov.adj.small.h.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj_small_h",j, ".rds"))
      res.cov.adj.med.h.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj_med_h",j, ".rds"))
      res.cov.adj.large.h.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj_large_h",j, ".rds"))

      res.cov.adj.no.reg.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj_no_reg",j, ".rds"))
      res.cov.adj.small.h.no.reg.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj_small_h_no_reg",j, ".rds"))
      res.cov.adj.med.h.no.reg.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj_med_h_no_reg",j, ".rds"))
      res.cov.adj.large.h.no.reg.tmp <- readRDS(paste0("data/05_pwl_pred_cov_adj_large_h_no_reg",j, ".rds"))

      res.z.score.tmp <- readRDS(paste0("data/05_pwl_pred_zscore",j, ".rds"))
    } else {
      res.tmp <- readRDS(paste0("data/05_pred",j, ".rds"))
      res.no.reg.tmp <- readRDS(paste0("data/05_pred_no_reg",j, ".rds"))

      res.cov.adj.tmp <- readRDS(paste0("data/05_pred_cov_adj",j, ".rds"))
      res.cov.adj.small.h.tmp <- readRDS(paste0("data/05_pred_cov_adj_small_h",j, ".rds"))
      res.cov.adj.med.h.tmp <- readRDS(paste0("data/05_pred_cov_adj_med_h",j, ".rds"))
      res.cov.adj.large.h.tmp <- readRDS(paste0("data/05_pred_cov_adj_large_h",j, ".rds"))

      res.cov.adj.no.reg.tmp <- readRDS(paste0("data/05_pred_cov_adj_no_reg",j, ".rds"))
      res.cov.adj.small.h.no.reg.tmp <- readRDS(paste0("data/05_pred_cov_adj_small_h_no_reg",j, ".rds"))
      res.cov.adj.med.h.no.reg.tmp <- readRDS(paste0("data/05_pred_cov_adj_med_h_no_reg",j, ".rds"))
      res.cov.adj.large.h.no.reg.tmp <- readRDS(paste0("data/05_pred_cov_adj_large_h_no_reg",j, ".rds"))

      res.z.score.tmp <- readRDS(paste0("data/05_pred_zscore",j, ".rds"))
    }


    res <- rbind(res, res.tmp)
    res.no.reg <- rbind(res.no.reg, res.no.reg.tmp)

    res.cov.adj <- rbind(res.cov.adj, res.cov.adj.tmp)
    res.cov.adj.small.h <- rbind(res.cov.adj.small.h, res.cov.adj.small.h.tmp)
    res.cov.adj.med.h <- rbind(res.cov.adj.med.h, res.cov.adj.med.h.tmp)
    res.cov.adj.large.h <- rbind(res.cov.adj.large.h, res.cov.adj.large.h.tmp)

    res.cov.adj.no.reg <- rbind(res.cov.adj.no.reg, res.cov.adj.no.reg.tmp)
    res.cov.adj.small.h.no.reg <- rbind(res.cov.adj.small.h.no.reg, res.cov.adj.small.h.no.reg.tmp)
    res.cov.adj.med.h.no.reg <- rbind(res.cov.adj.med.h.no.reg, res.cov.adj.med.h.no.reg.tmp)
    res.cov.adj.large.h.no.reg <- rbind(res.cov.adj.large.h.no.reg, res.cov.adj.large.h.no.reg.tmp)

    res.z.score <- rbind(res.z.score, res.z.score.tmp)

  }

  n.sims = nrow(res.z.score)
  n.set <- round(c(100,200,500,1000, 2000, 5000)/13)*13 # smoothed it over

  J = length(n.set)


  data.list = list(res,
                   res.no.reg,
                   res.cov.adj,
                   res.cov.adj.small.h,
                   res.cov.adj.med.h,
                   res.cov.adj.large.h,
                   res.cov.adj.no.reg,
                   res.cov.adj.small.h.no.reg,
                   res.cov.adj.med.h.no.reg,
                   res.cov.adj.large.h.no.reg,
                   res.z.score)

  mean.vec <- c()
  sd.vec <- c()

  for(i in seq(length(data.list))){
      data = data.list[[i]]
      mean.vec <- c(mean.vec, colMeans(data))
      sd.vec <- c(sd.vec, colSDs(data))
  }
  se.vec <- sd.vec/sqrt(n.sims)

  method.vec <- c(rep("No Covariates", J),
                  rep("No Covariates", J), #unreg
                  rep("Binomial", J),
                  rep("Gaussian h = 1", J),
                  rep("Gaussian h = 3", J),
                  rep("Gaussian h = 9", J),
                  rep("Binomial", J), #unreg
                  rep("Gaussian h = 1", J), #unreg
                  rep("Gaussian h = 3", J), #unreg
                  rep("Gaussian h = 9", J), #unreg
                  rep("Z score Matching", J)
                  )
  reg.vec <- c(rep("Yes", J),
               rep("No", J),
               rep("Yes", J),
               rep("Yes", J),
               rep("Yes", J),
               rep("Yes", J),
               rep("No", J),
               rep("No", J),
               rep("No", J),
               rep("No", J),
               rep("No", J)
               )

  sample.size.vec <- rep(n.set, times = length(data.list))




  res.data <- data.frame("method" = method.vec,
                         "regularization"= reg.vec,
                         "SampleSize" = sample.size.vec,
                         "ce" = -mean.vec,
                         "ce_se" = se.vec)

  res.data.no.z <- res.data %>% filter(method != "Z score Matching")
  res.data.z  <- res.data %>% filter(method == "Z score Matching")
  res.data.z$ce <- as.numeric(res.data.z$ce)
  #res.data.z$SampleSize


  plt.predictions <- ggplot(res.data, aes(x = SampleSize, y = ce,
                                          group = interaction(method, regularization),
                                          color = method, linetype = regularization)) +
    geom_line() +
    #geom_point() +
    geom_errorbar(aes(ymin = ce - 2*ce_se, ymax = ce + 2*ce_se)) +
    #geom_line(data = res.data.z, aes(x = SampleSize, y = ce, group = method), color = "black", linetype = 4) +
    #geom_line(data = res.data.z, aes(x = SampleSize, y = ce,linetype = 3, group = method, color = "black")) +
    #geom_point(data = res.data.z, aes(x = SampleSize, y = ce,group = method, color = "black")) +
    #geom_errorbar(data = res.data.z, aes(x = SampleSize, ymin = ce - 2*ce_se, ymax = ce + 2*ce_se)) +
    ggtitle("Cross Entropy Of Competing Methods") +
    xlab("Sample Size (n)") +
    ylab("Cross Entropy") +
    scale_color_manual(values = c("No Covariates" = "orange",
                                  "Binomial" = "olivedrab4",
                                  "Gaussian h = 1" = "deepskyblue1",
                                  "Gaussian h = 3" = "deepskyblue3",
                                  "Gaussian h = 9" = "deepskyblue4",
                                  "Z score Matching" = "black")) +
    coord_cartesian(
      xlim =c(0,5000),
      ylim = c(2.5,6.5)
    ) # can change the window for smoother plots


  #geom_line(aes(x = n, y = rmse, color = method)) #+
  #

  plt.predictions

  if(latent.pwl){
    png(filename = "plots/prediction_pwl_sim_cross_entropy_full.png",
        width = png.width, height = png.height, res = png.res)

    plt.predictions
    # Close the pdf file
    dev.off()
  } else {
    png(filename = "plots/prediction_sim_cross_entropy_full.png",
        width = png.width, height = png.height, res = png.res)

    plt.predictions
    # Close the pdf file
    dev.off()
  }




  plt.predictions <- ggplot(res.data, aes(x = SampleSize, y = ce,
                                          group = interaction(method, regularization),
                                          color = method, linetype = regularization)) +
    geom_line() +
    #geom_point() +
    geom_errorbar(aes(ymin = ce - 2*ce_se, ymax = ce + 2*ce_se)) +
    #geom_line(data = res.data.z, aes(x = SampleSize, y = ce, group = method), color = "black", linetype = 4) +
    #geom_line(data = res.data.z, aes(x = SampleSize, y = ce,linetype = 3, group = method, color = "black")) +
    #geom_point(data = res.data.z, aes(x = SampleSize, y = ce,group = method, color = "black")) +
    #geom_errorbar(data = res.data.z, aes(x = SampleSize, ymin = ce - 2*ce_se, ymax = ce + 2*ce_se)) +
    ggtitle("Cross Entropy Of Competing Methods") +
    xlab("Sample Size (n)") +
    ylab("Cross Entropy") +
    scale_color_manual(values = c("No Covariates" = "orange",
                                  "Binomial" = "olivedrab4",
                                  "Gaussian h = 1" = "deepskyblue1",
                                  "Gaussian h = 3" = "deepskyblue3",
                                  "Gaussian h = 9" = "deepskyblue4",
                                  "Z score Matching" = "black")) +
    coord_cartesian(
      xlim =c(0,5000),
      ylim = c(2.5,3.5)
    ) # can change the window for smoother plots


  #geom_line(aes(x = n, y = rmse, color = method)) #+
  #

  plt.predictions

  if(latent.pwl){
    png(filename = "plots/prediction_pwl_sim_cross_entropy_zoom.png",
        width = png.width, height = png.height, res = png.res)

    plt.predictions
    # Close the pdf file
    dev.off()
  } else {
    png(filename = "plots/prediction_sim_cross_entropy_zoom.png",
        width = png.width, height = png.height, res = png.res)

    plt.predictions
    # Close the pdf file
    dev.off()
  }



}







