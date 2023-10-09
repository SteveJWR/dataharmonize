

#In this script, we illustrate the consistency of the model selection procedure.
rm(list = ls())
library(dnoiseR)

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}

set.seed(id)
if(id %% 2 == 1){
  kernel="Gaussian"
  ker.true <- gaussian_kernel
} else {
  kernel="Exponential"
  ker.true <- exponential_kernel
}

R.bins = 1000

# True Conditional Model
h.true <- 2


N <- 30

cond.true <- conditional_mkm(N,ker.true, h.true)


two.obs.ratio = 1
n.seq = c(100,500,1000,5000)

J = length(n.seq)


# grid for the values of mu
L = 200
mu.set <- seq(0,2,length.out = L)

# folds for the cross validation
folds <- 5

# correct model selection fraction of sims.
n.sims <- 50
# dimensions of the results (sim, samplesize, tuning parameter)
res <- array(NA,c(n.sims,length(n.seq),L))
res.lat.dist <- array(NA,c(n.sims,length(n.seq),L))

#simulation parameters
alpha1 <- 1.2
alpha2 <- 3


x.grid = seq(0,1,length.out = R.bins)

latent.true = (1/3)*dbeta(x.grid,alpha1,alpha2) + (2/3)*dbeta(x.grid,alpha2,alpha1)

for(j in seq(J)){
  n = n.seq[j]
  for(sim in seq(n.sims)){
    cat(paste0("Simulation ", sim, "/",n.sims), end = "\r")
    gamma <- c(rbeta((1/3)*(n),alpha1,alpha2),rbeta((2/3)*(n),alpha2,alpha1))
    gamma <- gamma[sample(seq(length(gamma)))] # shuffle


    Y <- simulate_test_cond(obs.set = c(rep(1,length(gamma))),cond.true ,gamma)
    model <- mu_selection(mu.set,cond.true,Y,R.bins,folds, verbose = F, latent.true = latent.true)
    res[sim,j,] <- model$cv.lik
    res.lat.dist[sim,j,] <- model$cv.lat.dist
  }
}


saveRDS(res, paste0("data/regularization_selection_results_",kernel,ceiling(id/2), ".rds"))
saveRDS(res.lat.dist, paste0("data/regularization_latent_distance_results_",kernel,ceiling(id/2), ".rds"))






make.plots = F

if(make.plots){
  library(ggpubr)
  library(abind)

  kernel = "Gaussian" # "Gaussian" # "Exponential"
  png.width = 1200
  png.height = 1000
  png.res = 200

  # grid.parameters
  n.seq = c(100,500,1000,5000)
  J = length(n.seq)
  # grid for the values of mu
  L = 200
  mu.set <- seq(0,2,length.out = L)

  # TODO: update this section to concatenate the results
  res <- readRDS(paste0("data/regularization_selection_results_",kernel,1, ".rds"))

  for(j in seq(2,200)){
    file.name <- paste0("data/regularization_selection_results_",kernel,j, ".rds")
    if(file.exists(file.name)){
      res.tmp <- readRDS(file.name)
      res <- abind(res, res.tmp, along = 1)
    }
  }

  n.sims = dim(res)[1]

  # subdivide by sample size
  lik.mean = matrix(NA, nrow = L, ncol = J)
  lik.sd = matrix(NA, nrow = L, ncol = J)
  for(j in seq(J)){
    res.tmp = as.matrix(res[,j,])
    lik.mean[,j]= colMeans(res.tmp, na.rm = T)
    lik.sd[,j] = colSDs(res.tmp, na.rm = T)
  }
  lik.sd = lik.sd/sqrt(n.sims) # standard error of the mean function.




  lik.mean.vec = c(lik.mean[,1]/(-max(lik.mean[,1])),
                   lik.mean[,2]/(-max(lik.mean[,2])),
                   lik.mean[,3]/(-max(lik.mean[,3])),
                   lik.mean[,4])/(-max(lik.mean[,4]))
  block.1 <- (lik.mean[,1] - min(lik.mean[,1]))/(max(lik.mean[,1])-min(lik.mean[,1]))
  block.2 <- (lik.mean[,2] - min(lik.mean[,2]))/(max(lik.mean[,2])-min(lik.mean[,2]))
  block.3 <- (lik.mean[,3] - min(lik.mean[,3]))/(max(lik.mean[,3])-min(lik.mean[,3]))
  block.4 <- (lik.mean[,4] - min(lik.mean[,4]))/(max(lik.mean[,4])-min(lik.mean[,4]))

  lik.mean.scaled.vec = c(block.1,block.2,block.3,block.4)
  block.1 <- (lik.sd[,1])/(max(lik.mean[,1])-min(lik.mean[,1]))
  block.2 <- (lik.sd[,2])/(max(lik.mean[,2])-min(lik.mean[,2]))
  block.3 <- (lik.sd[,3])/(max(lik.mean[,3])-min(lik.mean[,3]))
  block.4 <- (lik.sd[,4])/(max(lik.mean[,4])-min(lik.mean[,4]))

  lik.sd.scaled.vec = c(block.1,block.2,block.3,block.4)

  res.data <- data.frame("SampleSize" = c(rep(paste0("n = ", n.seq[1]), L),
                                          rep(paste0("n = ", n.seq[2]), L),
                                          rep(paste0("n = ", n.seq[3]), L),
                                          rep(paste0("n = ", n.seq[4]), L)),
                         "mu" = c(mu.set,mu.set,mu.set,mu.set),
                         "lik.mean" = c(lik.mean[,1], lik.mean[,2], lik.mean[,3], lik.mean[,4]),
                         "lik.sd" = c(lik.sd[,1], lik.sd[,2], lik.sd[,3], lik.sd[,4]),
                         "lik.mean.scaled" = lik.mean.scaled.vec,
                         "lik.sd.scaled" = lik.sd.scaled.vec)

  res.data$SampleSize <- factor(res.data$SampleSize, levels = paste0("n = ", n.seq))

  plt.mu <- ggplot(res.data, aes(x = mu, y = lik.mean.scaled, color = SampleSize)) +
    geom_line() +
    geom_ribbon(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled, fill = SampleSize),alpha=0.3) +
    ggtitle(paste0(kernel, " Kernel Cross Validation")) +
    xlab("Mu") +
    ylab("Normalized CV Likelihood")
    #geom_errorbar(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled))

  plt.mu

  png(filename = paste0("plots/regularization_cv_sim",kernel,".png"),
      width = png.width, height = png.height, res = png.res)

  plt.mu
  # Close the pdf file
  dev.off()
}




if(make.plots){
  library(ggpubr)
  library(abind)

  kernel = "Gaussian" # "Gaussian" # "Exponential"
  png.width = 1200
  png.height = 1000
  png.res = 200

  # grid.parameters
  n.seq = c(100,500,1000,5000)
  J = length(n.seq)
  # grid for the values of mu
  L = 200
  mu.set <- seq(0,2,length.out = L)

  # TODO: update this section to concatenate the results
  res <- readRDS(paste0("data/regularization_latent_distance_results_",kernel,1, ".rds"))

  for(j in seq(2,200)){
    file.name <- paste0("data/regularization_latent_distance_results_",kernel,j, ".rds")
    if(file.exists(file.name)){
      res.tmp <- readRDS(file.name)
      res <- abind(res, res.tmp, along = 1)
    }
  }

  n.sims = dim(res)[1]

  # subdivide by sample size
  lat.dist.mean = matrix(NA, nrow = L, ncol = J)
  lat.dist.sd = matrix(NA, nrow = L, ncol = J)
  for(j in seq(J)){
    res.tmp = as.matrix(res[,j,])
    lat.dist.mean[,j]= colMeans(res.tmp, na.rm = T)
    lat.dist.sd[,j] = colSDs(res.tmp, na.rm = T)
  }
  lat.dist.sd = lat.dist.sd/sqrt(n.sims) # standard error of the mean function.




  lat.dist.mean.vec = c(lat.dist.mean[,1]/(-max(lat.dist.mean[,1])),
                        lat.dist.mean[,2]/(-max(lat.dist.mean[,2])),
                        lat.dist.mean[,3]/(-max(lat.dist.mean[,3])),
                        lat.dist.mean[,4])/(-max(lat.dist.mean[,4]))
  block.1 <- (lat.dist.mean[,1] - min(lat.dist.mean[,1]))/(max(lat.dist.mean[,1])-min(lat.dist.mean[,1]))
  block.2 <- (lat.dist.mean[,2] - min(lat.dist.mean[,2]))/(max(lat.dist.mean[,2])-min(lat.dist.mean[,2]))
  block.3 <- (lat.dist.mean[,3] - min(lat.dist.mean[,3]))/(max(lat.dist.mean[,3])-min(lat.dist.mean[,3]))
  block.4 <- (lat.dist.mean[,4] - min(lat.dist.mean[,4]))/(max(lat.dist.mean[,4])-min(lat.dist.mean[,4]))

  lat.dist.mean.scaled.vec = c(block.1,block.2,block.3,block.4)
  block.1 <- (lat.dist.sd[,1])/(max(lat.dist.mean[,1])-min(lat.dist.mean[,1]))
  block.2 <- (lat.dist.sd[,2])/(max(lat.dist.mean[,2])-min(lat.dist.mean[,2]))
  block.3 <- (lat.dist.sd[,3])/(max(lat.dist.mean[,3])-min(lat.dist.mean[,3]))
  block.4 <- (lat.dist.sd[,4])/(max(lat.dist.mean[,4])-min(lat.dist.mean[,4]))

  lat.dist.sd.scaled.vec = c(block.1,block.2,block.3,block.4)

  res.data <- data.frame("SampleSize" = c(rep(paste0("n = ", n.seq[1]), L),
                                          rep(paste0("n = ", n.seq[2]), L),
                                          rep(paste0("n = ", n.seq[3]), L),
                                          rep(paste0("n = ", n.seq[4]), L)),
                         "mu" = c(mu.set,mu.set,mu.set,mu.set),
                         "lat.dist.mean" = c(lat.dist.mean[,1], lat.dist.mean[,2], lat.dist.mean[,3], lat.dist.mean[,4]),
                         "lat.dist.sd" = c(lat.dist.sd[,1], lat.dist.sd[,2], lat.dist.sd[,3], lat.dist.sd[,4]),
                         "lat.dist.mean.scaled" = lat.dist.mean.scaled.vec,
                         "lat.dist.sd.scaled" = lat.dist.sd.scaled.vec)

  res.data$SampleSize <- factor(res.data$SampleSize, levels = paste0("n = ", n.seq))

  plt.mu <- ggplot(res.data, aes(x = mu, y = lat.dist.mean.scaled, color = SampleSize)) +
    geom_line() +
    geom_ribbon(aes(ymin = lat.dist.mean.scaled - 2*lat.dist.sd.scaled, ymax = lat.dist.mean.scaled + 2*lat.dist.sd.scaled, fill = SampleSize),alpha=0.3) +
    ggtitle(paste0(kernel, " Kernel Cross Validation")) +
    xlab("Mu") +
    ylab("Normalized CV Likelihood")
  #geom_errorbar(aes(ymin = lat.dist.mean.scaled - 2*lat.dist.sd.scaled, ymax = lat.dist.mean.scaled + 2*lat.dist.sd.scaled))

  plt.mu

  png(filename = paste0("plots/regularization_cv_sim",kernel,".png"),
      width = png.width, height = png.height, res = png.res)

  plt.mu
  # Close the pdf file
  dev.off()
}







