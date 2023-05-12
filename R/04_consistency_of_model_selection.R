

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




two.obs.ratio = 1
n1.seq = c(100,500,1000,5000)
n2.seq =  two.obs.ratio*n1.seq
J = length(n1.seq)


# grid for the values of h
h.set <- c(0.8,1,2,3,5,10)
H = length(h.set)

if(kernel == "Gaussian"){
  cond.set <- generate_mkm_list(N = N, ker = gaussian_kernel, h.set = h.set)
  cond.names <- paste0("Gaussian h = ",as.character(h.set))
} else {
  cond.set <- generate_mkm_list(N = N, ker = exponential_kernel, h.set = h.set)
  cond.names <- paste0("Exponential h = ",as.character(h.set))
}



# correct model selection fraction of sims.
n.sims <- 5
res <- matrix(NA, nrow= n.sims, ncol =length(n1.seq))
res <- array(NA, c(n.sims, length(n1.seq), H))

#simulation parameters
alpha1 <- 1.2
alpha2 <- 3
for(h in seq(H)){
  i.true = h
  h.true = h.set[i.true]
  cond.true <- conditional_mkm(N,ker.true, h.true)
  for(j in seq(J)){
    n1 = n1.seq[j]
    n2 = n2.seq[j]


    for(sim in seq(n.sims)){
      cat(paste0("Simulation ", sim, "/",n.sims), end = "\r")
      gamma <- c(rbeta((1/3)*(n1 + n2),alpha1,alpha2),rbeta((2/3)*(n1 + n2),alpha2,alpha1))
      gamma <- gamma[sample(seq(length(gamma)))] # shuffle


      Y <- simulate_test_cond(obs.set = c(rep(1,n1), rep(2,n2)),cond.true ,gamma)
      model <- error_model_selection_bivariate(cond.set,Y,R.bins,cond.names)
      i.max <- which.max(model$lik_vals)
      res[sim,j, h] <- TRUE*(i.max == i.true) + FALSE*(i.max != i.true)
    }
  }
}


saveRDS(res, paste0("data/model_selection_results_",kernel,ceiling(id/2), ".rds"))



make.plots = F

if(make.plots){
  library(ggpubr)
  library(abind)

  kernel = "Gaussian" # "Gaussian" # "Exponential"
  png.width = 1200
  png.height = 1000
  png.res = 200

  kernel = "Gaussian" #"Exponential", "Gaussian"
  # grid.parameters
  n.seq = c(100,500,1000,5000)
  J = length(n.seq)


  # TODO: update this section to concatenate the results
  res <- readRDS(paste0("data/model_selection_results_",kernel,1, ".rds"))

  for(j in seq(2,100)){
    file.name <- paste0("data/model_selection_results_",kernel,j, ".rds")
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











