

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
#h.true <- 2
N <- 30


two.obs.ratio = 1
one.obs.ratio = 1
n.seq = c(100,500,1000,5000,10000)
n1.seq = one.obs.ratio*n.seq
n2.seq =  two.obs.ratio*n.seq
J = length(n1.seq)


# grid for the values of h
h.set <- c(0.8,1,2,3,5,10)
h.set <- c(seq(1,6,length.out = 6))
H = length(h.set)

if(kernel == "Gaussian"){
  cond.set <- generate_mkm_list(N = N, ker = gaussian_kernel, h.set = h.set)
  cond.names <- paste0("Gaussian h = ",as.character(h.set))
} else {
  cond.set <- generate_mkm_list(N = N, ker = exponential_kernel, h.set = h.set)
  cond.names <- paste0("Exponential h = ",as.character(h.set))
}



# correct model selection fraction of sims.
n.sims <- 2
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
      res[sim,j, h] <- abs(h.set[i.max]  - h.true)
    }
  }
}


saveRDS(res, paste0("data/model_selection_results_",kernel,ceiling(id/2), ".rds"))



make.plots = F

if(make.plots){
  library(ggpubr)
  library(abind)

  png.width = 1200
  png.height = 1000
  png.res = 200

  kernel = "Exponential" #"Exponential", "Gaussian"
  # grid.parameters
  n.seq = c(100,500,1000,5000,10000)
  J = length(n.seq)


  # grid for the values of h
  h.set <- c(seq(1,6,length.out = 6))
  H = length(h.set)


  res <- readRDS(paste0("data/model_selection_results_",kernel,1, ".rds"))
  bad.ids <- c() #c(6,38)
  for(j in seq(2,500)){
    file.name <- paste0("data/model_selection_results_",kernel,j, ".rds")
    if(file.exists(file.name) & ! j %in% bad.ids){
      res.tmp <- readRDS(file.name)
      res <- abind(res, res.tmp, along = 1)
    }
  }

  n.sims = dim(res)[1]

  # fix the results
  correct.model.vec <- rep(paste(kernel, " Kernel h = ", h.set), each = J)
  # subdivide by sample size
  correct.mean.vec = c()
  correct.sd.vec = c()
  for(h in seq(H)){
    res.tmp = as.matrix(res[,,h])
    correct.mean.vec = c(correct.mean.vec, colMeans(res.tmp, na.rm = T))
    correct.sd.vec = c(correct.sd.vec, colSDs(res.tmp, na.rm = T))
  }
  #correct.sd.vec = sqrt(correct.mean.vec*(1 - correct.mean.vec))
  correct.se.vec = correct.sd.vec/sqrt(n.sims)
  sample.size.vec = rep(n.seq)

  res.data <- data.frame("SampleSize" = sample.size.vec,
                         "TrueModel" = correct.model.vec,
                         "ModelDev" = correct.mean.vec,
                         "ModelDev_sd" = correct.se.vec)


  plt.mod.sel <- ggplot(res.data, aes(x = log(SampleSize), y = ModelDev, group = TrueModel,color = TrueModel)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = ModelDev - 2*ModelDev_sd, ymax = ModelDev + 2*ModelDev_sd)) +
    ggtitle(paste0(kernel, " Kernel Model Selection")) +
    xlab("log-Sample Size") +
    ylab("Devation of (h) from Correct Model")
  #geom_errorbar(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled))

  plt.mod.sel

  png(filename = paste0("plots/model_selection",kernel,".png"),
      width = png.width, height = png.height, res = png.res)

  plt.mod.sel
  # Close the pdf file
  dev.off()
}











