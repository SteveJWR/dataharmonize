

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

#simulation parameters
alpha1 <- 1.2
alpha2 <- 3


for(j in seq(J)){
  n = n.seq[j]
  for(sim in seq(n.sims)){
    cat(paste0("Simulation ", sim, "/",n.sims), end = "\r")
    gamma <- c(rbeta((1/3)*(n),alpha1,alpha2),rbeta((2/3)*(n),alpha2,alpha1))
    gamma <- gamma[sample(seq(length(gamma)))] # shuffle


    Y <- simulate_test_cond(obs.set = c(rep(1,length(gamma))),cond.true ,gamma)
    model <- mu_selection(mu.set,cond.true,Y,R.bins,folds, verbose = F)
    res[sim,j,] <- model$cv.lik
  }
}

saveRDS(res, paste0("data/regularization_selection_results_",kernel,id, ".rds"))




