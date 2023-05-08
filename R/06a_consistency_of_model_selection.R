

#In this script, we illustrate the consistency of the model selection procedure.
rm(list = ls())


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}

if(id %% 2 == 1){
  kernel="Gaussian"
  ker.true <- gaussian_kernel
} else {
  kernel="Laplace"
  ker.true <- laplace_kernel
}

R.bins = 1000

# True Conditional Model
h.true <- 2


N <- 30

cond.true <- conditional_mkm(N,ker.true, h.true)


two.obs.ratio = 1
n1.seq = c(100,200,500,1000)
n2.seq =  two.obs.ratio*n1.seq
J = length(n1.seq)


# grid for the values of h
h.seq <- c(0.8,1,2,3,5,10)

i.true = 3
if(kernel == "Gaussian"){
  cond.set <- dnoiseR::generate_mkm_list(N = N, ker = gaussian_kernel, h.set = h.seq)
  cond.names <- paste0("Gaussian h = ",as.character(h.seq))
} else {
  cond.set <- dnoiseR::generate_mkm_list(N = N, ker = laplace_kernel, h.set = h.seq)
  cond.names <- paste0("Laplace h = ",as.character(h.seq))
}



# correct model selection fraction of sims.
n.sims <- 5
res <- matrix(NA, nrow= n.sims, ncol =length(n1.seq))

#simulation parameters
alpha1 <- 1.2
alpha2 <- 3

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
    res[sim,j] <- TRUE*(i.max == i.true) + FALSE*(i.max != i.true)
  }
}

saveRDS(res, paste0("data/model_selection_results_",kernel,id, ".rds"))


















