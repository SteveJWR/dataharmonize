

## Feasibility Tests
#In this script, we illustrate the feasibility tests
rm(list = ls())
library(dnoiseR)

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 2
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
n.seq = c(100,200,500,1000)

J = length(n.seq)


# grid for the values of h
h.set <- c(0.5,0.8,1,2,3,5,10,20)
H = length(h.set)

i.true = 3
if(kernel == "Gaussian"){
  cond.set <- generate_mkm_list(N = N, ker = gaussian_kernel, h.set = h.set)
  cond.names <- paste0("Gaussian h = ",as.character(h.set))
} else {
  cond.set <- generate_mkm_list(N = N, ker = exponential_kernel, h.set = h.set)
  cond.names <- paste0("Exponential h = ",as.character(h.set))
}



# correct model selection fraction of sims.
n.sims <- 5
res1 <- array(NA, c(n.sims,length(n.seq),H))
res2 <- res1

#simulation parameters
alpha1 <- 1.2
alpha2 <- 3

for(j in seq(J)){
  n = n.seq[j]

  for(i in seq(length(cond.set))){

    cond <- cond.set[[i]]
    A.matrix <- compute_A_matrix(R.bins, cond)
    A.tensor <- compute_A_tensor(R.bins, cond)
    for(sim in seq(n.sims)){

      cat(paste0("Simulation ", sim, "/",n.sims), end = "\r")
      gamma <- c(rbeta((1/3)*(n),alpha1,alpha2),rbeta((2/3)*(n),alpha2,alpha1))
      gamma <- gamma[sample(seq(length(gamma)))] # shuffle


      Y <- simulate_test_cond(obs.set = c(rep(2,length(gamma))),cond.true ,gamma)
      p.hat <- compute_edf(Y[,1], N)
      p.hat.2 <- compute_edf_2(Y,N)
      model1 <- estimate_mixing_npem(p.hat,A.matrix,mu = 0)
      model2 <- estimate_mixing_npem_2(Y,A.matrix,A.tensor,mu = 0)

      p.first.order <- first_order_feasibility_test(model1$observed, p.hat, length(gamma))
      p.second.order <- second_order_feasibility_test(model2$latent,A.tensor,p.hat.2, length(gamma))

      res1[sim,j,i] <- p.first.order
      res2[sim,j,i] <- p.second.order
    }
  }

}


saveRDS(res1, paste0("data/feasibility_test_order_1_results_",kernel,ceiling(id/2), ".rds"))
saveRDS(res2, paste0("data/feasibility_test_order_2_results_",kernel,ceiling(id/2), ".rds"))





















