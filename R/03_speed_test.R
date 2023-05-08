






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

#N.seq = c(100,500,1000,5000)
N.set <- c(5,10,20,50,100)
J = length(N.set)




# grid for the values of mu
L = 5
mu.set <- seq(0,2,length.out = L)



# dimensions of the results (sim, samplesize, tuning parameter)
res <- array(NA,c(J,length(n.seq),L))
time.npem <- res
time.numeric <- res
like.npem <- res
like.numeric <- res

# 4 different marginals
marg.1 <- function(N){return(rep(1/(N+1),(N + 1)))}
marg.2 <- function(N){return(-(seq(0,N) - N/2)**2 + 2*(N/2)^2)}
marg.3 <- function(N){return(exp(-(seq(0,N))**2) + exp(- (seq(0,N) - N)**2))}
marg.4 <- function(N){return(sin(seq(0,N)/N * 4* pi) + 1)}

marg.list <- list(marg.1, marg.2, marg.3, marg.4)

K <- length(marg.list)
uniform.latent <- rep(1/R.bins,R.bins)
for(j in seq(J)){
  N = N.set[j]
  cond.true <- conditional_mkm(N,ker.true, h.true)
  # model for the distribution
  A.matrix <- compute_A_matrix(R.bins, cond.true)


  for(k in seq(K)){
    marg <- marg.list[[k]]
    p.hat <- marg(N)
    p.hat <- p.hat/sum(p.hat)

    cat(paste0("Conditional ", k, "/",K), end = "\r")
    for(i in seq(length(mu.set))){
      mu = mu.set[i]

      time1 <- Sys.time()
      model.npem <- estimate_mixing_npem(p.hat,A.matrix, mu)
      time2 <- Sys.time()
      model.numeric <- estimate_mixing_numeric(p.hat,A.matrix, mu)
      time3 <- Sys.time()

      time.npem <- as.numeric(difftime(time2, time1, units = "secs"))
      time.numeric <- as.numeric(difftime(time3, time2, units = "secs"))
      if(mu == 0){
        like.npem <- -kl_divergence(p.hat,model.npem$observed)
        like.numeric <- -kl_divergence(p.hat,model.numeric$observed)

      } else {
        like.npem <- -kl_divergence(p.hat,model.npem$observed) - mu*kl_divergence(uniform.latent, model.npem$latent)
        like.numeric <- -kl_divergence(p.hat,model.numeric$observed) - mu*kl_divergence(uniform.latent, model.numeric$latent)
      }
      time.npem[j,k,i] <- time.npem
      time.numeric[j,k,i] <- time.numeric

      like.npem[j,k,i] <- like.npem
      like.numeric[j,k,i] <- like.numeric

    }
  }
}

saveRDS(time.npem, paste0("data/fitting_speed_npem_results_",kernel,id, ".rds"))
saveRDS(time.cvxr, paste0("data/fitting_speed_cvxr_results_",kernel,id, ".rds"))

saveRDS(like.npem, paste0("data/fitting_likelihood_npem_results_",kernel,id, ".rds"))
saveRDS(like.cvxr, paste0("data/fitting_likelihood_cvxr_results_",kernel,id, ".rds"))












