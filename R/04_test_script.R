

## test script for simulations 
## we want to make sure all of the basic functions can 
## run on the cluster

source("R/01_functions.R")


library(parallel)
library(tictoc)


# Preparation 

# picking a set of grid tuning parameters for simulations. 
h.set<- c(1,2) #c(0.25,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 
h.set.feasibility <- c(0.25,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 
ker.set <- list(gaussian_kernel, epanechnikov_kernel) #list(gaussian_kernel, exponential_kernel, epanechnikov_kernel)
mu.set <- c(0,0.1)#c(0,exp(seq(log(0.001),log(0.1), length.out = 3)))
mu.set.long <- c(0,exp(seq(log(0.001),log(.1), length.out = 15))) # For selection of Mu 
h.set.feasibility.test <- c(0.25,0.5,1,2,4,7,10,14,20,25,30)
h.set.speed.test <- c(0.25,0.5, 0.75,1.0,2.0 ,3.0, 5.0, 8.0) 

# number of individuals in the data set
dataset.size <- 100
# number of simulations 
n.sims <- 8
N = 30 
# set of possible matrices corresponding to kernel bandwidth pairs
A.matrix.set <- readRDS("Data/A_matrix_set.RDS")
if(!exists("A.matrix.set")){
  # Pre computing A matrices
  A.matrix.set <- list()
  for(j in 1:length(h.set)){
    h.tmp <- h.set[j]
    A.matrix.set.tmp <- list()
    for(k in 1:length(ker.set)){
      
      ker.tmp <- ker.set[[k]]
      A.mat <- A.matrix.compute(R_bins = 1000, N = 30, ker = ker.tmp, h = h.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.mat[is.nan(A.mat)] <- 0
      A.matrix.set.tmp[[k]] <- A.mat 
      
    }
    A.matrix.set[[j]] <- A.matrix.set.tmp
  }
}
#saveRDS(A.matrix.set, "Data/A_matrix_set.RDS")


# Simulation 1 Intrinsic Variability setup 
# ---------------------------


### true model parameters (model 1)
k.model1 <- gaussian_kernel
h.model1 <- 2

beta1.model1 <- 12 
beta2.model1 <- 5

cond.model1 <- conditional_mkm(N, ker = k.model1, h = h.model1)


results.array.intrinsic.1 <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))
# list of possible tuning parameters 
hyper.param.idx1 <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))

# compact support kernels must have h >= 1
idx1 <- hyper.param.idx1[,1]  %in% which(h.set <= 1)
idx2 <- hyper.param.idx1[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.idx1 <- hyper.param.idx1[!(idx1 & idx2),]




# Simulation 2 Intrinsic Variability setup 
# ---------------------------


### true model parameters (model 1)
k.model2 <- exponential_kernel
h.model2 <- 1

beta1.model2 <- 6 
beta2.model2 <- 6

cond.model2 <- conditional_mkm(N, ker = k.model2, h = h.model2)


results.array.intrinsic.2 <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))
# list of possible tuning parameters 
hyper.param.idx2 <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))

# compact support kernels must have h >= 1
idx1 <- hyper.param.idx2[,1]  %in% which(h.set <= 1)
idx2 <- hyper.param.idx2[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.idx2 <- hyper.param.idx2[!(idx1 & idx2),]











## ---- Simulation functions ---- 

# pass the following objects
simulation_intrinsic_variability_model1 <- function(sim.number){
  # loading the required functions 
  #source("R/DataHarmonizationFunctions.R")
  #source("R/02a_simulations_setup.R")
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  intrinsic.variability.each.model <- sapply(1:nrow(hyper.param.idx1), function(x){
    
    j = hyper.param.idx1[x,1]
    k = hyper.param.idx1[x,2]
    l = hyper.param.idx1[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- mu.set[l]
    A.matrix <- A.matrix.set[[j]][[k]]
    
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix, mu = mu.tmp)
    
    cond <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    model.observed.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
      model.observed.list[[m]] <- model.estimate$observed
    }
    # intrinsic variability results
    iv.res <- intrinsic_variability(pair.obs = sim.data, latent.mix.list = latent.mix.list, 
                                    model.observed.list = model.observed.list, n.samp = 25, cond = cond)
    
    return(iv.res)
  })
  return(intrinsic.variability.each.model)
}


# Settings


RUN_PARALLEL = TRUE


RNGkind("L'Ecuyer-CMRG")
set.seed(1)



k.model1 <- gaussian_kernel
h.model1 <- 2

beta1.model1 <- 12 
beta2.model1 <- 5

cond.model1 <- conditional_mkm(N, ker = k.model1, h = h.model1)


sim.data <- simulate_beta(n.ind = 100, n.obs.per.ind = 2, 
                          beta1.y = beta1.model1, beta2.y = beta2.model1, 
                          cond.y = cond.model1, pair.obs = F)


if (RUN_PARALLEL) {
  # Detect number of cores, use all but 1
  no_cores <- parallel::detectCores() - 1
  # Initiate cluster
  tictoc::tic()
  
  cl <- makeCluster(no_cores, type="FORK")

  
  # Run computation
  result = parLapply(cl = cl, X = 1:n.sims,
                               fun = simulation_intrinsic_variability_model1)
  # Stop cluster
  parallel::stopCluster(cl)
  tictoc::toc()
} else {
  tictoc::tic()
  result = lapply(X = 1:n.sims,
                  FUN = simulation_intrinsic_variability_model1)
  tictoc::toc()
}

saveRDS(result, file = "Data/results/test_script_results")



