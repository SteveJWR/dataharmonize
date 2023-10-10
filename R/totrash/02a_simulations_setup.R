
source("R/01_functions.R")

# picking a set of grid tuning parameters for simulations. 
h.set<-  c(0.25,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0)  #c(1,2) #c(0.25,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 

ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)# list(gaussian_kernel, epanechnikov_kernel) #list(gaussian_kernel, exponential_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 3))) #c(0,0.1) #c(0,exp(seq(log(0.001),log(0.1), length.out = 3)))
mu.set.long <- c(0,exp(seq(log(0.001),log(.1), length.out = 15))) # For selection of Mu 

# Find these to be optimal
mu.y <- 0.01930
mu.z <- 0.01390 


h.set.conversion <- c(0.25,0.5,0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5,8.0)  #c(1,2) #c(0.25,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 
#mu.set.conversion <- c(0,exp(seq(log(0.001),log(mu.z), length.out = 3)), exp(seq(log(mu.y),log(.3), length.out = 3)))

mu.set.conversion <- c(0,0.0037,0.0139,0.0193,0.0761,0.3)
ker.set.conversion <- list(gaussian_kernel, exponential_kernel)

h.set.feasibility <- c(0.25,0.5,1,2,4,7,10,14,20,25,30)

h.set.speed <- c(0.25,0.5, 0.75,1.0,2.0 ,3.0, 5.0, 8.0) 

# number of latent bins used in the speed test 
R_bins.speed <- 300

#number of bins used in the feasibility test 

# number of individuals in the data set
dataset.size <- 100
# number of simulations 
n.sims <- 100
n.sims.conversion <- 10 #TODO: REturn to 25  # Takes much longer than the others
n.sims.feasibility <- 300
n.sims.speed <- 25 # Takes much longer than the others

# list of possible tuning parameters for the intrinsic variability 
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))

# compact support kernels must have h >= 1
idx1 <- hyper.param.idx[,1]  %in% which(h.set <= 1)
idx2 <- hyper.param.idx[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.idx <- hyper.param.idx[!(idx1 & idx2),]


# list of possible tuning parameters for the conversion
hyper.param.conversion.idx <- expand.grid(1:length(h.set.conversion),1:length(ker.set.conversion), 1:length(mu.set.conversion))


# compact support kernels must have h > 1
idx1 <- hyper.param.conversion.idx[,1]  %in% which(h.set.conversion <= 1)
idx2 <- hyper.param.conversion.idx[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.conversion.idx <- hyper.param.conversion.idx[!(idx1 & idx2),]

# Indexing guide for hyper parameter with NPMLE i.e. mu = 0
hyper.param.conversion.ml.idx <- hyper.param.conversion.idx[hyper.param.conversion.idx[,3] == 1,]


N = 30 
Ny = N 
Nz = N
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
      cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
      A.mat <- compute_A_matrix(R_bins = 1000, cond = cond.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.mat[is.nan(A.mat)] <- 0
      A.matrix.set.tmp[[k]] <- A.mat 
      
    }
    A.matrix.set[[j]] <- A.matrix.set.tmp
  }
}
#saveRDS(A.matrix.set, "Data/A_matrix_set.RDS")


# set of possible matrices corresponding to kernel bandwidth pairs
A.matrix.set.conversion <- readRDS("Data/A_matrix_set_conversion.RDS")
if(!exists("A.matrix.set")){
  # Pre computing A matrices
  A.matrix.set.conversion <- list()
  for(j in 1:length(h.set.conversion)){
    h.tmp <- h.set.conversion[j]
    A.matrix.set.tmp <- list()
    for(k in 1:length(ker.set.conversion)){
      
      ker.tmp <- ker.set[[k]]
      cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
      A.mat <- compute_A_matrix(R_bins = 1000, cond = cond.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.mat[is.nan(A.mat)] <- 0
      A.matrix.set.tmp[[k]] <- A.mat 
      
    }
    A.matrix.set.conversion[[j]] <- A.matrix.set.tmp
  }
}
#saveRDS(A.matrix.set.conversion, "Data/A_matrix_set_conversion.RDS")



# Simulation 1 Intrinsic Variability setup 
# ---------------------------


### true model parameters (model 1)
k.model1 <- gaussian_kernel
h.model1 <- 2

beta1.model1 <- 12 
beta2.model1 <- 5

cond.model1 <- conditional_mkm(N, ker = k.model1, h = h.model1)
A.matrix.model1 <- compute_A_matrix(R_bins = 1000, cond = cond.model1, numeric.points = 400)
A.tensor.model1 <- compute_A_two_obs_tensor(R_bins = 1000, cond = cond.model1, numeric.points = 400)





# Simulation 2 Intrinsic Variability setup 
# ---------------------------


### true model parameters (model 1)
k.model2 <- exponential_kernel
h.model2 <- 1

beta1.model2 <- 6 
beta2.model2 <- 6

cond.model2 <- conditional_mkm(N, ker = k.model2, h = h.model2)
A.matrix.model2 <- compute_A_matrix(R_bins = 1000, cond = cond.model2, numeric.points = 400)
A.tensor.model2 <- compute_A_two_obs_tensor(R_bins = 1000, cond = cond.model2, numeric.points = 400)





# ---- Feasibility Test Grid options ---- 

# set of possible matrices corresponding to kernel bandwidth pairs
A.matrix.set.feasibility <- readRDS("Data/A_matrix_set_feasibility.RDS")
if(!exists("A.matrix.set.feasibility")){
  # Pre computing A matrices
  A.matrix.set.feasibility <- list()
  for(j in 1:length(h.set.feasibility)){
    h.tmp <- h.set.feasibility[j]
    A.matrix.set.feasibility.tmp <- list()
    for(k in 1:length(ker.set)){
      
      ker.tmp <- ker.set[[k]]
      cond.tmp <- conditional_mkm(N, ker = ker.tmp, h = h.tmp)
      A.mat <- compute_A_matrix(R_bins = 1000, cond = cond.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.mat[is.nan(A.mat)] <- 0
      A.matrix.set.feasibility.tmp[[k]] <- A.mat 
      
    }
    A.matrix.set.feasibility[[j]] <- A.matrix.set.feasibility.tmp
  }
}
#saveRDS(A.matrix.set.feasibility, "Data/A_matrix_set_feasibility.RDS")




# set of possible tensors corresponding to kernel bandwidth pairs
A.tensor.set.feasibility <- readRDS("Data/A_tensor_set_feasibility.RDS")
if(!exists("A.tensor.set.feasibility")){
  # Pre computing A tensors
  A.tensor.set.feasibility <- list()
  for(j in 1:length(h.set.feasibility)){
    h.tmp <- h.set.feasibility[j]
    A.tensor.set.feasibility.tmp <- list()
    for(k in 1:length(ker.set)){
      
      ker.tmp <- ker.set[[k]]
      cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
      A.tensor <- compute_A_two_obs_tensor(R_bins = 1000, cond = cond.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.tensor[is.nan(A.tensor)] <- 0
      A.tensor.set.feasibility.tmp[[k]] <- A.tensor 
      
    }
    A.tensor.set.feasibility[[j]] <- A.tensor.set.feasibility.tmp
  }
}
#saveRDS(A.tensor.set.feasibility, "Data/A_tensor_set_feasibility.RDS")



# list of possible tuning parameters for the feasibility test 
hyper.param.feasibility.idx <- expand.grid(1:length(h.set.feasibility),1:length(ker.set))

# compact support kernels must have h >= 1
idx1 <- hyper.param.feasibility.idx[,1]  %in% which(h.set.feasibility < 1)
idx2 <- hyper.param.feasibility.idx[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.feasibility.idx <- hyper.param.feasibility.idx[!(idx1 & idx2),]




### Pre computing speed test matrices

# ---- Feasibility Test Grid options ---- 

# set of possible matrices corresponding to kernel bandwidth pairs
A.matrix.set.speed <- readRDS("Data/A_matrix_set_speed.RDS")
if(!exists("A.matrix.set.speed")){
  # Pre computing A matrices
  A.matrix.set.speed <- list()
  for(j in 1:length(h.set.speed)){
    h.tmp <- h.set.speed[j]
    A.matrix.set.speed.tmp <- list()
    ker.tmp <- gaussian_kernel
    cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    A.mat <- compute_A_matrix(R_bins = R_bins.speed, cond = cond.tmp, numeric.points = 400)
    
    
    # replacing NAN with 0 due to numerical rounding error 
    A.mat[is.nan(A.mat)] <- 0
      
    
    A.matrix.set.speed[[j]] <- A.mat
  }
}
#saveRDS(A.matrix.set.speed, "Data/A_matrix_set_speed.RDS")







## ------------ Simulation functions ------------ 


simulation_intrinsic_variability_model1 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  intrinsic.variability.each.model <- sapply(1:nrow(hyper.param.idx), function(x){
    
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    l = hyper.param.idx[x,3]
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


simulation_intrinsic_variability_model2 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model2, beta2.y = beta2.model2, 
                            cond.y = cond.model2, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  intrinsic.variability.each.model <- sapply(1:nrow(hyper.param.idx), function(x){
    
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    l = hyper.param.idx[x,3]
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





## -------- mu selection --------
## --------    Model 1   --------


simulation_smoothing_selection_model1 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  neg.log.like.each.mu <- sapply(mu.set.long, function(x){
    mu.tmp <- x
    
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix.model1, mu = mu.tmp)
    
    cond <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
    }
    # intrinsic variability results
    log.like <- compute_two_obs_loglikelihood(sim.data, A.tensor.model1, latent.mix.list)
    
    
    neg.log.like <- -log.like
    return(neg.log.like)
  })
  return(neg.log.like.each.mu)
}


## --------    Model 2   --------

simulation_smoothing_selection_model2 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model2, beta2.y = beta2.model2, 
                            cond.y = cond.model2, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  neg.log.like.each.mu <- sapply(mu.set.long, function(x){
    mu.tmp <- x
    
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix.model2, mu = mu.tmp)
    
    cond <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
    }
    # intrinsic variability results
    log.like <- compute_two_obs_loglikelihood(sim.data, A.tensor.model2, latent.mix.list)
    
    
    neg.log.like <- -log.like
    return(neg.log.like)
  })
  return(neg.log.like.each.mu)
}



## --------   Conversion Cross Entropy   --------

### population test set 


test.grid <- expand.grid(0:Ny, 0:Nz)
# required for estimating the population cross entropy 
p.yz <- compute_joint_dist_beta(beta1.model.y = beta1.model1, beta2.model.y = beta2.model1,
                                beta1.model.z = beta1.model2, beta2.model.z = beta2.model2,
                                cond.y = cond.model1, cond.z = cond.model2, grid.size = 10000)




simulation_conversion_cross_entropy <- function(sim.number){
  set.seed(sim.number)
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1,
                            beta1.z = beta1.model2, beta2.z = beta2.model2, 
                            cond.z = cond.model2,
                            pair.obs = T)
  
  # only keep the single observations 
  sim.data <- sim.data[,c(1,3)]
  train.p.hat.y <- compute_edf(sim.data[,1], N)
  train.p.hat.z <- compute_edf(sim.data[,2], N)
  
  
  ### Fixed Model Estimate of y under correctly learned model 
  model.estimate.y <- estimate_mixing_numeric(p.hat = train.p.hat.y, A.matrix = A.matrix.model1, mu = mu.y)
  model.estimate.y.ml <- estimate_mixing_numeric(p.hat = train.p.hat.y, A.matrix = A.matrix.model1, mu = 0)
  
  
  ## Computes the intrinsic variability across all possible models in our set 
  conversion.each.model <- sapply(1:nrow(hyper.param.conversion.idx), function(x){
    
    j = hyper.param.conversion.idx[x,1]
    k = hyper.param.conversion.idx[x,2]
    l = hyper.param.conversion.idx[x,3]
    h.tmp <- h.set.conversion[j]
    ker.tmp <- ker.set.conversion[[k]]
    mu.tmp <- mu.set.conversion[l]
    A.matrix.z.tmp <- A.matrix.set.conversion[[j]][[k]]
    
    model.estimate.z <- estimate_mixing_numeric(p.hat = train.p.hat.z, 
                                                A.matrix = A.matrix.z.tmp, 
                                                mu = mu.tmp)
    
    cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    
    latent.mix.list.y <- list()
    model.observed.list.y <- list()
    latent.mix.list.z <- list()
    model.observed.list.z <- list()
    for (m in 1:nrow(test.grid)) {
      latent.mix.list.y[[m]] <- model.estimate.y$latent
      latent.mix.list.z[[m]] <- model.estimate.z$latent
    }
    
    # conversion population cross entropy 
    ce.pop <- convert_score_ce(test.pairs = test.grid, 
                               latent.mix.list.y = latent.mix.list.y, 
                               latent.mix.list.z = latent.mix.list.z, 
                               cond.y = cond.model1, 
                               cond.z = cond.tmp, 
                               joint.prob = p.yz,
                               grid.size = 1000) 
    
    return(ce.pop)
  })
  
  
  conversion.each.model.ml <- sapply(1:nrow(hyper.param.conversion.ml.idx), function(x){
    
    j = hyper.param.conversion.ml.idx[x,1]
    k = hyper.param.conversion.ml.idx[x,2]
    l = hyper.param.conversion.ml.idx[x,3]
    h.tmp <- h.set.conversion[j]
    ker.tmp <- ker.set.conversion[[k]]
    A.matrix.z.tmp <- A.matrix.set.conversion[[j]][[k]]
    cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    
    model.estimate.z.ml <- estimate_mixing_numeric(p.hat = train.p.hat.z, 
                                                A.matrix = A.matrix.z.tmp, 
                                                mu = 0)
    
    
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    
    latent.mix.list.y.ml <- list()
    latent.mix.list.z.ml <- list()
    for (m in 1:nrow(test.grid)) {
      latent.mix.list.y.ml[[m]] <- model.estimate.y.ml$latent
      latent.mix.list.z.ml[[m]] <- model.estimate.z.ml$latent
    }
    
    # conversion population cross entropy 
    ce.pop.ml <- convert_score_ce(test.pairs = test.grid, 
                                  latent.mix.list.y = latent.mix.list.y.ml, 
                                  latent.mix.list.z = latent.mix.list.z.ml, 
                                  cond.y = cond.model1, 
                                  cond.z = cond.tmp, 
                                  joint.prob = p.yz,
                                  grid.size = 1000) 
    
    return(ce.pop.ml)
  })
  
  results <- list("smoothed" = conversion.each.model,
                  "npmle" = conversion.each.model.ml)
  return(results)
}





## feasibility test simulations

simulation_feasibility_test_model1 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  train.bivariate.p.hat <- compute_bivariate_edf(sim.data, N)
    
  ## Computes the feasibility test for each model
  feasibility.test.each.model <- lapply(1:nrow(hyper.param.feasibility.idx), function(x){
    
    j = hyper.param.feasibility.idx[x,1]
    k = hyper.param.feasibility.idx[x,2]

    
    A.matrix <- A.matrix.set.feasibility[[j]][[k]]
    A.tensor <- A.tensor.set.feasibility[[j]][[k]]
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, 
                                              A.matrix = A.matrix, 
                                              mu = 0)
    
    
    # feasibility tests
    p.order1 <- test_feasibility_first_order(latent.mixture = model.estimate$latent, 
                                             A.matrix = A.matrix, 
                                             p.hat = train.p.hat, 
                                             sample.size = dataset.size)
    

    p.order2 <- test_feasibility_second_order(latent.mixture = model.estimate$latent, 
                                              A.two.sample.tensor = A.tensor, 
                                              p.hat = train.bivariate.p.hat, 
                                              sample.size = dataset.size)
      
    out.tests <- list("first_order" = p.order1, 
                      "second_order" = p.order2)
    
    return(out.tests)
  })
  return(feasibility.test.each.model)
}



### ---- Speed tests ---- 

simulation_speed_test <- function(sim.number){
  set.seed(sim.number)
  ## Computes the time it takes to fit each model
  # setup of parameters of the simulation
  n.mc.samp <- 1000
  em.tau <- seq(0,1,length.out = R_bins.speed)
  # array for computing the results for the latent distribution 
  results.time <- array(NA, dim = c(length(h.set.speed), 2, length(mu.set)))
  results.likelihood <- array(NA, dim = c(length(h.set.speed), 2, length(mu.set)))
  
  ## Computes the time for each algorithm to be applied to the model 
  for(l in 1:length(h.set.speed)){
    h.tmp <- h.set.speed[l]
    A.model <- A.matrix.set.speed[[l]]
    
    
    cond.model.speed <- conditional_mkm(N, ker = gaussian_kernel, h = h.tmp)
    
    sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                              beta1.y = beta1.model1, beta2.y = beta2.model1, 
                              cond.y = cond.model.speed, pair.obs = F)
    
    train.p.hat <- compute_edf(sim.data[,1], N)
    
    for(k in 1:length(mu.set)){
      mu.tmp <- mu.set[k]
      
      # time this
      em.time <- system.time(em.quantiles <- fit_nonpar_em(p.hat  = train.p.hat, 
                                                           cond = cond.model.speed,  
                                                           R_bins = R_bins.speed, 
                                                           mu = mu.tmp, 
                                                           n.mc.samp = n.mc.samp, 
                                                           verbose = F)) 
                               
      
      #time this 
      gp.time <- system.time(gp.model.estimate.binned <- estimate_mixing_numeric(p.hat = train.p.hat, 
                                                                                 A.matrix = A.model, 
                                                                                 mu = mu.tmp))
      
      
      em.p.ma <- compute_p_ma(tau = em.tau,
                              latent.trait.quantiles = em.quantiles,
                              cond = cond,
                              numeric.points = 100)
      
      
      em.lik <- compute_loglikelihood_from_latent(p.hat = train.p.hat, 
                                                  p.ma = em.p.ma, 
                                                  tau = em.tau, 
                                                  latent.trait.quantiles = em.quantiles,
                                                  mu = mu.tmp)
      
      
      # Defines the quantiles of the latent bins 
      gp.tau <- inv_quantiles_from_weights(gp.model.estimate.binned$latent)
      gp.quantiles <- seq(0,1, length.out = length(gp.tau))
      
      gp.p.ma <- as.vector(gp.model.estimate.binned$observed)
      
      gp.lik <- compute_loglikelihood_from_latent(p.hat = train.p.hat, 
                                                  p.ma = gp.p.ma, 
                                                  tau = gp.tau, 
                                                  latent.trait.quantiles = gp.quantiles,
                                                  mu = mu.tmp)
      
      
      results.time[l,1,k] <- em.time[3]
      results.time[l,2,k] <- gp.time[3]
      
      results.likelihood[l,1,k] <- em.lik
      results.likelihood[l,2,k] <- gp.lik
        
      }
    }
  
  out <- list("time" = results.time, "likelihood" = results.likelihood)
  return(out)
}



### ---- simulation framework function --- 

simulate_experiment <- function(RUN_PARALLEL, sim_function, n.sims, start.sim = 1){
  if (RUN_PARALLEL) {
    # Detect number of cores, use all but 1
    no_cores <- detectCores() - 1
    # Initiate cluster
    tictoc::tic()
    
    cl <- makeCluster(no_cores, type="FORK",timeout=timeout)
    
    
    # Run computation
    results = parLapply(cl = cl, X = 1:n.sims, fun = sim_function)
    # Stop cluster
    stopCluster(cl)
    tictoc::toc()
  } else {
    tictoc::tic()
    results = lapply(X = start.sim:n.sims, FUN = function(z){
      out <- sim_function(z)
      cat(paste0("Simulation: ", z, "/",n.sims), end = "\r")
      return(out)
      })
    tictoc::toc()
  }
  return(results)
}






