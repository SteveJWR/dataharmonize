
# TODO: DELETE 
# setwd("/Users/Owner/Documents/PhD/Latent Variable Modelling with Application to Data Harmonization/GitHub/Data-Harmonization-Nonparametric")
source("R/01_functions.R")
source("R/03a_application_setup.R")



#### MMSE Model Selection 
p.hat.list.y.val <- sort_conditionals(data = y.val, 
                                      x.design = x.design, 
                                      p.hat.list = p.hat.short.list.y)


### Applying this to three choices of mu for completeness
cond.y.model.0 <- select_cond(pair.obs = y.val,  
                              p.hat = p.hat.list.y.val, 
                              cond.list = cond.list, 
                              cond.names = cond.names, 
                              mu = 0, 
                              n.samp = 10, 
                              R_bins = R_bins)

cond.y.model.0.001 <- select_cond(pair.obs = y.val,  
                                  p.hat = p.hat.list.y.val, 
                                  cond.list = cond.list, 
                                  cond.names = cond.names, 
                                  mu = 0.001, 
                                  n.samp = 10, 
                                  R_bins = R_bins)

cond.y.model.0.1 <- select_cond(pair.obs = y.val,  
                                p.hat = p.hat.list.y.val, 
                                cond.list = cond.list, 
                                cond.names = cond.names, 
                                mu = 0.1, 
                                n.samp = 10, 
                                R_bins = R_bins)



saveRDS(cond.y.model.0, "Data/results/conditional_model_selection_mu0_mmse.rds")
saveRDS(cond.y.model.0.001, "Data/results/conditional_model_selection_mu0001_mmse.rds")
saveRDS(cond.y.model.0.1, "Data/results/conditional_model_selection_mu01_mmse.rds")








##### MOCA Model Selection 


p.hat.list.z.val <- sort_conditionals(data = z.val, 
                                      x.design = x.design, 
                                      p.hat.list = p.hat.short.list.z)


### Applying this to three choices of mu for completeness
cond.z.model.0 <- select_cond(pair.obs = z.val,  
                              p.hat = p.hat.list.z.val, 
                              cond.list = cond.list, 
                              cond.names = cond.names, 
                              mu = 0, 
                              n.samp = 10, 
                              R_bins = R_bins)

cond.z.model.0.001 <- select_cond(pair.obs = z.val,  
                                  p.hat = p.hat.list.z.val, 
                                  cond.list = cond.list, 
                                  cond.names = cond.names, 
                                  mu = 0.001, 
                                  n.samp = 10, 
                                  R_bins = R_bins)

cond.z.model.0.1 <- select_cond(pair.obs = z.val,  
                                p.hat = p.hat.list.z.val, 
                                cond.list = cond.list, 
                                cond.names = cond.names, 
                                mu = 0.1, 
                                n.samp = 10, 
                                R_bins = R_bins)



saveRDS(cond.z.model.0, "Data/results/conditional_model_selection_mu0_moca.rds")
saveRDS(cond.z.model.0.001, "Data/results/conditional_model_selection_mu0001_moca.rds")
saveRDS(cond.z.model.0.1, "Data/results/conditional_model_selection_mu01_moca.rds")



#Optimal Models 
cond.y.opt <- cond.z.model.0.001$optimal # Binomial|| TV dist: 0.0308
cond.z.opt <- cond.z.model.0.1$optimal # Exponential Kernel, h = 1.32|| TV dist: 0.072




#### Smoothing Selection MMSE 

y.smooth <- select_mu(pair.obs = y.val, 
                      p.hat = p.hat.list.y.val, 
                      cond = cond.y.opt, 
                      mu.vec = mu.set.smooth, 
                      R_bins = R_bins)


#### Smoothing Selection MOCA 

z.smooth <- select_mu(pair.obs = z.val, 
                      p.hat = p.hat.list.z.val, 
                      cond = cond.z.opt, 
                      mu.vec = mu.set.smooth, 
                      R_bins = R_bins)



saveRDS(y.smooth, "Data/results/smoothing_selection_mmse.rds")
saveRDS(z.smooth, "Data/results/smoothing_selection_moca.rds")


mu.y.opt <- y.smooth$optimal # 0.00379 index: 12
mu.z.opt <- z.smooth$optimal # 0.011 index: 15


### Conversion from Y to Z  

conversion.results.array <- array(NA, dim = c(length(cond.list),length(mu.set.conversion)))
conversion.results.array.ml <- array(NA, dim = c(length(cond.list)))
conversion.results.array.para <- conversion.results.array.ml


###

# optimal A matrix for Y 
A.mat.y.opt <- compute_A_matrix(R_bins, cond = cond.y.opt)

p.hat.list.y.cw <- sort_conditionals(data = cw.test, 
                                     x.design = x.design, 
                                     p.hat.list = p.hat.short.list.y)


p.hat.list.z.cw <- sort_conditionals(data = cw.test, 
                                     x.design = x.design, 
                                     p.hat.list = p.hat.short.list.z)

y.param <- logitnorm_model_fit(y.train.wide[,1], X = y.train.wide[,2:6], cond = cond.y.opt)

## list of latent y. 
p.hat.list.short <- unique(p.hat.list.y.cw)
latent.mix.list.short <- list()
for(k in 1:length(p.hat.list.short)){
  latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.y.opt, mu = mu.y.opt)
  latent.mix.list.short[[k]] <- latent.model$latent
}

latent.mix.y.list.cw <- list()
for(i in 1:length(p.hat.list.y.cw)){
  p.hat.tmp <- p.hat.list.y.cw[[i]]
  short.idx <- which(sapply(p.hat.list.short, function(s) all(s == p.hat.tmp)))
  latent.mix.y.list.cw[[i]] <- latent.mix.list.short[[short.idx]]
}


## list of latent y with no regularization
p.hat.list.short <- unique(p.hat.list.y.cw)
latent.mix.list.short.ml <- list()
for(k in 1:length(p.hat.list.short)){
  latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.y.opt, mu = 0)
  latent.mix.list.short.ml[[k]] <- latent.model$latent
}

latent.mix.y.list.cw.ml <- list()
for(i in 1:length(p.hat.list.y.cw)){
  p.hat.tmp <- p.hat.list.y.cw[[i]]
  short.idx <- which(sapply(p.hat.list.short, function(s) all(s == p.hat.tmp)))
  latent.mix.y.list.cw.ml[[i]] <- latent.mix.list.short.ml[[short.idx]]
}



# use all models for conversion: 
# 



for(j in 1:length(cond.list)){
  cond.tmp <- cond.list[[j]]
  A.mat.tmp <- compute_A_matrix(R_bins = R_bins, cond = cond.tmp)
  # replacing error with some kernels
  A.mat.tmp[is.na(A.mat.tmp)] <- 0
  ce.samp.mu.vec <- sapply(1:length(mu.set.conversion), function(l){
    mu.tmp <- mu.set.conversion[l]
    
    # here we introduce a shortcut for computing. We only used the unique elements of the list
    # so we don't have to re-estimate mixing distributions that will be the same. 
    p.hat.list.short <- unique(p.hat.list.z.cw)
    latent.mix.list.short <- list()
    for(k in 1:length(p.hat.list.short)){
      latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.tmp, mu = mu.tmp)
      latent.mix.list.short[[k]] <- latent.model$latent
    }
    
    latent.mix.z.list.cw <- list()
    for(i in 1:length(p.hat.list.z.cw)){
      p.hat.tmp <- p.hat.list.z.cw[[i]]
      short.idx <- which(sapply(p.hat.list.short, function(s) all(s == p.hat.tmp)))
      latent.mix.z.list.cw[[i]] <- latent.mix.list.short[[short.idx]]
    }
    
    
    # sample cross entropy 
    ce.samp  <- convert_score_metric(test.pairs = cw.test, 
                                     latent.mix.list.y = latent.mix.y.list.cw, 
                                     latent.mix.list.z = latent.mix.z.list.cw, 
                                     cond.y = cond.y.opt, 
                                     cond.z = cond.tmp, # using the current iteration of cond.z 
                                     grid.size = 1000) 
    
    
  })
  conversion.results.array[j,] <- ce.samp.mu.vec
  
  ### ML (unregularized estimate)
  
  p.hat.list.short <- unique(p.hat.list.z.cw)
  latent.mix.list.short.ml <- list()
  for(k in 1:length(p.hat.list.short)){
    latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.tmp, mu = 0)
    latent.mix.list.short.ml[[k]] <- latent.model$latent
  }

  latent.mix.z.list.cw.ml <- list()
  for(i in 1:length(p.hat.list.z.cw)){
    p.hat.tmp <- p.hat.list.z.cw[[i]]
    short.idx <- which(sapply(p.hat.list.short, function(s) all(s == p.hat.tmp)))
    latent.mix.z.list.cw.ml[[i]] <- latent.mix.list.short.ml[[short.idx]]
  }


  # sample cross entropy
  ce.samp.ml <-  convert_score_metric(test.pairs = cw.test,
                                      latent.mix.list.y = latent.mix.y.list.cw.ml,
                                      latent.mix.list.z = latent.mix.z.list.cw.ml,
                                      cond.y = cond.y.opt,
                                      cond.z = cond.tmp, # using the current iteration of cond.z
                                      grid.size = 1000) # TODO: Replace with larger value

  conversion.results.array.ml[j] <- ce.samp.ml


  ## Parametric Model
  z.param <- logitnorm_model_fit(z.train.wide[,1], X = z.train.wide[,2:6], cond = cond.tmp)

  ce.param <- compute_conversion_parametric_ce(test.pairs = cw.test.wide[,c(1,2)], X = cw.test.wide[,3:7],
                                               params.y = y.param, params.z = z.param,
                                               cond.y = cond.y.opt, cond.z = cond.tmp,
                                               grid.size = 1000)

  conversion.results.array.para[j] <- ce.param

  cat(paste0("Model: ", j, "/", length(cond.list), " complete" ), end = "\r")
  
}



## Naive Z-score method 



mu.y <- mean(y.train$y)
sd.y <- sd(y.train$y)
mu.z <- mean(z.train$z)
sd.z <- sd(z.train$z)



ce.zscore = 0

for(i in 1:nrow(cw.test)){
  z.prob = naive_conversion_prob(y = cw.test[i,1], muy = mu.y, sdy = sd.y,
                                 muz = mu.z, sdz = sd.z, Nz = Nz)
  ce.zscore = ce.zscore - log(z.prob[cw.test[i,2] + 1])
  cat(paste0("Test sample ", i, "/",nrow(cw.test)), end = "\r")
}
cat(end = "\n")


saveRDS(conversion.results.array, "Data/results/regularized_cross_entropy.rds")
saveRDS(conversion.results.array.ml, "Data/results/unregularized_cross_entropy.rds")
saveRDS(conversion.results.array.para, "Data/results/parametric_cross_entropy.rds")
saveRDS(ce.zscore, "Data/results/naive_cross_entropy.rds")


