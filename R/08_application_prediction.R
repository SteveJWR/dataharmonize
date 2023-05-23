

rm(list = ls())

# suppose that we want to study the impact on education
# and the aging decline across normal individuals
# Our initial data
library(dplyr)
library(gee)
library(dnoiseR)
#library(geepack)

png.width = 1200
png.height = 1000
png.res = 200

R.bins = 1000

Ny = 30
Nz = 30

# imputation application

y.train <- read.csv("data/NACCMMSE_training.csv")
z.train <- read.csv("data/MOCATOTS_training.csv")

y.val <- read.csv("data/NACCMMSE_validation.csv")
z.val <- read.csv("data/MOCATOTS_validation.csv")

y.val <- y.val[,c("y1", "y2", "age", "group")]
z.val <- z.val[,c("z1", "z2", "age", "group")]



### First we find the best Conditional Model for the population

Y.two.obs <- y.val[,c(1,2)]
Z.two.obs <- y.val[,c(1,2)]


h.set <- c(0.5,0.8,1,2,3,4,5,6,7,8,9,10)

ker1 <- gaussian_kernel
ker2 <- exponential_kernel

cond.y.set <- c(generate_mkm_list(Ny, ker1, h.set), generate_mkm_list(Ny, ker2, h.set), generate_cond_binomial(Ny))
cond.y.names <- c(paste0("Gaussain: h=", h.set), paste0("Exponential: h=", h.set), "Binomial")

cond.z.set <- c(generate_mkm_list(Nz, ker1, h.set), generate_mkm_list(Nz, ker2, h.set), generate_cond_binomial(Nz))
cond.z.names <- c(paste0("Gaussain: h=", h.set), paste0("Exponential: h=", h.set), "Binomial")


estimate.cond = F
if(estimate.cond){

  cond.model.y.select <- error_model_selection_bivariate(cond.y.set, Y.two.obs, R.bins, cond.y.names, verbose = T)
  print(cond.model.y.select$opt_model_name)
  plot(cond.model.y.select$lik_vals)


  cond.model.z.select <- error_model_selection_bivariate(cond.z.set, Z.two.obs, R.bins, cond.z.names, verbose = T)
  print(cond.model.z.select$opt_model_name)
  plot(cond.model.z.select$lik_vals)
}




y.train[,4] <- NA
z.train[,4] <- NA

colnames(y.train)[1] <- "y1"
colnames(y.train)[4] <- "y2"

colnames(z.train)[1] <- "z1"
colnames(z.train)[4] <- "z2"

y.train <- y.train[,c("y1", "y2", "age", "group")]
z.train <- z.train[,c("z1", "z2", "age", "group")]

Y.dat <- rbind(y.train, y.val)
Z.dat <- rbind(z.train, z.val)

Y <- Y.dat[,c(1,2)]
X.Y <- Y.dat[,c(3,4)]
Z <- Z.dat[,c(1,2)]
X.Z <- Z.dat[,c(3,4)]

### learn the measurement model

Ny = 30
Nz = 30
R.bins = 1000
#

cond.bin.y <- generate_cond_binomial(Ny)
cond.bin.z <- generate_cond_binomial(Nz)


cond.y <- cond.bin.y
cond.z <- cond.bin.z

# Cross Validation for the conditional Model



mu.set <- c(0,0.005,0.01,0.05,0.1,0.5,1,5,10,50)
# kernel set for the regression model

# function for cross validation
# ker.params is a list of values for

# Smoothing Params for the model
# 0.5 is used for categorized predictors to break them into blocks.
ker.params <- list(seq(1,5,length.out = 5), 0.5)
mu.params <- mu.set
# we use a grid search to jointly select the regression and the latent regularization tuning parameters
reg.kers <- list(gaussian_kernel, uniform_kernel)

eval.cv = F
if(eval.cv){
  y.model.cv <- cv_regression_gridsearch(y.train, cond.y,
                                         folds = 5, outcome = 1,
                                         X.cols = c(3,4),
                                         ker.params = ker.params,
                                         reg.kers = reg.kers,
                                         mu.params = mu.params,
                                         threshold = 10**(-5), max.iter = 50)


  z.model.cv <- cv_regression_gridsearch(z.train, cond.z,
                                         folds = 5, outcome = 1,
                                         X.cols = c(3,4),
                                         ker.params = ker.params,
                                         reg.kers = reg.kers,
                                         mu.params = mu.params,
                                         threshold = 10**(-5), max.iter = 50)


  save.cv = F
  if(save.cv){
    saveRDS(y.model.cv, "data/mmse_cv.rds")
    saveRDS(z.model.cv, "data/moca_cv.rds")
  }
}


# Add a cross validation which accounts for the tuning parameter of the regression piece as well

# select.mu.y <- mu_selection_regression(Y,X.Y,cond.y,mu.set,R.bins.small, ker.set)
# select.mu.z <- mu_selection_regression(Z,X.Z,cond.z,mu.set,R.bins.small, ker.set)
#



# select.mu.y <- mu_selection(mu.set,cond.y,
#                             y.train[,c(1,2)],y.val[,c(1,2)],
#                             R.bins.small)
# select.mu.z <- mu_selection(mu.set,cond.z,
#                             z.train[,c(1,2)],z.val[,c(1,2)],
#                             R.bins.small)
#
# plot(log(mu.set + mu.set[2]), select.mu.y$ce, type = "l")
# plot(log(mu.set + mu.set[2]), select.mu.z$ce, type = "l")
#
# mu.y <- select.mu.y$opt.mu
# mu.z <- select.mu.z$opt.mu

# if already selected these are the cross validated parameters
mu.y <- 0.01
mu.z <- 0.01

unif.4 <- scale_kernel(uniform_kernel,4)
unif.0.5 <- scale_kernel(uniform_kernel,0.5) # groups the categories in the other case
ker.set <- list(unif.4,unif.0.5)


# A.matrix.y <- compute_A_matrix_2(R.bins,cond.y)
# A.tensor.y <- compute_A_tensor_2(R.bins,cond.y)
#
# A.matrix.z <- compute_A_matrix_2(R.bins,cond.z)
# A.tensor.z <- compute_A_tensor_2(R.bins,cond.z)

# model.y <- estimate_mixing_numeric_2(y.train[,c(1,2)], A.matrix.y,
#                                      A.tensor.y, mu.y)
# plot(model.y$latent, type = "l")
#
# model.z <- estimate_mixing_numeric_2(z.train[,c(1,2)], A.matrix.z,
#                                      A.tensor.z, mu.z)
# plot(model.z$latent, type = "l")

##### best conversion to MOCA as a function of cross entropy

cw.data <- read.csv("data/NACCMMSE_to_MOCATOTS_test.csv")


grid.size = 50
mu1.set = seq(0,1, length.out = grid.size)
mu2.set = seq(0,1, length.out = grid.size)

mu.grid = expand.grid(mu1.set, mu2.set)
J = nrow(mu.grid)

Y.train <- y.train[,c("y1", "y2", "age", "group")]
Z.train <- z.train[,c("z1", "z2", "age", "group")]


ref.cols = c("age", "group")
mu.y = 1
mu.z = 1
res = rep(NA,J)
z.test = cw.data$z

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}


subset.idx = seq((id - 1)*25 + 1,(id)*25)

for(j in subset.idx){
  mu1 = mu.grid[j,1]
  mu2 = mu.grid[j,2]
  cat(paste0("smoothing combination: ", j, "/", J), end = "\r")
  pred.dist <- predictedDistributions(cw.data,cw.data$y,
                                      Y.train,Z.train,cond.y,cond.z,
                                      mu1,mu2,ref.cols, ker.set, R.bins = 1000, verbose = F)

  #TODO: fix the cross entropy term in the package then make this not negative
  ce = -empirical_cross_entropy(pred.dist,z.test)
  res[j] = ce
}


results = cbind(mu.grid, res)
results = as.data.frame(results)
colnames(results) = c("mu1", "mu2", "ce")
sim.idx = which(!is.na(results$ce))
results <- results[sim.idx,]
saveRDS(results, paste0("data/mmse_moca_conversion_mu_results",id,".rds"))



plot.results = F

if(plot.results){
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200


  results = readRDS(paste0("data/mmse_moca_conversion_mu_results",1,".rds"))
  for(i in seq(2,100)){
    res.tmp <-readRDS(paste0("data/mmse_moca_conversion_mu_results",i,".rds"))
    results <- rbind(results, res.tmp)
  }
  results <- results[!is.na(results$ce),]
  j.min = which.min(results$ce)

  plt <- ggplot(data = results, aes(mu1, mu2, fill= ce)) +
                geom_tile() +
                scale_fill_gradient(low="yellow", high="blue",
                                    limits=c(min(results$ce), max(results$ce))) +
    geom_point(aes(x=0.01,y=0.01),colour="red") +
    geom_point(aes(x=results[j.min,1],y=results[j.min,2]),colour="black") +
    theme_bw() +
    ggtitle("Smoothing Parameter Conversion Cross Entropy") +
    xlab("\u03bc 1") +
    ylab("\u03bc 2")


  png(filename = "plots/mmse_moca_conversion_prediction.png",
      width = png.width, height = png.height, res = png.res)

  plt
  # Close the pdf file
  dev.off()


  # naive conversion comparison
  compare.conversions = T

  if(compare.conversions){
    mean.y = mean(y.train$y, na.rm = T)
    mean.z = mean(z.train$z, na.rm = T)
    sd.y = sd(y.train$y, na.rm = T)
    sd.z = sd(z.train$z, na.rm = T)

    y.test <- cw.data$y
    z.test <- cw.data$z

    naive.pred <- normalScoreConversionProb(y.test,mean.y,sd.y,mean.z,sd.z, Nz)


    mu.y1 = 0.01
    mu.y2 = results[j.min,1]

    mu.z1 = 0.01
    mu.z2 = results[j.min,2]

    ref.cols = c("age", "group")
    unif.4 <- scale_kernel(uniform_kernel,4)
    unif.0.5 <- scale_kernel(uniform_kernel,0.5) # groups the categories in the other case
    ker.set <- list(unif.4,unif.0.5)
    cond.pred.cov.adj <- predictedDistributions(cw.data,y.test,
                                                y.train,z.train,cond.y,cond.z,
                                                mu.y1,mu.z1,ref.cols, ker.set, R.bins = 1000)

    # picking the optimal parameters from the crosswalk smoothed dataset.
    cond.pred.cov.adj.cw.opt <- predictedDistributions(cw.data,y.test,
                                                       y.train,z.train,cond.y,cond.z,
                                                       mu.y2,mu.z2,ref.cols, ker.set, R.bins = 1000)

    ce.naive = empirical_cross_entropy(naive.pred,z.test)
    ce.test = empirical_cross_entropy(cond.pred.cov.adj,z.test)
    ce.cw.opt = empirical_cross_entropy(cond.pred.cov.adj.cw.opt,z.test)
    print(c(ce.naive,ce.test,ce.cw.opt))
  }

}






# include the comparison naive conversion







