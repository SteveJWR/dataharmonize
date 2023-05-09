

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

saveRDS(res, paste0("data/regularization_selection_results_",kernel,ceiling(id/2), ".rds"))





make.plots = F

if(make.plots){
  library(ggpubr)
  library(abind)

  kernel = "Gaussian" # "Exponential"
  png.width = 1200
  png.height = 1000
  png.res = 200

  # grid.parameters
  n.seq = c(100,500,1000,5000)
  J = length(n.seq)
  # grid for the values of mu
  L = 200
  mu.set <- seq(0,2,length.out = L)

  # update this section to concatenate the results
  res <- readRDS(paste0("data/regularization_selection_results_",kernel,1, ".rds"))

  for(j in seq(2,200)){

    res.tmp <- readRDS(paste0("data/regularization_selection_results_",kernel,j, ".rds"))

    res <- rbind(cover.cc, cover.cc.tmp)

  }

  n.sims = nrow(cover.cc)
  n.set

  cc.mean.bias <- colMeans(dev.cc, na.rm = T)
  mean.bias <- colMeans(dev, na.rm = T)
  mean.bias.cov.adj <- colMeans(dev.cov.adj, na.rm = T)
  mean.bias.true.latent <- colMeans(dev.true.latent, na.rm = T)
  mean.bias.bootstrap <- colMeans(dev.bootstrap, na.rm = T)
  z.score.mean.bias <- colMeans(dev.z.score, na.rm = T)
  quantile.bias <- colMeans(dev.quantile, na.rm = T)

  cc.rmse <- sqrt(colMeans(abs(dev.cc)^2, na.rm = T))
  rmse <-sqrt( colMeans(abs(dev)^2, na.rm = T))
  rmse.cov.adj <-sqrt( colMeans(abs(dev.cov.adj)^2, na.rm = T))
  rmse.true.latent <-sqrt( colMeans(abs(dev.true.latent)^2, na.rm = T))
  rmse.bootstrap <-sqrt( colMeans(abs(dev.bootstrap)^2, na.rm = T))
  z.score.rmse <- sqrt(colMeans(abs(dev.z.score)^2, na.rm = T))
  quantile.rmse <-sqrt( colMeans(abs(dev.quantile)^2, na.rm = T))

  cc.rmse.sd <- colSDs(dev.cc, na.rm = T)/sqrt(n.sims)
  rmse.sd <- colSDs(dev, na.rm = T)/sqrt(n.sims)
  rmse.cov.adj.sd <- colSDs(dev.cov.adj, na.rm = T)/sqrt(n.sims)
  rmse.true.latent.sd <- colSDs(dev.true.latent, na.rm = T)/sqrt(n.sims)
  rmse.bootstrap.sd <- colSDs(dev.bootstrap, na.rm = T)/sqrt(n.sims)
  z.score.rmse.sd <- colSDs(dev.z.score, na.rm = T)/sqrt(n.sims)
  quantile.rmse.sd <- colSDs(dev.quantile, na.rm = T)/sqrt(n.sims)


  res.data <- data.frame("method" = c(rep("Complete Case", length(n.set)),
                                      rep("DNOISE", length(n.set)),
                                      rep("DNOISE (cov.adj.)", length(n.set)),
                                      rep("DNOISE (T.L.)", length(n.set)),
                                      rep("DNOISE (Bootstrap)", length(n.set)),
                                      rep("Z Score", length(n.set)),
                                      rep("Quantile", length(n.set))),
                         "n" = c(n.set,n.set,n.set,
                                 n.set,n.set,n.set, n.set),
                         "bias" = c(cc.mean.bias,mean.bias, mean.bias.cov.adj,
                                    mean.bias.true.latent, mean.bias.bootstrap, z.score.mean.bias,quantile.bias),
                         "rmse" = c(cc.rmse,rmse,rmse.cov.adj,
                                    rmse.true.latent, rmse.bootstrap,
                                    z.score.rmse,quantile.rmse),
                         "rmse_sd" = c(cc.rmse.sd,rmse.sd,rmse.cov.adj.sd,
                                       rmse.true.latent.sd,
                                       rmse.true.latent.sd,
                                       z.score.rmse.sd,quantile.rmse.sd))


  plt.bias <- ggplot(res.data, aes(x = log(n), y = bias, color = method)) +
    geom_line()  #+
  #geom_line(aes(x = n, y = rmse, color = method)) #+
  #geom_errorbar(aes(ymin = bias - 2*rmse, ymax = bias + 2*rmse))

  plt.bias

  png(filename = "plots/sim_binomial_bias.png",
      width = png.width, height = png.height, res = png.res)

  plt.bias
  # Close the pdf file
  dev.off()

  plt.rmse <- ggplot(res.data, aes(x = log(n), y = log(rmse), color = method)) +
    geom_line() +
    geom_errorbar(aes(ymin = log(rmse - 2*rmse_sd), ymax = log(rmse + 2*rmse_sd)))

  plt.rmse
  png(filename = "plots/sim_binomial_rmse.png",
      width = png.width, height = png.height, res = png.res)

  plt.rmse
  # Close the pdf file
  dev.off()



  cc.coverage <- colMeans(cover.cc, na.rm = T)
  coverage <- colMeans(cover, na.rm = T)
  coverage.cov.adj <- colMeans(cover.cov.adj, na.rm = T)
  coverage.true.latent <- colMeans(cover.true.latent, na.rm = T)
  coverage.bootstrap <- colMeans(cover.bootstrap, na.rm = T)
  z.score.coverage <- colMeans(cover.z.score, na.rm = T)
  quantile.coverage  <- colMeans(cover.quantile, na.rm = T)

  cov.vec <- c(cc.coverage,coverage,coverage.cov.adj,
               coverage.true.latent, coverage.bootstrap, z.score.coverage,quantile.coverage)
  cov.error <- sqrt(cov.vec*(1 - cov.vec)/sum(!is.na(cover[,1])))

  cov.data <- data.frame("method" = c(rep("Complete Case", length(n.set)),
                                      rep("DNOISE", length(n.set)),
                                      rep("DNOISE (cov. adj.)", length(n.set)),
                                      rep("DNOISE (T.L.)", length(n.set)),
                                      rep("DNOISE (Bootstrap)", length(n.set)),
                                      rep("Z Score", length(n.set)),
                                      rep("Quantile", length(n.set))),
                         "n" = c(n.set,n.set,n.set,n.set, n.set, n.set, n.set),
                         "coverage" = cov.vec,
                         "error" = cov.error)


  #cov.data <- cov.data %>% filter(n != 500)
  plt.coverage <- ggplot(cov.data, aes(x = log(n), y = coverage, color = method)) +
    geom_line(position=position_dodge(width=0.2)) +
    geom_errorbar(aes(ymin = coverage - 2*error, ymax = coverage + 2*error), width=0.5,
                  linewidth=0.5, position=position_dodge(width=0.2)) +
    geom_hline(yintercept=0.95, linetype='dotted', col = 'black')

  plt.coverage

  png(filename = "plots/sim_binomial_coverage.png",
      width = png.width, height = png.height, res = png.res)

  plt.coverage
  # Close the pdf file
  dev.off()
}


