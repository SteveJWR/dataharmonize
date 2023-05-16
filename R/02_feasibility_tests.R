

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
h.true <- 3

N <- 30

cond.true <- conditional_mkm(N,ker.true, h.true)


two.obs.ratio = 1
n.seq = c(100,500,1000,5000, 10000)

J = length(n.seq)


# grid for the values of h
h.set <- c(0.8,1,2,3,5,10,20)
H = 35
h.set <- exp(seq(log(0.8), log(20), length.out = H))
#H = length(h.set)


if(kernel == "Gaussian"){
  cond.set <- generate_mkm_list(N = N, ker = gaussian_kernel, h.set = h.set)
  cond.names <- paste0("Gaussian h = ",as.character(h.set))
} else {
  cond.set <- generate_mkm_list(N = N, ker = exponential_kernel, h.set = h.set)
  cond.names <- paste0("Exponential h = ",as.character(h.set))
}



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








make.plots = F

if(make.plots){
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200
  n.seq = c(100,500,1000,5000, 10000)

  alpha = 0.01 #
  h.true = 3
  # grid for the values of h
  H = 35
  h.set <- exp(seq(log(0.8), log(20), length.out = H))

  kernel <- "Gaussian" #"Exponential" #"Gaussian"
  # update this section to concatenate the results
  res1 <- readRDS(paste0("data/feasibility_test_order_1_results_",kernel,1, ".rds"))
  res2 <- readRDS(paste0("data/feasibility_test_order_2_results_",kernel,1, ".rds"))


  for(j in seq(2,100)){

    res1.tmp <- readRDS(paste0("data/feasibility_test_order_1_results_",kernel,j, ".rds"))
    res2.tmp <- readRDS(paste0("data/feasibility_test_order_2_results_",kernel,j, ".rds"))

    res1 <- abind(res1, res1.tmp, along = 1)
    res2 <- abind(res2, res2.tmp, along = 1)

  }


  n.sims = dim(res1)[1]
  num.samples = dim(res1)[2]
  H = dim(res1)[3]

  samp.size.vec = rep(paste0("n = ", n.seq), each = H)
  p.vec <- c() # rejection probabilities
  for(j in seq(num.samples)){
    p.vec <- c(p.vec, colMeans(as.matrix(res1[,j,] < alpha), na.rm = T))
  }

  se.vec <- sqrt(p.vec*(1 - p.vec)/n.sims)
  res.data <- data.frame("SampleSize" = samp.size.vec,
                         "Bandwidth" = rep(h.set, times = num.samples),
                         "rejectProb" = p.vec,
                         "rejectProb_sd" = se.vec)

  res.data$SampleSize <- factor(res.data$SampleSize, levels = paste0("n = ", n.seq))
  plt.1 <- ggplot(res.data, aes(x = log(Bandwidth), y = rejectProb, color = SampleSize)) +
    geom_line()  +
    geom_vline(aes(xintercept = log(h.true)), linetype = "dashed") +
    geom_ribbon(aes(ymin = rejectProb - 2*rejectProb_sd, ymax = rejectProb + 2*rejectProb_sd, fill = SampleSize),alpha=0.3) +
    ggtitle(paste0(kernel, " Kernel First Order Feasibility Test : \u03B1 = ", alpha )) +
    xlab("Bandwidth (h)") +
    ylab("Power")

  #geom_line(aes(x = n, y = rmse, color = method)) #+
  #geom_errorbar(aes(ymin = bias - 2*rmse, ymax = bias + 2*rmse))

  plt.1

  png(filename = paste0("plots/feasibility_test_first_order",kernel,".png"),
      width = png.width, height = png.height, res = png.res)

  plt.1
  # Close the pdf file
  dev.off()




  n.sims = dim(res2)[1]
  num.samples = dim(res2)[2]
  H = dim(res2)[3]

  samp.size.vec = rep(paste0("n = ", n.seq), each = H)
  p.vec <- c() # rejection probabilities
  for(j in seq(num.samples)){
    p.vec <- c(p.vec, colMeans(as.matrix(res2[,j,] < alpha), na.rm = T))
  }

  se.vec <- sqrt(p.vec*(1 - p.vec)/n.sims)
  res.data <- data.frame("SampleSize" = samp.size.vec,
                         "Bandwidth" = rep(h.set, times = num.samples),
                         "rejectProb" = p.vec,
                         "rejectProb_sd" = se.vec)

  res.data$SampleSize <- factor(res.data$SampleSize, levels = paste0("n = ", n.seq))
  plt.2 <- ggplot(res.data, aes(x = log(Bandwidth), y = rejectProb, color = SampleSize)) +
    geom_line()  +
    geom_vline(aes(xintercept = log(h.true)), linetype = "dashed") +
    geom_ribbon(aes(ymin = rejectProb - 2*rejectProb_sd, ymax = rejectProb + 2*rejectProb_sd, fill = SampleSize),alpha=0.3) +
    ggtitle(paste0(kernel, " Kernel Second Order Feasibility Test : \u03B1 = ", alpha )) +
    xlab("log-Bandwidth (h)") +
    ylab("Power")

  #geom_line(aes(x = n, y = rmse, color = method)) #+
  #geom_errorbar(aes(ymin = bias - 2*rmse, ymax = bias + 2*rmse))

  plt.2

  png(filename = paste0("plots/feasibility_test_second_order",kernel,".png"),
      width = png.width, height = png.height, res = png.res)

  plt.2
  # Close the pdf file
  dev.off()


}


















