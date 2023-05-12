
#In this script,compare the speeds of the corresponding methods
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
  ker <- gaussian_kernel
} else {
  kernel="Exponential"
  ker <- exponential_kernel
}

R.bins = 1000

# Conditional Models
h.set <- c(0.8,1,2,3,5,10)
H = length(h.set)


#N.seq = c(100,500,1000,5000)
N.set <- c(5,10,20,50,100)
J = length(N.set)




# grid for the values of mu
L = 5
mu.set <- seq(0,2,length.out = L)




# 4 different marginals
marg.1 <- function(N){return(rep(1/(N+1),(N + 1)))}
marg.2 <- function(N){return(-(seq(0,N) - N/2)**2 + 2*(N/2)^2)}
marg.3 <- function(N){return(exp(-(seq(0,N))**2) + exp(- (seq(0,N) - N)**2))}
marg.4 <- function(N){return(sin(seq(0,N)/N * 4* pi) + 1)}

marg.list <- list(marg.1, marg.2, marg.3, marg.4)
K <- length(marg.list)

# dimensions of the results (numquestions,marginal, regularization values, bandwidth)
res <- array(NA,c(J,K,L,H))
time.npem <- res
time.cvxr <- res
like.npem <- res
like.cvxr <- res


uniform.latent <- rep(1/R.bins,R.bins)
# TODO: Also do this as a function of h
for(j in seq(J)){
  N = N.set[j]

  for(h in seq(length(h.set))){
    h.bandwidth = h.set[h]
    cond <- conditional_mkm(N,ker, h.bandwidth)
    # model for the distribution
    A.matrix <- compute_A_matrix(R.bins, cond)

    for(k in seq(K)){
      marg <- marg.list[[k]]
      p.hat <- marg(N)
      p.hat <- p.hat/sum(p.hat)

      cat(paste0("Marginal ", k, "/",K), end = "\r")
      for(i in seq(length(mu.set))){
        mu = mu.set[i]

        time1 <- Sys.time()
        model.npem <- estimate_mixing_npem(p.hat,A.matrix, mu)
        time2 <- Sys.time()
        model.numeric <- estimate_mixing_numeric(p.hat,A.matrix, mu)
        time3 <- Sys.time()

        if(mu == 0){
          like.npem.tmp <- -kl_divergence(p.hat,model.npem$observed)
          like.numeric.tmp <- -kl_divergence(p.hat,model.numeric$observed)

        } else {
          like.npem.tmp <- -kl_divergence(p.hat,model.npem$observed) - mu*kl_divergence(uniform.latent, model.npem$latent)
          like.numeric.tmp <- -kl_divergence(p.hat,model.numeric$observed) - mu*kl_divergence(uniform.latent, model.numeric$latent)
        }
        time.npem[j,k,i,h] <- as.numeric(difftime(time2, time1, units = "secs"))
        time.cvxr[j,k,i,h] <- as.numeric(difftime(time3, time2, units = "secs"))

        like.npem[j,k,i,h] <- like.npem.tmp
        like.cvxr[j,k,i,h] <- like.numeric.tmp

      }
    }
  }
}



saveRDS(time.npem, paste0("data/fitting_speed_npem_results_",kernel,ceiling(id/2), ".rds"))
saveRDS(time.cvxr, paste0("data/fitting_speed_cvxr_results_",kernel,ceiling(id/2), ".rds"))

saveRDS(like.npem, paste0("data/fitting_likelihood_npem_results_",kernel,ceiling(id/2), ".rds"))
saveRDS(like.cvxr, paste0("data/fitting_likelihood_cvxr_results_",kernel,ceiling(id/2), ".rds"))





make.plots = F

if(make.plots){
  library(ggpubr)
  library(abind)


  png.width = 1200
  png.height = 1000
  png.res = 200

  kernel = "Gaussian" # "Gaussian" # "Exponential"

  # grid.parameters
  # Conditional Models
  h.set <- c(0.8,1,2,3,5,10)
  H = length(h.set)


  #N.seq = c(100,500,1000,5000)
  N.set <- c(5,10,20,50,100)
  J = length(N.set)




  # grid for the values of mu
  L = 5
  mu.set <- seq(0,2,length.out = L)

  # two groups, one with mu = 0
  # one with mu = 2


  # TODO: update this section to concatenate the results
  time.npem <- readRDS(paste0("data/fitting_speed_npem_results_",kernel,1, ".rds"))
  time.cvxr <- readRDS(paste0("data/fitting_speed_cvxr_results_",kernel,1, ".rds"))

  like.npem <- readRDS(paste0("data/fitting_likelihood_npem_results_",kernel,1, ".rds"))
  like.cvxr <- readRDS(paste0("data/fitting_likelihood_cvxr_results_",kernel,1, ".rds"))

  dim(time.npem)
  plot.data <- data.frame("N" = c(),
                          "Dist" = c(),
                          "Bandwidth" = c(),
                          "mu" = c(),
                          "Time" = c(),
                          "KL" = c(),
                          "Method" = c())

  # K is the number of the different distributions to fit
  K = dim(time.npem)[2]
  for(j in seq(J)){
      for(k in seq(K)){
        for(l in seq(L)){
          for(h in seq(H)){

            new.row.npem <- c(N.set[j],k,h.set[h],mu.set[l],time.npem[j,k,l,h] ,like.npem[j,k,l,h] ,"NPEM")
            new.row.cvxr <- c(N.set[j],k,h.set[h],mu.set[l],time.cvxr[j,k,l,h] ,like.cvxr[j,k,l,h] ,"CVXR")
            plot.data <- rbind(plot.data,new.row.npem)
            plot.data <- rbind(plot.data,new.row.cvxr)
          }
        }
      }
  }
  colnames(plot.data) = c("N","Dist","Bandwidth",
                        "mu","Time","KL","Method")
  plot.data$N <- factor(plot.data$N, levels = N.set)
  plot.data$Time = as.numeric(plot.data$Time)
  plot.data$KL = as.numeric(plot.data$KL)

  for(k in seq(K)){
    for(l in seq(L)){
        cond.idx = k
        mu.tmp = mu.set[l]
        plot.data.tmp <- plot.data %>% filter(mu == mu.tmp, Dist == cond.idx)

        title <- paste0("Time Comparison ", kernel, " Kernel : \u03bc = ", mu.tmp," --- Conditional: ",cond.idx )
        plot.time <-  ggplot(plot.data.tmp, aes(Bandwidth, log(Time), group = interaction(N, Method),
                       color = N, linetype = Method)) +
          geom_line() + geom_point() +
          ggtitle(title) +
          xlab("Bandwidth (h)")
        plot.time

        png(filename = paste0("plots/time_comparison_mu_",mu.tmp,"_marg_",k,kernel,".png"),
            width = png.width, height = png.height, res = png.res)

        plot.time
        # Close the pdf file
        dev.off()
        title <- paste0("Fit Comparison ", kernel, " Kernel : \u03bc = ", mu.tmp," --- Conditional: ",cond.idx )

        plot.lik <-  ggplot(plot.data.tmp, aes(x = Bandwidth, y = KL, group = interaction(N, Method),
                                                color = N, linetype = Method)) +
          geom_line() + geom_point() +
          ggtitle(title) +
          xlab("Bandwidth (h)")
        plot.lik


        png(filename = paste0("plots/fit_comparison_mu_",mu.tmp,"_marg_",k,kernel,".png"),
            width = png.width, height = png.height, res = png.res)

        plot.lik
        # Close the pdf file
        dev.off()
      }
  }
}








