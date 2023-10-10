

# model 1
# model 2

png.width = 1200
png.height = 1000
png.res = 200 

n.set <- c(100,200,500,1000, 2000, 5000)
rho = 0.1
#theta = 0.5
n.true <- 10000000



n.sims <- 500

Ny <- 15
Nz <- 35 


cond.y <- generate_cond_binomial(Ny)
cond.z <- generate_cond_binomial(Nz)
# h1 <- 2
# h2 <- 0.5
# cond.y <- conditional_mkm(Ny, gaussian_kernel, h1)
# cond.z <- conditional_mkm(Nz, gaussian_kernel, h2)

# uniform age sample
X <- round(runif(n.true, min = 54.5, max = 80.5))
lat.mu <- -1/20*(X - 50) +  2.5

U <- runif(n.true)
gamma.true.seq <- logistic(lat.mu + (U <= 1/2)*(2 + rnorm(length(X),sd = 1)) -(U > 1/2)*(2 + rnorm(length(X),sd = 1))) 
#gamma.true.seq <- logistic(lat.mu +rnorm(length(X),sd = 1))
# plot(hist(gamma.true.seq))
# plot(hist(sqrt(gamma.true.seq)))
#z <- simulate_test_cond(rep(c(2,2),n.true/2),cond.z,gamma.true.seq)
z <- rbinom(length(gamma.true.seq),Nz, gamma.true.seq)
#z <- simulate_test_cond_uni(cond.z ,gamma.true.seq)
dat.pop <- data.frame("Z" = z, "age" = X)

mod.full.pop <- glm(Z ~ age, data = dat.pop)
beta.true <- mod.full.pop$coefficients[2]
beta.true
plot(hist(z,breaks = 36))


#TODO: add a coverage result 
#TODO: Add also naive z-score and quantile matching 


mu.y = 0.1 #0.3
mu.z = 0.1 #0.3
n.impute = 50 
fmla <- formula(outcome ~ age )
# we want to plot coverage and bias of conversion
# TODO: Try again in case gamma2 is bimodal
n.sims = 500
beta.set <- rep(NA,n.sims)
cover.set <- rep(NA,n.sims)
for(sim in 1:n.sims){
  cat(paste0("Number of sims: ", sim, "/",n.sims), end = "\n")

  n = 5000
  X <- round(runif(n, min = 54.5, max = 80.5))

  lat.mu2 <- -1/20*(X - 50) +  2.5
  U <- runif(n)
  gamma2 <- logistic(lat.mu2 +  (U <= 1/2)*(2 + rnorm(length(X),sd = 1)) -(U > 1/2)*(2 + rnorm(length(X),sd = 1)))

  z1 <- rbinom(length(gamma2),Nz, gamma2)
  
  cc.sim.data <- data.frame("Z" = z1, "age" = X)
  
  mod.full.pop <- glm(Z ~ age, data = cc.sim.data)
  beta.hat.cc <- mod.full.pop$coefficients[2]
  sd.cc <- sqrt(sandwich(mod.full.pop)[2,2])
  #quantile matching 
  print(beta.hat.cc)
  print(sd.cc)
  beta.set[sim] = beta.hat.cc 
  cover.set[sim] = beta.hat.cc - 2.0*sd.cc  <= beta.true & beta.hat.cc + 2.0*sd.cc >= beta.true
}

plot(hist(beta.set))
mean(cover.set)
vline()
#########

# dev.cc <- readRDS("data/06c_mkm_dev_cc.rds")
# dev <- readRDS("data/06c_mkm_dev.rds")
# dev.z.score <- readRDS("data/06c_mkm_dev_zscore.rds")
# dev.quantile <- readRDS("data/06c_mkm_dev_quantile.rds")
# 
# cover.cc <- readRDS("data/06c_mkm_cover_cc.rds")
# cover <- readRDS("data/06c_mkm_cover.rds")
# cover.z.score <- readRDS("data/06c_mkm_cover_zscore.rds")
# cover.quantile <- readRDS("data/06c_mkm_cover_quantile.rds")
# 


cc.mean.bias <- colMeans(dev.cc, na.rm = T)
mean.bias <- colMeans(dev, na.rm = T)
z.score.mean.bias <- colMeans(dev.z.score, na.rm = T)
quantile.bias <- colMeans(dev.quantile, na.rm = T)

cc.rmse <- sqrt(colMeans(abs(dev.cc)^2, na.rm = T))
rmse <-sqrt( colMeans(abs(dev)^2, na.rm = T))
z.score.rmse <- sqrt(colMeans(abs(dev.z.score)^2, na.rm = T))
quantile.rmse <-sqrt( colMeans(abs(dev.quantile)^2, na.rm = T))

res.data <- data.frame("method" = c(rep("Complete Case", length(n.set)), 
                                    rep("DNOIS", length(n.set)), 
                                    rep("Z Score", length(n.set)), 
                                    rep("Quantile", length(n.set))), 
                       "n" = c(n.set,n.set,n.set,n.set),
                       "bias" = c(cc.mean.bias,mean.bias,z.score.mean.bias,quantile.bias),
                       "rmse" = c(cc.rmse,rmse,z.score.rmse,quantile.rmse))


plt.bias <- ggplot(res.data, aes(x = log(n), y = bias, color = method)) + 
  geom_line() #+ 
#geom_line(aes(x = n, y = rmse, color = method)) #+ 
#geom_errorbar(aes(ymin = bias - 2*rmse, ymax = bias + 2*rmse))


png(filename = "Plots/sim_binomial_bias.png",
    width = png.width, height = png.height, res = png.res)

plt.bias
# Close the pdf file
dev.off() 

plt.rmse <- ggplot(res.data, aes(x = log(n), y = log(rmse), color = method)) + 
  geom_line()

png(filename = "Plots/sim_binomial_rmse.png",
    width = png.width, height = png.height, res = png.res)

plt.rmse
# Close the pdf file
dev.off() 



cc.coverage <- colMeans(cover.cc, na.rm = T)
coverage <- colMeans(cover, na.rm = T)
z.score.coverage <- colMeans(cover.z.score, na.rm = T)
quantile.coverage  <- colMeans(cover.quantile, na.rm = T)

cov.vec <- c(cc.coverage,coverage,z.score.coverage,quantile.coverage)
cov.error <- sqrt(cov.vec*(1 - cov.vec)/sum(!is.na(cover[,1])))

cov.data <- data.frame("method" = c(rep("Complete Case", length(n.set)), 
                                    rep("DNOIS", length(n.set)), 
                                    rep("Z Score", length(n.set)), 
                                    rep("Quantile", length(n.set))), 
                       "n" = c(n.set,n.set,n.set,n.set),
                       "coverage" = cov.vec,
                       "error" = cov.error)


plt.coverage <- ggplot(cov.data, aes(x = log(n), y = coverage, color = method)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = coverage - 2*error, ymax = coverage + 2*error), width=0.1,
                size=0.5) +
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black')

png(filename = "Plots/sim_binomial_coverage.png",
    width = png.width, height = png.height, res = png.res)

plt.coverage
# Close the pdf file
dev.off() 



`### Current setup is pretty good, 
#TODO: Add this idea of using quantile matching 
# or a generalized optimal transport method. 

# colSDs <- function(X, ...){
#   X.means <- colMeans(X, ...)
#   v.out <- sqrt(colMeans((X.means - X)^2, ...))
#   return(v.out)
# }





# expectile curve, 




Ny = 200
Nz = 200
n = 500000
X <- round(runif(n, min = 54.5, max = 80.5))

lat.mu2 <- -1/20*(X - 50) +  2.5
#lat.mu2 <- (-1/5*(X - 50) +  2.5)*sin(X)



U <- runif(n)
gamma2 <- logistic(lat.mu2 +  (U <= 1/2)*(2 + rnorm(length(X),sd = 1)) -(U > 1/2)*(2 + rnorm(length(X),sd = 1)))
gamma2 <- logistic(lat.mu2 +  (U <= 1/2)*( rnorm(length(X),sd = 1)) -(U > 1/2)*( rnorm(length(X),sd = 1)))
#gamma2 <- logistic(lat.mu2 +rnorm(length(X),sd = 1)) 
# one to one mapping of latent traits
gamma1 <- gamma2**2
#gamma1 <- gamma2*(1/5) + 4/5
# true.mix.y <- hist(gamma1, breaks = seq(0,1,length.out = 1001))$density
# true.mix.z <- hist(gamma2, breaks = seq(0,1,length.out = 1001))$density
# true.mix.y <- true.mix.y/sum(true.mix.y)
# true.mix.z <- true.mix.z/sum(true.mix.z)
# plot(hist(gamma1, breaks = 50))
# plot(hist(gamma2, breaks = 50))
y1 <- rbinom(length(gamma1),Ny, gamma1)
z1 <- rbinom(length(gamma2),Nz, gamma2)




quantile.vec <- rep(NA,Ny + 1)
expectile.vec <- rep(NA,Ny + 1)




for(y.idx in seq(0,Ny)){
  print(y.idx)
  which.idx <- which(y1 == y.idx)
  
  cond.mean <- mean(z1[which.idx], na.rm = T)
  expectile.vec[y.idx + 1] = cond.mean
  
  p <- mean(y1 <= y.idx)
  #q <- ceiling(quantile(z1, p, type = 3))
  q <- quantile(z1, p, na.rm = T)
  
  quantile.vec[y.idx + 1] = q
}


plot(seq(0,Ny), quantile.vec, type = "l", col = "blue")
lines(seq(0,Ny), expectile.vec, type = "l", col = "green")
legend(60, 60, legend=c("Quantile (OT)", "Conditional Expectation"),
       col=c("blue", "green"), lty=c(1,1), cex=0.8)







library(mvtnorm)
theta.seq <- c(-1,-0.5,0,.3,.5,.8,.9,1)


n = 1000000
Ny = 60
Nz = 60
conditional.expect <- matrix(NA, nrow = length(theta.seq), ncol = Ny + 1)

for(i in seq(length(theta.seq))){
  theta = theta.seq[i]
  sigma <- matrix(c(1,theta,theta,1), 2,2)
  X <- rmvnorm(n, sigma = sigma)
  p1 <- pnorm(X[,1])
  p2 <- pnorm(X[,2])
  y <- qbinom(p1,Ny, .5)
  z <- qbinom(p2,Nz, .8)
  expectile.vec <- rep(NA,Ny + 1)
  for(y.idx in seq(0,Ny)){
    print(y.idx)
    which.idx <- which(y == y.idx)
    
    cond.mean <- mean(z[which.idx], na.rm = T)
    expectile.vec[y.idx + 1] = cond.mean
  }
  conditional.expect[i,] = expectile.vec
}




quantile.vec <- rep(NA,Ny + 1)
for(y.idx in seq(0,Ny)){
  print(y.idx)
  which.idx <- which(y == y.idx)
  p <- mean(y <= y.idx)
  q <- quantile(z, p, na.rm = T)
  quantile.vec[y.idx + 1] = q
}



lty=2

plot(seq(0,Ny), quantile.vec, type = "l", col = "blue",xlim = c(0,Ny*4/2), ylim = c(30,Nz), xlab = "y", ylab = "z")
lines(seq(0,Ny), conditional.expect[1,], type = "l", col = "red", lty = 1)
lines(seq(0,Ny), conditional.expect[2,], type = "l", col = "red", lty = 2)
lines(seq(0,Ny), conditional.expect[3,], type = "l", col = "purple", lty = 1)
lines(seq(0,Ny), conditional.expect[4,], type = "l", col = "green", lty = 1)
lines(seq(0,Ny), conditional.expect[5,], type = "l", col = "green", lty = 2)
lines(seq(0,Ny), conditional.expect[6,], type = "l", col = "green", lty = 3)
lines(seq(0,Ny), conditional.expect[7,], type = "l", col = "green", lty = 4)
lines(seq(0,Ny), conditional.expect[8,], type = "l", col = "darkgreen", lty = 5)
legend(Ny, Nz, legend=c("Quantile (OT)", "Conditional (theta = -1)",
                        "Conditional (theta = -.5)",
                        "Conditional (theta = 0)",
                        "Conditional (theta = 0.3)",
                        "Conditional (theta = 0.5)",
                        "Conditional (theta = 0.8)",
                        "Conditional (theta = 0.9)",
                        "Conditional (theta = 1)"),
       col=c("blue","red", "red", "purple", 
             rep("green", 4), "darkgreen"), lty=c(1,1,2,1,1:6), cex=0.8)





