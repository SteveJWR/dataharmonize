

### true model

#n2 number of paired observations 
#n1 number of single observations
#mixture a betas latent model. 
x <- seq(0,1,length.out = R.bins)
mixture.bimodal <- (2/3)*dbeta(x,35,15) + (1/3)*dbeta(x,15,35)
mixture.bimodal <- mixture.bimodal/sum(mixture.bimodal)

plot(x,mixture.bimodal)
plot(seq(0,R.bins -1)/R.bins,mixture.bimodal, type = "l")

h.true <- 2
ker.true <- gaussian_kernel
N <- 30
cond.true <- conditional_mkm(N,ker.true, h.true)



n1 = 600
n2 = 200 
Y.sim <- simulate_test(n2,n1,cond.true,mixture.bimodal)



# beta1 <- 6
# beta2 <- 6
# 
# 
# sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
#                           beta1.y = beta1.model2, beta2.y = beta2.model2, 
#                           cond.y = cond.model2, pair.obs = F)


h.set <- c(0.5,0.8,1,1.5,2,2.5,3,3.5,4,5,10,15,20)


# the list of models to choose
gaussian.cond.list <- generate_mkm_list(N = N, ker = gaussian_kernel, h.set = h.set)
gaussian.cond.names <- paste0("Gaussian Kernel, h = ", round(h.set,2))

exponential.cond.list <- generate_mkm_list(N = N, ker = exponential_kernel, h.set = h.set)
exponential.cond.names <- paste0("Exponential Kernel, h = ", round(h.set,2))

triangle.cond.list <- generate_mkm_list(N = N, ker = triangle_kernel_smooth, h.set = h.set)
triangle.cond.names <- paste0("Triangle Kernel, h = ", round(h.set,2))

epanechnikov.cond.list <- generate_mkm_list(N = N, ker = epanechnikov_kernel_smooth, h.set = h.set)
epanechnikov.cond.names <- paste0("Epanechnikov Kernel, h = ", round(h.set,2))

err.gaussian <- error_model_selection(gaussian.cond.list,Y.sim,R.bins)
err.exponential <- error_model_selection(exponential.cond.list,Y.sim,R.bins)
err.triangle <- error_model_selection(triangle.cond.list,Y.sim,R.bins)
err.epanechnikov <- error_model_selection(epanechnikov.cond.list,Y.sim,R.bins)

plot(log(h.set), err.gaussian$lik_vals, type = "l")
lines(log(h.set), err.exponential$lik_vals, type = "l", col = "red")
lines(log(h.set), err.triangle$lik_vals, type = "l", col = "blue")
lines(log(h.set), err.epanechnikov$lik_vals, type = "l", col = "green")



n1 = 800
n2 = 100
Y.sim <- simulate_test(n2,n1,cond.true,mixture.bimodal)



mu.set <- c(0,0.00001,0.0000001,0.0001,0.001,0.01,0.1,1,10, 1000)
mu.set <- c(0,10**seq(-5,0,length.out = 12))
Y.val <-  simulate_test(n2 + n1,0,cond.true,mixture.bimodal)

cond <- cond.true
lik.mu <- mu_selection(mu.set, cond, Y.sim,Y.val, R.bins)
plot(log(mu.set + 10**(-5)),lik.mu$ce, type = "l")

# This tends to work 

n1 = 700
n2 = 300
set.seed(2)
##### introduce covariates

X.sim.age <- sample(seq(50,80), (n1 + n2), replace = T)
X.sim.educ <- sample(seq(0,1), (n1 + n2), prob = c(3/5,2/5), replace = T)
X.sim.gen <- sample(seq(0,1), (n1 + n2), prob = c(3/5,2/5), replace = T)

X.sim <- matrix(c(X.sim.age,X.sim.educ,X.sim.gen), ncol = 3)

##### introduce covariates
X.val.age <- sample(seq(50,80), (n2), replace = T)
#X.val.age <- sample(c(50,55,60,65,70,75,80), (n2), replace = T)
X.val.educ <- sample(seq(0,1), (n2), prob = c(3/5,2/5), replace = T)
X.val.gen <- sample(seq(0,1), (n2), prob = c(3/5,2/5), replace = T)

X.val <- matrix(c(X.val.age,X.val.educ,X.val.gen), ncol = 3)


# 
mu.lat <- 2 + -(1/2)*(X.val[,1] - 50)^(1/2) + 0.2*(X.val[,2]) + 0.2*(X.val[,3])
sd.lat <- 0.3
gamma.val <- logistic(rnorm(length(mu.lat),mu.lat, sd.lat))
plot(X.val[,1],gamma.val)

z <- seq(-10,10,length.out = 10000)
plot(logistic(z),dnorm(z,mu.lat[2], sd.lat)*(1/(logistic(z)*(1-logistic(z)))))

#plot(gamma.val)

mu.lat <- 2 + -(1/2)*(X.sim[,1] - 50)^(1/2) + 0.2*(X.sim[,2]) + 0.2*(X.sim[,3])
sd.lat <- 0.3
gamma.sim <- logistic(rnorm(length(mu.lat),mu.lat, sd.lat))
plot(X.sim[,1],gamma.sim)

Y.sim <- sapply(gamma.sim,function(x){
    p.vec <- cond(x)
    y <- sample(seq(0,length(p.vec) - 1), 2, prob = p.vec, replace = T)
    return(y)
    })
Y.sim <- t(Y.sim)

Y.val <- sapply(gamma.val,function(x){
  p.vec <- cond(x)
  y <- sample(seq(0,length(p.vec) - 1), 2, prob = p.vec, replace = T)
  return(y)
})
Y.val <- t(Y.val)

unif.3 <- scale_kernel(uniform_kernel,3)
unif.0.5 <- scale_kernel(uniform_kernel,0.5)
ker.set <- list(unif.3,unif.0.5,unif.0.5)



R.bins = 1000

R.bins = 500
# fitting the model for every mu at each 
# term might in fact be too slow 



mu.set <- c(0,10**seq(-3,3,length.out = 31))

cond <- cond.true
lik.mu <- mu_selection(mu.set,cond,
                       Y.sim,Y.val,
                       R.bins)

plot(log(mu.set + 10**(-3)),lik.mu$ce, type = "l")



R.bins = 300
mu.reg.set <- c(0,10**seq(-3,1,length.out = 13))

cond <- cond.true
lik.mu.reg <- mu_selection_regression(Y.sim, X.sim, 
                                      Y.val, X.val, 
                                      cond, mu.reg.set, R.bins, ker.set)

plot(log(mu.reg.set + 10**(-3)),lik.mu.reg$ce, type = "l")






