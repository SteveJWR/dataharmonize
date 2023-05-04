


# the new selection criteria


library(CVXR)
#library(tensorr)
source("R/01_functions.R")

# N = 30
# n = 600 
# R.bins <- 1000
# cond <- conditional_mkm(N, gaussian_kernel, 2)
# 
# cond2 <- conditional_mkm(N, exponential_kernel, 5)
# cond3 <- generate_cond_binomial(N)
# 
# cond.set <- list(cond, cond2, cond3)
# 
# Y <- matrix(NA, nrow = n, ncol = 2)
# Y[,1] <- rbinom(n,N, 0.5)
# Y[seq(n/2),2] <- rbinom(n/2,N, 0.5)
# Y.nm <- Y[seq(n/2),]

tensor_prod <- function(A.tensor, mixture){
  out <- apply(A.tensor, c(1,2), function(Z){
    t(Z) %*% mixture
  })
  return(out)
}


uniform_ce <- function(mixture){
  R.bins <- length(mixture)
  return(sum(log(R.bins*mixture)/R.bins))
}

# weights are in general a function of x
population_bivariate_likelihood <- function(Y, A.matrix, A.tensor, mixture, mu, weights){
  if(missing(weights)){
    weights <- rep(1, nrow(Y))
  }
  
  if(missing(mu)){
    mu = 0
  }
  count.vector <- rep(0,N + 1)
  count.matrix <- matrix(0,N+1, N+1)
  
  #NA's can be allowed in second column only
  for(i in seq(nrow(Y))){
    if(!any(is.na(Y[i,]))){
      count.matrix[Y[i,1] + 1,Y[i,2] + 1] <- count.matrix[Y[i,1] + 1,Y[i,2] + 1] + weights[i]
    } else {
      count.vector[Y[i,1] + 1] <- count.vector[Y[i,1] + 1] + weights[i]
    }
  }
  
  p.ma <- as.vector(A.matrix %*% mixture)
  p.ma2 <- tensor_prod(A.tensor,mixture)
  
  logl <-  sum(t(count.matrix) %*% log(p.ma)) + sum(diag( t(count.matrix) %*% log(p.ma2))) + mu*uniform_ce(mixture)
  
  return(logl)
 
  
}


error_model_selection <- function(cond.set, Y, R.bins, cond.names){
  if(missing(cond.names)){
    cond.names = NULL
  }
  res <- rep(NA, length(cond.set))
  for(i in seq(length(cond.set))){
    cat(paste("Model",i,"/",length(cond.set)), end = "\r")
    cond <- cond.set[[i]]
    A.matrix <- compute_A_matrix_2(R.bins, cond)
    A.tensor <- compute_A_tensor_2(R.bins, cond)
    lat.list <- estimate_mixing_numeric_2(Y,A.matrix,A.tensor, mu = 0)
    mixture <- lat.list$latent
    res[i] <- population_bivariate_likelihood(Y,A.matrix,A.tensor,mixture, mu = 0)
  }
  i.max <- which.max(res)
  out.list <- list("opt_model" = cond.set[[i.max]], "lik_vals" = res, "opt_model_name" = cond.names[i.max])
  return(out.list)
} 


mu_selection <- function(mu.set, cond, Y, R.bins){
  res <- rep(NA, length(mu.set))
  A.matrix <- compute_A_matrix_2(R.bins, cond)
  A.tensor <- compute_A_tensor_2(R.bins, cond)
  for(i in seq(length(mu.set))){
    cat(paste("Model",i,"/",length(mu.set)), end = "\r")
    mu <- mu.set[[i]]
    lat.list <- estimate_mixing_numeric_2(Y,A.matrix,A.tensor, mu = mu)
    mixture <- as.vector(lat.list$latent)
    res[i] <- conversion_cross_entropy(Y,A.matrix,A.tensor,mixture)
  }
  i.min <- which.min(res)
  out.list <- list("opt.mu" = mu.set[[i.min]], "ce" = res)
  return(out.list)
} 

mu_selection <- function(mu.set, cond, Y.train, Y.val, R.bins){
  res <- rep(NA, length(mu.set))
  A.matrix <- compute_A_matrix_2(R.bins, cond)
  A.tensor <- compute_A_tensor_2(R.bins, cond)
  for(i in seq(length(mu.set))){
    cat(paste("Model",i,"/",length(mu.set)), end = "\r")
    mu <- mu.set[[i]]
    lat.list <- estimate_mixing_numeric_2(Y.train,A.matrix,A.tensor, mu = mu)
    mixture <- as.vector(lat.list$latent)
    res[i] <- conversion_cross_entropy(Y.val,A.matrix,A.tensor,mixture)
  }
  i.min <- which.min(res)
  out.list <- list("opt.mu" = mu.set[[i.min]], "ce" = res)
  return(out.list)
} 




# includes regression component implicitly through the mixture. 
conversion_cross_entropy <- function(Y, A.matrix, A.tensor, mixture){
  if(any(is.na(Y))){
    stop("Y cannot have any missing values")
  }
  if(is.vector(mixture)){
    p.ma <- as.vector(A.matrix %*% mixture)
    p.ma2 <- tensor_prod(A.tensor, mixture)
    # first observation is row, second observation is column 
    p.ma.cond <- p.ma2/p.ma
    out <- 0
    for(i in seq(nrow(Y))){
      # if(length(log(p.ma.cond[Y[i,1] + 1,Y[i,2] + 1])) == 0 ){
      #   print(i)
      # }
      out <- out - log(p.ma.cond[Y[i,1] + 1,Y[i,2] + 1])
    }
  } else if(is.matrix(mixture)){
    if(nrow(mixture) != nrow(Y)){
      stop("mixture matrix must have matching number of rows to Y")
    } else if(ncol(mixture) != ncol(A.matrix)){
      stop("mixture matrix must have matching number of columns to A (i.e. bin number does not match)")
    }
    out <- 0
    for(i in seq(nrow(Y))){
      mixture.row <- as.vector(mixture[i,])
      p.ma <- as.vector(A.matrix %*% mixture.row)
      p.ma2 <- tensor_prod(A.tensor, mixture.row)
      # first observation is row, second observation is column 
      p.ma.cond <- p.ma2/p.ma
      out <- out - log(p.ma.cond[Y[i,1] + 1,Y[i,2] + 1])
    }
  }
  return(out)
}


conversion_cross_entropy_2 <- function(Y, A.matrix, A.tensor, mixture){
  if(any(is.na(Y))){
    stop("Y cannot have any missing values")
  }
  if(is.vector(mixture)){
    p.ma <- as.vector(A.matrix %*% mixture)
    p.ma2 <- tensor_prod(A.tensor, mixture)
    # first observation is row, second observation is column 
    p.ma.cond <- p.ma2/p.ma
    out <- 0
    for(i in seq(nrow(Y))){
      # if(length(log(p.ma.cond[Y[i,1] + 1,Y[i,2] + 1])) == 0 ){
      #   print(i)
      # }
      out <- out - log(p.ma2[Y[i,1] + 1,Y[i,2] + 1])
    }
  } else if(is.matrix(mixture)){
    if(nrow(mixture) != nrow(Y)){
      stop("mixture matrix must have matching number of rows to Y")
    } else if(ncol(mixture) != ncol(A.matrix)){
      stop("mixture matrix must have matching number of columns to A (i.e. bin number does not match)")
    }
    out <- 0
    for(i in seq(nrow(Y))){
      mixture.row <- as.vector(mixture[i,])
      p.ma <- as.vector(A.matrix %*% mixture.row)
      p.ma2 <- tensor_prod(A.tensor, mixture.row)
      # first observation is row, second observation is column 
      p.ma.cond <- p.ma2/p.ma
      out <- out - log(p.ma2[Y[i,1] + 1,Y[i,2] + 1])
    }
  }
  return(out)
}



# measuring the difference in likelihood value to the true 
# i.e. we want the true maximizer to be no more than 

select_R <- function(Y, cond, mu = 0.001, R.max = 30000, threshold = 10**(-4), stepsize = 5){
  
  # generally 100 is quite smooth
  R.start <- 100
  R.prev <- R.start
  ratio <- Inf
  A.matrix <- compute_A_matrix_2(R.max, cond)
  A.tensor <- compute_A_tensor_2(R.max, cond)
  lat.list <- estimate_mixing_numeric_2(Y,A.matrix,A.tensor, mu = mu)
  mixture <- lat.list$latent
  exact.lik <- population_bivariate_likelihood(Y,A.matrix,A.tensor,mixture, mu = mu, weights = NULL)
  
  
  while(ratio > threshold){
    A.matrix <- compute_A_matrix_2(R.prev, cond)
    A.tensor <- compute_A_tensor_2(R.prev, cond)
    lat.list <- estimate_mixing_numeric_2(Y,A.matrix,A.tensor, mu = mu)
    mixture <- lat.list$latent
    approx.lik <- population_bivariate_likelihood(Y,A.matrix,A.tensor,mixture, mu = mu, weights = NULL)
    ratio <- (exact.lik - approx.lik)/(abs(exact.lik))
    
    R.step <- round(exp((stepsize*log(R.prev) + log(R.max))/(stepsize + 1)))
    R.prev <- R.prev + R.step 
    cat(paste("Relative error ratio:", ratio, "Num Bins: ", R.prev), end = "\r")
  }
  R.prev <- R.prev - R.step
  paste("Selected number of bins:",R.prev)
  return(R.prev)
}

# occationally, the optimization can be better due to approximation error


# faster than the previous version

compute_A_matrix_2 <- function(R.bins, cond){
  grid <- seq(0,R.bins -1)/(R.bins - 1)
  A.out <- sapply(grid,function(z){
    out <- cond(z)
    out <- out/sum(out)
    return(out)
  })
  return(A.out)
} 


compute_A_tensor_2 <- function(R.bins, cond){
  grid <- seq(0,R.bins -1)/(R.bins - 1)
  #grid <- array(grid,c(length(grid),1,1))
  N <- length(cond(0.5)) - 1
  A.out <- vapply(grid,function(z){
    out <- cond(z)
    out.mat <- outer(out,out,"*")
    out.mat <- out.mat/sum(out.mat)
    return(out.mat)
  }, FUN.VALUE = array(0,c(N + 1,N + 1)))
  return(A.out)
}
# R.bins = 300
# A.matrix <- compute_A_matrix_2(R.bins, cond)
# A.tensor <- compute_A_tensor_2(R.bins, cond)

estimate_mixing_numeric_2 <- function(Y, A.matrix, A.tensor, mu, weights){
  
  
  if(ncol(Y) != 2){
    stop("Y must have 2 columns" )
  }
  if(any(is.na(Y[,1]))){
    stop("Y cannot have missingness in first column")
  }
  if(missing(weights)){
    weights <- rep(1,nrow(Y))
  }
  if(missing(mu)){
    mu = 0
  }
  n.obs <- nrow(Y)
  N <- nrow(A.matrix) - 1 # number of scores 
  R.bins <- ncol(A.matrix)
  count.vector <- rep(0,N + 1)
  count.matrix <- matrix(0,N+1, N+1)
  
  
  A.ten.mat <- A.tensor
  # converting the dimension of the tensor to a matrix
  
  dim(A.ten.mat) <- c((N + 1)^2,R.bins)
  
  #NA's can be allowed in second column only
  for(i in seq(nrow(Y))){
    if(!any(is.na(Y[i,]))){
      count.matrix[Y[i,1] + 1,Y[i,2] + 1] <- count.matrix[Y[i,1] + 1,Y[i,2] + 1] + weights[i]
    } else {
      count.vector[Y[i,1] + 1] <- count.vector[Y[i,1] + 1] + weights[i]
    }
    
  }
  # need to vectorize the matrix counts 
  count.mat.vec <- count.matrix
  dim(count.mat.vec) <- (N + 1)^2
  
  theta <- Variable(R.bins, name = "latent discretized distribution") # values of the weight vector 

  data.obj1 <-  t(count.vector) %*% log(A.matrix %*% theta)
  data.obj2 <-  t(count.mat.vec) %*% log(A.ten.mat %*% theta)
  pen.obj <- n.obs*(mu/R.bins)*t(rep(1, R.bins)) %*% log(theta) 
  
  constraints <- list(
    #obs.dist == A.matrix %*% theta,
    sum(theta) <= 1#,
    #mu/(R.bins*(1 + mu)) <= theta
  )
  
  obj.arg <- data.obj1 + data.obj2 + pen.obj
  obj <- Maximize(obj.arg)
  prob <- Problem(obj, constraints)
  
  
  #value(theta) <- rep(1/R.bins, R.bins) # initial guess of a uniform distribution
  result <- solve(prob, solver = "MOSEK")
  #result <- solve(p, verbose = TRUE)
  
  p.m <- result$getValue(theta)
  p.ma <- A.matrix %*% p.m
  out.list <- list("latent" = p.m, "observed" = p.ma)
  return(out.list)
}


scale_kernel <- function(ker,scale){
  force(scale)
  func <- function(x){
    (1/scale)*ker(x/scale)
  }
  return(func)
}

weight_vec <- function(x,X, ker.set){
  if(length(ker.set) != length(x)){
    stop("ker.set and x must have the same length")
  }
  if(ncol(X) != length(x)){
    stop("x and ncol(X) must have the same length")
  }
  w.mat <- matrix(NA,nrow = nrow(X), ncol = ncol(X))
  for(i in seq(length(ker.set))){
    ker.tmp <- ker.set[[i]]
    w.tmp <- ker.tmp(x[i] - X[,i])
    w.mat[,i] <- w.tmp
  }
  w.vec <- exp(rowSums(log(w.mat)))
  w.vec <- nrow(X)*w.vec/sum(w.vec)
  return(w.vec)
}

unique_X <- function(X){
  X.unique <- unique(X, MARGIN = 1)
  return(X.unique)
}


estimate_mixing_numeric_regression <- function(Y,A.matrix,A.tensor, mu, X, ker.set, X.out){
  # efficient estimation here due to the regression approach
  if(missing(X.out)){
    X.out <- X
  }
  X.un <- unique_X(X.out)
  R.bins <- ncol(A.matrix)
  mixture.block <- matrix(NA, nrow =  nrow(X.out), ncol = R.bins)

  for(i in seq(nrow(X.un))){
    cat(paste("Unique mixture estimate:",i,"/",nrow(X.un)), end = "\r")
    x <- as.numeric(X.un[i,])
    weights <- weight_vec(x,X,ker.set)
    mix <- estimate_mixing_numeric_2(Y, A.matrix, A.tensor, mu, weights)
    match.idx <- which(apply(X.out, 1, function(z) return(all(z == x))))
    for(j in match.idx){
      mixture.block[j,] <- as.vector(mix$latent)
    }
    
  }
  return(mixture.block)
}



# TOD0: speed up conversion_cross_entropy ? ()
# TODO: Add an automatic splitting function to 
# create a train and test set in the mu.selection and regression function. 
mu_selection_regression <- function(Y, X, cond, mu.set, R.bins, ker.set){
  res <- rep(NA, length(mu.set)) 
  A.matrix <- compute_A_matrix_2(R.bins, cond)
  A.tensor <- compute_A_tensor_2(R.bins, cond)
  complete.rows <- !(is.na(Y[,1])|is.na(Y[,2]))
  Y.complete <- Y[complete.rows,]
  for(i in seq(length(mu.set))){
    cat(paste("Model",i,"/",length(mu.set)), end = "\r")
    mu <- mu.set[[i]]
    mixture.block <- estimate_mixing_numeric_regression(Y,A.matrix,A.tensor, mu = mu,X, ker.set)
    mixture.block.complete <- mixture.block[complete.rows,]
    res[i] <- conversion_cross_entropy(Y.complete,A.matrix,A.tensor,mixture.block.complete)
  }
  i.max <- which.max(res)
  out.list <- list("opt.mu" = mu.set[[i.max]], "ce" = res)
  return(out.list)
} 

mu_selection_regression <- function(Y.train, X.train, Y.val, X.val, cond, mu.set, R.bins, ker.set){
  res <- rep(NA, length(mu.set)) 
  A.matrix <- compute_A_matrix_2(R.bins, cond)
  A.tensor <- compute_A_tensor_2(R.bins, cond)
  #complete.rows <- !(is.na(Y[,1])|is.na(Y[,2]))
  #Y.complete <- Y[complete.rows,]

  for(i in seq(length(mu.set))){
    cat(paste("Model",i,"/",length(mu.set)), end = "\r")
    mu <- mu.set[[i]]
    mixture.block <- estimate_mixing_numeric_regression(Y.train, A.matrix, 
                                                        A.tensor, mu = mu, 
                                                        X.train, ker.set, 
                                                        X.val)
    
    
    
    res[i] <- conversion_cross_entropy(Y.val,A.matrix,A.tensor,mixture.block)
    #res[i] <- conversion_cross_entropy_2(Y.val,A.matrix,A.tensor,mixture.block)
    
  }
  i.min <- which.min(res)
  out.list <- list("opt.mu" = mu.set[[i.min]], "ce" = res)
  return(out.list)
} 




quantile_map <- function(mixture,n.q){
  q.seq <- seq(0,1,length.out = n.q)
  cdf <- cumsum(mixture)
  cdf <- cdf/max(cdf)
  qf <- rep(0,n.q)
  R.bins <- length(mixture)
  grid <- seq(0,1,length.out = R.bins)
  for(j in seq(length(q.seq))){
    idx <- min(which(cdf > q.seq[j]))
    qf[j] <- grid[idx]
  }
  qf[is.na(qf)] <- max(qf, na.rm = T) 
  return(qf)
}


compute_joint_tensor <- function(R.bins,cond.y,cond.z){
  grid <- seq(0,R.bins -1)/(R.bins - 1)
  #grid <- array(grid,c(length(grid),1,1))
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  
  A.out <- vapply(grid,function(x){
    out.y <- cond.y(x)
    out.z <- cond.z(x)
    out.mat <- outer(out.y,out.z,"*")
    out.mat <- out.mat/sum(out.mat)
    return(out.mat)
  }, FUN.VALUE = array(0,c(Ny + 1,Nz + 1)))
  return(A.out)
}



joint_scores_quantile_tensor <- function(cond.y,cond.z,mixture.y,mixture.z,n.q, R.bins){
  Ny <- length(cond.y(.5)) - 1
  Nz <- length(cond.z(.5)) - 1
  qf.y <- quantile_map(mixture.y,n.q)
  qf.z <- quantile_map(mixture.z,n.q)
  grid <- seq(0,R.bins -1)/(R.bins - 1)
  A.quantile.tens <- vapply(seq(length(grid)),function(x){
    out.y <- cond.y(qf.y[x])
    out.z <- cond.z(qf.z[x])
    out.mat <- outer(out.y,out.z,"*")
    out.mat <- out.mat/sum(out.mat)
    return(out.mat)
  }, FUN.VALUE = array(0,c(Ny + 1,Nz + 1)))
  return(A.quantile.tens)
}


conditional_imputation_model <- function(cond.y,cond.z,mixture.y,mixture.z,n.q,R.bins){
  A.joint.tensor <- joint_scores_quantile_tensor(cond.y,cond.z,mixture.y,mixture.z,n.q, R.bins)
  # corresponds to weights placed at each quantile
  q.s <- rep(1/R.bins,length.out = R.bins)
  p.yz <- tensor_prod(A.joint.tensor, q.s)
  
  #A.mat.y <- compute_A_matrix_2(R.bins,cond.y)
  
  p.y <- rowSums(p.yz)
  p.z.cond.y <- p.yz/p.y
  return(p.z.cond.y)
}

impute_scores <- function(n.impute,p.z.cond.y,y){
  p.vec <- p.z.cond.y[y + 1,]
  N.z <- ncol(p.z.cond.y) - 1
  Z.samp <- sample(seq(0,N.z),n.impute,prob = p.vec)
  return(Z.samp)
}

impute_dataset <- function(n.impute,y.vec,cond.y,cond.z,latent.block.y,latent.block.z,R.bins){
  n <- length(y.vec)
  Z.imp <- matrix(0,n,n.impute)
  for(i in 1:n){
    cat(paste("Imputed",i,"/",n),end = "\r")
    mixture.y <- as.vector(latent.block.y[i,])
    mixture.z <- as.vector(latent.block.z[i,])
    p.z.cond.y <- conditional_imputation_model(cond.y,cond.z,mixture.y,mixture.z,n.q,R.bins)
    z.row <- impute_scores(n.impute,p.z.cond.y,y.vec[i])
    #return(z.row)
    Z.imp[i,] <- z.row
  }
  return(Z.imp)
}

reshape_imputed_dataset <- function(X,Z.imp){
  v = ncol(Z.imp)
  X.long <- X[rep(1:nrow(X), times = v), ]
  Z.long <- as.vector(t(Z.imp))
  out <- list("X" = X.long, "Z" = Z.long)
  return(out)
}


reference_regression <- function(X.ref, Y.train, X.train,  A.matrix,A.tensor, mu, ker.set){
  # efficient estimation here due to the regression approach
  X.un <- unique_X(X.ref)
  R.bins <- ncol(A.matrix)
  mixture.block <- matrix(NA, nrow =  nrow(X.ref), ncol = R.bins)
  for(i in seq(nrow(X.un))){
    cat(paste("Unique mixture estimate:",i,"/",nrow(X.un)), end = "\r")
    x <- as.numeric(X.un[i,])
    weights <- weight_vec(x,X.train,ker.set)
    mix <- estimate_mixing_numeric_2(Y.train, A.matrix, A.tensor, mu, weights)
    match.idx <- which(apply(X.ref, 1, function(z) return(all(z == x))))
    mixture.block[match.idx,] <- as.vector(mix$latent)
  }
  return(mixture.block)
}


impute_dataset_2 <- function(X.ref,y.ref,n.impute,Y.train,Z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q){
  X.un <- unique_X(X.ref)
  
  Ny <- length(cond.y(.5)) - 1
  Nz <- length(cond.z(.5)) - 1
  
  A.matrix.y <- compute_A_matrix_2(R.bins,cond.y)
  A.tensor.y <- compute_A_tensor_2(R.bins,cond.y)
  
  A.matrix.z <- compute_A_matrix_2(R.bins,cond.z)
  A.tensor.z <- compute_A_tensor_2(R.bins,cond.z)
  
  ## Think about how to generalize this
  # These should always be the same
  y.tr <- Y.train[,c("y1","y2")]
  z.tr <- Z.train[,c("z1","z2")]
  
  # 
  ref.cols <- colnames(X.ref)
  X.train.Y <- Y.train[,ref.cols]
  X.train.Z <- Z.train[,ref.cols]
  
  cond.block <- array(NA, c(nrow(X.ref),Ny + 1,Nz + 1))
  for(i in seq(nrow(X.un))){
    cat(paste("Unique mixture estimate:",i,"/",nrow(X.un)), end = "\r")
    x <- as.numeric(X.un[i,])
    weights.y <- weight_vec(x,X.train.Y,ker.set)
    weights.z <- weight_vec(x,X.train.Z,ker.set)
    
    mix.y <- estimate_mixing_numeric_2(y.tr, A.matrix, A.tensor, mu.y, weights.y)
    mix.z <- estimate_mixing_numeric_2(z.tr, A.matrix, A.tensor, mu.z, weights.z)
    mixture.y <- mix.y$latent
    mixture.z <- mix.z$latent
    p.z.cond.y.slice <- conditional_imputation_model(cond.y,cond.z,mixture.y,mixture.z,n.q,R.bins)
    if(sum(abs(rowSums(p.z.cond.y.slice) - 1)) > 0.01){
      break
    }
    match.idx <- which(apply(X.ref, 1, function(z) return(all(z == x))))
    # this piece is incorrect
    # Fixed this error in the assignment of the 
    
    # for loop and only assigning 
    # is pretty simple and fast
    for(id in match.idx){
      cond.block[id,,] <- p.z.cond.y.slice
    }
    
  }
  
  cond.mat <- matrix(NA,nrow(X.ref),Nz + 1)
  #relatively fast
  for(j in seq(length(y.ref))){
    # the offset of the block 
    # when y.ref = 0 strange results occur 
    if(!is.na(y.ref[j])){
      cond.mat[j,] <- cond.block[j,y.ref[j] + 1,]
    } else {
      cond.mat[j,] <- rep(1/(Nz + 1), Nz + 1)
    }
    
  }
  
  # much faster than before
  Z.imp <- matrix(NA,length(y.ref),n.impute)
  for(k in 1:length(y.ref)){
    if(!is.na(y.ref[k])){
      p.vec <- as.numeric(cond.mat[k,])
      Z.imp.row <- sample(seq(0,Nz), n.impute, prob = p.vec, replace = F)
      Z.imp[k,] <- Z.imp.row
    } else {
      Z.imp[k,] <- rep(NA,n.impute)
    }
    
  }
  return(Z.imp)
}

impute_regression <- function(formula,Z.full, X.Y.id,Z.impute){
  n.impute <- ncol(Z.impute)
  n.vars <- length(formula) + 1
  bin.naive.fit <- gee(formula,
                       id = id, data = Z.full,corstr = "exchangeable")
  n.vars <- length(bin.naive.fit$coefficients)
  
  Vg <- matrix(0,n.vars,n.vars)
  impute.betas <- matrix(0,n.vars,n.impute)
  for(j in 1:n.impute){
    impute.block <- cbind(as.numeric(Z.impute[,j]),X.Y.id)
    
    colnames(impute.block)[1] = "z"
    
    impute.block$bin_z <- 1*(impute.block$z >= 26)
    impute.data <- rbind(impute.block, Z.full)
    fit <- gee(formula,
               id = id, data = impute.data,corstr = "exchangeable")
    Vg <- Vg + fit$robust.variance
    impute.betas[,j] <- fit$coefficients
  }
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(fit$coefficients)
  
  Vb <- matrix(0,n.vars,n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  #colnames(mean.centered.betas) <- names(fit$coefficients)
  
  rownames(mean.centered.betas) <- names(fit$coefficients)
  colnames(Vb) <- names(fit$coefficients)
  rownames(Vb) <- names(fit$coefficients)
  colnames(Vg) <- names(fit$coefficients)
  rownames(Vg) <- names(fit$coefficients)
  for(j in 1:n.impute){
    Vb <- Vb + outer(mean.centered.betas[,j], mean.centered.betas[,j], "*")
  }
  Vb <- (1/(n.impute - 1))*Vb
  
  rubin.var <- Vg + (1 + 1/n.impute)*Vb
  
  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  cc.z.scores <- bin.naive.fit$coefficients/sqrt(diag(bin.naive.fit$robust.variance))
  out.list <- list("coefficients" = beta.impute.mean, 
                   "variance" = rubin.var, 
                   "z-scores" = z.scores,
                   "cc-coefficients" = bin.naive.fit$coefficients, 
                   "cc-variance" = bin.naive.fit$robust.variance, 
                   "cc-z-scores" = cc.z.scores)
  
  return(out.list)
}



impute_logistic_regression <- function(formula,Z.full, X.Y.id,Z.impute){
  n.impute <- ncol(Z.impute)
  n.vars <- length(formula) + 1
  bin.naive.fit <- gee(formula,
                       id = id, data = Z.full,corstr = "exchangeable", family = binomial(link = "logit"))
  n.vars <- length(bin.naive.fit$coefficients)
  
  Vg <- matrix(0,n.vars,n.vars)
  impute.betas <- matrix(0,n.vars,n.impute)
  for(j in 1:n.impute){
    impute.block <- cbind(as.numeric(Z.impute[,j]),X.Y.id)
    
    colnames(impute.block)[1] = "z"
    
    impute.block$bin_z <- 1*(impute.block$z >= 26)
    impute.data <- rbind(impute.block, Z.full)
    fit <- gee(formula,
               id = id, data = impute.data,corstr = "exchangeable", family = binomial(link = "logit"))
    Vg <- Vg + fit$robust.variance
    impute.betas[,j] <- fit$coefficients
  }
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(fit$coefficients)
  
  Vb <- matrix(0,n.vars,n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  #colnames(mean.centered.betas) <- names(fit$coefficients)
  
  rownames(mean.centered.betas) <- names(fit$coefficients)
  colnames(Vb) <- names(fit$coefficients)
  rownames(Vb) <- names(fit$coefficients)
  colnames(Vg) <- names(fit$coefficients)
  rownames(Vg) <- names(fit$coefficients)
  for(j in 1:n.impute){
    Vb <- Vb + outer(mean.centered.betas[,j], mean.centered.betas[,j], "*")
  }
  Vb <- (1/(n.impute - 1))*Vb
  
  rubin.var <- Vg + (1 + 1/n.impute)*Vb
  
  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  cc.z.scores <- bin.naive.fit$coefficients/sqrt(diag(bin.naive.fit$robust.variance))
  out.list <- list("coefficients" = beta.impute.mean, 
                   "variance" = rubin.var, 
                   "z-scores" = z.scores,
                   "cc-coefficients" = bin.naive.fit$coefficients, 
                   "cc-variance" = bin.naive.fit$robust.variance, 
                   "cc-z-scores" = cc.z.scores)
  
  return(out.list)
}

# requires specific names on the data frame
# 
lag_pairs_impute_regression <- function(formula,Z.full, X.Y.id,Z.impute){
  Z.1 <- Z.full %>% filter(visit == 1) 
  Z.2 <- Z.full %>% filter(visit == 2) 
  
  z.diff <- Z.1$z - Z.2$z
   
  ######### TODO: finish the editing of this function 
  
  
  
  n.impute <- ncol(Z.impute)
  n.vars <- length(formula) + 1
  bin.naive.fit <- gee(formula,
                       id = id, data = Z.full,corstr = "exchangeable")
  n.vars <- length(bin.naive.fit$coefficients)
  
  Vg <- matrix(0,n.vars,n.vars)
  impute.betas <- matrix(0,n.vars,n.impute)
  for(j in 1:n.impute){
    impute.block <- cbind(as.numeric(Z.impute[,j]),X.Y.id)
    
    colnames(impute.block)[1] = "z"
    
    impute.block$bin_z <- 1*(impute.block$z >= 26)
    impute.data <- rbind(impute.block, Z.full)
    fit <- gee(formula,
               id = id, data = impute.data,corstr = "exchangeable")
    Vg <- Vg + fit$robust.variance
    impute.betas[,j] <- fit$coefficients
  }
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(fit$coefficients)
  
  Vb <- matrix(0,n.vars,n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  #colnames(mean.centered.betas) <- names(fit$coefficients)
  
  rownames(mean.centered.betas) <- names(fit$coefficients)
  colnames(Vb) <- names(fit$coefficients)
  rownames(Vb) <- names(fit$coefficients)
  colnames(Vg) <- names(fit$coefficients)
  rownames(Vg) <- names(fit$coefficients)
  for(j in 1:n.impute){
    Vb <- Vb + outer(mean.centered.betas[,j], mean.centered.betas[,j], "*")
  }
  Vb <- (1/(n.impute - 1))*Vb
  
  rubin.var <- Vg + (1 + 1/n.impute)*Vb
  
  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  cc.z.scores <- bin.naive.fit$coefficients/sqrt(diag(bin.naive.fit$robust.variance))
  out.list <- list("coefficients" = beta.impute.mean, 
                   "variance" = rubin.var, 
                   "z-scores" = z.scores,
                   "cc-coefficients" = bin.naive.fit$coefficients, 
                   "cc-variance" = bin.naive.fit$robust.variance, 
                   "cc-z-scores" = cc.z.scores)
  
  return(out.list)
}




imputed_5ydiff <- function(formula,Z.impute,X.ref){
  n.impute <- ncol(Z.impute)
  Z.full = cbind(as.numeric(Z.impute[,1]),X.ref)
  Z.full <- Z.full[Z.full$complete == 1,]
  colnames(Z.full)[1] = "z"
  if(nrow(Z.full)== 0){
    naive.fit <- list()
    naive.fit$coefficients = NULL
    naive.fit$robust.variance = NULL
    impute.data <- cbind(as.numeric(Z.impute[,1]),X.ref)
    
    colnames(impute.data)[1] = "z"
    
    false.fit <- gee(formula, id = id, data = impute.data)
    n.vars <- length(false.fit$coefficients)
  } else {
    naive.fit <- gee(formula, id = id, data = Z.full)
    
    n.vars <- length(naive.fit$coefficients)
  }
  Vg <- matrix(0,n.vars,n.vars)
  impute.betas <- matrix(0,n.vars,n.impute)
  for(j in 1:n.impute){
    impute.data <- cbind(as.numeric(Z.impute[,j]),X.ref)
    
    colnames(impute.data)[1] = "z"
    
    fit <- gee(formula, id = id, data = impute.data)
    Vg <- Vg + fit$robust.variance
    impute.betas[,j] <- fit$coefficients
  }
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(fit$coefficients)
  
  Vb <- matrix(0,n.vars,n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  #colnames(mean.centered.betas) <- names(fit$coefficients)
  
  rownames(mean.centered.betas) <- names(fit$coefficients)
  colnames(Vb) <- names(fit$coefficients)
  rownames(Vb) <- names(fit$coefficients)
  colnames(Vg) <- names(fit$coefficients)
  rownames(Vg) <- names(fit$coefficients)
  for(j in 1:n.impute){
    Vb <- Vb + outer(mean.centered.betas[,j], mean.centered.betas[,j], "*")
  }
  Vb <- (1/(n.impute - 1))*Vb
  
  rubin.var <- Vg + (1 + 1/n.impute)*Vb
  
  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  cc.z.scores <- naive.fit$coefficients/sqrt(diag(naive.fit$robust.variance))
  out.list <- list("coefficients" = beta.impute.mean, 
                   "variance" = rubin.var, 
                   "z-scores" = z.scores,
                   "cc-coefficients" = naive.fit$coefficients, 
                   "cc-variance" = naive.fit$robust.variance, 
                   "cc-z-scores" = cc.z.scores)
  
  return(out.list)
}


impute_dataset_3 <- function(joined.tests,y.ref,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q){
  X.ref1 <- joined.tests[,c("age", "group")]
  X.ref2 <- joined.tests[,c("age_2", "group_2")]
  colnames(X.ref2) <- c("age", "group")
  y.ref1 <- joined.tests$y
  y.ref2 <- joined.tests$y_2
  
  Z1.impute <- impute_dataset_2(X.ref1,y.ref1,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)
  Z2.impute <- impute_dataset_2(X.ref2,y.ref2,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)
  idx1 <- which(!is.na(joined.tests$z))
  idx2 <- which(!is.na(joined.tests$z_2))
  Z1.impute[idx1,] = joined.tests$z[idx1]
  Z2.impute[idx2,] = joined.tests$z_2[idx2]
  z.diff <- Z2.impute - Z1.impute
  return(z.diff)
}

fit_to_table <- function(fit){
  n = length(fit$coefficients)
  out.table <- matrix(NA,nrow = n, ncol = 6)
  rownames(out.table) = names(fit$coefficients)
  colnames(out.table) = c("C.C. coef", 
                          "C.C. sd",
                          "C.C. Zscore",
                          "Imputation coef", 
                          "Imputation sd",
                          "Imputation Zscore")
  
  out.table[,1] <-  fit$`cc-coefficients`
  out.table[,2] <-  sqrt(diag(fit$`cc-variance`))
  out.table[,3] <-  fit$`cc-z-scores`
  out.table[,4] <-  fit$`coefficients`
  out.table[,5] <-  sqrt(diag(fit$`variance`))
  out.table[,6] <-  fit$`z-scores`
  return(out.table)
}

# I think that the generalized imputation idea is a bit too rough. 
# 



# 
time.lag <- as.Date("2006/01/01") - as.Date("2000/01/01")
# assume X is the cleaned test, this requires the knowledge of outcomes
# default X is the time window 
lag_pair_visit_label <- function(X,time.lag, time.window = as.Date("2001/01/01") - as.Date("2000/01/01")){
  
  X.out <- X %>% 
    group_by(id) %>% 
    arrange(date, by_group = T) %>% 
    mutate(date_diff = date - min(date))  
  
  # only take first day or time lagged day
  X.out <- X.out %>% 
    filter((date_diff >= time.lag & date_diff <= time.lag + time.window) | date_diff == 0)
  
  X.out <- X.out %>% 
   group_by(id) %>% 
   arrange(date, by_group = T) %>% 
   mutate(visit = row_number()) 
  
  X.out <- X.out %>% 
    filter(visit <= 2)
  
  X.out <- X.out %>% 
    group_by(id) %>% 
    mutate(num_visits = max(visit)) 
  X.out <- X.out %>% filter(num_visits == 2)
  return(X.out)
}


