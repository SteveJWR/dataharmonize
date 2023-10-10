

# This will consist of the minimal set of functions needed for the package
# DNOISe (Thoughts on the name, kind of a pun as we often denoise latent variables)
# DataHarmonization via 
# Nonparametric Outcome ImputationS (e)

# TODO: Check CRAN to see if this exists 

library(gee)
library(sandwich)


#devtools::install_github("richardkwo/multChernoff")
library(multChernoff)
source("R/01_functions.R")

# Use the training and test set to 
# establish the imputation latent model
# sim.data,sim.data$Y,sim.data$Z,n.impute,
# Y,Z,cond.y,cond.z,
# mu.y,mu.z,R.bins = 1000,nq = 1000

X.ref = sim.data
y.ref = sim.data$Y
z.ref = sim.data$Z
n.impute = 50


ImputeOutcomes <- function(X.ref,y.ref,z.ref,n.impute,
                           Y.train,Z.train,cond.y,cond.z,
                           mu.y,mu.z,ref.cols, ker.set,R.bins = 1000,nq = 1000){
  # when missing reference columns, we no longer adjust for covariates. 
  if(missing(ref.cols) | missing(ker.set)){

    
    if(missing(n.impute)){
      stop("must indicate number of imputations")
    }
    if(missing(cond.y) | missing(cond.y)){
      stop("must specify conditional distributions")
    }
    if(missing(mu.y) | missing(mu.z)){
      stop("must specify regularization parameters")
    }
    if(missing(Y.train) | missing(Z.train)){
      stop("must specify training data")
    }
    
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
    

    
    cond.block <- array(NA, c(nrow(X.ref),Ny + 1,Nz + 1))
    
    mix.y <- estimate_mixing_numeric_2(y.tr, A.matrix.y, A.tensor.y, mu.y)
    mix.z <- estimate_mixing_numeric_2(z.tr, A.matrix.z, A.tensor.z, mu.z)
    mixture.y <- mix.y$latent
    mixture.z <- mix.z$latent
    #TODO: verify if the conditional imputation model is correct
    p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,cond.z,mixture.y,mixture.z,nq = 10000)
    for(id in 1:nrow(X.ref)){
      cond.block[id,,] <- p.z.cond.y.slice
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
        Z.imp.row <- sample(seq(0,Nz), n.impute, prob = p.vec, replace = T)
        Z.imp[k,] <- Z.imp.row
      } else {
        Z.imp[k,] <- rep(NA,n.impute)
      }
      
    }
    
    # placing the true outcomes in the NA 
    # blocks for ease of imputation. 
    na.idx <- which(is.na(Z.imp[,1]))
    for(i in na.idx){
      Z.imp[i,] <- z.ref[i]
    }
    
  } else {
    if(length(ref.cols) == 1){
      X.un <- data.frame("tmp" = unique((X.ref[,ref.cols])))
    } else {
      X.un <-unique_X(X.ref[,ref.cols])
    }
    colnames(X.un) = ref.cols

    if(length(ref.cols) != length(ker.set)){
      stop("ref.cols and ker.set must be the same length")
    }
    if(missing(n.impute)){
      stop("must indicate number of imputations")
    }
    if(missing(cond.y) | missing(cond.y)){
      stop("must specify conditional distributions")
    }
    if(missing(mu.y) | missing(mu.z)){
      stop("must specify regularization parameters")
    }
    if(missing(Y.train) | missing(Z.train)){
      stop("must specify training data")
    }

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
    #ref.cols <- colnames(X.ref)
    X.train.Y <- as.data.frame(Y.train[,ref.cols])
    X.train.Z <- as.data.frame(Z.train[,ref.cols])

    colnames(X.train.Y) = ref.cols
    colnames(X.train.Z) = ref.cols
    cond.block <- array(NA, c(nrow(X.ref),Ny + 1,Nz + 1))
    for(i in seq(nrow(X.un))){
      cat(paste("Unique mixture estimate:",i,"/",nrow(X.un)), end = "\r")

      x <- as.numeric(X.un[i,])

      weights.y <- weight_vec(x,X.train.Y,ker.set)
      weights.z <- weight_vec(x,X.train.Z,ker.set)

      mix.y <- estimate_mixing_numeric_2(y.tr, A.matrix.y, A.tensor.y, mu.y, weights.y)
      mix.z <- estimate_mixing_numeric_2(z.tr, A.matrix.z, A.tensor.z, mu.z, weights.z)
      mixture.y <- mix.y$latent
      mixture.z <- mix.z$latent
      #TODO: verify if the conditional imputation model is correct
      p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,cond.z,mixture.y,mixture.z,nq = 10000)
      #TODO: turn this into an error
      if(sum(abs(rowSums(p.z.cond.y.slice) - 1)) > 0.01){
        break
      }

      if(length(ref.cols) > 1){
        match.idx <- which(apply(X.ref[,ref.cols], 1, function(z) return(all(z == x))))
      } else {
        match.idx <- which(X.ref[,ref.cols] == x)
      }

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
        Z.imp.row <- sample(seq(0,Nz), n.impute, prob = p.vec, replace = T)
        Z.imp[k,] <- Z.imp.row
      } else {
        Z.imp[k,] <- rep(NA,n.impute)
      }

    }

    # placing the true outcomes in the NA
    # blocks for ease of imputation.
    na.idx <- which(is.na(Z.imp[,1]))
    for(i in na.idx){
      Z.imp[i,] <- z.ref[i]
    }
  }
  
  
  return(Z.imp)
}


ImputeOutcomesKnownLatent <- function(X.ref,n.impute,
                                      cond.y,cond.z,
                                      true.mix.y, true.mix.z, 
                                      R.bins = 500){
  # when missing reference columns, we no longer adjust for covariates. 
  

  
  Ny <- length(cond.y(.5)) - 1
  Nz <- length(cond.z(.5)) - 1
  
  A.matrix.y <- compute_A_matrix_2(R.bins,cond.y)
  A.tensor.y <- compute_A_tensor_2(R.bins,cond.y)
  
  A.matrix.z <- compute_A_matrix_2(R.bins,cond.z)
  A.tensor.z <- compute_A_tensor_2(R.bins,cond.z)
  
  ## Think about how to generalize this
  # These should always be the same
  y.ref <- sim.data$Y
  z.ref <- sim.data$Z
  
  # 
  #ref.cols <- colnames(X.ref)
  
  cond.block <- array(NA, c(nrow(X.ref),Ny + 1,Nz + 1))

  mixture.y <- true.mix.y
  mixture.z <- true.mix.z
  #TODO: verify if the conditional imputation model is correct
  p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,cond.z,mixture.y,mixture.z,nq = 10000)
  for(id in 1:nrow(X.ref)){
    cond.block[id,,] <- p.z.cond.y.slice
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
      Z.imp.row <- sample(seq(0,Nz), n.impute, prob = p.vec, replace = T)
      Z.imp[k,] <- Z.imp.row
    } else {
      Z.imp[k,] <- rep(NA,n.impute)
    }
    
  }
  
  # placing the true outcomes in the NA 
  # blocks for ease of imputation. 
  na.idx <- which(is.na(Z.imp[,1]))
  for(i in na.idx){
    Z.imp[i,] <- z.ref[i]
  }
    
  
  
  
  return(Z.imp)
}

conditional_imputation_model_2 <- function(cond.y,cond.z,mixture.y,mixture.z,nq = 10000){
  
  R.bins <- length(mixture.y)
  A.matrix.y <- compute_A_matrix_2(R.bins,cond.y)
  A.matrix.z <- compute_A_matrix_2(R.bins,cond.z)
  Ny <- length(cond.y(.5)) - 1
  Nz <- length(cond.z(.5)) - 1
  
  nearest.qy <- cumsum(mixture.y)
  nearest.qy <- nearest.qy/max(nearest.qy)
  
  nearest.qz <- cumsum(mixture.z)
  nearest.qz <- nearest.qz/max(nearest.qz)
  
  U <- runif(nq)
  gamma.seq <- seq(nq)
  zeta.seq <- seq(nq)
  p.yz <- matrix(0, nrow = Ny + 1, ncol = Nz + 1)
  for(i in seq(nq)){
    u = U[i]
    gamma.idx = which.min(abs(nearest.qy - u))
    zeta.idx = which.min(abs(nearest.qz - u))
    p.y <- A.matrix.y[,gamma.idx]
    p.z <- A.matrix.z[,zeta.idx]
    p.yz <- p.yz + outer(p.y, p.z,"*")
  }
  p.yz <- p.yz/sum(p.yz)
  
  p.marg.y <- rowSums(p.yz)
  p.z.cond.y <- p.yz/p.marg.y
  p.z.cond.y <-  matrix(0, nrow = Ny + 1, ncol = Nz + 1)
  for(k in seq(Ny + 1)){
    p.z.cond.y[k,] <- p.yz[k,]/p.marg.y[k]
  }
  
  return(p.z.cond.y)
}



# requires an indicator of sequential observations. 
# requires #id and #visit variables 
# only for a specific application 
# X.ref <- clean.tests.lag.3y
# y.ref <- y.ref.3y
# z.ref <- z.ref.3y

ImputeOutcomeDifferences <- function(X.ref,y.ref,z.ref,n.impute,
                                     Y.train,Z.train,cond.y,cond.z,
                                     mu.y,mu.z,ref.cols,ker.set,R.bins = 1000,nq = 1000){
  
  Z.imp <- ImputeOutcomes(X.ref,y.ref,z.ref,n.impute,
                          Y.train,Z.train,cond.y,cond.z,
                          mu.y,mu.z,ref.cols,ker.set,R.bins,nq)
  
  id.set <- unique(X.ref$id)
  X.out <- data.frame(matrix(NA, nrow = 0, ncol = ncol(X.ref)))
  Z.diff <- matrix(NA, nrow = 0, ncol = ncol(Z.imp))
  colnames(X.out) <- colnames(X.ref)
  k = 1
  cat(end = "\n")
  complete.vec <- c()
  for(id in id.set){
    
    if(k %% 100 == 0 ){
      cat(paste0("Computing Score Differences: ",k,"/",length(id.set)), end = "\r")
    }
    
    k = k + 1
    idx <- which(X.ref$id == id & X.ref$visit == 1)
    idx2 <- which(X.ref$id == id & X.ref$visit == 2)
    if(length(idx) != 1){
      break 
    } 
    if(length(idx2) != 1){
      break 
    }
    # fixed difference vector 
    Z.diff <- rbind(Z.diff,  Z.imp[idx2,] - Z.imp[idx,])
    X.out <- rbind(X.out, X.ref[idx,])
    z.obs1 <- !is.na(X.ref$z[idx])
    z.obs2 <- !is.na(X.ref$z[idx2])
    if(z.obs1 & z.obs1){
      complete.vec <- c(complete.vec, 1)
    }else {
      complete.vec <- c(complete.vec, 0)
    }
    
  }
  X.out$complete <- complete.vec
  return(list("X" = X.out, "Z" = Z.diff))
}



ZScoreMatchDifferences <- function(X.ref,y.train,z.train, Ny = 30,Nz = 30) { 
  mu.y <- mean(y.train$y1, na.rm = T) 
  sd.y <- sd(y.train$y1, na.rm = T) 
  
  mu.z <- mean(z.train$z1, na.rm = T) 
  sd.z <- sd(z.train$z1, na.rm = T) 
  
  
  #quantile.vec <- rep(NA, Ny + 1)
  conversion.vec <- rep(NA, Ny + 1)
  
  for(y in seq(0,Ny)){
    q <- round((sd.z/sd.y)*(y - mu.y) + mu.z)
    conversion.vec[y + 1] = max(min(q,Nz),0)
  }
  
  X.ref.new <- X.ref
  X.ref.new <- X.ref.new %>% arrange(id,visit)
  idx.set.1 <- 2*(1:(nrow(X.ref.new)/2) ) - 1 
  idx.set.2 <- 2*(1:(nrow(X.ref.new)/2) )
  
  y.out <- X.ref.new$y
  z.out <- X.ref.new$z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]
  
  z1.vec <- z.out[idx.set.1]
  z2.vec <- z.out[idx.set.2]
  z.diff <- z2.vec - z1.vec
  
  X.out <- X.ref.new[idx.set.1,]
  return(list("X" = X.out, "Z" = z.diff))
}

ZScoreConversion <- function(X.ref,y.train,z.train, Ny = 30,Nz = 30) { 
  mu.y <- mean(y.train$y1, na.rm = T) 
  sd.y <- sd(y.train$y1, na.rm = T) 
  
  mu.z <- mean(z.train$z1, na.rm = T) 
  sd.z <- sd(z.train$z1, na.rm = T) 
  
  
  #quantile.vec <- rep(NA, Ny + 1)
  conversion.vec <- rep(NA, Ny + 1)
  
  for(y in seq(0,Ny)){
    q <- round((sd.z/sd.y)*(y - mu.y) + mu.z)
    conversion.vec[y + 1] = max(min(q,Nz),0)
  }
  

  
  y.out <- X.ref$Y
  z.out <- X.ref$Z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]
  
  X.out <- X.ref
  return(list("X" = X.out, "Z" = z.out))
}

QuantileMatchDifferences <- function(X.ref,y.train,z.train, Ny = 30,Nz = 30) { 
  
  #quantile.vec <- rep(NA, Ny + 1)
  conversion.vec <- rep(NA, Ny + 1)
  
  for(y in seq(0,Ny)){
    p <- mean(y.train$y1 <= y)
    q <- round(quantile(z.train$z1, p))
    conversion.vec[y + 1] = q
  }
  
  X.ref.new <- X.ref
  X.ref.new <- X.ref.new %>% arrange(id,visit)
  idx.set.1 <- 2*(1:(nrow(X.ref.new)/2) ) - 1 
  idx.set.2 <- 2*(1:(nrow(X.ref.new)/2) )
  
  y.out <- X.ref.new$y
  z.out <- X.ref.new$z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]
  
  z1.vec <- z.out[idx.set.1]
  z2.vec <- z.out[idx.set.2]
  z.diff <- z2.vec - z1.vec
  
  X.out <- X.ref.new[idx.set.1,]
  return(list("X" = X.out, "Z" = z.diff))
}


QuantileConversion <- function(X.ref,y.train,z.train, Ny = 30,Nz = 30) { 
  conversion.vec <- rep(NA, Ny + 1)
  y.out <- X.ref$Y
  z.out <- X.ref$Z
  for(y in seq(0,Ny)){
    p <- mean(y.out <= y, na.rm = T)
    q <- ceiling(quantile(z.out, p, na.rm = T))
    conversion.vec[y + 1] = q
  }
  
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]
  
  X.out <- X.ref
  return(list("X" = X.out, "Z" = z.out))
}


# should be able to pass in all the info of GEE 
# Rely on GEE package
# pass on all information to the GLM function 

ImputationRegressionGLM <- function(formula, X, Z.impute, fit.cc = T, ...){
  n.impute <- ncol(Z.impute)
  if(!("complete" %in% colnames(X))){
    stop("Must indicate which rows have complete data") 
  }
  n.cc <- nrow(X[X$complete == 1,])
  if(fit.cc){
    data.full = cbind(as.numeric(Z.impute[,1]),X)
   
    data.full <- data.full[data.full$complete == 1,]
    n.cc <- nrow(data.full)
    colnames(data.full)[1] <- "outcome"
    naive.fit <- tryCatch(                       # Applying tryCatch
      expr = {                      # Specifying expression
        glm(formula, data = data.full, ...)
        #message("Everything was fine.")
      },
      
      error = function(e){          # Specifying error message
        message("Error in Complete Cases Fit")
      },
      
      warning = function(w){        # Specifying warning message
        message("Warning in Complete Cases Fit")
      },
      
      finally = {                   # Specifying final message
        #message("tryCatch is finished.")
      }
    )
  } else {
    #warning("No complete cases, cc refers to first imputation")
    data.full = cbind(as.numeric(Z.impute[,1]),X)
    if(!("complete" %in% colnames(data.full))){
      stop("Must indicate which rows have complete data") 
    }
    n.cc <- nrow(data.full)
    colnames(data.full)[1] <- "outcome"
    naive.fit <- tryCatch(                       # Applying tryCatch
      expr = {                      # Specifying expression
        glm(formula, data = data.full, ...)
        #message("Everything was fine.")
      },
      
      error = function(e){          # Specifying error message
        message("Error in Complete Cases Fit")
      },
      
      warning = function(w){        # Specifying warning message
        message("Warning in Complete Cases Fit")
      },
      
      finally = {                   # Specifying final message
        #message("tryCatch is finished.")
      }
    )
  }
  
  
  
  n.vars <- length(naive.fit$coefficients)
  

  
  Vg <- matrix(0,n.vars,n.vars)
  impute.betas <- matrix(0,n.vars,n.impute)
  for(j in 1:n.impute){
    data.imp = cbind(as.numeric(Z.impute[,j]),X)
    colnames(data.imp)[1] <- "outcome"
    imp.fit <- glm(formula, data = data.imp, ...)
    #imp.fit <- glm(formula, data = data.imp)
    Vg <- Vg + sandwich(imp.fit)
    impute.betas[,j] <- imp.fit$coefficients
  }
  # normalize the Rubin Variance
  Vg <- Vg/n.impute
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(imp.fit$coefficients)
  
  Vb <- matrix(0,n.vars,n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  #colnames(mean.centered.betas) <- names(imp.fit$coefficients)
  
  rownames(mean.centered.betas) <- names(imp.fit$coefficients)
  colnames(Vb) <- names(imp.fit$coefficients)
  rownames(Vb) <- names(imp.fit$coefficients)
  colnames(Vg) <- names(imp.fit$coefficients)
  rownames(Vg) <- names(imp.fit$coefficients)
  for(j in 1:n.impute){
    Vb <- Vb + outer(mean.centered.betas[,j], mean.centered.betas[,j], "*")
  }
  Vb <- (1/(n.impute - 1))*Vb
  
  rubin.var <- Vg + (1 + 1/n.impute)*Vb
  impute.varfrac <- max(diag((((1/n.impute)*Vb)/rubin.var)))
  print(paste0("Residual Fraction of Variance due to imputation: ", round(impute.varfrac, 3)))
  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  p.vals.imp <- 2*pnorm(-abs(z.scores))
  
  if(fit.cc){
    cc.var <- sandwich(naive.fit)
    cc.beta <- naive.fit$coefficients
    cc.z.scores <- cc.beta/sqrt(diag(cc.var))
    cc.p.vals <- 2*pnorm(-abs(cc.z.scores))
  } else {
    cc.var <- NA
    cc.beta <- NA
    cc.z.scores <- NA
    cc.p.vals <- NA
    n.cc <- NA
  }
  
  
  out.list <- list("coefficients" = beta.impute.mean, 
                   "variance" = rubin.var, 
                   "z-scores" = z.scores,
                   "p-values" = p.vals.imp,
                   "cc-coefficients" = cc.beta, 
                   "cc-variance" = cc.var, 
                   "cc-z-scores" = cc.z.scores,
                   "cc-p-values" = cc.p.vals, 
                   "cc-n" = n.cc,
                   "impute.variance.fraction" = impute.varfrac)
  
  return(out.list)
}


fit_to_table <- function(fit){
  n = length(fit$coefficients)
  out.table <- matrix(NA,nrow = n, ncol = 8)
  rownames(out.table) = names(fit$coefficients)
  colnames(out.table) = c("C.C. coef", 
                          "C.C. sd",
                          "C.C. Zscore",
                          "C.C. p-values",
                          "Imputation coef", 
                          "Imputation sd",
                          "Imputation Zscore",
                          "Imputation p-values",)
  
  out.table[,1] <-  fit$`cc-coefficients`
  out.table[,2] <-  sqrt(diag(fit$`cc-variance`))
  out.table[,3] <-  fit$`cc-z-scores`
  out.table[,4] <-  fit$`cc-p-values`

  
  out.table[,5] <-  fit$`coefficients`
  out.table[,6] <-  sqrt(diag(fit$`variance`))
  out.table[,7] <-  fit$`z-scores`
  out.table[,8] <-  fit$`p-values`
  return(out.table)
}



ImputationRegressionGEE <- function(formula, X, Z.impute, ...){
  n.impute <- ncol(Z.impute)
  data.full = cbind(as.numeric(Z.impute[,1]),X)
  if(!("complete" %in% colnames(data.full))){
    stop("Must indicate which rows have complete data") 
  }
  data.full <- data.full[data.full$complete == 1,]
  colnames(data.full)[1] <- "z"
  naive.fit <- tryCatch(                       # Applying tryCatch
    expr = {                      # Specifying expression
      gee(formula, data = data.full, ...)
      #message("Everything was fine.")
    },
    
    error = function(e){          # Specifying error message
      message("Error in Complete Cases Fit")
    },
    
    warning = function(w){        # Specifying warning message
      message("Warning in Complete Cases Fit")
    },
    
    finally = {                   # Specifying final message
      #message("tryCatch is finished.")
    }
  )
  
  
  n.vars <- length(naive.fit$coefficients)
  
  
  
  Vg <- matrix(0,n.vars,n.vars)
  impute.betas <- matrix(0,n.vars,n.impute)
  for(j in 1:n.impute){
    data.imp = cbind(as.numeric(Z.impute[,j]),X)
    colnames(data.imp)[1] <- "z"
    imp.fit <- gee(formula, data = data.imp, ...)
    Vg <- Vg + imp.fit$robust.variance
    impute.betas[,j] <- imp.fit$coefficients
  }
  Vg <- Vg/n.impute
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(imp.fit$coefficients)
  
  Vb <- matrix(0,n.vars,n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  #colnames(mean.centered.betas) <- names(imp.fit$coefficients)
  
  rownames(mean.centered.betas) <- names(imp.fit$coefficients)
  colnames(Vb) <- names(imp.fit$coefficients)
  rownames(Vb) <- names(imp.fit$coefficients)
  colnames(Vg) <- names(imp.fit$coefficients)
  rownames(Vg) <- names(imp.fit$coefficients)
  for(j in 1:n.impute){
    Vb <- Vb + outer(mean.centered.betas[,j], mean.centered.betas[,j], "*")
  }
  Vb <- (1/(n.impute - 1))*Vb
  
  rubin.var <- Vg + (1 + 1/n.impute)*Vb
  impute.varfrac <- max(diag((((1/n.impute)*Vb)/rubin.var)))
  print(paste0("Residual Fraction of Variance due to imputation: ", round(impute.varfrac, 3)))
  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  p.vals.imp <- 2*pnorm(-abs(z.scores))
  cc.var <- naive.fit$robust.variance
  cc.beta <- naive.fit$coefficients
  cc.z.scores <- cc.beta/sqrt(diag(cc.var))
  cc.p.vals <- 2*pnorm(-abs(cc.z.scores))
  out.list <- list("coefficients" = beta.impute.mean, 
                   "variance" = rubin.var, 
                   "z-scores" = z.scores,
                   "p-values" = p.vals.imp,
                   "cc-coefficients" = cc.beta, 
                   "cc-variance" = cc.var, 
                   "cc-z-scores" = cc.z.scores,
                   "cc-p-values" = cc.p.vals, 
                   "impute.variance.fraction" = impute.varfrac)
  
  return(out.list)
}



FeasibilityTest <- function(y,cond){
  
}


# requires a pre-fit on the univariate data
BivariateFeasibilityTest <- function(yz.pairs, cond.y, cond.z, 
                                     latent.y, latent.z, method = "Guo"){
  qy <- quantiles_from_weights(latent.y)
  qz <- quantiles_from_weights(latent.z)
  
  if(length(qy) != length(qz)){
    stop("latent distributions must have the same grid")
  }
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  
  nq <- length(qy)
  
  p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
  
  for(i in seq(nq)){
    gam <- qy[i]
    zet <- qz[i]
    py <- cond.y(gam)
    pz <- cond.z(zet)
    p.yz.slice <- outer(py,pz,"*")
    p.yz <- p.yz + p.yz.slice
    
  }
  #p.yz <- p.yz/nq
  #normalization
  p.yz <- p.yz/sum(p.yz)
  
  emp.p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
  
  for(j in seq(nrow(yz.pairs))){
    emp.p.yz[yz.pairs[j,1] + 1, yz.pairs[j,2] + 1] = emp.p.yz[yz.pairs[j,1] + 1, yz.pairs[j,2] + 1] + 1
  }
  emp.p.yz <- emp.p.yz/sum(emp.p.yz)
  
  kldiv <- kl_divergence(emp.p.yz, emp.p.yz)
  dim.k <- length(emp.p.yz) - 1
  n.obs <- nrow(yz.pairs)
  
  if(method == "Guo"){
    p.cons <- tryCatch(                       # Applying tryCatch
      expr = {                      # Specifying expression
        tailProbBound(kldiv,dim.k, n.obs)
        #message("Everything was fine.")
      },
      
      error = function(e){          # Specifying error message
        message("Error in Guo Method")
      },
      
      warning = function(w){        # Specifying warning message
        message("Warning in Guo Method")
      },
      
      finally = {                   # Specifying final message
        #message("tryCatch is finished.")
      }
    )
    
  } else {
    p.cons <- mardia_tail_bound(kldiv,dim.k, n.obs)
  }
  return(min(p.cons,1))
}

















