
# New application regarding the 
# imputation of an outcome. 





# suppose that we want to study the impact on education 
# and the aging decline across normal individuals
# Our initial data 
library(dplyr)
library(gee)
#library(geepack)


n.q <- 1000 #number of quantiles

cond.y <- conditional_mkm(Ny, gaussian_kernel, 2)
cond.z <- conditional_mkm(Nz, gaussian_kernel, 2)






# imputation application

y.train <- read.csv("Data/NACCMMSE_training.csv")
z.train <- read.csv("Data/MOCATOTS_training.csv")

y.val <- read.csv("Data/NACCMMSE_validation.csv")
z.val <- read.csv("Data/MOCATOTS_validation.csv")

y.val <- y.val[,c("y1", "y2", "age", "group")]
z.val <- z.val[,c("z1", "z2", "age", "group")]

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
h.set <- c(0.5,0.8,1,1.5,2,2.5,3,4,5,7.5,10)
cond.bin.y <- generate_cond_binomial(Ny)
cond.bin.z <- generate_cond_binomial(Nz)

cond.gaus.ker.set.y <- generate_mkm_list(N = Ny, ker = gaussian_kernel, h.set = h.set)
cond.gaus.ker.set.z <- generate_mkm_list(N = Nz, ker = gaussian_kernel, h.set = h.set)


cond.list.y <- append(cond.gaus.ker.set.y,cond.bin.y)
cond.list.z <- append(cond.gaus.ker.set.z,cond.bin.z)

gaussian.cond.names <- paste0("Gaussian Kernel, h = ", round(h.set,2))
cond.names.y <- c(gaussian.cond.names, "Binomial")
cond.names.z <- c(gaussian.cond.names, "Binomial")

select.y <- error_model_selection(cond.list.y, Y, R.bins, cond.names.y)
select.z <- error_model_selection(cond.list.z, Z, R.bins, cond.names.z)

# optimal models in both cases are the binomial 
cond.y <- select.y$opt_model
cond.z <- select.z$opt_model

cond.y <- cond.bin.y
cond.z <- cond.bin.z

# group must match, age must be within 3 years 
# TODO: Fix the cross validation of the mu.set

mu.set <- c(0,0.005,0.01,0.05,0.1,0.5,1,5,10,50)

mu.set <- c(0,10**seq(-3,3,length.out = 31))

unif.3 <- scale_kernel(uniform_kernel,3)
unif.0.5 <- scale_kernel(uniform_kernel,0.5)
ker.set <- list(unif.3,unif.0.5)

R.bins.small = 1000
#TODO: make this faster if possible 

select.mu.y <- mu_selection_regression(Y,X.Y,cond.y,mu.set,R.bins.small, ker.set)
select.mu.z <- mu_selection_regression(Z,X.Z,cond.z,mu.set,R.bins.small, ker.set)




select.mu.y <- mu_selection(mu.set,cond.y,
                            y.train[,c(1,2)],y.val[,c(1,2)],
                            R.bins.small)
select.mu.z <- mu_selection(mu.set,cond.z,
                            z.train[,c(1,2)],z.val[,c(1,2)],
                            R.bins.small)

lines(log(mu.set + mu.set[2]), select.mu.y$ce)
lines(log(mu.set + mu.set[2]), select.mu.z$ce)

mu.y <- select.mu.y$opt.mu
mu.z <- select.mu.z$opt.mu

# if already selected  
mu.y <- 0.03981072
mu.z <- 0.06309573

A.matrix.y <- compute_A_matrix_2(R.bins,cond.y)
A.tensor.y <- compute_A_tensor_2(R.bins,cond.y)

A.matrix.z <- compute_A_matrix_2(R.bins,cond.z)
A.tensor.z <- compute_A_tensor_2(R.bins,cond.z)










################ new data cleaning ################


all_tests <- read.csv("Local_Data/investigator_nacc47.csv")

# indicates cognitively normal individuals 
#w_normal = all_tests$CDRGLOB==0

# TODO: can also include the different cdr scores and cdr/age interaction in the model 
# discretize it into 3 groups, >= 1, 0.5, 0 

# Maxima of possible tests scores
Ny <- 30
Nz <- 30

# Between the age of 59 and 85
w_age = (all_tests$NACCAGE>59) & (all_tests$NACCAGE <= 85)

all_tests = all_tests[w_age,]

n = nrow(all_tests)
col.subset <- c("NACCMMSE", "MOCATOTS", "NACCAGE", "SEX", "EDUC", "NACCID", "CDRGLOB", "NACCNE4S" ,"date")

# create date column 

textdate <- paste0(all_tests$VISITYR,"/",all_tests$VISITMO,"/",all_tests$VISITDAY)

all_tests$date <- as.Date(textdate)

clean.tests <- all_tests[,col.subset]
colnames(clean.tests) <- c("y","z", "age", "sex", "educ", "id", "cdr", "ne4s", "date")
clean.tests$y <- ifelse(clean.tests$y >= 0 & clean.tests$y <= Ny,  clean.tests$y, NA)
clean.tests$z <- ifelse(clean.tests$z >= 0 & clean.tests$z <= Nz,  clean.tests$z, NA)
clean.tests$id <- as.factor(clean.tests$id)
clean.tests$sex <- clean.tests$sex - 1

missing.both <- is.na(clean.tests$y) & is.na(clean.tests$z) & is.na(clean.tests$age)

clean.tests <- clean.tests[!missing.both,]
clean.tests <- clean.tests[clean.tests$educ < 50,] # error variable here, make sure this does not include the 99 error code 
clean.tests <- categorize(clean.tests)
clean.tests$educ_binary <- 1*(clean.tests$educ >= 16)
clean.tests$cdr_group <- case_when(
  clean.tests$cdr == 0 ~ "0",
  clean.tests$cdr == 0.5 ~ "0.5",
  clean.tests$cdr == 1 ~ ">=1"
)
clean.tests$cdr_group <- as.factor(clean.tests$cdr_group)
clean.tests$ne4s[clean.tests$ne4s == 9] = NA
clean.tests$ne4s_group <- as.factor(clean.tests$ne4s)








##

































































