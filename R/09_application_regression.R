
# New application regarding the 
# imputation of an outcome. 





# suppose that we want to study the impact on education 
# and the aging decline across normal individuals
# Our initial data 
library(dplyr)
library(gee)
#library(geepack)

png.width = 1200
png.height = 1000
png.res = 200 

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

# select.mu.y <- mu_selection_regression(Y,X.Y,cond.y,mu.set,R.bins.small, ker.set)
# select.mu.z <- mu_selection_regression(Z,X.Z,cond.z,mu.set,R.bins.small, ker.set)
# 



select.mu.y <- mu_selection(mu.set,cond.y,
                            y.train[,c(1,2)],y.val[,c(1,2)],
                            R.bins.small)
select.mu.z <- mu_selection(mu.set,cond.z,
                            z.train[,c(1,2)],z.val[,c(1,2)],
                            R.bins.small)

plot(log(mu.set + mu.set[2]), select.mu.y$ce, type = "l")
plot(log(mu.set + mu.set[2]), select.mu.z$ce, type = "l")

mu.y <- select.mu.y$opt.mu
mu.z <- select.mu.z$opt.mu

# if already selected  
mu.y <- 2.511886 #0.03981072  #2.511886
mu.z <- 1.584893 #0.06309573  #1.584893

A.matrix.y <- compute_A_matrix_2(R.bins,cond.y)
A.tensor.y <- compute_A_tensor_2(R.bins,cond.y)

A.matrix.z <- compute_A_matrix_2(R.bins,cond.z)
A.tensor.z <- compute_A_tensor_2(R.bins,cond.z)

model.y <- estimate_mixing_numeric_2(y.train[,c(1,2)], A.matrix.y, 
                          A.tensor.y, mu.y)
plot(model.y$latent, type = "l")

model.z <- estimate_mixing_numeric_2(z.train[,c(1,2)], A.matrix.z, 
                                     A.tensor.z, mu.z)
plot(model.z$latent, type = "l")

##### feasibility of the model through the crosswalk dataset 

cw.data <- read.csv("Data/NACCMMSE_to_MOCATOTS_test.csv")

yz.pairs <- cw.data[,c(1,2)]

model.y <- estimate_mixing_numeric_2(Y.train[,c(1,2)], 
                                      A.matrix.y, A.tensor.y ,
                                      0)
lat.mix.y <- as.numeric(model.y$latent)

plot(lat.mix.y, type = "l")

model.z <- estimate_mixing_numeric_2(Z.train[,c(1,2)], 
                                      A.matrix.z, A.tensor.z ,
                                      0)
lat.mix.z <- as.numeric(model.z$latent)
plot(lat.mix.z, type = "l")

p.feasibility <- BivariateFeasibilityTest(yz.pairs, cond.y, cond.z,
                                          lat.mix.y, lat.mix.z, method = "Mardia")


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


no.tests <- is.na(clean.tests$y) & is.na(clean.tests$z)
clean.tests <- clean.tests[!no.tests,]
#clean.tests.cdr0 <- clean.tests %>% filter(cdr_group == '0')



##
time.lag.3y <- as.Date("2003/01/01") - as.Date("2000/01/01")
clean.tests.lag.3y <- lag_pair_visit_label(clean.tests, time.lag.3y)

y.ref.3y <- clean.tests.lag.3y$y
z.ref.3y <- clean.tests.lag.3y$z

n.impute <- 50
list.3y <- ImputeOutcomeDifferences(clean.tests.lag.3y,
                                    y.ref.3y,z.ref.3y,n.impute,
                                    y.train,z.train,cond.y,cond.z,
                                    mu.y,mu.z,ref.cols,ker.set,
                                    R.bins = 1000,n.q = 1000)

#### Impute Without Covariates
list.3y.nocov <- ImputeOutcomeDifferences(clean.tests.lag.3y,
                                          y.ref.3y,z.ref.3y,n.impute,
                                          y.train,z.train,cond.y,cond.z,
                                          mu.y,mu.z,
                                          R.bins = 1000,n.q = 1000)


list.3y.zscore <- ZScoreMatchDifferences(clean.tests.lag.3y,
                                         y.train,z.train)
list.3y.quantile <- QuantileMatchDifferences(clean.tests.lag.3y,
                                             y.train,z.train)

z.score.dat.3y <- cbind(list.3y.zscore$Z, list.3y.zscore$X)
colnames(z.score.dat.3y)[1] = "outcome"

quantile.dat.3y <- cbind(list.3y.quantile$Z, list.3y.quantile$X)
colnames(quantile.dat.3y)[1] = "outcome"

fmla <- formula(outcome ~ age + sex + educ_binary + ne4s_group)
X <- list.3y$X
Z.imp <- list.3y$Z

X.nocov <- list.3y.nocov$X
Z.imp.nocov <- list.3y.nocov$Z

fmla.1 <- formula(outcome ~ age + sex +  educ_binary *ne4s_group + cdr_group)
#fmla.1 <- formula(outcome ~ age)

imp.reg.results.3y <- ImputationRegressionGLM(fmla.1, X, Z.imp)
imp.reg.results.3y.nocov <- ImputationRegressionGLM(fmla.1, X.nocov, Z.imp.nocov)


imp.reg.results.3y$coefficients
imp.reg.results.3y$`cc-coefficients`

imp.reg.results.3y$`p-values`
imp.reg.results.3y$`cc-p-values`


imp.reg.results.3y.nocov$coefficients
imp.reg.results.3y.nocov$`cc-coefficients`

imp.reg.results.3y.nocov$`p-values`
imp.reg.results.3y.nocov$`cc-p-values`



#naive z matching 
fit.z.score <- glm(fmla.1, data = z.score.dat.3y)
z.score.coefs <- fit.z.score$coefficients
z.score.coefs

#quantile matching 
fit.quantile <- glm(fmla.1, data = quantile.dat.3y)
quantile.match.coefs <- fit.quantile$coefficients
quantile.match.coefs

#####


term <- names(quantile.match.coefs)

# complete case 
m1 <- data.frame("term" = term, 
                 "estimate" = imp.reg.results.3y$`cc-coefficients`, 
                 "std.error" = sqrt(diag(imp.reg.results.3y$`cc-variance`)))
m1$model = "Complete Case"
m2 <- data.frame("term" = term, 
                 "estimate" = imp.reg.results.3y.nocov$`coefficients`, 
                 "std.error" = sqrt(diag(imp.reg.results.3y.nocov$variance)))
m2$model = "DNOISe"
m3 <- data.frame("term" = term, 
                 "estimate" = imp.reg.results.3y$`coefficients`, 
                 "std.error" = sqrt(diag(imp.reg.results.3y$variance)))
m3$model = "DNOISe (cov. smoothed)"
m4 <- data.frame("term" = term, 
                 "estimate" = z.score.coefs, 
                 "std.error" = sqrt(diag(sandwich(fit.z.score))))
m4$model = "Z Score Matching"
m5 <- data.frame("term" = term, 
                 "estimate" = quantile.match.coefs, 
                 "std.error" = sqrt(diag(sandwich(fit.quantile))))
m5$model = "Quantile Matching"

model.coefs <- rbind(m1,m2,m3,m4,m5)

### dot and whisker plot
library(dotwhisker)


p.3y <- dwplot(model.coefs, show_intercept = FALSE)

png(filename = "Plots/3_year_impute_e4_education_interaction.png",
    width = png.width, height = png.height, res = png.res)
p.3y


# Close the pdf file
dev.off() 






tab.reg.3y.1 <- fit_to_table(imp.reg.results.3y)

fmla.2 <- formula(outcome ~ ne4s_group)
imp.reg.results.3y.2 <- ImputationRegressionGLM(fmla.2, X, Z.imp)
tab.reg.3y.2 <- fit_to_table(imp.reg.results.3y.2)

fmla.3 <- formula(outcome ~ educ_binary)
imp.reg.results.3y.3 <- ImputationRegressionGLM(fmla.3, X, Z.imp)
tab.reg.3y.3 <- fit_to_table(imp.reg.results.3y.3)

fmla.4 <- formula(outcome ~ sex)
imp.reg.results.3y.4 <- ImputationRegressionGLM(fmla.4, X, Z.imp)
tab.reg.3y.4 <- fit_to_table(imp.reg.results.3y.4)

fmla.5 <- formula(outcome ~ age)
imp.reg.results.3y.5 <- ImputationRegressionGLM(fmla.5, X, Z.imp)
tab.reg.3y.5 <- fit_to_table(imp.reg.results.3y.5)












#### 6 year trends
time.lag.6y <- as.Date("2006/01/01") - as.Date("2000/01/01")
clean.tests.lag.6y <- lag_pair_visit_label(clean.tests, time.lag.6y)
y.ref.6y <- clean.tests.lag.6y$y
z.ref.6y <- clean.tests.lag.6y$z

n.impute <- 50
list.6y <- ImputeOutcomeDifferences(clean.tests.lag.6y,
                                    y.ref.6y,z.ref.6y,n.impute,
                                    y.train,z.train,cond.y,cond.z,
                                    mu.y,mu.z,ref.cols,ker.set,
                                    R.bins = 1000,n.q = 1000)

fmla <- formula(outcome ~ age + sex + educ_binary + ne4s_group)
X <- list.6y$X
Z.imp <- list.6y$Z
fmla.1 <- formula(outcome ~ age + sex +  educ_binary *ne4s_group )

#TODO: fix conformable arrays error
imp.reg.results.6y <- ImputationRegressionGLM(fmla.1, X, Z.imp, fit.cc = F)


imp.reg.results.6y$coefficients
imp.reg.results.6y$`cc-coefficients`

imp.reg.results.6y$`p-values`
imp.reg.results.6y$`cc-p-values`


list.6y.zscore <- ZScoreMatchDifferences(clean.tests.lag.6y,
                                         y.train,z.train)
list.6y.quantile <- QuantileMatchDifferences(clean.tests.lag.6y,
                                             y.train,z.train)




z.score.dat.6y <- cbind(list.6y.zscore$Z, list.6y.zscore$X)
colnames(z.score.dat.6y)[1] = "outcome"

quantile.dat.6y <- cbind(list.6y.quantile$Z, list.6y.quantile$X)
colnames(quantile.dat.6y)[1] = "outcome"



#naive z matching 
fit.z.score <- glm(fmla.1, data = z.score.dat.6y)
z.score.coefs <- fit.z.score$coefficients
z.score.coefs

#quantile matching 
fit.quantile <- glm(fmla.1, data = quantile.dat.6y)
quantile.match.coefs <- fit.quantile$coefficients
quantile.match.coefs

#####


term <- names(quantile.match.coefs)

# complete case 
m1 <- data.frame("term" = term, 
                 "estimate" = imp.reg.results.6y$`cc-coefficients`, 
                 "std.error" = sqrt(diag(imp.reg.results.6y$`cc-variance`)))
m1$model = "Complete Case"
m2 <- data.frame("term" = term, 
                 "estimate" = imp.reg.results.6y$`coefficients`, 
                 "std.error" = sqrt(diag(imp.reg.results.6y$variance)))
m2$model = "DNOISe"
m3 <- data.frame("term" = term, 
                 "estimate" = z.score.coefs, #
                 "std.error" = sqrt(diag(sandwich(fit.z.score))))
m3$model = "Z Score Matching"
m4 <- data.frame("term" = term, 
                 "estimate" = quantile.match.coefs, 
                 "std.error" = sqrt(diag(sandwich(fit.quantile))))
m4$model = "Quantile Matching"

model.coefs <- rbind(m1,m2,m3,m4)
model.coefs <- rbind(m2,m3,m4)
### dot and whisker plot
library(dotwhisker)
p.6y <- dwplot(model.coefs, show_intercept = FALSE)

png(filename = "Plots/6_year_impute_e4_education_interaction.png",
    width = png.width, height = png.height, res = png.res)
p.6y


# Close the pdf file
dev.off() 




tab.reg.6y.1 <- fit_to_table(imp.reg.results.6y)

fmla.2 <- formula(outcome ~ ne4s_group)
imp.reg.results.6y.2 <- ImputationRegressionGLM(fmla.2, X, Z.imp)
tab.reg.6y.2 <- fit_to_table(imp.reg.results.6y.2)

fmla.3 <- formula(outcome ~ educ_binary)
imp.reg.results.6y.3 <- ImputationRegressionGLM(fmla.3, X, Z.imp)
tab.reg.6y.3 <- fit_to_table(imp.reg.results.6y.3)

fmla.4 <- formula(outcome ~ sex)
imp.reg.results.6y.4 <- ImputationRegressionGLM(fmla.4, X, Z.imp)
tab.reg.6y.4 <- fit_to_table(imp.reg.results.6y.4)

fmla.5 <- formula(outcome ~ age)
imp.reg.results.6y.5 <- ImputationRegressionGLM(fmla.5, X, Z.imp)
tab.reg.6y.5 <- fit_to_table(imp.reg.results.6y.5)










#### 8 year trends
time.lag.8y <- as.Date("2008/01/01") - as.Date("2000/01/01")
time.window <- as.Date("2002/01/01") - as.Date("2000/01/01")
clean.tests.lag.8y <- lag_pair_visit_label(clean.tests, time.lag.8y, time.window = time.window)
y.ref.8y <- clean.tests.lag.8y$y
z.ref.8y <- clean.tests.lag.8y$z

n.impute <- 50
list.8y <- ImputeOutcomeDifferences(clean.tests.lag.8y,
                                    y.ref.8y,z.ref.8y,n.impute,
                                    y.train,z.train,cond.y,cond.z,
                                    mu.y,mu.z,ref.cols,ker.set,
                                    R.bins = 1000,n.q = 1000)

#### Impute Without Covariates
list.8y.nocov <- ImputeOutcomeDifferences(clean.tests.lag.8y,
                                          y.ref.8y,z.ref.8y,n.impute,
                                          y.train,z.train,cond.y,cond.z,
                                          mu.y,mu.z,
                                          R.bins = 1000,n.q = 1000)

fmla <- formula(outcome ~ age + sex + educ_binary + ne4s_group)
X <- list.8y$X
Z.imp <- list.8y$Z

X.nocov <- list.8y.nocov$X
Z.imp.nocov <- list.8y.nocov$Z


fmla.1 <- formula(outcome ~ age + sex +  educ_binary *ne4s_group + cdr_group)
#fmla.1 <- formula(outcome ~ age + sex +  educ_binary *ne4s_group + cdr_group)
#fmla.1 <- formula(outcome ~ age)


#fmla.1 <- formula(outcome ~ age + sex +  educ_binary  + ne4s_group + cdr_group)

#TODO: fix conformable arrays error
imp.reg.results.8y <- ImputationRegressionGLM(fmla.1, X, Z.imp, fit.cc = F)
imp.reg.results.8y.nocov <- ImputationRegressionGLM(fmla.1, X.nocov, Z.imp.nocov, fit.cc = F)


imp.reg.results.8y$coefficients
imp.reg.results.8y$`cc-coefficients`

imp.reg.results.8y$`p-values`
imp.reg.results.8y$`cc-p-values`


imp.reg.results.8y.nocov$coefficients
imp.reg.results.8y.nocov$`cc-coefficients`

imp.reg.results.8y.nocov$`p-values`
imp.reg.results.8y.nocov$`cc-p-values`




list.8y.zscore <- ZScoreMatchDifferences(clean.tests.lag.8y,
                                         y.train,z.train)
list.8y.quantile <- QuantileMatchDifferences(clean.tests.lag.8y,
                                             y.train,z.train)




z.score.dat.8y <- cbind(list.8y.zscore$Z, list.8y.zscore$X)
colnames(z.score.dat.8y)[1] = "outcome"

quantile.dat.8y <- cbind(list.8y.quantile$Z, list.8y.quantile$X)
colnames(quantile.dat.8y)[1] = "outcome"



#naive z matching 
fit.z.score <- glm(fmla.1, data = z.score.dat.8y)
z.score.coefs <- fit.z.score$coefficients
z.score.coefs

#quantile matching 
fit.quantile <- glm(fmla.1, data = quantile.dat.8y)
quantile.match.coefs <- fit.quantile$coefficients
quantile.match.coefs

#####


term <- names(quantile.match.coefs)

# complete case 
m1 <- data.frame("term" = term, 
                 "estimate" = imp.reg.results.8y.nocov$`coefficients`, 
                 "std.error" = sqrt(diag(imp.reg.results.8y.nocov$variance)))
m1$model = "DNOISe"
m2 <- data.frame("term" = term, 
                 "estimate" = imp.reg.results.8y$`coefficients`, 
                 "std.error" = sqrt(diag(imp.reg.results.8y$variance)))
m2$model = "DNOISe (cov.smooothed)"
m3 <- data.frame("term" = term, 
                 "estimate" = z.score.coefs, #
                 "std.error" = sqrt(diag(sandwich(fit.z.score))))
m3$model = "Z Score Matching"
m4 <- data.frame("term" = term, 
                 "estimate" = quantile.match.coefs, 
                 "std.error" = sqrt(diag(sandwich(fit.quantile))))
m4$model = "Quantile Matching"

model.coefs <- rbind(m1,m2,m3,m4)
#model.coefs <- rbind(m2,m3,m4)
### dot and whisker plot
library(dotwhisker)

p.8y <- dwplot(model.coefs, show_intercept = FALSE)
p.8y


png(filename = "Plots/8_year_impute_e4_education_interaction.png",
    width = png.width, height = png.height, res = png.res)
p.8y


# Close the pdf file
dev.off() 


imp.reg.results.8y$`p-values`

imp.reg.results.8y.nocov$`p-values`





