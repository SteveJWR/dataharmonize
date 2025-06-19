rm(list = ls())

# suppose that we want to study the impact on education
# and the aging decline across normal individuals
# Our initial data
library(dplyr)
library(gee)
library(dnoiseR)
library(sandwich)
source("R/difference_regression_functions.R")

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




### Estimated Model Form Application Prediction Script

cond.y <- generate_cond_binomial(Ny)
cond.z <- generate_cond_binomial(Nz)

mu.y1 <- 0.01
mu.z1 <- 0.01

mu.y2 <- 0.08163265 # optimized via conversion from cross entropy
mu.z2 <- 0.3061224

unif.4 <- scale_kernel(uniform_kernel,4)
unif.0.5 <- scale_kernel(uniform_kernel,0.5) # groups the categories in the other case
ker.set <- list(unif.4,unif.0.5)



################ regression data cleaning ################



all_tests <- read.csv("Local_Data/investigator_nacc47.csv")

# indicates cognitively normal individuals
#w_normal = all_tests$CDRGLOB==0

# can also include the different cdr scores and cdr/age interaction in the model
# discretize it into 3 groups, >= 1, 0.5, 0

# Maxima of possible tests scores
Ny <- 30
Nz <- 30

# Between the age of 60 and 85
w_age = (all_tests$NACCAGE>59) & (all_tests$NACCAGE <= 85)

all_tests = all_tests[w_age,]

n = nrow(all_tests)

col.subset <- c("NACCMMSE", "MOCATOTS", "NACCAGE", "SEX", "EDUC", "NACCID", "CDRGLOB", "NACCNE4S" ,"date")

# create date column

textdate <- paste0(all_tests$VISITYR,"/",all_tests$VISITMO,"/",all_tests$VISITDAY)

all_tests$date <- as.Date(textdate)

clean.tests <- all_tests[,col.subset]
# MMSE, MOCA, age, sex, education (years), id, cdr global score, number of e4 alleles, date
colnames(clean.tests) <- c("y","z", "age", "sex", "educ", "id", "cdr", "ne4s", "date")
clean.tests$y <- ifelse(clean.tests$y >= 0 & clean.tests$y <= Ny,  clean.tests$y, NA)
clean.tests$z <- ifelse(clean.tests$z >= 0 & clean.tests$z <= Nz,  clean.tests$z, NA)
clean.tests$id <- as.factor(clean.tests$id)


missing.both <- is.na(clean.tests$y) & is.na(clean.tests$z) & is.na(clean.tests$age)

clean.tests <- clean.tests[!missing.both,]
clean.tests <- clean.tests[clean.tests$educ < 50,] # error variable here, make sure this does not include the 99 error code
clean.tests <- categorize(clean.tests)



table(clean.tests$group)

clean.tests$educ_binary <- 1*(clean.tests$educ >= 16)
clean.tests$cdr_group <- case_when(
  clean.tests$cdr == 0 ~ "0",
  clean.tests$cdr == 0.5 ~ "0.5",
  clean.tests$cdr == 1 ~ ">=1"
)

clean.tests$sex <- case_when(
  clean.tests$sex == 1 ~ "F",
  clean.tests$sex == 2 ~ "M",
)

clean.tests$cdr_group <- factor(clean.tests$cdr_group, levels = c('0', '0.5', '>=1'))
clean.tests$ne4s[clean.tests$ne4s == 9] = NA #missing
clean.tests$ne4s_group <- as.factor(clean.tests$ne4s)


no.tests <- is.na(clean.tests$y) & is.na(clean.tests$z)
clean.tests <- clean.tests[!no.tests,]
#clean.tests.cdr0 <- clean.tests %>% filter(cdr_group == '0')



## First application for which there is a 3 year time delay for the subsequent visits.
time.lag.3y <- as.Date("2003/01/01") - as.Date("2000/01/01")

clean.tests.lag.3y <- lag_pair_visit_label(clean.tests, time.lag.3y)

summary(clean.tests.lag.3y)

mean(clean.tests.lag.3y$age)
sd(clean.tests.lag.3y$age)

mean(clean.tests.lag.3y$sex == 2)
sd(clean.tests.lag.3y$sex == 2)


mean(clean.tests.lag.3y$educ_binary)
sd(clean.tests.lag.3y$educ_binary)

table(clean.tests.lag.3y$cdr_group)/nrow(clean.tests.lag.3y)

table(clean.tests.lag.3y$ne4s_group)/nrow(clean.tests.lag.3y)


summary(clean.tests.lag.3y)

y.ref.3y <- clean.tests.lag.3y$y
z.ref.3y <- clean.tests.lag.3y$z

clean.tests.lag.3y <-  categorize(clean.tests.lag.3y)

n.impute <- 50


ref.cols = c("age", "group")


list.3y <- ImputeOutcomeDifferences(clean.tests.lag.3y,
                                    y.ref.3y,z.ref.3y,n.impute,
                                    y.train,z.train,cond.y,cond.z,
                                    mu.y1,mu.z1,ref.cols,ker.set,
                                    R.bins = 1000, verbose = T)



list.3y.cw.opt <- ImputeOutcomeDifferences(clean.tests.lag.3y,
                                           y.ref.3y,z.ref.3y,n.impute,
                                           y.train,z.train,cond.y,cond.z,
                                           mu.y2,mu.z2,ref.cols,ker.set,
                                           R.bins = 1000,  verbose = T)

#### Impute Without Covariates
list.3y.nocov <- ImputeOutcomeDifferences(clean.tests.lag.3y,
                                          y.ref.3y,z.ref.3y,n.impute,
                                          y.train,z.train,cond.y,cond.z,
                                          mu.y1,mu.z1,
                                          R.bins = 1000, verbose = T)


list.3y.zscore <- ZScoreMatchDifferences(clean.tests.lag.3y,
                                         y.train,z.train)



list.3y.quantile <- QuantileMatchDifferences(clean.tests.lag.3y,
                                             y.train,z.train)



#outcome is the difference of scores
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
imp.reg.results.3y.cwopt <- ImputationRegressionGLM(fmla.1, list.3y.cw.opt$X, list.3y.cw.opt$Z)


imp.reg.results.3y$coefficients
imp.reg.results.3y$`cc-coefficients`

imp.reg.results.3y$`p-values`
imp.reg.results.3y$`cc-p-values`


imp.reg.results.3y.nocov$coefficients
imp.reg.results.3y.nocov$`cc-coefficients`

imp.reg.results.3y.nocov$`p-values`
imp.reg.results.3y.nocov$`cc-p-values`



#naive z matching
fit.z.score.3y <- glm(fmla.1, data = z.score.dat.3y)
z.score.coefs.3y <- fit.z.score.3y$coefficients
z.score.coefs.3y

#quantile matching
fit.quantile.3y <- glm(fmla.1, data = quantile.dat.3y)
quantile.match.coefs.3y <- fit.quantile.3y$coefficients
quantile.match.coefs.3y

#####





run.boot = F
if(run.boot){
  boot.model.3y <-  ImputationRegressionDifferencesGLMBootstrap(fmla.1, clean.tests.lag.3y, y.ref.3y, z.ref.3y,
                                                                n.impute, y.train, z.train,
                                                                cond.y, cond.z, mu.y1, mu.z1,
                                                                ref.cols, ker.set, R.bins = 1000,
                                                                threshold = 5 * 10^(-5), max.iter = 3,
                                                                B.boot = 200, verbose = T)
}

save.boot = F
if(save.boot){
  saveRDS(boot.model.3y, "data/bootstrap_regression_3y.rds")
}

load.boot = T
if(load.boot){
  boot.model.3y <- readRDS("data/bootstrap_regression_3y.rds")
}








term <- names(quantile.match.coefs.3y)

# complete case
m1 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.3y$`cc-coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.3y$`cc-variance`)))
m1$model = "Only MOCA"

m2 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.3y.nocov$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.3y.nocov$variance)))
m2$model = "DNOISe (no cov)"

m3 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.3y$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.3y$variance)))
m3$model = "DNOISe"

m4 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.3y.cwopt$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.3y.cwopt$variance)))
m4$model = "DNOISe (CW OPT)"

m5 <- data.frame("term" = term,
                 "estimate" = boot.model.3y$coefficients,
                 "std.error" = sqrt(diag(boot.model.3y$variance)))
m5$model = "DNOISe Bootstrap"


m6 <- data.frame("term" = term,
                 "estimate" = z.score.coefs.3y,
                 "std.error" = sqrt(diag(sandwich(fit.z.score.3y))))
m6$model = "Z Score Matching"

m7 <- data.frame("term" = term,
                 "estimate" = quantile.match.coefs.3y,
                 "std.error" = sqrt(diag(sandwich(fit.quantile.3y))))
m7$model = "Quantile Matching"

model.coefs.3y <- rbind(m1,m5,m6,m7)


imp.reg.results.3y$`cc-coefficients`

quantile.match.coefs.3y

### dot and whisker plot
plot.results = T
if(plot.results){
  library(dotwhisker)
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200

  p.3y <- dwplot(model.coefs.3y, show_intercept = FALSE)

  png(filename = "plots/3_year_impute_e4_education_interaction.png",
      width = png.width, height = png.height, res = png.res)
  p.3y


  # Close the pdf file
  dev.off()




  # tab.reg.3y.1 <- fit_to_table(imp.reg.results.3y)
  #
  # fmla.2 <- formula(outcome ~ ne4s_group)
  # imp.reg.results.3y.2 <- ImputationRegressionGLM(fmla.2, X, Z.imp)
  # tab.reg.3y.2 <- fit_to_table(imp.reg.results.3y.2)
  #
  # fmla.3 <- formula(outcome ~ educ_binary)
  # imp.reg.results.3y.3 <- ImputationRegressionGLM(fmla.3, X, Z.imp)
  # tab.reg.3y.3 <- fit_to_table(imp.reg.results.3y.3)
  #
  # fmla.4 <- formula(outcome ~ sex)
  # imp.reg.results.3y.4 <- ImputationRegressionGLM(fmla.4, X, Z.imp)
  # tab.reg.3y.4 <- fit_to_table(imp.reg.results.3y.4)
  #
  # fmla.5 <- formula(outcome ~ age)
  # imp.reg.results.3y.5 <- ImputationRegressionGLM(fmla.5, X, Z.imp)
  # tab.reg.3y.5 <- fit_to_table(imp.reg.results.3y.5)
}













#### 6 year trends
time.lag.6y <- as.Date("2006/01/01") - as.Date("2000/01/01")
time.window <- as.Date("2001/01/01") - as.Date("2000/01/01") # error of time window between times is 2 year
clean.tests.lag.6y <- lag_pair_visit_label(clean.tests, time.lag.6y, time.window = time.window)
y.ref.6y <- clean.tests.lag.6y$y
z.ref.6y <- clean.tests.lag.6y$z

clean.tests.lag.6y <-  categorize(clean.tests.lag.6y)

n.impute <- 50


ref.cols = c("age", "group")


list.6y <- ImputeOutcomeDifferences(clean.tests.lag.6y,
                                    y.ref.6y,z.ref.6y,n.impute,
                                    y.train,z.train,cond.y,cond.z,
                                    mu.y1,mu.z1,ref.cols,ker.set,
                                    R.bins = 1000, verbose = T)



list.6y.cw.opt <- ImputeOutcomeDifferences(clean.tests.lag.6y,
                                           y.ref.6y,z.ref.6y,n.impute,
                                           y.train,z.train,cond.y,cond.z,
                                           mu.y2,mu.z2,ref.cols,ker.set,
                                           R.bins = 1000,  verbose = T)

#### Impute Without Covariates
list.6y.nocov <- ImputeOutcomeDifferences(clean.tests.lag.6y,
                                          y.ref.6y,z.ref.6y,n.impute,
                                          y.train,z.train,cond.y,cond.z,
                                          mu.y1,mu.z1,
                                          R.bins = 1000, verbose = T)


list.6y.zscore <- ZScoreMatchDifferences(clean.tests.lag.6y,
                                         y.train,z.train)



list.6y.quantile <- QuantileMatchDifferences(clean.tests.lag.6y,
                                             y.train,z.train)



z.score.dat.6y <- cbind(list.6y.zscore$Z, list.6y.zscore$X)
colnames(z.score.dat.6y)[1] = "outcome"

quantile.dat.6y <- cbind(list.6y.quantile$Z, list.6y.quantile$X)
colnames(quantile.dat.6y)[1] = "outcome"

fmla <- formula(outcome ~ age + sex + educ_binary + ne4s_group)
X <- list.6y$X
Z.imp <- list.6y$Z

X.nocov <- list.6y.nocov$X
Z.imp.nocov <- list.6y.nocov$Z

fmla.1 <- formula(outcome ~ age + sex +  educ_binary *ne4s_group + cdr_group)
#fmla.1 <- formula(outcome ~ age)

# Here there are no individuals with scores over an 8 year period which measured MOCA at each time
imp.reg.results.6y <- ImputationRegressionGLM(fmla.1, X, Z.imp, fit.cc = F)
imp.reg.results.6y.nocov <- ImputationRegressionGLM(fmla.1, X.nocov, Z.imp.nocov, fit.cc = F)
imp.reg.results.6y.cwopt <- ImputationRegressionGLM(fmla.1, list.6y.cw.opt$X, list.6y.cw.opt$Z, fit.cc = F)



imp.reg.results.6y$coefficients
imp.reg.results.6y$`cc-coefficients`

imp.reg.results.6y$`p-values`
imp.reg.results.6y$`cc-p-values`


imp.reg.results.6y.nocov$coefficients
imp.reg.results.6y.nocov$`cc-coefficients`

imp.reg.results.6y.nocov$`p-values`
imp.reg.results.6y.nocov$`cc-p-values`



#naive z matching
fit.z.score.6y <- glm(fmla.1, data = z.score.dat.6y)
z.score.coefs.6y <- fit.z.score.6y$coefficients
z.score.coefs.6y

#quantile matching
fit.quantile.6y <- glm(fmla.1, data = quantile.dat.6y)
quantile.match.coefs.6y <- fit.quantile.6y$coefficients
quantile.match.coefs.6y

#####





run.boot = F
if(run.boot){
  boot.model.6y <-  ImputationRegressionDifferencesGLMBootstrap(fmla.1, clean.tests.lag.6y, y.ref.6y, z.ref.6y,
                                                                n.impute, y.train, z.train,
                                                                cond.y, cond.z, mu.y1, mu.z1,
                                                                ref.cols, ker.set, R.bins = 1000,
                                                                threshold = 5 * 10^(-5), max.iter = 3,
                                                                B.boot = 200, verbose = T)
}

load.boot = T
if(load.boot){
  boot.model.6y <- readRDS("data/bootstrap_regression_6y.rds")
}


boot.model.6y$coefficients - imp.reg.results.6y$coefficients
boot.model.6y$coefficients - imp.reg.results.6y.nocov$coefficients
boot.model.6y$coefficients - imp.reg.results.6y.cwopt$coefficients

save.boot = F
if(save.boot){
  saveRDS(boot.model.6y, "data/bootstrap_regression_6y.rds")
}







term <- names(quantile.match.coefs.6y)

# complete case
# m1 <- data.frame("term" = term,
#                  "estimate" = imp.reg.results.6y$`cc-coefficients`,
#                  "std.error" = sqrt(diag(imp.reg.results.6y$`cc-variance`)))
# m1$model = "Only MOCA"

m2 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.6y.nocov$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.6y.nocov$variance)))
m2$model = "DNOISe (no cov)"

m3 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.6y$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.6y$variance)))
m3$model = "DNOISe"

#TODO: Fix these with the corresponding models.
m4 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.6y.cwopt$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.6y.cwopt$variance)))
m4$model = "DNOISe (CW OPT)"

m5 <- data.frame("term" = term,
                 "estimate" = boot.model.6y$coefficients,
                 "std.error" = sqrt(diag(boot.model.6y$variance)))
m5$model = "DNOISe Bootstrap"


m6 <- data.frame("term" = term,
                 "estimate" = z.score.coefs.6y,
                 "std.error" = sqrt(diag(sandwich(fit.z.score.6y))))
m6$model = "Z Score Matching"

m7 <- data.frame("term" = term,
                 "estimate" = quantile.match.coefs.6y,
                 "std.error" = sqrt(diag(sandwich(fit.quantile.6y))))
m7$model = "Quantile Matching"

# no complete cases exist
model.coefs.6y <- rbind(m5,m6,m7)


imp.reg.results.6y$`cc-coefficients`

quantile.match.coefs.6y

### dot and whisker plot
plot.results = T
if(plot.results){
  library(dotwhisker)
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200

  p.6y <- dwplot(model.coefs.6y, show_intercept = FALSE)

  png(filename = "plots/6_year_impute_e4_education_interaction.png",
      width = png.width, height = png.height, res = png.res)
  p.6y


  # Close the pdf file
  dev.off()




  # tab.reg.6y.1 <- fit_to_table(imp.reg.results.6y)
  #
  # fmla.2 <- formula(outcome ~ ne4s_group)
  # imp.reg.results.6y.2 <- ImputationRegressionGLM(fmla.2, X, Z.imp)
  # tab.reg.6y.2 <- fit_to_table(imp.reg.results.6y.2)
  #
  # fmla.3 <- formula(outcome ~ educ_binary)
  # imp.reg.results.6y.3 <- ImputationRegressionGLM(fmla.3, X, Z.imp)
  # tab.reg.6y.3 <- fit_to_table(imp.reg.results.6y.3)
  #
  # fmla.4 <- formula(outcome ~ sex)
  # imp.reg.results.6y.4 <- ImputationRegressionGLM(fmla.4, X, Z.imp)
  # tab.reg.6y.4 <- fit_to_table(imp.reg.results.6y.4)
  #
  # fmla.5 <- formula(outcome ~ age)
  # imp.reg.results.6y.5 <- ImputationRegressionGLM(fmla.5, X, Z.imp)
  # tab.reg.6y.5 <- fit_to_table(imp.reg.results.6y.5)
}













#### 8-10 year trends
time.lag.8y <- as.Date("2008/01/01") - as.Date("2000/01/01")
time.window <- as.Date("2002/01/01") - as.Date("2000/01/01") # error of time window between times is 2 year
clean.tests.lag.8y <- lag_pair_visit_label(clean.tests, time.lag.8y, time.window = time.window)

summary(clean.tests.lag.8y)

mean(clean.tests.lag.8y$age)
sd(clean.tests.lag.8y$age)

mean(clean.tests.lag.8y$sex == 2)
sd(clean.tests.lag.8y$sex == 2)


mean(clean.tests.lag.8y$educ_binary)
sd(clean.tests.lag.8y$educ_binary)

table(clean.tests.lag.8y$cdr_group)/nrow(clean.tests.lag.8y)

table(clean.tests.lag.8y$ne4s_group)/nrow(clean.tests.lag.8y)



#colMeans(as.data.frame(clean.tests.lag.8y))


y.ref.8y <- clean.tests.lag.8y$y
z.ref.8y <- clean.tests.lag.8y$z

clean.tests.lag.8y <-  categorize(clean.tests.lag.8y)

n.impute <- 50


ref.cols = c("age", "group")


list.8y <- ImputeOutcomeDifferences(clean.tests.lag.8y,
                                    y.ref.8y,z.ref.8y,n.impute,
                                    y.train,z.train,cond.y,cond.z,
                                    mu.y1,mu.z1,ref.cols,ker.set,
                                    R.bins = 1000, verbose = T)



list.8y.cw.opt <- ImputeOutcomeDifferences(clean.tests.lag.8y,
                                           y.ref.8y,z.ref.8y,n.impute,
                                           y.train,z.train,cond.y,cond.z,
                                           mu.y2,mu.z2,ref.cols,ker.set,
                                           R.bins = 1000,  verbose = T)

#### Impute Without Covariates
list.8y.nocov <- ImputeOutcomeDifferences(clean.tests.lag.8y,
                                          y.ref.8y,z.ref.8y,n.impute,
                                          y.train,z.train,cond.y,cond.z,
                                          mu.y1,mu.z1,
                                          R.bins = 1000, verbose = T)


list.8y.zscore <- ZScoreMatchDifferences(clean.tests.lag.8y,
                                         y.train,z.train)



list.8y.quantile <- QuantileMatchDifferences(clean.tests.lag.8y,
                                             y.train,z.train)



z.score.dat.8y <- cbind(list.8y.zscore$Z, list.8y.zscore$X)
colnames(z.score.dat.8y)[1] = "outcome"

quantile.dat.8y <- cbind(list.8y.quantile$Z, list.8y.quantile$X)
colnames(quantile.dat.8y)[1] = "outcome"

fmla <- formula(outcome ~ age + sex + educ_binary + ne4s_group)
X <- list.8y$X
Z.imp <- list.8y$Z

X.nocov <- list.8y.nocov$X
Z.imp.nocov <- list.8y.nocov$Z

fmla.1 <- formula(outcome ~ age + sex +  educ_binary *ne4s_group + cdr_group)
#fmla.1 <- formula(outcome ~ age)

# Here there are no individuals with scores over an 8 year period which measured MOCA at each time
imp.reg.results.8y <- ImputationRegressionGLM(fmla.1, X, Z.imp, fit.cc = F)
imp.reg.results.8y.nocov <- ImputationRegressionGLM(fmla.1, X.nocov, Z.imp.nocov, fit.cc = F)
imp.reg.results.8y.cwopt <- ImputationRegressionGLM(fmla.1, list.8y.cw.opt$X, list.8y.cw.opt$Z, fit.cc = F)



imp.reg.results.8y$coefficients
imp.reg.results.8y$`cc-coefficients`

imp.reg.results.8y$`p-values`
imp.reg.results.8y$`cc-p-values`


imp.reg.results.8y.nocov$coefficients
imp.reg.results.8y.nocov$`cc-coefficients`

imp.reg.results.8y.nocov$`p-values`
imp.reg.results.8y.nocov$`cc-p-values`



#naive z matching
fit.z.score.8y <- glm(fmla.1, data = z.score.dat.8y)
z.score.coefs.8y <- fit.z.score.8y$coefficients
z.score.coefs.8y

#quantile matching
fit.quantile.8y <- glm(fmla.1, data = quantile.dat.8y)
quantile.match.coefs.8y <- fit.quantile.8y$coefficients
quantile.match.coefs.8y

#####




run.boot = F
if(run.boot){
  boot.model.8y <-  ImputationRegressionDifferencesGLMBootstrap(fmla.1, clean.tests.lag.8y, y.ref.8y, z.ref.8y,
                                                                n.impute, y.train, z.train,
                                                                cond.y, cond.z, mu.y1, mu.z1,
                                                                ref.cols, ker.set, R.bins = 1000,
                                                                threshold = 5 * 10^(-5), max.iter = 3,
                                                                B.boot = 200, verbose = T)
}
load.boot = T
if(load.boot){
  boot.model.8y <- readRDS("data/bootstrap_regression_8y.rds")
}

boot.model.8y$coefficients - imp.reg.results.8y$coefficients
boot.model.8y$coefficients - imp.reg.results.8y.nocov$coefficients
boot.model.8y$coefficients - imp.reg.results.8y.cwopt$coefficients

save.boot = F
if(save.boot){
  saveRDS(boot.model.8y, "data/bootstrap_regression_8y.rds")
}








term <- names(quantile.match.coefs.8y)

# complete case
# m1 <- data.frame("term" = term,
#                  "estimate" = imp.reg.results.8y$`cc-coefficients`,
#                  "std.error" = sqrt(diag(imp.reg.results.8y$`cc-variance`)))
# m1$model = "Only MOCA"

m2 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.8y.nocov$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.8y.nocov$variance)))
m2$model = "DNOISe (no cov)"

m3 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.8y$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.8y$variance)))
m3$model = "DNOISe"

#TODO: Fix these with the corresponding models.
m4 <- data.frame("term" = term,
                 "estimate" = imp.reg.results.8y.cwopt$`coefficients`,
                 "std.error" = sqrt(diag(imp.reg.results.8y.cwopt$variance)))
m4$model = "DNOISe (CW OPT)"

m5 <- data.frame("term" = term,
                 "estimate" = boot.model.8y$coefficients,
                 "std.error" = sqrt(diag(boot.model.8y$variance)))
m5$model = "DNOISe Bootstrap"


m6 <- data.frame("term" = term,
                 "estimate" = z.score.coefs.8y,
                 "std.error" = sqrt(diag(sandwich(fit.z.score.8y))))
m6$model = "Z Score Matching"

m7 <- data.frame("term" = term,
                 "estimate" = quantile.match.coefs.8y,
                 "std.error" = sqrt(diag(sandwich(fit.quantile.8y))))
m7$model = "Quantile Matching"


# no complete cases exist
model.coefs.8y <- rbind(m5,m6,m7)


imp.reg.results.8y$`cc-coefficients`

quantile.match.coefs.8y

### dot and whisker plot
plot.results = F
if(plot.results){
  library(dotwhisker)
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200

  p.8y <- dwplot(model.coefs.8y, show_intercept = FALSE)

  png(filename = "plots/8_year_impute_e4_education_interaction.png",
      width = png.width, height = png.height, res = png.res)
  p.8y


  # Close the pdf file
  dev.off()


  # tab.reg.8y.1 <- fit_to_table(imp.reg.results.8y)
  #
  # fmla.2 <- formula(outcome ~ ne4s_group)
  # imp.reg.results.8y.2 <- ImputationRegressionGLM(fmla.2, X, Z.imp)
  # tab.reg.8y.2 <- fit_to_table(imp.reg.results.8y.2)
  #
  # fmla.3 <- formula(outcome ~ educ_binary)
  # imp.reg.results.8y.3 <- ImputationRegressionGLM(fmla.3, X, Z.imp)
  # tab.reg.8y.3 <- fit_to_table(imp.reg.results.8y.3)
  #
  # fmla.4 <- formula(outcome ~ sex)
  # imp.reg.results.8y.4 <- ImputationRegressionGLM(fmla.4, X, Z.imp)
  # tab.reg.8y.4 <- fit_to_table(imp.reg.results.8y.4)
  #
  # fmla.5 <- formula(outcome ~ age)
  # imp.reg.results.8y.5 <- ImputationRegressionGLM(fmla.5, X, Z.imp)
  # tab.reg.8y.5 <- fit_to_table(imp.reg.results.8y.5)
}





plot.results = F
if(plot.results){
  library(dotwhisker)
  library(ggpubr)
  png.width = 1200
  png.height = 1000
  png.res = 200


  coef.vec <- c(imp.reg.results.3y$`cc-coefficients`,
                boot.model.3y$coefficients,
                z.score.coefs.3y,
                quantile.match.coefs.3y,
                boot.model.8y$coefficients,
                z.score.coefs.8y,
                quantile.match.coefs.8y)
  sd.vec <- c(sqrt(diag(imp.reg.results.3y$`cc-variance`)),
              sqrt(diag(boot.model.3y$variance)),
              sqrt(diag(sandwich(fit.z.score.3y))),
              sqrt(diag(sandwich(fit.quantile.3y))),
              sqrt(diag(boot.model.8y$variance)),
              sqrt(diag(sandwich(fit.z.score.8y))),
              sqrt(diag(sandwich(fit.quantile.8y))))
  term <- names(quantile.match.coefs.3y)

  coef.names <- rep(term, times = 7)

  n.coef <- length(term)
  coef.x <- rep(seq(1,n.coef), times = 7)
  coef.x <- coef.x + rep(seq(-3,3)/(14), each = n.coef)
  method.names <- c(rep(c("Only MOCA","DNOISe Bootstrap","Z Score Matching","Quantile Matching"), each = n.coef),
                   rep(c("DNOISe Bootstrap","Z Score Matching","Quantile Matching"), each = n.coef))

  dataset <- c(rep("3-4 year", 4*n.coef),
               rep("8-10 year", 3*n.coef))

  plot.data <- data.frame("Coefficients" = coef.vec,
                          "Xcoef" = coef.x,
                          "sd" = sd.vec,
                          "Parameters" = coef.names,
                          "Method" = method.names,
                          "Dataset" = dataset
  )
  plot.data <- plot.data %>% filter(Method %in% c("Only MOCA","DNOISe Bootstrap","Z Score Matching"))
  plot.data$Method <- factor(plot.data$Method, levels = c("Only MOCA","DNOISe Bootstrap","Z Score Matching","Quantile Matching"))
  title = "Cognitive Reserve Regression Coefficients"

  plt <- ggplot(plot.data,aes(x = Xcoef, y = Coefficients, group = interaction(Dataset, Method),
                color = Method, linetype = Dataset)) +
      geom_point(aes(shape=Dataset), size = 2) +
      geom_errorbar(aes(ymin = Coefficients - 2*sd, ymax = Coefficients + 2*sd), width = 0.2) +
      ggtitle(title) +
      xlab("Parameter") +
      scale_linetype_manual(values=c("twodash", "solid")) +
      scale_x_continuous(labels=term, breaks=seq(1,n.coef))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # add levels to coef.names
  png(filename = paste0("plots/impute_e4_education_interaction.png"),
      width = (1.8)*png.width, height = png.height, res = png.res)

  plt
  # Close the pdf file
  dev.off()

  plot.data.subset <- plot.data %>% filter(Parameters %in% c("age","sexM","educ_binary"))
  plt <- ggplot(plot.data.subset,aes(x = Xcoef, y = Coefficients, group = interaction(Dataset, Method),
                              color = Method, linetype = Dataset)) +
    geom_point(aes(shape=Dataset), size = 2) +
    geom_errorbar(aes(ymin = Coefficients - 2*sd, ymax = Coefficients + 2*sd), width = 0.2) +
    ggtitle(paste0(title," Subset")) +
    xlab("Parameter") +
    scale_linetype_manual(values=c("twodash", "solid")) +
    scale_x_continuous(labels=term, breaks=seq(1,n.coef))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # add levels to coef.names
  png(filename = paste0("plots/impute_e4_education_interaction_subset.png"),
      width = (1.3)*png.width, height = png.height, res = png.res)

  plt
  # Close the pdf file
  dev.off()



     }







