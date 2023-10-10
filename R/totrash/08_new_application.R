

# suppose that we want to study the impact on education 
# and the aging decline across normal individuals
# Our initial data 
library(dplyr)
library(gee)
#library(geepack)
#library(multgee)

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

unif.3 <- scale_kernel(uniform_kernel,3)
unif.0.5 <- scale_kernel(uniform_kernel,0.5)
ker.set <- list(unif.3,unif.0.5)

R.bins.small = 1000
#TODO: make this faster if possible 
select.mu.y <- mu_selection_regression(Y,X.Y,cond.y,mu.set,R.bins.small, ker.set)
select.mu.z <- mu_selection_regression(Z,X.Z,cond.z,mu.set,R.bins.small, ker.set)


mu.y <- select.mu.y$opt.mu
mu.z <- select.mu.z$opt.mu
# 
mu.y <- 1
mu.z <- 1

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





#######

# 3 year gap 
time.gap.3y <- as.Date("2003/01/01") - as.Date("2000/01/01")

# 8 year gap 
time.gap.8y <- as.Date("2008/01/01") - as.Date("2000/01/01")


clean.tests <- clean.tests %>% 
  group_by(id) %>% 
  arrange(date) %>% 
  mutate(time_since_first = date - min(date)) 


tests.3y <- clean.tests %>% filter(time_since_first == 0 | time_since_first > time.gap.3y)

tests.3y <- tests.3y %>% 
  group_by(id) %>% 
  mutate(num_visits = n())
  
tests.3y <- tests.3y %>% 
  group_by(id) %>% 
  arrange(date) %>% 
  mutate(visit = row_number())

tests.3y <- tests.3y %>% 
  filter(visit %in% c(1,2)) %>% 
  filter(num_visits > 1)


tests.8y <- clean.tests %>% filter(time_since_first == 0 | time_since_first > time.gap.8y)

tests.8y <- tests.8y %>% 
  group_by(id) %>% 
  mutate(num_visits = n())

tests.8y <- tests.8y %>% 
  group_by(id) %>% 
  arrange(date) %>% 
  mutate(visit = row_number())

tests.8y <- tests.8y %>% 
  filter(visit %in% c(1,2)) %>% 
  filter(num_visits > 1)


design.cols <- c("age", "sex","educ","id","cdr","ne4s",
                 "date","group","educ_binary","cdr_group",       
                 "ne4s_group","time_since_first","num_visits","visit")
X.ref.3y <- tests.3y[,c("age", "group")]
y.ref.3y <- tests.3y[,"y"]
z.ref.3y <- tests.3y[,"z"]

z.impute.3y <- impute_dataset_2(X.ref.3y,y.ref.3y,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)


design.cols <- c("age", "sex","educ","id","cdr","ne4s",
                 "date","group","educ_binary","cdr_group",       
                 "ne4s_group","time_since_first","num_visits","visit")
X.ref.8y <- tests.8y[,c("age", "group")]
y.ref.8y <- tests.8y[,"y"]

z.impute.8y <- impute_dataset_2(X.ref.8y,y.ref.8y,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)






#######









clean.tests2 <- clean.tests
colnames(clean.tests2) <- paste0(colnames(clean.tests),"_2")
colnames(clean.tests2)[6] <- "id"

first.visits <- clean.tests %>% 
  group_by(id) %>% 
  arrange(date, by_group = T) %>% 
  mutate(idx = row_number())
first.visits <- first.visits %>% filter(idx == 1)

# 3 year gap 
time.gap <- as.Date("2003/01/01") - as.Date("2000/01/01")

joined.tests <- inner_join(first.visits, clean.tests2, by = "id")
joined.tests <- joined.tests %>% mutate(datediff = date_2 - date)
joined.tests <- joined.tests %>% filter(datediff > time.gap)
joined.tests <- joined.tests %>% 
  group_by(id) %>% 
  arrange(date_2, by_group = T) %>% 
  mutate(idx = row_number())
joined.tests <- joined.tests %>% filter(idx == 1)
joined.tests <- joined.tests %>% mutate(complete = 1*(!is.na(z) & !is.na(z_2)))

idx.missing.score1 <- is.na(joined.tests$y) & is.na(joined.tests$z)
idx.missing.score2 <- is.na(joined.tests$y_2) & is.na(joined.tests$z_2)
joined.tests <- joined.tests[!idx.missing.score1 & ! idx.missing.score2,]
mean(joined.tests$complete)

sum(joined.tests$complete)


joined.tests <- joined.tests %>% mutate(zdiff = z_2 - z)

complete.joined.tests <- joined.tests %>% filter(complete == 1)
n.impute = 5



# training and testing data

Z.impute <- impute_dataset_2(X.ref,y.ref,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)
X.ref

y.ref <- joined.tests$y
z.diff <- impute_dataset_3(joined.tests,y.ref,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)


X.ref <- joined.tests[,c("age", "sex", "id","educ_binary","cdr_group","ne4s_group", "complete")]

form.3y.1 <- formula(z ~ age + sex + educ_binary)
form.3y.2 <- formula(z ~ age + sex + educ_binary + educ_binary:age)
form.3y.3 <- formula(z ~ age + sex + educ_binary + educ_binary:age + cdr_group + ne4s_group)
form.3y.4 <- formula(z ~ age + sex + educ_binary + ne4s_group )



# TODO: look at the 10 year change 
# highlight that we could not do this before

#TODO: Appendix of univariate expressions. 

form.3y.univ.1 <- formula(z ~ age)
form.3y.univ.2 <- formula(z ~ sex)
form.3y.univ.3 <- formula(z ~ educ_binary)
form.3y.univ.4 <- formula(z ~ cdr_group)
form.3y.univ.5 <- formula(z ~ ne4s_group)


fit.3y.1 <- imputed_5ydiff(form.3y.1,z.diff,X.ref)
fit.3y.2 <- imputed_5ydiff(form.3y.2,z.diff,X.ref)
fit.3y.3 <- imputed_5ydiff(form.3y.3,z.diff,X.ref)
fit.3y.4 <- imputed_5ydiff(form.3y.4,z.diff,X.ref)

fit.3y.univ.1 <- imputed_5ydiff(form.3y.univ.1,z.diff,X.ref)
fit.3y.univ.2 <- imputed_5ydiff(form.3y.univ.2,z.diff,X.ref)
fit.3y.univ.3 <- imputed_5ydiff(form.3y.univ.3,z.diff,X.ref)
fit.3y.univ.4 <- imputed_5ydiff(form.3y.univ.4,z.diff,X.ref)
fit.3y.univ.5 <- imputed_5ydiff(form.3y.univ.5,z.diff,X.ref)

fit.3y.1$coefficients
fit.3y.1$`cc-coefficients`
sqrt(diag(fit.3y.1$variance))
sqrt(diag(fit.3y.1$`cc-variance`))
fit.3y.1$`z-scores`
fit.3y.1$`cc-z-scores`


#### Results Set 1
fit.3y.4$coefficients
fit.3y.4$`cc-coefficients`
sqrt(diag(fit.3y.4$variance))
sqrt(diag(fit.3y.4$`cc-variance`))
fit.3y.4$`z-scores`
fit.3y.4$`cc-z-scores`




# 6 year gap 
time.gap <- as.Date("2006/01/01") - as.Date("2000/01/01")

joined.tests <- inner_join(first.visits, clean.tests2, by = "id")
joined.tests <- joined.tests %>% mutate(datediff = date_2 - date)
joined.tests <- joined.tests %>% filter(datediff > time.gap)
joined.tests <- joined.tests %>% 
  group_by(id) %>% 
  arrange(date_2, by_group = T) %>% 
  mutate(idx = row_number())
joined.tests <- joined.tests %>% filter(idx == 1)
joined.tests <- joined.tests %>% mutate(complete = 1*(!is.na(z) & !is.na(z_2)))

idx.missing.score1 <- is.na(joined.tests$y) & is.na(joined.tests$z)
idx.missing.score2 <- is.na(joined.tests$y_2) & is.na(joined.tests$z_2)
joined.tests <- joined.tests[!idx.missing.score1 & ! idx.missing.score2,]
mean(joined.tests$complete)

sum(joined.tests$complete)


joined.tests <- joined.tests %>% mutate(zdiff = z_2 - z)

complete.joined.tests <- joined.tests %>% filter(complete == 1)
n.impute = 5



# training and testing data

#Z.impute <- impute_dataset_2(X.ref,y.ref,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)
X.ref



z.diff <- impute_dataset_3(joined.tests,y.ref,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)


X.ref <- joined.tests[,c("age", "sex", "id","educ_binary","cdr_group","ne4s_group", "complete")]

form.6y.1 <- formula(z ~ age + sex + educ_binary)
form.6y.2 <- formula(z ~ age + sex + educ_binary + educ_binary:age)
form.6y.3 <- formula(z ~ age + sex + educ_binary + cdr_group + ne4s_group)
form.6y.4 <- formula(z ~ age + sex + educ_binary + educ_binary:age + cdr_group + ne4s_group + ne4s_group:age)



# TODO: look at the 10 year change 
# highlight that we could not do this before

form.6y.univ.1 <- formula(z ~ age)
form.6y.univ.2 <- formula(z ~ sex)
form.6y.univ.3 <- formula(z ~ educ_binary)
form.6y.univ.4 <- formula(z ~ cdr_group)
form.6y.univ.5 <- formula(z ~ ne4s_group)


fit.6y.1 <- imputed_5ydiff(form.6y.1,z.diff,X.ref)
fit.6y.2 <- imputed_5ydiff(form.6y.2,z.diff,X.ref)
fit.6y.3 <- imputed_5ydiff(form.6y.3,z.diff,X.ref)
fit.6y.4 <- imputed_5ydiff(form.6y.4,z.diff,X.ref)

fit.6y.univ.1 <- imputed_5ydiff(form.6y.univ.1,z.diff,X.ref)
fit.6y.univ.2 <- imputed_5ydiff(form.6y.univ.2,z.diff,X.ref)
fit.6y.univ.3 <- imputed_5ydiff(form.6y.univ.3,z.diff,X.ref)
fit.6y.univ.4 <- imputed_5ydiff(form.6y.univ.4,z.diff,X.ref)
fit.6y.univ.5 <- imputed_5ydiff(form.6y.univ.5,z.diff,X.ref)

fit.6y.1$coefficients
fit.6y.1$`cc-coefficients`
sqrt(diag(fit.6y.1$variance))
sqrt(diag(fit.6y.1$`cc-variance`))
fit.6y.1$`z-scores`
fit.6y.1$`cc-z-scores`



#### Results Set 2
fit.6y.4$coefficients
fit.6y.4$`cc-coefficients`
sqrt(diag(fit.6y.4$variance))
sqrt(diag(fit.6y.4$`cc-variance`))
fit.6y.4$`z-scores`
fit.6y.4$`cc-z-scores`


#### Results Set 2
fit.6y.3$coefficients
fit.6y.3$`cc-coefficients`
sqrt(diag(fit.6y.3$variance))
sqrt(diag(fit.6y.3$`cc-variance`))
fit.6y.3$`z-scores`
fit.6y.3$`cc-z-scores`

#TODO: Model 3 is better, no need for additional interaction effects with age 


# TODO: try to make the imputation difference function as a single dataset
# so that someone can do a complete data analysis 







  
y.obs.tests <- clean.tests[!is.na(clean.tests$y),]
z.obs.tests <- clean.tests[!is.na(clean.tests$z),]




X.ref <- y.obs.tests[,c("age","group")]

X.id.ref <- y.obs.tests[,c("age","group", "id")]
sum(is.na(X.ref$group))
sum(is.na(X.id.ref$group))
# TODO: memory issues 
# may occur here. But we cannot make it sparse. 
# it may be possible to make it approximately sparse. 
# i.e. decompose the distribution to a flat amount, then 
# the rest. Will update if it becomes a problem

#latent.block.y <- reference_regression(X.ref, Y, X.Y,  A.matrix.y,A.tensor.y, mu.y, ker.set)
# reference block trained on the X.Y set! 
#latent.block.z <- reference_regression(X.ref, Z, X.Z,  A.matrix.z,A.tensor.z, mu.z, ker.set)


y.ref <- y.obs.tests$y
n.impute = 5



# training and testing data

Z.impute <- impute_dataset_2(X.ref,y.ref,n.impute,y.train,z.train,cond.y,cond.z,mu.y,mu.z,ker.set,R.bins,n.q)


Z.imp.long <- reshape_imputed_dataset(X.id.ref,Z.impute)

Z.full <- z.obs.tests[,c("z","age", "sex", "educ_binary", "cdr_group", "ne4s_group","id")]
X.Y.id <- y.obs.tests[,c("age", "sex", "educ_binary", "cdr_group", "ne4s_group","id")]

Z.full$bin_z <- 1*(Z.full$z >= 26)


# computes the standard error of the coefficients 
# using a simple Rubin's rules approach 


# TODO: also can just treat the problem as a continuous linear regression. 
# TODO: Consider APOE4 variable, make it into category. 
# At least one category can be included. (6 are defined, any can be in the category)
# TODO: 
# Z.full$sex <- Z.full$sex - 1

# Z.full$z <- ordered(Z.full$z, levels = seq(0,Ny))
# 
# naive.fit <- ordgee(z ~ age + sex + educ_binary,id = id, data = Z.full,corstr = "exchangeable", mean.link = "logit")

# tmp.fit <- ordLORgee(y ~ factor(trt) + factor(time) + factor(baseline),
#                      data = arthritis, id = id, repeated = time)


### ordinal logistic regression packages seem to not be able to 
# fit to the data, for now, we use a discretization of the output and study this 
# as a binary threshold for normal/abnormal 
# MOCA >= 26 is normal 




# Z.full$num_id = NA
# X.Y.id$num_id = NA
# for(j in 1:length(unique(clean.tests$id))){
#   idx <- which(unique(clean.tests$id)[j] == Z.full$id)
#   Z.full$num_id[idx] = j
#   
#   idx <- which(unique(clean.tests$id)[j] == X.Y.id$id)
#   X.Y.id$num_id[idx] = j
# }



# TODO: What is the 5 - 7 year trend in the MOCA score
# in the complete case, the sample size is small. 
# first measurement vs. 5-7th year measurement
# not a gee, but just a simple linear/ logistic regression
# we want to highlight that with a complete case, taking this difference 
# is actually going to help



form.bin.1 <- formula(bin_z ~ age + sex + educ_binary)
form.bin.2 <- formula(bin_z ~ age + sex + educ_binary + educ_binary:age)
form.bin.3 <- formula(bin_z ~ age + sex + educ_binary + educ_binary:age + cdr_group + ne4s_group)
form.bin.4 <- formula(bin_z ~ age + sex + educ_binary + educ_binary:age + cdr_group + ne4s_group + ne4s_group:age)


form.univ.bin.1 <- formula(bin_z ~ age)
form.univ.bin.2 <- formula(bin_z ~ sex)
form.univ.bin.3 <- formula(bin_z ~ educ_binary)
form.univ.bin.4 <- formula(bin_z ~ cdr_group)
form.univ.bin.5 <- formula(bin_z ~ ne4s_group)

bin.naive.fit <- gee(form.bin.4,
                     id = id, data = Z.full,corstr = "exchangeable", family = binomial(link = "logit"))

summary(bin.naive.fit)

impute.bin.fit.1 <- imputed_gee_binary(form.bin.1, Z.full, X.Y.id,Z.impute)
impute.bin.fit.2 <- imputed_gee_binary(form.bin.2, Z.full, X.Y.id,Z.impute)
impute.bin.fit.3 <- imputed_gee_binary(form.bin.3, Z.full, X.Y.id,Z.impute)
impute.bin.fit.4 <- imputed_gee_binary(form.bin.4, Z.full, X.Y.id,Z.impute)

impute.bin.fit.4$coefficients
impute.bin.fit.4$`cc-coefficients`
impute.bin.fit.4$variance
impute.bin.fit.4$`cc-variance`
impute.bin.fit.4$`z-scores`
impute.bin.fit.4$`cc-z-scores`


impute.univ.bin.fit.1 <- imputed_gee_binary(form.univ.bin.1, Z.full, X.Y.id,Z.impute)
impute.univ.bin.fit.2 <- imputed_gee_binary(form.univ.bin.2, Z.full, X.Y.id,Z.impute)
impute.univ.bin.fit.3 <- imputed_gee_binary(form.univ.bin.3, Z.full, X.Y.id,Z.impute)
impute.univ.bin.fit.4 <- imputed_gee_binary(form.univ.bin.4, Z.full, X.Y.id,Z.impute)
impute.univ.bin.fit.5 <- imputed_gee_binary(form.univ.bin.5, Z.full, X.Y.id,Z.impute)





form.1 <- formula(z ~ age + sex + educ_binary)
form.2 <- formula(z ~ age + sex + educ_binary + educ_binary:age)
form.3 <- formula(z ~ age + sex + educ_binary + educ_binary:age + cdr_group + ne4s_group)
form.4 <- formula(z ~ age + sex + educ_binary + educ_binary:age + cdr_group + ne4s_group + ne4s_group:age)

bin.naive.fit <- gee(form.4,
                     id = id, data = Z.full,corstr = "exchangeable")

summary(bin.naive.fit)

form.univ.1 <- formula(z ~ age)
form.univ.2 <- formula(z ~ sex)
form.univ.3 <- formula(z ~ educ_binary)
form.univ.4 <- formula(z ~ cdr_group)
form.univ.5 <- formula(z ~ ne4s_group)


impute.fit.1 <- imputed_gee_continuous(form.1, Z.full, X.Y.id,Z.impute)
impute.fit.2 <- imputed_gee_continuous(form.2, Z.full, X.Y.id,Z.impute)
impute.fit.3 <- imputed_gee_continuous(form.3, Z.full, X.Y.id,Z.impute)
impute.fit.4 <- imputed_gee_continuous(form.4, Z.full, X.Y.id,Z.impute)

impute.univ.fit.1 <- imputed_gee_continuous(form.univ.1, Z.full, X.Y.id,Z.impute)
impute.univ.fit.2 <- imputed_gee_continuous(form.univ.2, Z.full, X.Y.id,Z.impute)
impute.univ.fit.3 <- imputed_gee_continuous(form.univ.3, Z.full, X.Y.id,Z.impute)
impute.univ.fit.4 <- imputed_gee_continuous(form.univ.4, Z.full, X.Y.id,Z.impute)
impute.univ.fit.5 <- imputed_gee_continuous(form.univ.5, Z.full, X.Y.id,Z.impute)



impute.fit.4$coefficients
impute.fit.4$`cc-coefficients`
impute.fit.4$variance
impute.fit.4$`cc-variance`
impute.fit.4$`z-scores`
impute.fit.4$`cc-z-scores`








# can take a look at 

#only major difference seems to be the effect of education 
#this seems to result in a lower estimate of the effect of education 

# when imputing, the education seems to be 
# TODO: can also do regression for each variable separately and 
# include these as well. 
X.Y.id$id <- as.factor(X.Y.id$id)
Z.full$id <- as.factor(Z.full$id)

results <- imputed_gee_binary(formula,Z.full, X.Y.id,Z.impute)


results$coefficients



# The 









##### application to look at 5 year change as the outcome























library(stargazer)
stargazer(fit_to_table(impute.bin.fit.1), summary = F)
stargazer(fit_to_table(impute.bin.fit.2), summary = F)
stargazer(fit_to_table(impute.bin.fit.3), summary = F)
stargazer(fit_to_table(impute.bin.fit.4), summary = F)

stargazer(fit_to_table(impute.univ.bin.fit.1), summary = F)
stargazer(fit_to_table(impute.univ.bin.fit.2), summary = F)
stargazer(fit_to_table(impute.univ.bin.fit.3), summary = F)
stargazer(fit_to_table(impute.univ.bin.fit.4), summary = F)
stargazer(fit_to_table(impute.univ.bin.fit.5), summary = F)


stargazer(fit_to_table(impute.fit.1), summary = F)
stargazer(fit_to_table(impute.fit.2), summary = F)
stargazer(fit_to_table(impute.fit.3), summary = F)
stargazer(fit_to_table(impute.fit.4), summary = F)

stargazer(fit_to_table(impute.univ.fit.1), summary = F)
stargazer(fit_to_table(impute.univ.fit.2), summary = F)
stargazer(fit_to_table(impute.univ.fit.3), summary = F)
stargazer(fit_to_table(impute.univ.fit.4), summary = F)
stargazer(fit_to_table(impute.univ.fit.5), summary = F)
























