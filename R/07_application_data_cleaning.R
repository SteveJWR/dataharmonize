# Here we begin to clean the data
# We will convert the MMSE and the MOCA scores

source("R/01_functions.R")



library(dplyr)
library(gridExtra)

# dataset of all cognitive tests and covariates
# TODO: update this full description, request to NACC, to replicate this you may get a different number.

all_tests <- read.csv("Local_Data/investigator_nacc47.csv")

# indicates individuals with a known CDR score above 0
w_cog = all_tests$CDRGLOB>=0

# Between the age of 59 and 85
w_age = (all_tests$NACCAGE>59) & (all_tests$NACCAGE <= 85)

all_tests = all_tests[w_cog&w_age,]
n = nrow(all_tests)

score1 <- "NACCMMSE" #
score2 <- "MOCATOTS" #

# Maxima of each test scores
Ny <- 30
Nz <- 30

first_score1 <- paste0("first_",score1)
first_score2 <- paste0("first_",score2)


### Filtering Missing codes so all are just listed as NA
all_tests = all_tests %>% mutate(NACCMMSE = ifelse(((NACCMMSE >=0)&(NACCMMSE <= Ny)),NACCMMSE, NA))
all_tests = all_tests %>% mutate(MOCATOTS = ifelse(((MOCATOTS >=0)&(MOCATOTS <= Nz)),MOCATOTS, NA))


# Numbering by ID
all_tests <- all_tests %>% group_by(NACCID) %>%
  arrange(NACCID,VISITYR,VISITMO,VISITDAY) %>%
  mutate(row.idx = row_number())


#branching off y variable tests
col.sub.y <- c("NACCMMSE", "NACCID", "VISITYR","VISITMO","VISITDAY", "row.idx", "NACCAGE", "SEX", "EDUC")
tests.y <- all_tests[,col.sub.y]

# including only the covariates of age and sex
col.newnames.y <- c("y", "id", "year","mo","day", "row.idx","age", "sex", "educ")
colnames(tests.y) <- col.newnames.y

# only including those who actually took test y
tests.y <- tests.y[!is.na(tests.y$y),]

# filtering to only include a single first test
tests.y1 <- tests.y[tests.y$row.idx == 1,c("y", "age", "sex", "educ")]
# joining all pairs of visits
tests.y2 <- tests.y %>% left_join(tests.y[,c("y", "id", "year","mo","day", "row.idx")], by = 'id')
colnames(tests.y2) <- c("y1", "id", "year1","mo1","day1", "row.idx1","age", "sex", "educ",
                        "y2", "year2","mo2","day2", "row.idx2")

tests.y2 <- tests.y2 %>% filter(row.idx1 + 1 == row.idx2)
tests.y2 <- tests.y2 %>% mutate(date.diff = as.Date(paste0(year2, '/', mo2, '/',day2)) -
                                  as.Date(paste0(year1, '/', mo1, '/',day1)))

# including only differences of at most 500 days
tests.y2 <- tests.y2 %>% filter(date.diff <= 500)
tests.y2 <- tests.y2 %>% group_by(id)  %>%  mutate(row.pair.idx = row_number())
tests.y2 <- tests.y2[tests.y2$row.pair.idx == 1, ]



col.sub.z <- c("MOCATOTS", "NACCID", "VISITYR","VISITMO","VISITDAY", "row.idx", "NACCAGE", "SEX", "EDUC")
tests.z <- all_tests[,col.sub.z]
col.newnames.z <- c("z", "id", "year","mo","day", "row.idx","age", "sex", "educ")
colnames(tests.z) <- col.newnames.z


tests.z <- tests.z[!is.na(tests.z$z),]

# filtering to only include a single first test
tests.z1 <- tests.z[tests.z$row.idx == 1,c("z", "age", "sex", "educ")]

tests.z2 <- tests.z %>% left_join(tests.z[,c("z", "id", "year","mo","day", "row.idx")], by = 'id')
colnames(tests.z2) <- c("z1", "id", "year1","mo1","day1", "row.idx1","age", "sex", "educ",
                        "z2", "year2","mo2","day2", "row.idx2")

# filtering subsequent visits
tests.z2 <- tests.z2 %>% filter(row.idx1 + 1 == row.idx2)
tests.z2 <- tests.z2 %>% mutate(date.diff = as.Date(paste0(year2, '/', mo2, '/',day2)) -
                                  as.Date(paste0(year1, '/', mo1, '/',day1)))

# including only differences of at most 500 days
tests.z2 <- tests.z2 %>% filter(date.diff <= 500)
tests.z2 <- tests.z2 %>% group_by(id)  %>%  mutate(row.pair.idx = row_number())
tests.z2 <- tests.z2[tests.z2$row.pair.idx == 1, ]

# function to make a category out of the 4 possible groups in the
tests.y2 <- categorize(tests.y2)
tests.z2 <- categorize(tests.z2)

tests.y2.val <- tests.y2[,c('y1','age','y2','group')]
tests.z2.val <- tests.z2[,c('z1','age','z2','group')]

tests.y1.train <- categorize(tests.y1)
tests.z1.train <- categorize(tests.z1)


colnames(tests.y1.train)[1] <- "y"
colnames(tests.z1.train)[1] <- "z"


write.csv(tests.y2.val,"Data/NACCMMSE_validation.csv", row.names = F)
write.csv(tests.z2.val,"Data/MOCATOTS_validation.csv", row.names = F)

write.csv(tests.y1.train,"Data/NACCMMSE_training.csv", row.names = F)
write.csv(tests.z1.train,"Data/MOCATOTS_training.csv", row.names = F)






# Similarly Constructing the test data

# crosswalk Data set
cw.full <- read.csv(file = "Local_Data/crosswalk_wang08012018.csv")



# indicates cognitively normal individuals with 0
# or 0.5
w_cog = cw.full$CDRGLOB>=0

# Between the age of 59 and 85
w_age = (cw.full$NACCAGE>59) & (cw.full$NACCAGE <= 85)

cw.clean <- cw.full[w_cog & w_age,]

score1 <- "NACCMMSE"
score2 <- "MOCATOTS"

# Maxima of possible tests scores
Ny <- 30
Nz <- 30

### Filtering Missing codes so all are just listed as NA
cw.clean = cw.clean %>% mutate(NACCMMSE = ifelse(((NACCMMSE >=0)&(NACCMMSE <= Ny)),NACCMMSE, NA))
cw.clean = cw.clean %>% mutate(C1WMOCATOTS = ifelse(((C1WMOCATOTS >=0)&(C1WMOCATOTS <= Nz)),C1WMOCATOTS, NA))
cw.clean <- cw.clean[!is.na(cw.clean$NACCMMSE) &!is.na(cw.clean$C1WMOCATOTS) ,]

cw.columns <- c("NACCMMSE", "C1WMOCATOTS", "NACCAGE", "SEX", "EDUC")
cw.clean <- cw.clean[,cw.columns]
colnames(cw.clean) <- c("y","z","age","sex","educ")
cw.clean <- categorize(cw.clean)

write.csv(cw.clean,"Data/NACCMMSE_to_MOCATOTS_test.csv", row.names = F)

























##

































































