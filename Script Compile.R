#############################################

# Compiling NST data from SPARK 1b #

#############################################

#This Version: Version 1: 2024 10-11

#Note: will first recode into long format and then re-widen.

######## LOAD PACKAGES ###########

# Data Processing:
library(psych) #psych stats package, includes factor analysis
library(bestNormalize) #yeojohnson and other transformations
library(mice) #multiple imputation
#library(Matching) #propensity score matching
#library(MatchIt) #simpler matching
library(robustbase) #adjusted box plot
library(performance) #check_outliers for mahalanobis distance, influence
#library(NbClust) #Clustering - how many
library(Rtsne) #t-SNE dimensionality reduction

# General Stats:
library(effsize) #cohen d, etc.
library(lm.beta) #standardized beta
library(DescTools) #Fisher r to z tranformation
library(nlme) #mixed effects
#library(irr) # interrater correlation
library(RVAideMemoire) #Pairwise comparisons and other stats
#library(lme4) #expands on nlme with generalized linear regressions
#library(lavaan) #cFA, mediation analysis
library(moments) #skewness and kurtosis

# Graphing
library(gridExtra) #arranges ggplot grobs
library(Hmisc) #works with ggplot and others
library(corrplot) # correlation plots
library(ggpubr) #pairwise and statistics for ggplot
#library(hexbin) #Binning and plotting functions for hexagonal bins - network
#library(networkD3) #Creates D3 JavaScript network, tree, dendrogram graphs from R.
#library(tidytext) #Text mining
#library(plotly) #Interactive graphs
#library(repr)#String and binary representations - Network
library(RColorBrewer)
library(ggthemes) #Preset themes for ggplot
library(ggnewscale) #Can switch color scale
library(patchwork) #Compile ggplots

# Standard Packages (all projects)
library(plyr) #data manipulations
library(naniar) #na functions
library(arsenal) #tableby
library(tidyverse) #dplyr, ggplot, stringr, et al


set.seed(123)


######## LOAD DATA ###########

raw <- read.csv("NST1b_DATA_2024-09-24_1108-NOPHI.csv")
raw <- filter(raw, visit_1_t1t2_clinician_complete == 2)




######## Demo ###########

### Select

demo_raw <- select(raw, record_id, screen_sciddx, starts_with("demo"))


### Dx
demo_raw$demo_dx <- factor(demo_raw$screen_sciddx, 
                           levels=1:2, labels=c("schizophrenia", "schizoaffective"))
table(demo_raw$demo_dx)

### Sex
demo_raw$demo_sex <- factor(demo_raw$demo_sex, 
                           levels=1:2, labels=c("Male", "Female"))
table(demo_raw$demo_sex)

### gender
demo_raw$demo_gender <- factor(demo_raw$demo_gender, 
                            levels=1:3, labels=c("Man", "Woman", "Non-Binary"))
table(demo_raw$demo_gender)

### ethnicity
demo_raw$demo_ethnicity <- factor(demo_raw$demo_ethnicity, 
                               levels=1:2, labels=c("Not Hispanic", "Hispanic"))
table(demo_raw$demo_ethnicity)

### race

# 1-White, 2-Black, 3-NativeAm, 4-Asian, 5-Pacific, 6-Other, 7-Multiple, 8-Unknown
demo_raw$demo_race <- rep(NA, nrow(demo_raw))

for (i in 1:nrow(demo_raw)) {
      if (rowSums(demo_raw[i,7:12]) > 1) { #UPDATE columns!
            demo_raw$demo_race[i] <- 7
      }
      else if (demo_raw$demo_race___7[i] == 1) {
            demo_raw$demo_race[i] <- 7 
      }
      else if (demo_raw$demo_race___6[i] == 1) {
            demo_raw$demo_race[i] <- 6
      }
      else if (demo_raw$demo_race___5[i] == 1) {
            demo_raw$demo_race[i] <- 3
      }
      else if (demo_raw$demo_race___4[i] == 1) {
            demo_raw$demo_race[i] <- 5
      }
      else if (demo_raw$demo_race___3[i] == 1) {
            demo_raw$demo_race[i] <- 4
      }
      else if (demo_raw$demo_race___2[i] == 1) {
            demo_raw$demo_race[i] <- 2
      }
      else if (demo_raw$demo_race___1[i] == 1) {
            demo_raw$demo_race[i] <- 1
      }
      else {
            demo_raw$demo_race[i] <- 8
      }
}

demo_raw$demo_race <- factor(demo_raw$demo_race, levels = c(1,2,4,6,7), 
                             labels = c("White", "Black", "Asian", "Other",
                                        "Multiple"))
table(demo_raw$demo_race)


### Marital status
demo_raw$demo_marital <- factor(demo_raw$demo_marital, 
                            levels=c(1,2,4), labels=c("Never married", "Married", "Divorced"))
table(demo_raw$demo_marital)

### employment status
demo_raw$demo_employment <- factor(demo_raw$demo_employment, 
                                levels=c(1,2,4,5), labels=c("FT Employed", "PT Employed", "Unemployed", "Unemployed"))
table(demo_raw$demo_employment)

### insurance status
demo_raw$demo_insurance <- factor(demo_raw$demo_insurance, 
                                levels=c(0:2, 4,5), labels=c("Uninsured", "Medicare", "Medicaid", "Private", "Other"))
table(demo_raw$demo_insurance)

### treatment location
demo_raw$demo_txtype <- factor(demo_raw$demo_txtype, 
                                  levels=c(1:5), labels=c("Private Clinic", "Academic Program", "Psychosis Specialty Care",
                                                          "Community MH Center", "Other"))
table(demo_raw$demo_txtype)



### Checking others
demo_raw$demo_age
demo_raw$demo_education
demo_raw$demo_edumother
demo_raw$demo_edumother[which(demo_raw$demo_edumother==99)] <- NA
demo_raw$demo_edufather
demo_raw$demo_edumother[which(demo_raw$demo_edufather==99)] <- NA
table(demo_raw$demo_med_yn)
demo_raw$demo_meds
demo_raw$demo_hospitalizations
table(demo_raw$demo_insurance)


# Select Final Features
demo <- select(demo_raw, record_id, demo_dx, demo_age, demo_sex, demo_gender, demo_race, demo_ethnicity, 
               demo_education, demo_edumother, demo_edufather, 
               demo_marital, demo_employment, 
               demo_med_yn, demo_meds, demo_hospitalizations,
               demo_insurance, demo_txtype)





######## NST-ATT Attitudes and Progress toward goals ###########

### Set up empty df
nstatt_raw1 <- as.data.frame(matrix(nrow = 20, ncol = 12))
colnames(nstatt_raw1) <- c("record_id", "timepoint", "nstatt_goaleasy", "nstatt_goalimp", "nstatt_intenteasy", "nstatt_intentimp", 
                      "nstatt_confidenteasy", "nstatt_confidentimp", 
                      "nstatt_progresseasy", "nstatt_progressimp", "nstatt_progdesceasy", "nstatt_progdescimp")

nstatt_raw2 <- as.data.frame(matrix(nrow = 20, ncol = 12))
colnames(nstatt_raw2) <- c("record_id", "timepoint", "nstatt_goaleasy", "nstatt_goalimp", "nstatt_intenteasy", "nstatt_intentimp", 
                           "nstatt_confidenteasy", "nstatt_confidentimp", 
                           "nstatt_progresseasy", "nstatt_progressimp", "nstatt_progdesceasy", "nstatt_progdescimp")

nstatt_raw3 <- as.data.frame(matrix(nrow = 20, ncol = 12))
colnames(nstatt_raw3) <- c("record_id", "timepoint", "nstatt_goaleasy", "nstatt_goalimp", "nstatt_intenteasy", "nstatt_intentimp", 
                           "nstatt_confidenteasy", "nstatt_confidentimp", 
                           "nstatt_progresseasy", "nstatt_progressimp", "nstatt_progdesceasy", "nstatt_progdescimp")

### Id
nstatt_raw1$record_id <- raw$record_id
nstatt_raw2$record_id <- raw$record_id
nstatt_raw3$record_id <- raw$record_id

### timepoint
nstatt_raw1$timepoint <- rep(1, 20)
nstatt_raw2$timepoint <- rep(2, 20)
nstatt_raw3$timepoint <- rep(3, 20)

### Goals
nstatt_raw1$nstatt_goaleasy <- raw$nstatt_02taska_t1
nstatt_raw2$nstatt_goaleasy <- raw$nstatt_02taska_t1
nstatt_raw3$nstatt_goaleasy <- raw$nstatt_02taska_t1

nstatt_raw1$nstatt_goalimp <- raw$nstatt_03taskb_t1
nstatt_raw2$nstatt_goalimp <- raw$nstatt_03taskb_t1
nstatt_raw3$nstatt_goalimp <- raw$nstatt_03taskb_t1

### Intent
nstatt_raw1$nstatt_intenteasy <- raw$nstatt_04intent_t1
nstatt_raw1$nstatt_intentimp <- raw$nstatt_05intent_t1

nstatt_raw2$nstatt_intenteasy <- raw$nstatt_01intent_t2
nstatt_raw2$nstatt_intentimp <- raw$nstatt_03intent_t2

### Confident
nstatt_raw1$nstatt_confidenteasy <- raw$nstatt_04confident_t1
nstatt_raw1$nstatt_confidentimp <- raw$nstatt_05confidentt_t1

nstatt_raw2$nstatt_confidenteasy <- raw$nstatt_01confident_t2
nstatt_raw2$nstatt_confidentimp <- raw$nstatt_03confident_t2

### Progress
nstatt_raw3$nstatt_progresseasy[which(raw$nstatt_01taska_t3 == 1)] <- 1
nstatt_raw3$nstatt_progresseasy[which(raw$nstatt_01taska_t3 == 2)] <- 0
nstatt_raw3$nstatt_progresseasy
raw$nstatt_01taska_t3

nstatt_raw3$nstatt_progressimp[which(raw$nstatt_02taskb_t3 == 1)] <- 1
nstatt_raw3$nstatt_progressimp[which(raw$nstatt_02taskb_t3 == 2)] <- 0

nstatt_raw3$nstatt_progdesceasy <- raw$nstatt_01taska_desc_t3
nstatt_raw3$nstatt_progdescimp <- raw$nstatt_02taskb_desc_t3

# Combine 3 separate bits
nstatt_raw1 %>% head()
nstatt_raw2 %>% head()
nstatt_raw3 %>% head()

nstatt <- rbind(nstatt_raw1, nstatt_raw2)
nstatt <- rbind(nstatt, nstatt_raw3)

#Fix cell
nstatt$nstatt_progressimp[which(nstatt$record_id==35 & nstatt$timepoint==3)] <- 1

# means
nstatt$nstatt_intentmean <- (nstatt$nstatt_intenteasy + nstatt$nstatt_intentimp)/2
nstatt$nstatt_confidentmean <- (nstatt$nstatt_confidenteasy + nstatt$nstatt_confidentimp)/2

# any progress
nstatt$nstatt_progressany <- nstatt$nstatt_progresseasy + nstatt$nstatt_progressimp

# factors
nstatt$nstatt_progresseasy <- as.factor(nstatt$nstatt_progresseasy)
nstatt$nstatt_progressimp <- as.factor(nstatt$nstatt_progressimp)
nstatt$nstatt_progressany <- as.factor(nstatt$nstatt_progressany)

nstatt

######## Semi-structured feedback ###########

### Set up empty df
semi_raw2 <- as.data.frame(matrix(nrow = 20, ncol = 12))
colnames(semi_raw2) <- c("record_id", "timepoint", "semi_first", "semi_most", "semi_least", "semi_negative", 
                         "semi_able", "semi_prescribe", "semi_otherfun", "semi_other", "semi_use", "semi_tell")

semi_raw3 <- as.data.frame(matrix(nrow = 20, ncol = 12))
colnames(semi_raw3) <- c("record_id", "timepoint", "semi_first", "semi_most", "semi_least", "semi_negative", 
                         "semi_able", "semi_prescribe", "semi_otherfun", "semi_other", "semi_use", "semi_tell")

### Id
semi_raw2$record_id <- raw$record_id
semi_raw3$record_id <- raw$record_id

### timepoint
semi_raw2$timepoint <- rep(2, 20)
semi_raw3$timepoint <- rep(3, 20)


### Fill In
semi_raw2$semi_first <- raw$semi_01impression_t2
semi_raw2$semi_most <- raw$semi_02most_t2
semi_raw2$semi_least <- raw$semi_03least_t2
semi_raw2$semi_negative <- raw$semi_04effects_t2
semi_raw2$semi_able <- raw$semi_05use_t2
semi_raw2$semi_prescribe <- raw$semi_06doctor_t2
semi_raw2$semi_otherfun <- raw$semi_07fun_t2
semi_raw2$semi_other <- raw$semi_08feedback_t2

semi_raw3$semi_use <- raw$semi_01use_t3
semi_raw3$semi_tell <- raw$semi_02share_t3
semi_raw3$semi_negative <- raw$semi_03effects_t3
semi_raw3$semi_prescribe <- raw$semi_04doctor_t3
semi_raw3$semi_otherfun <- raw$semi_05fun_t3
semi_raw3$semi_other <- raw$semi_06feedback_t3

head(semi_raw2)
head(semi_raw3)

### Combine
semi <- rbind(semi_raw2, semi_raw3)


######## defeatist beliefs scale ###########

### Timepoint 1
dbs_raw1 <- select(raw, dbs_01respect_t1:dbs_05failure_t1)
dbs_raw1 <- cbind(rep(1, 20), dbs_raw1)
dbs_raw1 <- cbind(demo$record_id, dbs_raw1)

colnames(dbs_raw1) <- c( "record_id", "timepoint", "dbs_respect", "dbs_admire", "dbs_inferior", "dbs_nopoint", "dbs_failure")

dbs_raw1


### Timepoint 2
dbs_raw2 <- select(raw, dbs_01respect_t2:dbs_05failure_t2)
dbs_raw2 <- cbind(rep(2, 20), dbs_raw2)
dbs_raw2 <- cbind(demo$record_id, dbs_raw2)

colnames(dbs_raw2) <- c( "record_id", "timepoint", "dbs_respect", "dbs_admire", "dbs_inferior", "dbs_nopoint", "dbs_failure")

dbs_raw2


### Timepoint 3
dbs_raw3 <- select(raw, dbs_01respect_t3:dbs_05failure_t3)
dbs_raw3 <- cbind(rep(3, 20), dbs_raw3)
dbs_raw3 <- cbind(demo$record_id, dbs_raw3)

colnames(dbs_raw3) <- c( "record_id", "timepoint", "dbs_respect", "dbs_admire", "dbs_inferior", "dbs_nopoint", "dbs_failure")

dbs_raw3

### combine
dbs <- rbind(dbs_raw1, dbs_raw2)
dbs <- rbind(dbs, dbs_raw3)


### Calculate total score
dbs$dbs_total <- rowSums(dbs[,3:7])

dbs




######## self-esteem ###########

### Timepoint 1
bses_raw1 <- select(raw, bses_01successful_t1, bses_02attractive_t1, bses_03popular_t1, bses_04independent_t1,
                    bses_05honest_t1, bses_06desire_t1, bses_07strong_t1, bses_08smart_t1, bses_09power_t1, bses_10lovable_t1,
                    bses_11pleasant_t1, bses_12efficient_t1, bses_13responsible_t1, bses_14generous_t1, bses_15worthwhile_t1,
                    bses_16interesting_t1, bses_17knowledge_t1, bses_18good_t1)
bses_raw1 <- cbind(rep(1, 20), bses_raw1)
bses_raw1 <- cbind(demo$record_id, bses_raw1)

colnames(bses_raw1) <- c( "record_id", "timepoint", "bses_successful", "bses_attractive", "bses_popular", "bses_independent",
                          "bses_honest", "bses_desire", "bses_strong", "bses_smart", "bses_power", "bses_lovable",
                          "bses_pleasant", "bses_efficient", "bses_responsible", "bses_generous", "bses_worthwhile",
                          "bses_interesting", "bses_knowledge", "bses_good")

bses_raw1


### Timepoint 2
bses_raw2 <- select(raw, bses_01successful_t2, bses_02attractive_t2, bses_03popular_t2, bses_04independent_t2,
                    bses_05honest_t2, bses_06desire_t2, bses_07strong_t2, bses_08smart_t2, bses_09power_t2, bses_10lovable_t2,
                    bses_11pleasant_t2, bses_12efficient_t2, bses_13responsible_t2, bses_14generous_t2, bses_15worthwhile_t2,
                    bses_16interesting_t2, bses_17knowledge_t2, bses_18good_t2)
bses_raw2 <- cbind(rep(2, 20), bses_raw2)
bses_raw2 <- cbind(demo$record_id, bses_raw2)

colnames(bses_raw2) <- c("record_id", "timepoint", "bses_successful", "bses_attractive", "bses_popular", "bses_independent",
                          "bses_honest", "bses_desire", "bses_strong", "bses_smart", "bses_power", "bses_lovable",
                          "bses_pleasant", "bses_efficient", "bses_responsible", "bses_generous", "bses_worthwhile",
                          "bses_interesting", "bses_knowledge", "bses_good")

bses_raw2

### Timepoint 3
bses_raw3 <- select(raw, bses_01successful_t3, bses_02attractive_t3, bses_03popular_t3, bses_04independent_t3,
                    bses_05honest_t3, bses_06desire_t3, bses_07strong_t3, bses_08smart_t3, bses_09power_t3, bses_10lovable_t3,
                    bses_11pleasant_t3, bses_12efficient_t3, bses_13responsible_t3, bses_14generous_t3, bses_15worthwhile_t3,
                    bses_16interesting_t3, bses_17knowledge_t3, bses_18good_t3)
bses_raw3 <- cbind(rep(3, 20), bses_raw3)
bses_raw3 <- cbind(demo$record_id, bses_raw3)

colnames(bses_raw3) <- c("record_id", "timepoint", "bses_successful", "bses_attractive", "bses_popular", "bses_independent",
                         "bses_honest", "bses_desire", "bses_strong", "bses_smart", "bses_power", "bses_lovable",
                         "bses_pleasant", "bses_efficient", "bses_responsible", "bses_generous", "bses_worthwhile",
                         "bses_interesting", "bses_knowledge", "bses_good")

bses_raw3

### combine
bses <- rbind(bses_raw1, bses_raw2)
bses <- rbind(bses, bses_raw3)


### Calculate total score
colnames(bses)
bses$bses_total <- rowSums(bses[,3:20])

bses



######## AIM ###########

### Timepoint 2
aim_raw2 <- select(raw, aim_01approval_t2:aim_04welcome_t2)
aim_raw2 <- cbind(rep(2, 20), aim_raw2)
aim_raw2 <- cbind(demo$record_id, aim_raw2)

colnames(aim_raw2) <- c( "record_id", "timepoint", "aim_approval", "aim_appealing", "aim_like", "aim_welcome")

aim_raw2


### Timepoint 3
aim_raw3 <- select(raw, aim_01approval_t3:aim_04welcome_t3)
aim_raw3 <- cbind(rep(3, 20), aim_raw3)
aim_raw3 <- cbind(demo$record_id, aim_raw3)

colnames(aim_raw3) <- c( "record_id", "timepoint", "aim_approval", "aim_appealing", "aim_like", "aim_welcome")

aim_raw3

### combine
aim <- rbind(aim_raw2, aim_raw3)


### Calculate total score
colnames(aim)
aim$aim_mean <- rowSums(aim[,3:6]) / 4

aim



######## FIM ###########

### Timepoint 2
fim_raw2 <- select(raw, fim_01implement_t2:fim_04easy_t2)
fim_raw2 <- cbind(rep(2, 20), fim_raw2)
fim_raw2 <- cbind(demo$record_id, fim_raw2)

colnames(fim_raw2) <- c( "record_id", "timepoint", "fim_implementable", "fim_possible", "fim_doable", "fim_easy")

fim_raw2


### Timepoint 3
fim_raw3 <- select(raw, fim_01implement_t3:fim_04easy_t3)
fim_raw3 <- cbind(rep(3, 20), fim_raw3)
fim_raw3 <- cbind(demo$record_id, fim_raw3)

colnames(fim_raw3) <- c( "record_id", "timepoint", "fim_implementable", "fim_possible", "fim_doable", "fim_easy")

fim_raw3

### combine
fim <- rbind(fim_raw2, fim_raw3)


### Calculate total score
colnames(fim)
fim$fim_mean <- rowSums(fim[,3:6]) / 4

fim



######## Knit together & export ###########

### All data - Long
#demo, dbs, bses, aim, fim, nstatt, semi

long <- merge(demo, dbs, by = c("record_id"), all = TRUE)
long <- merge(long, bses, by = c("record_id", "timepoint"), all = TRUE)
long <- merge(long, aim, by = c("record_id", "timepoint"), all = TRUE)
long <- merge(long, fim, by = c("record_id", "timepoint"), all = TRUE)
long <- merge(long, nstatt, by = c("record_id", "timepoint"), all = TRUE)
long <- merge(long, semi, by = c("record_id", "timepoint"), all = TRUE)

View(long)


### Make wide format
colnames(long)
vars_topivot <- names(long)[18:72] # Check cols - should stuff thats tied to time, not participant level demo
wide <- pivot_wider(long, id_cols = record_id, names_from = timepoint, values_from = all_of(vars_topivot))

wide <- merge(demo, wide, by="record_id")


View(wide)

### Export all data
save(long, file = "long.R")
write.csv(long, file = "long.csv")

save(wide, file = "wide.R")
write.csv(wide, file = "wide.csv")


### Just feedback
justfeedback <- merge(demo, nstatt, by = c("record_id"), all = TRUE)
justfeedback <- merge(justfeedback, semi, by = c("record_id", "timepoint"), all = TRUE)
justfeedback <- justfeedback %>% relocate(timepoint, .after = record_id)

#View(justfeedback)
write.csv(justfeedback, "SPARK1b JustFeedback_V1original.csv")
