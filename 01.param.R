
# Load packages and data --------------------------------------------------
library(readstata13)
library(tidyverse)
library(MASS)

load("data/baseline_inj.Rdata")

ego <- readRDS("data/ego_2014_2020.RData")


# Parameters --------------------------------------------------------------
pop <- n <-  10000 #Size of hypothetical population for simulation
n_years <- 11 #simulate 11 years with the first year as burn-in period 
time_steps <- n_years*12 #each time step is a month 

## Detail and source of parameters can be found in Supplementary Table 1

params <- list()

#Demographic 
params$hcv_prev <- 0.466 #HCV prevalence within ALIVE cohort in 2018
params$hiv_prev <- 0.229 #HIV prevalence within ALIVE cohort in 2018
params$coinf_prev <- 0.087 #Prevalence of co-infection in ALIVE cohort 
params$hcv_prev_general <- 0.00880 #HCV prevalence in Maryland general population
params$hiv_prev_general <- 0.00643 #HIV prevalence in Maryland general population
params$m0f1 <- 0.324 #proportion of female in population
params$age_dist <- c(9.7, 0.2) #age of population at beginning of simulation Gamma(9.7, 0.2)
params$age_1st_inj_dist <- c(10.8, 0.5) #age of an individual's first injection Gamma(10.8, 0.5)
params$mortality <- c(rep(0,18), rep(25.7, 50-18), rep(36.6, 60-50), rep(54.1, 100-60)) #mortality rate per 1000PY 
params$mortality_covid <- c(rep(0,18), rep(26.6, 50-18), rep(46.6, 60-50), rep(48.5, 100-60)) #mortality rate during pandemic
params$mort_hiv <- 1.94 #mortality multiplier for people living with HIV
params$mort_hcv <- 1.50 #mortality multiplier for people living with HCV 
params$ltf <- c(0.33, 0.31, 0.27, 0.28) #differential probability of being lost to follow-up during pandemic by cessation trajectories

#Network
params$degree_pre <- 3.37 #average node degree (average number of drug-use ties each person has)
params$edges <- params$degree_pre*pop/2 #total number of edges
params$edge_duration <- 63.1 #average duration of edges (in months)
params$concurrent <- (179-9-38)/179*pop #Nodes with 2+ ties 
params$nodematch <- params$edges * (0.72) #Assortative mixing by age group 

#Injecting behavior
params$inj30d_dist <- 0.032 #monthly injection frequency Exponential(0.032)
params$prob_inj <- 0.79 #probability of injecting with drug-use ties 
params$inj_share <- 0.80 #proportion of monthly injections done with other PWID
params$syringe_share <- 0.33 #base probability of syringe sharing for each injection done with others

#Cessation trajectory
params$cess_dist <- c(42,13,26,19) #frequency of 4 cessation trajectories in population
params$cess_begin_dist <- c(2,240) #number of months before cessation trajectory begins Uniform(2, 240)

#Sexual behavior 
params$sexual_male <- 0.18 #probability of male having an IDU sexual tie
params$sexual_female <- 0.18*1.3 #probability of female having an IDU sexual tie
params$sex30d_dist <- 12 #average number of unprotected sexual acts per month Poisson(12)

#HIV transmission
params$hiv_blood <- 0.0063 #HIV transmission rate per shared injection 
params$hiv_sexual <- mean(c(0.04,0.08))/100 #HIV transmission rate per sexual act 
params$art_effect_s <- 0.94 #decrease in transmission with viral suppression
params$art_effect_f <- 0.73 #decrease in transmission with ART but no viral suppression
params$art_effect_given_hcv <- 0.846 #decrease in ART effect when HCV+
params$hiv_acute <- 5.3 #multiplier of HIV transmission rate when in acute stage 

#HIV treatment 
params$hiv_test <- 0.056 #monthly HIV testing probability 
params$hiv_test_ever <- 0.93 #proportion of PWID ever tested for HIV 
params$art_prescribe <- 0.83 #ART prescription probability once diagnosed
params$achieve_supp <- 0.915 #probability of achieving viral suppression when on ART 
params$art_adhere <- 0.85 #monthly ART adherence probability

#HCV transmission
params$hcv_blood <- 0.025 #HCV transmission rate per shared injection 
params$hcv_sexual <- 1/380000 #HCV transmission rate per sexual act
params$hcv_post_treat_dist <- c(1,9) #Reduction in HCV transmission after ever receiving treatment Beta(1,9)
params$hcv_given_hiv <- 6 #multiplier of HCV infection probability if HIV+
params$hcv_acute <- 2.7 #multiplier of HCV transmission rate when in acute stage 

#HCV treat
params$hcv_spont_clear_dist <- c(1.5, 2.67) #HCV spontaneous clearance probability Beta(1.5, 2.67)
params$spont_clear_given_hiv <- 0.33 #multiplier of HCV spontaneous clearance prob if HIV+
params$hcv_acute_length_dist <- c(2,8) #length of HCV acute stage Uniform(2,8)
params$hcv_treat <- 0.0088 #monthly HCV treatment probability 
params$hcv_treat_success <- 0.9 #probability of SVR once treated with DAA
params$hcv_ever_treated <- 0.365 #proportion of population ever treated for HCV at beginning of simulation

#Harm-reduction services 
params$moud_prob <- c(0.66, 0.51, 0.69, 0.72) #MOUD access prob by cessation trajectories
params$moud_effect_dist <- c(11.1, 9.5) #multiplier of injection events when on MOUD Beta(11.1, 9.5)
params$ssp_prob <- c(0.279, 0.043, 0.16, 0.50) #SSP access prob by cessation trajectories
params$ssp_effect_dist <- c(80.5, 5.1) #multiplier of syringe sharing prob when on SSP Beta(80.5, 5.1)



#Injection/cessation trajectories
stop_list <- data.frame(year_vec = seq(0,20,by = 1), 
                        early_line = c(1.0, 0.95, 0.7, 0.52, 0.4, 0.33, 0.25, 0.18, 0.13, 0.11, 0.09, rep(0.06, 10)), 
                        delay_line = c(1.0, 1.0, 0.97, 0.96, 0.95, 0.93, 0.90, 0.8, 0.7, 0.6, 0.41, 0.3, 0.24, 0.2, 0.17, 0.15, 0.13, 0.12, 0.11, 0.1, 0.1),
                        persistent_line = c(seq(1.0, 0.95, length = 15), seq(0.94,0.75,length=21-15)))

start <- list(beta = 2)
exp_model_early <- nls(early_line ~ exp(beta * year_vec), data = stop_list, start = start) 
exp_model_delay <- nls(delay_line ~ exp(beta * year_vec), data = stop_list, start = start)
exp_model_persistent <- nls(persistent_line ~ exp(beta * year_vec), data = stop_list, start = start)

### 1 = early cessation; 2 = delayed cessation; 3 = relapse; 4 = persistent 
test_years <- seq(0,20, length = 20*12)
exp_model_early_predict <- predict(exp_model_early, list(year_vec = test_years))
exp_model_delay_predict <- predict(exp_model_delay, list(year_vec = test_years))
exp_model_random_predict_prob <- 0.5
exp_model_persistent_predict <- predict(exp_model_persistent, list(year_vec = test_years))

cess_prob_values <- matrix(NA,length(exp_model_early_predict), 4)
cess_prob_values[,1] = exp_model_early_predict
cess_prob_values[,2] = exp_model_delay_predict
cess_prob_values[,3] = 0.5
cess_prob_values[,4] = exp_model_persistent_predict

params$cess_prob_values <- cess_prob_values #probability of injecting by each cessation trajectory 

params$cess_change_probs <- data.frame(early = c(0.607, 0.107, 0.143, 0.143),
                                       delay = c(0.902, 0.098, 0, 0),
                                       relapse = c(0.54, 0.069, 0.31, 0.08),
                                       persistent = c(0.316, 0.175, 0.105, 0.404)) #probabilities of switching between trajectories during pandemic




# Initial conditions ------------------------------------------------------

t0 <- list()

t0$age <- rgamma(pop*2, shape = params$age_dist[1], rate = params$age_dist[2])
t0$age <- t0$age[sample(which(t0$age > 22 & t0$age < 74), pop, replace = F)]
t0$age_group <- cut(t0$age, 
                    breaks = c(18, 24, 34, 44, 75), 
                    labels = c(1:4), 
                    right = TRUE, include.lowest = TRUE)

t0$age_1st_inj <- rgamma(pop*2, shape = params$age_1st_inj_dist[1], rate = params$age_1st_inj_dist[2])
t0$age_1st_inj <- t0$age_1st_inj[sample(which(t0$age_1st_inj >= 18 & t0$age_1st_inj <= 60), pop, replace = F)]

t0$gender <- rbinom(pop, 1, prob = params$m0f1)

t0$hiv_vec <- t0$hcv_vec <-  rep(0,pop) #set initial infected compartments based on known prevalence
t0$hiv_vec[sample(1:pop, pop*params$hiv_prev, replace=F)] <- 1 #HIV mono-infection
t0$hcv_vec[sample(which(t0$hiv_vec==1), pop*params$coinf_prev, replace = F)] <- 1 #co-infection
t0$hcv_vec[sample(which(t0$hiv_vec==0), pop*(params$hcv_prev - params$coinf_prev), replace = F)] <- 1 #HCV mono-infection
t0$status_vec <- rep("s", pop)
t0$status_vec[which(t0$hcv_vec==1 | t0$hiv_vec==1)] <- "i" #a disease status vector of either "s" or "i"

t0$inj30d <- rexp(pop*2, rate = params$inj30d_dist) %>% round()
t0$inj30d <- t0$inj30d[sample(which(t0$inj30d <= 200 & t0$inj30d >= 1), pop, replace=F)]

t0$sex30d <- rpois(pop, 12)

t0$cess_group <- sample(1:4, pop, replace=T, prob=params$cess_dist)
t0$cess_begin <- runif(pop, min = params$cess_begin_dist[1], max = params$cess_begin_dist[2])

t0$hcv_spont_clear_prob <- rbeta(pop, params$hcv_spont_clear_dist[1], params$hcv_spont_clear_dist[2])
t0$hcv_acute_length <- runif(pop, params$hcv_acute_length_dist[1], params$hcv_acute_length_dist[2])
t0$hcv_ever_treated <- rbinom(pop, 1, params$hcv_ever_treated)

#t0$moud <- sapply(cess_group, function(x) rbinom(1,size = 1, prob = params$moud_prob[x]))
#t0$ssp <- sapply(cess_group, function(x) rbinom(1,size = 1, prob = params$ssp_prob[x]))


## Set vectors for service availability
test_rate <- test_rate_covid  <- rep(params$hiv_test, time_steps)
test_rate_covid[37:(37+2)] <- params$hiv_test/2

art_prescribe <- art_prescribe_covid  <- rep(params$art_prescribe, time_steps)
art_prescribe_covid[37:(37+2)] <- params$art_prescribe/2

art_adhere <- art_adhere_covid <- rep(params$art_adhere, time_steps)

hcv_treat <- hcv_treat_covid <- rep(params$hcv_treat, time_steps)
hcv_treat_covid[37:(37+2)] <- (params$hcv_treat) * 0.5
hcv_treat_covid[(37+3):(37+11)] <- (params$hcv_treat) * 0.7

moud_on <- ssp_on <- moud_on_covid <- ssp_on_covid <- rep(1, time_steps)
moud_on_covid[37:(37+2)] <- ssp_on_covid[37:(37+2)] <- 0
moud_on_covid[(37+3):(37+11)] <- ssp_on_covid[(37+3):(37+11)] <- 0.5

##Set vector for on/off of pandemic-induced behavior changes
behav_change <- behav_change_covid <- rep(0, time_steps)
behav_change_covid[37:(37+40)] <- 1 



