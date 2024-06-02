
# Load packages and data --------------------------------------------------
library(readstata13)
library(tidyverse)
library(MASS)

load("data/baseline_inj.Rdata")

ego <- readRDS("data/ego_2014_2020.RData")


# Parameters --------------------------------------------------------------
pop <- n <-  10000
n_years <- 11 #simulate 11 years with the first year as burn-in period 
time_steps <- n_years*12 #each time step is a month 

params <- list()

params$degree_pre <- 3 #average node degree 
params$edges <- params$degree_pre*pop/2 #total number of edges 
params$edge_duration <- 63.55
params$concurrent <- (179- 9-38)/179*pop #Nodes with 2+ ties 
params$nodematch <- params$edges * (0.656) #Assortative mixing by age group 


## monthly injection frequency
ego$inject[ego$inject>=995] <- NA
ego$inj30d <- round(ego$inject/6) #convert 6m freq to 1m freq
ego <- ego[which(ego$inj30d > 0), ]
inj_30d <- spread(ego[,c("id","visquarter", "inj30d")], visquarter, inj30d)
inj_30d <- rowMeans(inj_30d[,18:21], na.rm=T) #monthly inj freq in 2018 
inj_30d <- inj_30d[!is.na(inj_30d)]
inj30d_fit <- fitdistr(inj_30d, "exponential") %>% coefficients()


## Age dist from BESURE
age.group <- sample(1:4, pop, replace=T, prob=c(2,21,16,61))
age <- sapply(age.group, function(x) if(x==1){sample(18:24,1)}
              else if (x==2){sample(25:34,1)}
              else if (x==3){sample(35:44,1)}
              else {sample(45:75,1)})
age_fit <- fitdistr(age, "gamma") %>% coefficients()


#Unborn pop age
age_1st_inj_freq <- as.vector(table(bdf$binjage1))
age_1stinj_fit <- fitdistr(bdf$binjage1[!is.na(bdf$binjage1)], "gamma") %>% coefficients()


#Injection risk trajectories
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

cess_change_probs <- data.frame(early = c(0.607, 0.107, 0.143, 0.143),
                                delay = c(0.902, 0.098, 0, 0),
                                relapse = c(0.54, 0.069, 0.31, 0.08),
                                persistent = c(0.316, 0.175, 0.105, 0.404))

mortality_rates <- c(rep(0,18), rep(25.7, 50-18), rep(36.6, 60-50), rep(54.1, 100-60))
mortality_rates_covid <- c(rep(0,18), rep(26.6, 50-18), rep(46.6, 60-50), rep(48.5, 100-60))

ever_treated_hcv <- 0.311*0.15 + 0.357*0.33 + 0.385*0.52


hcv_sexual <- 1/380000
hiv_sexual <- mean(c(0.04,0.08))/100


