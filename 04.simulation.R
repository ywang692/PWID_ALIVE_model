source("01.param.R")
source("02.net.R")
source("03.API_extension.R")

library(EpiModel)

## Set vectors for service availability
test_rate <- test_rate_covid <- test_rate_12m <- test_rate_24m <- rep(0.055, time_steps)
test_rate_covid[37:(37+2)] <- 0.055/2
test_rate_12m[37:(37+11)] <- 0
test_rate_24m[37:(37+23)] <- 0

art_prescribe <- art_prescribe_covid <-  art_prescribe_12m <- art_prescribe_24m <- rep(0.787, time_steps)
art_prescribe_covid[37:(37+2)] <- 0.787/2
art_prescribe_12m[37:(37+11)] <- 0
art_prescribe_24m[37:(37+23)] <- 0

art_adhere <- art_adhere_covid <- rep(0.85, time_steps)

hcv_treat <- hcv_treat_covid <- hcv_treat_12m <- hcv_treat_24m <- rep(0.1/12, time_steps)
hcv_treat_covid[37:(37+2)] <- (0.1/12) * 0.5
hcv_treat_covid[(37+3):(37+11)] <- (0.1/12) * 0.7
hcv_treat_12m[37:(37+11)] <- 0
hcv_treat_24m[37:(37+23)] <- 0

moud_on <- ssp_on <- moud_on_covid <- ssp_on_covid <- moud_on_12m <- ssp_on_12m <- moud_on_24m <- ssp_on_24m <- rep(1, time_steps)
moud_on_covid[37:(37+2)] <- ssp_on_covid[37:(37+2)] <- 0
moud_on_covid[(37+3):(37+11)] <- ssp_on_covid[(37+3):(37+11)] <- 0.5
moud_on_12m[37:(37+11)] <- ssp_on_12m[37:(37+11)] <- 0
moud_on_24m[37:(37+23)] <- ssp_on_24m[37:(37+23)] <- 0

##Set vector for on/off of pandemic-induced behavior changes
behav_change <- behav_change_covid <- behav_change_4y <- behav_change_6y <- behav_change_8y <- rep(0, time_steps)
behav_change_covid[37:(37+40)] <- 1 
behav_change_4y[37:(37+4*12)] <- 1 
behav_change_6y[37:(37+6*12)] <- 1 
behav_change_8y[37:132] <- 1 


hiv_vec <- hcv_vec <-  rep(0,pop) #set initial infected compartments based on known prevalence 
hiv_vec[sample(1:pop, pop*0.294, replace=F)] <- 1 #HIV mono-infection
hcv_vec[sample(which(hiv_vec==1), pop*0.117, replace = F)] <- 1 #co-infection
hcv_vec[sample(which(hiv_vec==0), pop*(0.409 - 0.117), replace = F)] <- 1 #HCV mono-infection

status_vec <- rep("s", pop)
status_vec[which(hcv_vec==1 | hiv_vec==1)] <- "i"

init <- init.net(status.vector = status_vec)


module.order <- c("resim_nets.FUN",
                  "setup.FUN",
                  "test.FUN",
                  "treatment.FUN",
                  "infection.FUN",
                  "hcvClearance.FUN",
                  "departures.FUN",
                  "aging.FUN",
                  "cess.FUN",
                  "arrivals.FUN",
                  "nwupdate.FUN",
                  "prevalence.FUN")

ncores <- 1
nsims <- 5

control <- control.net(type=NULL, nsims=nsims, nsteps=time_steps,
                       ncores = ncores,
                       setup.FUN = setup,
                       test.FUN = test,
                       treatment.FUN = treatment,
                       infection.FUN = infect,
                       hcvClearance.FUN = hcv.spont.clear,
                       departures.FUN = dfunc,
                       aging.FUN = aging,
                       cess.FUN = cess,
                       arrivals.FUN = bfunc,
                       nwupdate.FUN = nwup,
                       verbose.FUN = vb,
                       resimulate.network=T,
                       module.order = module.order,
                       verbose.int = 20,
                       save.trans=F)

param <- param.net(test_rate = test_rate,
                   art_prescribe = art_prescribe,
                   art_adhere = art_adhere,
                   moud_on = moud_on,
                   ssp_on = ssp_on,
                   hcv_treat = hcv_treat,
                   behav_change = behav_change)

sim <- netsim(est, param, init, control)




