setwd("for_repo")
source("01.param.R")
source("02.net.R")
source("03.API_extension.R")

init <- init.net(status.vector = t0$status_vec)

module.order <- c("resim_nets.FUN",
                  "setup.FUN",
                  "test.FUN",
                  "treatment.FUN",
                  "infection.FUN",
                  "hcvClearance.FUN",
                  "aging.FUN",
                  "departures.FUN",
                  "cess.FUN",
                  "arrivals.FUN",
                  "prevalence.FUN")

control <- control.net(type=NULL, 
                       nsims= 1, 
                       ncores = 1,
                       nsteps=time_steps,
                       setup.FUN = setup,
                       test.FUN = test,
                       treatment.FUN = treatment,
                       infection.FUN = infect,
                       hcvClearance.FUN = hcv.spont.clear,
                       aging.FUN = aging,
                       departures.FUN = dfunc,
                       cess.FUN = cess,
                       arrivals.FUN = bfunc,
                       verbose.FUN = vb,
                       resimulate.network=T,
                       module.order = module.order,
                       verbose.int = 1,
                       save.trans=F)

param <- param.net(test_rate = test_rate, #pandemic-induced service disruptions 
                   art_prescribe = art_prescribe,
                   art_adhere = art_adhere,
                   moud_on = moud_on,
                   ssp_on = ssp_on,
                   hcv_treat = hcv_treat,
                   behav_change = behav_change, #pandemic-induced behavioral changes
                   lost_to_followup = "HIGH")
#Risk level for individuals lost to follow-ups during the pandemic: 
#"HIGH", "MODERATE", "LOW" corresponds to "REDUCED RISK" "CONSISTENT RISK" and "INCREASED RISK" in the paper

sim <- netsim(est, param, init, control)


