
setup <- function(dat,at) {
  
  #Set up initial attributes for each individual
  if (at == 2){
    dat <- set_attr(dat, "active", rep(1,pop))
    dat <- set_attr(dat, "age", t0$age)
    dat <- set_attr(dat, "gender", t0$gender)
    dat <- set_attr(dat, "status", t0$status_vec) #s or i
    dat <- set_attr(dat, "hiv_status", t0$hiv_vec)
    dat <- set_attr(dat, "hcv_status", t0$hcv_vec)
    dat <- set_attr(dat, "coinfection", as.numeric(t0$hiv_vec + t0$hcv_vec == 2))
    dat <- set_attr(dat, "infTime_hiv", t0$hiv_vec) #Set infection to 1 for those who were infected at the beginning of simulations
    dat <- set_attr(dat, "infTime_hcv", t0$hcv_vec)
    dat <- set_attr(dat, "inj_freq", t0$inj30d)
    ## 1 = early cessation; 2 = delayed cessation; 3 = relapse; 4 = persistent 
    dat <- set_attr(dat, "cess_group", t0$cess_group) ## risk trajectory each individual is initially in
    dat <- set_attr(dat, "cess_begin", pmin(t0$cess_begin, time_steps)) #time when one's risk trajectory begins
    dat <- set_attr(dat, "cess_at", rep(0,pop)) #progress of one's trajectory (x axis of probability curve)
    dat <- set_attr(dat, "hcv_spont_clear_prob", t0$hcv_spont_clear_prob) 
    dat <- set_attr(dat, "hcv_acute_length", t0$hcv_acute_length) ## length of HCV acute infection for each individual
    dat <- set_attr(dat, "hcv_ever_treated", t0$hcv_ever_treated) ## whether an individual was ever treated for HCV
    dat <- set_attr(dat, "behav_change", rep(1,pop)) ##whether an individual would be affected by pandemic-induced behavior change
  }
  
  ## at time step 37 (March 2020), a fraction of individuals are assumed to be 'lost to follow-up' and are differently affected by pandemic-induced behavior changes
  if (at == 37){
    alive <- is.na(get_attr(dat,"exitTime")) %>% as.numeric()
    traj <- dat$attr$cess_group
    early_vec <- sample(which(traj == 1 & alive == 1), sum(traj == 1 & alive == 1)*params$ltf[1], replace = F)
    delay_vec <- sample(which(traj == 2 & alive == 1), sum(traj == 2 & alive == 1)*params$ltf[2], replace = F)
    relapse_vec <- sample(which(traj == 3 & alive == 1), sum(traj == 3 & alive == 1)*params$ltf[3], replace = F)
    persistent_vec <- sample(which(traj == 4 & alive == 1), sum(traj == 4 & alive == 1)*params$ltf[4], replace = F)
    lost <- c(early_vec, delay_vec, relapse_vec, persistent_vec)
    dat$attr$behav_change[lost] <- 0
  }
  
  return (dat)
}


test <- function(dat, at) {
  
  active <- get_attr(dat, "active")
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  
  ## HIV testing
  if (at == 2) {
    diag.status <- rep(0, length(active))
    infected <- which(get_attr(dat,"hiv_status") == 1)
    tested <- sample(infected, length(infected)*params$hiv_test_ever, replace=FALSE) #Proportion of individuals ever tested for HIV
    diag.status[tested] <- 1
    nTest <- length(tested)
  } else {
    diag.status <- get_attr(dat, "diag.status")
    hiv_status <- get_attr(dat,"hiv_status")
    infected <- which(hiv_status == 1)
    
    #Eligible individuals are active injectors who are susceptible or HIV+ but not diagnosed yet
    idsElig <- c(which(alive == 1 & hiv_status == 0), which(alive ==1 & hiv_status == 1 & diag.status == 0))
    nElig <- length(idsElig)

    #Draw ids to become tested at this time step 
    vecTest <- which(rbinom(nElig, 1, dat$param$test_rate[at]) == 1)
    idsTest <- idsElig[vecTest]
    diag.status[idsTest] <- 1
    diag.status[which(hiv_status==0)] <- 0 #only HIV+ individuals have their diagnosis status changed to 1
    nTest <- length(intersect(idsTest, infected))
  }
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_epi(dat, "nTest", at, ifelse(at==2, nTest ,nTest)) 

  return(dat)
}


treatment <- function(dat, at) {
  
  active <- get_attr(dat, "active")
  test.status <- get_attr(dat, "diag.status")
  hiv.status <- get_attr(dat,"hiv_status")
  hcv.status <- get_attr(dat,"hcv_status")
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  
  ## HIV treatment (ART)
  if (at == 2) {
    art.status <- rep(0, length(active))
    idsElig <- which(alive==1 & hiv.status ==1 & test.status==1)
    art.status[sample(idsElig, length(idsElig) * dat$param$art_prescribe[at], replace=F)] <- 1
    dat <- set_attr(dat, "art.status", art.status)
    idArt <- which(dat$attr$art.status == 1)
    nArt <- length(idArt)
    dat <- set_attr(dat, "viral.supp", rep(NA,length(active)))
    if (length(idArt) > 0){
      dat$attr$viral.supp[idArt] <- 0
      idSupp <- sample(idArt, length(idArt)*params$achieve_supp)
      dat$attr$viral.supp[idSupp] <- 1
    }
  } else {
    art.status <- get_attr(dat, "art.status")
    ids.art.on <- which(art.status==1)
    #At each time step individuals could lose their viral suppressions status due to ART non-adherence 
    nonadhere <- sample(ids.art.on, length(ids.art.on) * (1 - dat$param$art_adhere[at]))
    dat$attr$art.status[nonadhere] <- 0
    dat$attr$viral.supp[nonadhere] <- 0
    
    #Eligible individuals are those who are HIV+, diagnosed, and not on ART
    idsElig <- which(alive==1 & hiv.status ==1 & test.status==1 & art.status==0)
    nElig <- length(idsElig)

    idsArt <- idsElig[rbinom(nElig, 1, dat$param$art_prescribe[at]) == 1]
    nArt <- length(idsArt)
    art.status[idsArt] <- 1
    dat$attr$art.status <- art.status
    
    #Individuals who are HIV+, tested, treated but not yet suppressed, may become virally suppressed
    vecElig <- alive==1 & hiv.status ==1 & test.status==1 & art.status==1
    vecElig[which(dat$attr$viral.supp == 1)] <- F
    idsElig <- which(vecElig == T)
    viral.supp <- get_attr(dat, "viral.supp")
    viral.supp[idsArt] <- 0 #those who just initiated ART are not suppressed yet
    viral.supp[sample(idsArt, length(idsArt)*params$achieve_supp, replace=FALSE)] <- 1
    dat$attr$viral.supp <- viral.supp
  }
  dat <- set_epi(dat, "nArt", at, ifelse(at==2, nArt, nArt))
  
  
  ## HCV treatment 
  idsElig <- which(alive == 1 & hcv.status == 1)
  nElig <- length(idsElig)
  
  treat_vec <- rbinom(nElig, 1, prob = dat$param$hcv_treat[at])
  treated <- idsElig[as.logical(treat_vec)]
  if (length(treated) > 0){
    success_vec <- rbinom(length(treated), 1, prob=params$hcv_treat_success) #achieve SVR
    nTreated <- sum(success_vec)
    treat_success <- treated[as.logical(success_vec)]
    dat$attr$hcv_status[treat_success] <- 0 #change HCV status when successfully treated 
    dat$attr$coinfection[treat_success] <- 0 #no longer have dual infections
    if (length(treat_success) > 0) {
      hiv_sus <- which(dat$attr$hiv_status == 0)
      idsElig <- treat_success[treat_success %in% hiv_sus]
      dat$attr$status[idsElig] <- "s" #for individuals newly treated for HCV, if they are HIV-, then they change back to 'susceptible' status
    }
  } else {nTreated <- 0}
  dat$attr$hcv_ever_treated[treated] <- 1
  dat <- set_epi(dat, "nTreat_hcv", at, nTreated)
  
  return(dat)
}


infect <- function(dat, at) {
  
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  alive <- is.na(get_attr(dat,"exitTime")) %>% as.numeric()
  
  #Eligible individuals for transmission are active injectors who are infected by HCV and/or HIV
  idsInf <- which(alive == 1 & active == 1 & status == "i")
  nElig <- length(idsInf)
  
  if (nElig > 0) {
    
    nw <- data.frame(V1 = unlist(lapply(dat$nw[[1]]$mel, function(x) x$inl)),
                    V2 = unlist(lapply(dat$nw[[1]]$mel, function(x) x$outl))) #convert overall network to an edgelist 
    #compile the start&end point to active period of all edges
    nw$nrow <- unlist(lapply(dat$nw[[1]]$mel, function(x) nrow(x$atl$active)))
    nw$start <- nw$end <- NA
    for (r in which(nw$nrow > 0)){
      nrow <- nw$nrow[r]
      nw$start[r] <- dat$nw[[1]]$mel[[r]]$atl$active[nrow, 1]
      nw$end[r] <- dat$nw[[1]]$mel[[r]]$atl$active[nrow, 2]
    }
    
    nw <- nw[which(at >= nw$start & at < nw$end), ] #only keep edges that are active at the current step
  
    degree <- nrow(nw)*2/sum(active) #record mean degree of active drug use network 
    if (at == 2){
      dat$epi$mean_degree <- c(0,degree)
    } else {
      dat$epi$mean_degree[at] <- degree
    }
    
    # get info of the infected individuals
    gender <- get_attr(dat, "gender")[idsInf]
    behav <- get_attr(dat, "behav_change")[idsInf]
    ssp_1 <- sapply(dat$attr$cess_group[idsInf], function(x) rbinom(1,size = 1, prob = params$ssp_prob[x]))
    moud_1 <- sapply(dat$attr$cess_group[idsInf], function(x) rbinom(1,size = 1, prob = params$moud_prob[x]))
    ssp <- ssp_1 * rbeta(nElig, params$ssp_effect_dist[1], params$ssp_effect_dist[2]) * dat$param$ssp_on[at]
    moud <- moud_1 * rbeta(nElig, params$moud_effect_dist[1], params$moud_effect_dist[2]) * dat$param$moud_on[at]
    inj_freq <- get_attr(dat,"inj_freq")[idsInf]
    inj_freq <- inj_freq * (1-moud)
    if (dat$param$behav_change[at] == 1){
      if (dat$param$lost_to_followup == "LOW"){
        inj_freq <- inj_freq * ifelse(behav == 1, 1.44, 1)
      } else if (dat$param$lost_to_followup == "HIGH"){
        inj_freq <- inj_freq * ifelse(behav == 1, 1.44, 1.76)
      } else if (dat$param$lost_to_followup == "MODERATE") {
        inj_freq <- inj_freq * 1.44
      }
    }
    hiv_status <- get_attr(dat, "hiv_status")[idsInf]
    hcv_status <- get_attr(dat, "hcv_status")[idsInf]
    
    if (at == 2){
      acute_hiv <- rep(0, nElig)
    } else {
      acute_hiv <- ifelse(get_attr(dat,"infTime_hiv")[idsInf] == at-2, 1, 0) #determine whether in acute stage of HIV infection
    }
    acute_hiv_increase <- ifelse(acute_hiv == 1, params$hiv_acute, 1) #fold increase for acute HIV transmission
    diagnosed <- get_attr(dat,"diag.status")[idsInf]
    art <- get_attr(dat,"art.status")[idsInf]
    supp <- get_attr(dat,"viral.supp")[idsInf]
    viral_supp_scale <- ifelse(diagnosed + art == 2, params$art_effect_f, 0)
    viral_supp_scale[which(diagnosed + art + supp == 3)] <- params$art_effect_s #different scale of transmission reduction based on HIV viral suppresion
    viral_supp_scale[which(hcv_status == 1)] <- viral_supp_scale[which(hcv_status == 1)] * (1-params$art_effect_given_hcv)
    
    prob_syringe_share <- pmax((1 - ssp)* params$syringe_share, 0)
    if (dat$param$behav_change[at] == 1){
      if (dat$param$lost_to_followup == "LOW"){
        prob_syringe_share <- prob_syringe_share * 1/3
      } else if (dat$param$lost_to_followup == "HIGH"){
        prob_syringe_share <- prob_syringe_share * ifelse(behav == 1, 1/3, 5/3)
      } else if (dat$param$lost_to_followup == "MODERATE") {
        prob_syringe_share <- prob_syringe_share * ifelse(behav == 1, 1/3, 1)
      }
    }
    prob_hiv_per_share <- pmax(params$hiv_blood * acute_hiv_increase * (1 - viral_supp_scale), 0)
    prob_hiv_per_share[which(hiv_status == 0)] <- 0 #HIV- individuals have 0 transmission probability 
    prob_hiv_per_sex <- pmax(params$hiv_sexual * acute_hiv_increase *(1 - viral_supp_scale), 0)
    prob_hiv_per_sex[which(hiv_status == 0)] <- 0
    
    if (at == 2){
      acute_hcv <- rep(0, nElig)
    } else {
      acute_hcv <- ifelse(get_attr(dat,"hcv_acute_length")[idsInf] - (at - get_attr(dat,"infTime_hcv")[idsInf]) >= 0, 1, 0)
    } #determine whether in acute stage of HCV infection
    acute_hcv_increase <- ifelse(acute_hcv == 1, params$hcv_acute, 1)
    post_treatment_reduc <- rbeta(nElig, params$hcv_post_treat_dist[1], params$hcv_post_treat_dist[2])
    post_treatment_reduc[get_attr(dat,"hcv_ever_treated")[idsInf] != 1] <- 0 #reduced HCV transmission prob if ever treated
    prob_hcv_per_share <- pmax(params$hcv_blood * acute_hcv_increase * (1-post_treatment_reduc), 0)
    prob_hcv_per_share[which(hcv_status == 0)] <- 0 #HCV- individuals have 0 transmission probability
    prob_hcv_per_sex <- pmax(params$hcv_sexual * acute_hcv_increase * (1-post_treatment_reduc), 0)
    prob_hcv_per_sex[which(hcv_status == 0)] <- 0
  
    ## Start looping transmission simulation for each infected individual 
    for (ii in 1:nElig){
      index <- idsInf[ii]
      
      ## find all of index's drug-use partners who are also active at current time step 
      el <- unique(c(nw[which(nw$V1 == index), "V2"], nw[which(nw$V2 == index), "V1"]))
      el <- intersect(el, which(active == 1))

      if (length(el)>0){
        
        if (is.na(behav[ii])) browser()

        if (dat$param$behav_change[at]==1 & behav[ii] == 1){
          el <- el[as.logical(rbinom(length(el), 1, prob = 1/3))] ##during covid people cut 2/3 of their drug contacts
          if (length(el>0)) {
              inj_partners <- el[as.logical(rbinom(length(el), 1, prob = params$prob_inj*0.726))] ## during covid people are less likely to inject with drug contacts
            } else {
              inj_partners <- c()
            } 
        } else {
          inj_partners <- el[as.logical(rbinom(length(el), 1, prob = params$prob_inj))]
        } ## prob of drug-use partners inject with ego 

        if (length(inj_partners) > 0){
          n_inj_partners <- length(inj_partners)
          distribution <- get_attr(dat,"inj_freq")[inj_partners]
          inj_acts_per_partner <- sample(1:n_inj_partners, size=round(inj_freq[ii]*params$inj_share), replace=T, prob=distribution) 
            ## shared injections with alters are distributed based on alters' injection frequencies
          inj_acts_per_partner <- sapply(1:n_inj_partners, function(x) sum(inj_acts_per_partner == x))
          share_syringe_events_per_partner <- sapply(inj_acts_per_partner, function(x) sum(rbinom(x, size = 1, prob = prob_syringe_share[ii])))
            ## HIV+ individuals are more likely to get HCV
          prob_hcv <- rep(prob_hcv_per_share[ii], n_inj_partners)
          elig_partner <- sapply(inj_partners, function(x) ifelse(dat$attr$hiv_status[x]==1, 1, 0))
          prob_hcv[as.logical(elig_partner)] <- pmin(prob_hcv_per_share[ii] * params$hcv_given_hiv, 1) %>% pmin(.,1)
          new_hiv <- sapply(1:n_inj_partners, function(x) rbinom(1, size = share_syringe_events_per_partner[x], prob = prob_hiv_per_share[ii]))
          new_hcv <- sapply(1:n_inj_partners, function(x) rbinom(1, size = share_syringe_events_per_partner[x], prob = prob_hcv[x]))

          idsNewInf_hiv <- inj_partners[which(new_hiv > 0)]
          idsNewInf_hcv <- inj_partners[which(new_hcv > 0)]

          ## Sexual sub-network 
          sexual_prob <- ifelse(gender[ii]==0, params$sexual_male, params$sexual_female)
          if (dat$param$behav_change[at]==1 & behav[ii] == 1){
            sexual_prob <- sexual_prob * ifelse(gender[ii] == 0, 1.66, 1.83) ##PWID increased probability of IDU sexual tie during pandemic
          }
          sex_partners <- el[as.logical(rbinom(n_inj_partners, 1, prob=sexual_prob))] #prob of also being sexual contact given injection partnership
          n_sex_partners <- length(sex_partners)
          if (n_sex_partners > 0){
            sex_freq <- rpois(n_sex_partners, params$sex30d_dist) #frequency of unprotected sex
            new_hiv <-  sapply(1:n_sex_partners, function(x) rbinom(1, size = sex_freq[x], prob = prob_hiv_per_sex[ii]))
            idsNewInf_hiv <- unique(append(idsNewInf_hiv, sex_partners[which(new_hiv>0)]))
            prob_hcv <- rep(prob_hcv_per_sex[ii], n_sex_partners)
            elig_partner <- sapply(sex_partners, function(x) ifelse(dat$attr$hiv_status[x]==1 & dat$attr$hcv_status[x]==0, 1, 0))
            prob_hcv[as.logical(elig_partner)] <- pmin(prob_hcv_per_sex[ii] * params$hcv_given_hiv, 1)
            new_hcv <- rbinom(n_sex_partners, 1, prob = prob_hcv)
            idsNewInf_hcv <- unique(append(idsNewInf_hcv, sex_partners[as.logical(new_hcv)]))
            }

          dat$attr$infTime_hiv[idsNewInf_hiv[dat$attr$hiv_status[idsNewInf_hiv] == 0]] <- at #set time step of incident infection
          dat$attr$infTime_hcv[idsNewInf_hcv[dat$attr$hcv_status[idsNewInf_hcv] == 0]] <- at
          dat$attr$hiv_status[idsNewInf_hiv] <- 1
          dat$attr$hcv_status[idsNewInf_hcv] <- 1
          idsNewInf <- unique(append(idsNewInf_hiv, idsNewInf_hcv))
          dat$attr$status[idsNewInf] <- "i"

        }
      }
    }
  }
  coinf <- as.numeric(dat$attr$hiv_status + dat$attr$hcv_status == 2) # co-infected population
  coinf.flow <- length(which(coinf - dat$attr$coinfection == 1))
  dat$attr$coinfection <- coinf
  nInf_hiv <- length(which(dat$attr$hiv_status ==1 & dat$attr$infTime_hiv == at))
  nInf_hcv <- length(which(dat$attr$hcv_status ==1 & dat$attr$infTime_hcv == at))
  
  ## Set summary stats 
  if (at == 2) {
    dat$epi$hiv.si.flow <- c(0, nInf_hiv)
    dat$epi$hiv.i.num <- c(0, sum(active == 1 & dat$attr$hiv_status ==1))
    dat$epi$hcv.si.flow <- c(0, nInf_hcv)
    dat$epi$hcv.i.num <- c(0, sum(active == 1 & dat$attr$hcv_status ==1))
    dat$epi$s.num <- c(0, sum(active == 1 & dat$attr$status == "s"))
    dat$epi$coinf.flow <- c(0, coinf.flow)
    dat$epi$coinf.num <- c(0, sum(active == 1 & dat$attr$coinfection ==1))
  }
  else {
    dat$epi$hiv.si.flow[at] <- nInf_hiv
    dat$epi$hiv.i.num[at] <- sum(active == 1 & dat$attr$hiv_status ==1)
    dat$epi$hcv.si.flow[at] <- nInf_hcv
    dat$epi$hcv.i.num[at] <- sum(active == 1 & dat$attr$hcv_status ==1)
    dat$epi$s.num[at] <- sum(active == 1 & dat$attr$status == "s")
    dat$epi$coinf.flow[at] <- coinf.flow
    dat$epi$coinf.num[at] <- sum(active == 1 & dat$attr$coinfection ==1)
  }
  return(dat)
}


hcv.spont.clear <- function(dat, at){
  
  active <- get_attr(dat, "active")
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  hcv_status <- get_attr(dat, "hcv_status")
  hiv_status <- get_attr(dat, "hiv_status")
  infTime <- get_attr(dat, "infTime_hcv")
  acute_length <- get_attr(dat, "hcv_acute_length")
  prob_clear <- get_attr(dat, "hcv_spont_clear_prob") / acute_length #monthly probability of spontaneous clearance
  
  idsElig <- which(alive == 1 & hcv_status ==1 & infTime > 2) 
  
  if (at > 2){
    months_passed <- at - infTime[idsElig]
    idsElig <- idsElig[which(acute_length[idsElig] - months_passed >= 0)] #whether individuals are still in acute stage of HCV infection
    nElig <- length(idsElig)
    dat <- set_epi(dat, "n.hcv.clear", at, 0)
    if (length(idsElig)>0){
      hiv_scale <- sapply(hiv_status[idsElig], function(x) ifelse(x==1, params$spont_clear_given_hiv, 1)) #HIV+ indiv less likely to clear HCV infection
      new_hcv_vec <- rbinom(nElig, 1, prob = (1 - prob_clear[idsElig] * hiv_scale))
      dat$attr$hcv_status[idsElig] <- new_hcv_vec #Update individuals with chronic HCV
      new_sus <- idsElig[which(new_hcv_vec == 0)] 
      if (length(new_sus) > 0) {
        dat$attr$coinfection[new_sus] <- 0
        hiv_sus <- which(dat$attr$hiv_status == 0)
        idsElig <- new_sus[new_sus %in% hiv_sus]
        dat$attr$status[idsElig] <- "s" #If also HIV-, then go back to 'susceptible' compartment 
        dat <- set_epi(dat, "n.hcv.clear", at, length(new_sus))
      }
    }
  }
  return(dat)
}


aging <- function(dat,at){
  
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  dat$attr$age[as.logical(alive)] <- dat$attr$age[as.logical(alive)] + 1/12 #all alive individuals age
  
  #record mean age of active drug users at the current step 
  active <- get_attr(dat, "active")
  if (at == 2){
    dat$epi$meanAge <- c(NA_real_, mean(dat$attr$age[active == 1], na.rm=T))
  } else {
    dat$epi$meanAge[at] <- mean(dat$attr$age[active == 1], na.rm=T)
  }
  return(dat)
}


dfunc <- function(dat,at){
  
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  idsElig <- which(alive ==1)
  nElig <- length(idsElig)
  nDeaths <- 0
  if (nElig > 0) {
    ages <- dat$attr$age[idsElig]
    hiv <- as.logical(dat$attr$hiv_status[idsElig] == 1)
    hcv <- as.logical(dat$attr$hcv_status[idsElig] == 1)
    no_supp <- which(dat$attr$viral.supp == 0)
    if (dat$param$moud_on[at] == 1){
      death.rates <- sapply(ages, function(x) params$mortality[x]/1000/12)
    } else { #apply different death rates for pre-/inter-pandemic 
      death.rates <- sapply(ages, function(x) params$mortality_covid[x]/1000/12)
    }
    death.rates <- death.rates / (params$hiv_prev * params$mort_hiv + (1 - params$hiv_prev)) #first set everyone to have HIV- mortality
    death.rates[hiv] <- death.rates[hiv] * params$mort_hiv #then multiply mortality of HIV+ indiv
    death.rates <- death.rates / (params$hcv_prev * params$mort_hcv + (1 - params$hcv_prev)) #first set everyone to have HCV- mortality
    death.rates[hcv] <- death.rates[hcv] * params$mort_hcv #then multiply mortality of HCV+ indiv
    
    #Draw ids who passed at current time step
    vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
    idsDeaths <- idsElig[unique(c(vecDeaths, which(ages > 90)))] #individuals older than 90 exit automatically
    nDeaths <- length(idsDeaths)
    if (nDeaths > 0){
      dat$attr$active[idsDeaths] <- 0
      dat$attr$exitTime[idsDeaths] <- at #set exit time 
      deactivate.vertices(dat$nw[[1]], v = idsDeaths, onset = at, terminus = Inf, deactivate.edges = TRUE) #permanently deactivate dead nodes and associated edges
    }
  }
  if (at == 2) {
    dat$epi$d.flow <- c(0, nDeaths)
  } else {
    dat$epi$d.flow[at] <- nDeaths
  }
  return(dat)
}


cess <- function(dat,at){
  
  active <- which(dat$attr$active ==1)
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  change_indiv <- which(alive==1 & dat$attr$behav_change == 1)
  lost_indiv <- which(alive == 1 & dat$attr$behav_change == 0)
  
  ## At the beginning of pandemic, some individuals would change their risk trajectories 
  if (dat$param$behav_change[at]-dat$param$behav_change[at-1] == 1){
    dat$attr$cess_group[change_indiv] <- dat$attr$cess_group[change_indiv] %>%
      sapply(function(x) sample(1:4, 1, prob = params$cess_change_probs[ ,x], replace=T))
    if (dat$param$lost_to_followup == "LOW"){
      dat$attr$cess_group[lost_indiv] <- dat$attr$cess_group[lost_indiv] %>%
        sapply(function(x) sample(1:4, 1, prob = params$cess_change_probs[ ,x], replace=T))
    } else if (dat$param$lost_to_followup == "HIGH"){
      dat$attr$cess_group[lost_indiv] <- 4
    } 
  } 
  
  ## Eligible individuals are those who are alive and whose risk trajectories have begun before the current time step
  idsElig <- which(dat$attr$cess_begin < at & alive == 1)
  nElig <- length(idsElig)
  
  vecCess <- c()
  if (nElig > 0){
    #reset all eligible indiv as active injectors at each time step 
    dat$attr$active[idsElig] <- rep(1, nElig)
    
    inj_prob <- rep(NA,nElig)
    dat$attr$cess_at[idsElig] <- dat$attr$cess_at[idsElig] + 1 #add 1 time step on the progress of one's risk trajectory 
    
    inj_prob <- sapply(idsElig, function(x) params$cess_prob_values[dat$attr$cess_at[x] , dat$attr$cess_group[x]])
    vecCess <- which(rbinom(nElig, 1, prob=inj_prob) == 0) #draw individuals who are in cessation at current time step 
    
    if (length(vecCess) > 0){
      idsCess <- idsElig[vecCess]
      dat$attr$active[idsCess] <- 0
      #deactivate nodes in cessation (for one time step)
      dat$nw[[1]] <- deactivate.vertices(dat$nw[[1]], onset = at, terminus = at+1, v=idsCess, deactivate.edges=TRUE)
    }
  }
  ## Update summary stats (total individuals in cessation at current time step)
  if (at == 2){
    dat$epi$tot_cess <- c(0, length(vecCess))
  } else {
    dat$epi$tot_cess[at] <- length(vecCess)
  }
  return(dat)
}


bfunc <- function(dat, at){
  
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  active <- get_attr(dat, "active")
  active_pop <- sum(active==1 & alive==1)
  
  numNeeded <- pop - active_pop #add new injectors to keep population constant 
  
  if (numNeeded > 0){
    nBirths <- rpois(1, numNeeded)
  } else {
    nBirths <- 0
  }
  if (nBirths > 0){
    #assign new ids to new individuals 
    last_pid <- dat$`_last_unique_id`
    new_pid <- seq(1,nBirths) + last_pid
    dat$nw[[1]] <- add.vertices.networkDynamic(dat$nw[[1]], nv=nBirths, vertex.pid = as.character(new_pid))
    dat[["_last_unique_id"]] <- tail(new_pid, 1)
    newNodes <- (n+1) : (n+nBirths)
    dat$nw[[1]] <- activate.vertices(dat$nw[[1]], onset = at, terminus = Inf, v = newNodes) #activate new nodes
    dat$attr$unique_id <- c(dat$attr$unique_id, new_pid)
    #assign individual attributes to new nodes 
    dat$attr$active <- c(dat$attr$active, rep(1,nBirths))
      hiv_status <- rbinom(nBirths, 1, prob=params$hiv_prev_general)
      hcv_status <- rbinom(nBirths, 1, prob=params$hcv_prev_general)
    dat$attr$hiv_status <- c(dat$attr$hiv_status, hiv_status)
    dat$attr$hcv_status <- c(dat$attr$hcv_status, hcv_status)
    dat$attr$coinfection <- c(dat$attr$coinfection, as.numeric(hiv_status + hcv_status == 2))
    dat$attr$infTime_hiv <- c(dat$attr$infTime_hiv, rep(0, nBirths))
    dat$attr$infTime_hcv <- c(dat$attr$infTime_hcv, rep(0, nBirths))
      status <- rep("s", nBirths)
      status[which(hiv_status == 1 | hcv_status == 1)] <- "i"
    dat$attr$status <- c(dat$attr$status, status)
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
      unborn_age <- rgamma(nBirths, shape = params$age_1st_inj_dist[1], rate = params$age_1st_inj_dist[2]) %>%
        pmin(.,80) %>%
        pmax(.,18)
    dat$attr$age <- c(dat$attr$age, unborn_age)
    dat$attr$gender <- c(dat$attr$gender, rbinom(nBirths, 1, params$m0f1))
      new_freq <- rexp(nBirths, rate = params$inj30d_dist) %>% 
        round() %>% 
        pmax(.,1) %>% 
        pmin(.,200)
    dat$attr$inj_freq <- c(dat$attr$inj_freq, new_freq)
    dat$attr$cess_group <- c(dat$attr$cess_group, sample(1:4, nBirths, replace=T, prob=params$cess_dist))
      cess_begin <- runif(nBirths, min=params$cess_begin_dist[1], max=params$cess_begin_dist[2]) %>% round()
      cess_begin <- cess_begin + at 
      cess_begin <- pmin(cess_begin, time_steps)
    dat$attr$cess_at <- c(dat$attr$cess_at, rep(0, nBirths))
    dat$attr$cess_begin <- c(dat$attr$cess_begin, cess_begin)
    dat$attr$diag.status <- c(dat$attr$diag.status, rep(0, nBirths))
    dat$attr$art.status <- c(dat$attr$art.status, rep(0, nBirths))
    dat$attr$viral.supp <- c(dat$attr$viral.supp, rep(NA, nBirths))
    new.age.group <- cut(dat$attr$age, 
                         breaks = c(18, 24, 34, 44, 100), 
                         labels = c(1:4), 
                         right = TRUE, include.lowest = TRUE) %>% as.vector()
    dat$nw[[1]] <- network::set.vertex.attribute(dat$nw[[1]], "age.group", new.age.group)
    dat$attr$hcv_spont_clear_prob <- c(dat$attr$hcv_spont_clear_prob, rbeta(nBirths, params$hcv_spont_clear_dist[1], params$hcv_spont_clear_dist[2]))
    dat$attr$hcv_acute_length <- c(dat$attr$hcv_acute_length, runif(nBirths, params$hcv_acute_length_dist[1], params$hcv_acute_length_dist[2]))
    dat$attr$hcv_ever_treated <- c(dat$attr$hcv_ever_treated, rbinom(nBirths, 1, params$hcv_ever_treated))
    dat$attr$behav_change <- c(dat$attr$behav_change, rep(1, nBirths))
    
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
  }
  ## Set summary stats (total births, total active)
  if (at == 2){
    dat$epi$b.flow <- c(0,nBirths)
  } else {
    dat$epi$b.flow[at] <- nBirths
  }
  
  nActive <- length(which(dat$attr$active == 1))
  if (at == 2){
    dat$epi$nActive <- c(0,nActive)
  } else {
    dat$epi$nActive[at] <- nActive
  }
  
  alive <- is.na(get_attr(dat,"exitTime")) %>% as.numeric()
  
  return(dat)
}


vb <- function (x, type, s = 1, at = 2) {
  
  if (type == "startup") {
    if (x$verbose == TRUE) {
      cat("\nStarting Network Simulation...")
    }
  }
  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (x$control$verbose.int == 0 && at == x$control$nsteps) {
        cat("\nSim = ", s, "/", x$control$nsims, sep = "")
      }
      if (x$control$verbose.int > 0 && (at%%x$control$verbose.int == 
                                        0)) {
        cat("\f")
        cat("\nEpidemic Simulation")
        cat("\n----------------------------")
        cat("\nSimulation: ", s, "/", x$control$nsims, 
            sep = "")
        cat("\nTimestep: ", at, "/", x$control$nsteps, 
            sep = "")
        active <- x$attr$active
        cat("\nPopulation Size:", sum(active == 1))
        cat("\nHCV Prevalence:", round(x$epi$hcv.i.num[at]/sum(active == 1), 3))
        cat("\nHIV Prevalence:", round(x$epi$hiv.i.num[at]/sum(active == 1), 3))
        cat("\n----------------------------")
      }
    }
  }
}


