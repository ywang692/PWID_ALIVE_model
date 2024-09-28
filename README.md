# Impact of pandemic-induced service disruptions and behavioral changes on HCV and HIV transmission amongst people who inject drugs: a modeling study

This repository holds the source code to reproduce the simulations in our research article 'Impact of pandemic-induced service disruptions and behavioral changes on HCV and HIV transmission amongst people who inject drugs: a modeling study'. 

The model was written and executed in the R statistical software language and largely utilized the [EpiModel](http://epimodel.org/) package. 

Code to compile all parameters and generate initial conditions for simulations is in `01.param.R`\
Network simulations with ERGM are done in `02.net.R`\
Customized modules of the tranmission model are in `03.API_extension.R`. These codes are for the application programming interface (API) of the EpiModel package. 

For questions, please reach out to awesolowski@jhu.edu or ywang692@jhmi.edu



## Abstract
The COVID-19 pandemic may have disproportionally impacted vulnerable groups such as people who inject drugs (PWID) through reduced healthcare services as well as social changes from pandemic mitigation measures. Understanding how the COVID-19 pandemic and associated mitigation strategies subsequently changed the trajectory of hepatitis C virus (HCV) and HIV transmission is critical to estimating disease burdens, identifying outbreak risk, and developing informed intervention strategies. Using behavioral data from the AIDS Linked to the IntraVenous Experience (ALIVE) study, an ongoing community-based cohort of PWID in Baltimore, USA, and an individual-based network model, we explored the impacts of service disruptions combined with changes in social networks and injecting behaviors of PWID on HCV and HIV transmission. Analyses of ALIVE data showed that during the pandemic, there was an acceleration in injection cessation trajectories overall, but those who continued injecting increased the frequency of injection; at the same time, individual drug use networks became smaller and the probability of injecting with others decreased. Simulation results demonstrated that HCV and HIV prevalence increased from service disruptions alone, but these effects were mitigated when including observed behavior changes in addition. Model results combined with rich individual behavioral data indicated that pandemic-induced behavioral changes of PWID that lasted longer than service disruptions could have offset the increasing disease burden caused by disrupted service access during the pandemic.
 
