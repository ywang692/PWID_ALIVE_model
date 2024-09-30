library(EpiModel)

nw <- network::network.initialize(n=pop,directed=FALSE)
nw <- network::set.vertex.attribute(nw, "age.group", as.vector(t0$age_group))

formation <- ~edges + nodematch("age.group") + concurrent
target.stats <- c(params$edges, params$nodematch, params$concurrent)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = params$edge_duration)

est <- netest(nw, formation, target.stats, coef.diss, set.control.ergm = control.ergm(MCMLE.maxit = 100))



