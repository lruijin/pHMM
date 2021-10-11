library(foreach)
library(doParallel)
source('MHMM_v3.R')
source('evaluation_measure.R')
cl <- registerDoParallel(cores=30)
nrun=500

server.sim <- function(no.run, seed = 22562){
  set.seed(seed + no.run);
  daf <- genSim();
  x <- simulated(y = daf$y, inits.=c(daf$pars), nsim.= 4000, ksamp.= 1,
                 N.=daf$N, ni.=rep(6,daf$N), rx. = daf$rx, fitRx = daf$fitRx, report1.=1000,id=rep(1:daf$N,6),run. = no.run)
  
  beta.sim <- x$betaa[-(1:2000),1:length(daf$pars)]
  z.pred <- apply(x$Zmax[-(1:2000),],2,function(x) names(which.max(table(x))))
  z.true <- daf$Z
  res <- list(beta.m = apply(beta.sim,2,mean),
       beta.sd = apply(beta.sim,2,sd),
       beta.cover = (daf$pars > apply(beta.sim,2,quantile,0.025)) & (daf$pars < apply(beta.sim,2,quantile,0.975)),
       precision = pre(z.pred,z.true), recall = rec(z.pred,z.true), Fscore = Fscore(z.pred,z.true),
       accuracy = accuracy(z.pred,z.true), err = errRate(z.pred, z.true), precision_m = pre(z.pred,z.true,"micro"))
  return(res)
}

res <- foreach(no.run = 1:nrun,.packages = "MCMCpack",.combine = 'rbind',.multicombine = T,.errorhandling = 'pass') %dopar% {
  source("MHMM_v3.R",local=T);
  source('evaluation_measure.R',local = T);  
  server.sim(no.run)
}
save("res", file = paste0('../result/mhmm1_',nrun,'.Rdata'))
