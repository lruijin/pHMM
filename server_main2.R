library(foreach)
library(doParallel)
#library(doMC)
source('MHMM_v3.R')
source('evaluation_measure.R')
cl <- registerDoParallel(30)
nrun=500

server.sim <- function(no.run, seed = 2256){
  set.seed(seed + no.run);
  daf <- genSim(N = 354, K=3, ni=6, hs=3, tau=matrix(c(1.2,1.2,1.1,1.2,1.1,1.2,1.1,1.1,1.2),3,3),
                 residinvvar=c(10,10,10),reSigma=matrix(c(0.1,0.06,0.04,0.06,0.1,0.04,0.04,0.04,0.1),nr=3),
                 P0=matrix(c(0.1,0.2,0.1,0.85,0.6,0.1,0.05,0.2,0.8),nr=3),Pi0=c(0.6,0.25,0.15),
                 P1=matrix(c(0.85,0.6,0.1,0.1,0.2,0.1,0.05,0.2,0.8),nr=3),Pi1=c(0.6,0.25,0.15),
                 rx.=rep(c(0,1),each = 354/2),fitRx=c(FALSE,TRUE),
                 w=list(matrix(c(1,0,0),nr=1),matrix(c(0,1,0),nr=1),matrix(c(0,0,1),nr=1)));
  x <- simulated(y = daf$y, inits.=c(daf$pars), nsim.= 4000, ksamp.= 1,
                 N.=daf$N, K = 3, ni.=rep(6,daf$N), rx. = daf$rx, fitRx = daf$fitRx, report1.=100,id=rep(1:daf$N,6),
                 hyperpar=c(4,1,1,1,.001,.0002), run. = no.run)
  z.pred <- apply(x$Zmax[-(1:2000),],2,function(x) names(which.max(table(x))))
  z.true <- daf$Z 
  beta.sim <- x$betaa[-(1:2000),1:length(daf$pars)]
  res <- list(beta.m = apply(beta.sim,2,mean),
       beta.sd = apply(beta.sim,2,sd),
       beta.cover = (daf$pars > apply(beta.sim,2,quantile,0.025)) & (daf$pars < apply(beta.sim,2,quantile,0.975)),
       precision = pre(z.pred,z.true), recall = rec(z.pred, z.true), Fscore = Fscore(z.pred, z.true), accuracy = accuracy(z.pred, z.true),
       err = errRate(z.pred, z.true), precision_m = pre(z.pred, z.true, "micro"))
  return(res)
}

res <- foreach(no.run = 1:nrun,.packages = "MCMCpack",.combine = 'rbind',.multicombine = T, .errorhandling = 'pass') %dopar% {
  source("MHMM_v3.R",local=T);
  source('evaluation_measure.R', local = T);
  server.sim(no.run)
}
save("res", file = paste0('../result/mhmm2_',nrun,'.Rdata'))

