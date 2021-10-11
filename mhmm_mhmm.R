# generate data
params <- list(N = 390, K = 4, ni = 5, hs = 3,
               tau = matrix(c(34.102,-0.007,-2.636,3.791,-0.597,-0.04,
                              37.467,-0.54,-3.63,3.299,-0.05,-0.03),nr=3),
               residinvvar = c(.154,60.897,.207,52.557),
               reSigma = matrix(c(13.062,.032,9.663,-0.023,
                                  .032,0.026,.019,.01,
                                  9.663,.019,13.259,.098,
                                  -.023,.01,098,.029),4,4),
               P0 = matrix(c(.372,.017,.006,.551,.713,.017,.077,.27,.913),nr=3),
               pi0 = c(0.139,0.815,0.046),
               P1 = matrix(c(.185,.014,.006,.712,.747,.012,.103,.239,.982),nr=3),
               pi1 = c(0.139,0.815,0.046),
               w=list(matrix(c(1,0,0,0),nr=1),matrix(c(0,1,0,0),nr=1),matrix(c(0,0,1,0),nr=1),matrix(c(0,0,0,1),nr=1)))
sim_each <- function(no.run, seed = 2256, params = params,nsim.=8000,burnin=6000){
  source("MHMM_v3.R");
  source("evaluation_measure.R")
  set.seed(seed + no.run)
  daf = genSim_mhmm(N=params$N,K=params$K,ni=params$ni,hs=params$hs,tau=params$tau,
                    residinvvar=params$residinvvar,reSigma=params$reSigma,w=params$w,
                    P0=params$P0,Pi0=params$pi0,P1=params$P1,Pi1=params$pi1,fitRx=c(F,T))
  basepars = length(daf$pars)
  #Run Preliminary MCMC
  x <- simulated_mhmm(y = daf$y, inits.=c(daf$pars), hs = daf$hs, nsim.= nsim., ksamp.= 1, burnin = burnin,
                      N.=daf$N, K = params$K, ni.=daf$ni, rx. = daf$rx, fitRx = daf$fitRx, report1. = 1000,id=rep(1:daf$N,5),
                      run. = no.run,hyperpar=c(5,1,1,1,1,.001,.0002))
  z.pred <- apply(x$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  z.true <- daf$Z
  beta.sim <- x$betaa[-(1:burnin),1:length(daf$pars)]
  res <- list(beta.mu = apply(beta.sim,2,mean),
              beta.m = apply(beta.sim,2,median),
              beta.sd = apply(beta.sim,2,sd),
              beta.q1 = apply(beta.sim,2,qunatile,0.025),
              beta.q3 = apply(beta.sim,2,qunatile,0.975),
              beta.true = daf$pars, z.pred = z.pred, z.true=z.true,
              precision = pre(z.pred,z.true), recall = rec(z.pred, z.true), Fscore = Fscore(z.pred, z.true), accuracy = accuracy(z.pred, z.true),
              err = errRate(z.pred, z.true), precision_m = pre(z.pred, z.true, "micro"),WAIC = x$WAIC)
  return(res)
}

library(foreach)
library(doMPI)
nrun=1000
cl <- doMPI::startMPIcluster()
registerDoMPI(cl)
out <- foreach(i=1:nrun, .combine = "rbind", .multicombine = T,
               errorhandling="pass")%dopar%{
                sim_each(i) 
               }
save('out',file="mhmm_f.Rdata")
stopCluster(cl)
mpi.exit()
