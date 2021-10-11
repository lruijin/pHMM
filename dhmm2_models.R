library(doMPI)
library(foreach)

nrun=1000
sim_each <- function(no.run,nsim.= 5000,burnin=3000){
  source("MHMM_v3.R",local=T)
  source("DHMM_v1.R",local=T)
  # source the functions to get the initial values and evaluation methods
  source("utils.R",local = T)
  set.seed(2257 + no.run)
  daf<- genSim(N = 390, 
               rx. = rep(c(0,1),each = 390/2),fitRx = c(FALSE,TRUE,TRUE))
  x_dhmm2 <- simulated(y = daf$y, inits. = c(daf$pars), nsim.= nsim., burnin = burnin, ksamp.= 1, Km = c(2,2), hs = c(3,3,3),
                      N.=daf$N, ni.=daf$ni, rx. = daf$rx, fitRx = daf$fitRx, report1.=1000,id=rep(1:daf$N,5),run. = no.run)
  inits.pc <- get_inits_pc(daf$pars,x_dhmm2$pidx)
  x_mhmm_pc <- simulated_mhmm(y = daf$y, inits.= inits.pc, nsim.= nsim., ksamp.= 1,burnin = burnin,hs=3,
                              N.=daf$N, K = 4, ni.=rep(5,daf$N), rx. = daf$rx, fitRx = daf$fitRx, report1.=1000,id=rep(1:daf$N,5),
                              hyperpar=c(5,1,1,1,1,.001,.0002), run. = no.run)
  
  x_mhmm_c <- simulated_mhmm(y = list(daf$y[[1]],daf$y[[2]]), inits.= get_inits_c(daf$pars,x_dhmm2$pidx), nsim.= nsim., ksamp.= 1,
                             burnin = burnin,hs=3,N.=daf$N, K = 2, ni.=rep(5,daf$N), rx. = daf$rx, fitRx = c(F,T), report1.=1000,
                             id=rep(1:daf$N,5), hyperpar=c(3,1,1,.001,.0002), run. = no.run)
  
  x_mhmm_p <- simulated_mhmm(y = list(daf$y[[3]],daf$y[[4]]), inits.= get_inits_p(daf$pars,x_dhmm2$pidx), nsim.= nsim., ksamp.= 1,
                             burnin = burnin,hs=3,N.=daf$N, K = 2, ni.=rep(5,daf$N), rx. = daf$rx, fitRx = c(F,T), report1.=1000,
                             id=rep(1:daf$N,5), hyperpar=c(3,1,1,.001,.0002), run. = no.run)
  
  z.pred <- apply(x_dhmm2$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  z.true <- daf$Z
  confusion <- get_confusion(z.pred,z.true)
  for(m in 1:2){
    zm.true <- daf$Zm[[m]]
    zm.pred <- apply(x_dhmm2$Zmmax[[m]][-(1:burnin),],2,function(x) names(which.max(table(x))))
    confusion <- rbind(confusion,get_confusion(zm.pred,zm.true))
  }
  z.pc <- apply(x_mhmm_pc$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  confusion <- rbind(confusion,get_confusion(z.pc,z.true))
  
  z.c <- apply(x_mhmm_c$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  z.p <- apply(x_mhmm_p$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  Pm.est <- c(pm_est(z.pc,z.c),pm_est(z.pc,z.p))
  
  x_mhmm_c_fix <- simulated_mhmm(y = list(daf$y[[1]],daf$y[[2]]), inits.= get_fix_inits(x_mhmm_pc,burnin,inits.pc,3,4,1:2), 
                                 nsim.= nsim., ksamp.= 1,burnin = burnin,hs=3,N.=daf$N, K = 2, ni.=rep(5,daf$N), 
                                 rx. = daf$rx, fitRx = daf$fitRx, report1.=1000,id=rep(1:daf$N,5),
                                 hyperpar=c(3,1,1,.001,.0002), run. = no.run,emission.fix = T)
  x_mhmm_p_fix <- simulated_mhmm(y = list(daf$y[[3]],daf$y[[4]]), inits.= get_fix_inits(x_mhmm_pc,burnin,inits.pc,3,4,3:4), 
                                 nsim.= nsim., ksamp.= 1,burnin = burnin,hs=3, N.=daf$N, K = 2, ni.=rep(5,daf$N), 
                                 rx. = daf$rx, fitRx = daf$fitRx, report1.=1000,id=rep(1:daf$N,5),
                                 hyperpar=c(3,1,1,.001,.0002), run. = no.run,emission.fix = T)
  z.c <- apply(x_mhmm_c_fix$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  z.p <- apply(x_mhmm_p_fix$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  Pm.est_fix <- c(pm_est(z.pc,z.c),pm_est(z.pc,z.p))
  x_dhmm <- simulated(y = daf$y, inits. = c(daf$pars[-x_dhmm2$pidx$Pm1]), nsim.= nsim., burnin = burnin, ksamp.= 1, 
                       Km = c(2,2), hs = c(3,3,3), N.=daf$N, ni.=daf$ni, rx. = daf$rx, fitRx = c(F,T,F), report1.=1000,
                       id=rep(1:daf$N,5),run. = no.run)
  z.pred <- apply(x_dhmm$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
  confusion <- rbind(confusion,get_confusion(z.pred,z.true))
  for(m in 1:2){
    zm.true <- daf$Zm[[m]]
    zm.pred <- apply(x_dhmm$Zmmax[[m]][-(1:burnin),],2,function(x) names(which.max(table(x))))
    confusion <- rbind(confusion,get_confusion(zm.pred,zm.true))
  }
  res <- list(confusion = confusion, true.param = daf$pars, 
              WAIC = c(x_dhmm$WAIC,x_mhmm_pc$WAIC,x_mhmm_c$WAIC,x_mhmm_p$WAIC,
                       x_mhmm_c_fix$WAIC,x_mhmm_p_fix$WAIC,x_dhmm2$WAIC),
              Pm.est = Pm.est, Pm.est_fix = Pm.est_fix,
              dhmm2.m = colMeans(x_dhmm2$betaa[-(1:burnin),1:max(x_dhmm2$pidx$Pm1)]),
              dhmm2.sd = apply(x_dhmm2$betaa[-(1:burnin),1:max(x_dhmm2$pidx$Pm1)],2,sd),
              dhmm2.med = apply(x_dhmm2$betaa[-(1:burnin),1:max(x_dhmm2$pidx$Pm1)],2,median),
              dhmm2.q1 = apply(x_dhmm2$betaa[-(1:burnin),1:max(x_dhmm2$pidx$Pm1)],2,quantile,0.025),
              dhmm2.q2 = apply(x_dhmm2$betaa[-(1:burnin),1:max(x_dhmm2$pidx$Pm1)],2,quantile,0.975),
              dhmm2.cover = (daf$pars > apply(x_dhmm2$betaa[-(1:burnin),1:max(x_dhmm2$pidx$Pm1)],2,quantile,0.025))&
                  (daf$pars < apply(x_dhmm2$betaa[-(1:burnin),1:max(x_dhmm2$pidx$Pm1)],2,quantile,0.975)),
              dhmm.m = colMeans(x_dhmm$betaa[-(1:burnin),1:max(x_dhmm$pidx$P1)]),
              dhmm.sd = apply(x_dhmm$betaa[-(1:burnin),1:max(x_dhmm$pidx$P1)],2,sd),
              dhmm.med = apply(x_dhmm$betaa[-(1:burnin),1:max(x_dhmm$pidx$P1)],2,median),
              dhmm.q1 = apply(x_dhmm$betaa[-(1:burnin),1:max(x_dhmm$pidx$P1)],2,quantile,0.025),
              dhmm.q2 = apply(x_dhmm$betaa[-(1:burnin),1:max(x_dhmm$pidx$P1)],2,quantile,0.975),
              mhmm.m = colMeans(x_mhmm_pc$betaa[-(1:burnin),1:max(x_mhmm_pc$pidx$P1)]),
              mhmm.sd = apply(x_mhmm_pc$betaa[-(1:burnin),1:max(x_mhmm_pc$pidx$P1)],2,sd),
              mhmm.med = apply(x_mhmm_pc$betaa[-(1:burnin),1:max(x_mhmm_pc$pidx$P1)],2,median),
              mhmm.q1 = apply(x_mhmm_pc$betaa[-(1:burnin),1:max(x_mhmm_pc$pidx$P1)],2,quantile,0.025),
              mhmm.q2 = apply(x_mhmm_pc$betaa[-(1:burnin),1:max(x_mhmm_pc$pidx$P1)],2,quantile,0.975),
              mhmm_c.m = colMeans(x_mhmm_c$betaa[-(1:burnin),1:max(x_mhmm_c$pidx$P1)]),
              mhmm_c.sd = apply(x_mhmm_c$betaa[-(1:burnin),1:max(x_mhmm_c$pidx$P1)],2,sd),
              mhmm_c.med = apply(x_mhmm_c$betaa[-(1:burnin),1:max(x_mhmm_c$pidx$P1)],2,median),
              mhmm_c.q1 = apply(x_mhmm_c$betaa[-(1:burnin),1:max(x_mhmm_c$pidx$P1)],2,quantile,0.025),
              mhmm_c.q2 = apply(x_mhmm_c$betaa[-(1:burnin),1:max(x_mhmm_c$pidx$P1)],2,quantile,0.975),
              mhmm_p.m = colMeans(x_mhmm_p$betaa[-(1:burnin),1:max(x_mhmm_p$pidx$P1)]),
              mhmm_p.sd = apply(x_mhmm_p$betaa[-(1:burnin),1:max(x_mhmm_p$pidx$P1)],2,sd),
              mhmm_p.med = apply(x_mhmm_p$betaa[-(1:burnin),1:max(x_mhmm_p$pidx$P1)],2,median),
              mhmm_p.q1 = apply(x_mhmm_p$betaa[-(1:burnin),1:max(x_mhmm_p$pidx$P1)],2,quantile,0.025),
              mhmm_p.q2 = apply(x_mhmm_p$betaa[-(1:burnin),1:max(x_mhmm_p$pidx$P1)],2,quantile,0.975))
  return(res)
}

cl <- doMPI::startMPIcluster()
registerDoMPI(cl)
out <- foreach(i = 1: nrun,.combine="rbind",.multicombine=T,.errorhandling="pass")%dopar%{
  sim_each(i)
}
save('out',file="dhmm2_model.Rdata")
stopCluster(cl)
mpi.exit()