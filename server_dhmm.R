suppressMessages({
  library(foreach)
  library(doMPI)
})
#library(doMC)
source("DHMM.R")
source('evaluation_measure.R')
cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

nrun=10

server.sim <- function(no.run, seed = 2256,burnin = 7000, nsim. = 10000, hs=c(3,3,3)){
  set.seed(seed + no.run);
  daf <- genSim(N = 354, rx. = rep(c(0,1),each = 354/2),fitRx = c(FALSE,TRUE,FALSE));
  x <- simulated(y = daf$y, inits.=c(daf$pars), nsim.= nsim., burnin = burnin, ksamp.= 1, Km = c(2,2), hs = hs,
                 N.=daf$N, ni.=daf$ni, rx. = daf$rx, fitRx = daf$fitRx, report1.=1000,id=rep(1:daf$N,5),run. = no.run)
  
  z.pred <- apply(x$Zmax[-(1:burnin),], 2, function(x) names(which.max(table(x))))
  z.true <- daf$Z
  zm.true <- daf$Zm
  M <- length(zm.true)
  precision <- recall <- Fs <- acc <- err <- precision_m <- rep(NA, M+1)
  precision[1] <- pre(z.pred,z.true)
  recall[1] <- rec(z.pred,z.true)
  Fs[1] <- Fscore(z.pred,z.true)
  acc[1] <- accuracy(z.pred, z.true)
  err[1] <- errRate(z.pred, z.true)
  precision_m[1] <- pre(z.pred, z.true, "micro")
  for(m in 1:M){
    zm.true <- daf$Zm[[m]]
    zm.pred <- apply(x$Zmmax[[m]][-(1:burnin),],2,function(x) names(which.max(table(x))))
    precision[m+1] <- pre(zm.pred, zm.true)
    recall[m+1] <- rec(zm.pred, zm.true)
    Fs[m+1] <- Fscore(zm.pred, zm.true)
    acc[m+1] <- accuracy(zm.pred, zm.true)
    err[m+1] <- errRate(zm.pred, zm.true)
    precision_m[m+1] <- pre(zm.pred, zm.true, "micro")
  }
  beta.sim <- x$betaa[-(1:burnin),1:length(daf$pars)]
  res <- list(beta.m = apply(beta.sim,2,mean),
              beta.sd = apply(beta.sim,2,sd),
              beta.cover = (daf$pars > apply(beta.sim,2,quantile,0.025)) & (daf$pars < apply(beta.sim,2,quantile,0.975)),
              precision = precision, recall = recall, Fscore = Fs, accuracy = acc, errRate = err, precision_m = precision_m,
              WAIC = x$WAIC)
  return(res)
}

res <- foreach(no.run = 1:nrun,.packages = "MCMCpack",.combine = 'rbind',.multicombine = T, .errorhandling = 'pass') %dopar% {
  source("DHMM.R",local=T);
  source('evaluation_measure.R', local = T);
  server.sim(no.run, hs = c(4,4,4))
}
save("res", file = paste0('../result/dhmm_hs4_',nrun,'.Rdata'))

closeCluster(cl)
mpi.quit()
