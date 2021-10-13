####--------Functions to obtain the initial values-----######
# for the mhmm model from dhmm model
## pars: the parameters' values
## pidx: the index list for the parameters
get_inits_pc <- function(pars,pidx){
  return(pars[c(pidx$pi0,pidx$P0,pidx$tau,pidx$Sigma,pidx$invsigma,pidx$P1)])
}
# for the initial values from the family's model
get_inits_c <- function(pars,pidx){
  hs <- length(pidx$pi0) + 1
  tauidx <- pidx$tau[1:(2*hs)]
  Sigma <- matrix(pars[pidx$Sigma],4,4)
  Sigma <- Sigma[1:2,1:2]
  return(c(pars[c(pidx$pi0,pidx$P0,tauidx)],Sigma,pars[c(pidx$invsigma[1:(2*hs)],pidx$P1)]))
}
# for the parent's initial value from the family's model
get_inits_p <- function(pars,pidx){
  hs <- length(pidx$pi0) + 1
  tauidx <- pidx$tau[(2*hs+1):(4*hs)]
  Sigma <- matrix(pars[pidx$Sigma],4,4)
  Sigma <- Sigma[3:4,3:4]
  return(c(pars[c(pidx$pi0,pidx$P0,tauidx)],Sigma,pars[c(pidx$invsigma[(2*hs+1):(4*hs)],pidx$P1)]))
}
#add perception matrix to the dhmm model's initial

get_inits_dhmm <- function(pars,pidx,fitPm = F){
  Pm <- diag(length(pidx$pi0) + 1)
  if(!fitPm){
    return(c(pars[c(pidx$pi0, pidx$P0)],c(Pm[,-ncol(Pm)],Pm[,-ncol(Pm)]),pars[c(pidx$tau,pidx$Sigma,pidx$invsigma,pidx$P1)]))
  }else{
    return(c(pars[c(pidx$pi0, pidx$P0)],c(Pm[,-ncol(Pm)],Pm[,-ncol(Pm)]),
             pars[c(pidx$tau,pidx$Sigma,pidx$invsigma,pidx$P1)],
             c(Pm[,-ncol(Pm)],Pm[,-ncol(Pm)])))
  }
}
# get the confusion matrix
## z.pred : the true predicted labels
## z.true: the true labels
get_confusion <- function(z.pred,z.true){
  return(data.frame(precision = pre(z.pred,z.true), recall = rec(z.pred, z.true), Fscore = Fscore(z.pred, z.true), accuracy = accuracy(z.pred, z.true),
                    err = errRate(z.pred, z.true), precision_m = pre(z.pred, z.true, "micro")))
}
# obtain the perception matrix
## z.pc: the family labels
## z.m: the member's label
pm_est <-function(z.pc,z.m){
  tab <- table(z.pc,z.m)
  tab <- tab / matrix(rep(rowSums(tab),ncol(tab)),nr = nrow(tab))
  return(c(tab[,-ncol(tab)]))
}

# for the constrained model
## sample: the MCMC object from the mhmm model
## burnin: the number of samples to burn-in
## inits: the initial values for sample
## hs: the number of hidden states
## K: the number of variables in the overall model(the sample)
## Km: for which variables are searching for the fixed initial values

get_fix_inits <- function(sample,burnin,inits,hs,K,Km){
  pi <- inits[sample$pidx$pi0]
  P0 <- inits[sample$pidx$P0]
  tau <- colMeans(sample$betaa[-(1:burnin),sample$pidx$tau[((Km[1]-1)*hs+1):(Km[length(Km)]*hs)]])
  Sigma <- colMeans(sample$betaa[-(1:burnin),sample$pidx$Sigma])
  Sigma <- matrix(Sigma,K,K)
  Sigma <- as.numeric(Sigma[Km,Km])
  invsigma <- colMeans(sample$betaa[-(1:burnin),sample$pidx$invsigma[((Km[1]-1)*hs+1):(Km[length(Km)]*hs)]])
  P1 <- inits[sample$pidx$P1]
  return(c(pi,P0,tau,Sigma,invsigma,P1))
}

pre <- function(x.pred,x.true,average = "macro"){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K == 2){
    TP <- conf[1,1]
    FP <- conf[1,2]
    precision <- TP / (TP + FP)
  }else if(K > 2){
    if(average == "micro"){
      precision <-  sum(diag(conf)) / sum(conf)
    }
    else if(average == "macro"){
      prei = 0
      for(k in 1:K){
        TP = conf[k,k]
        FP = sum(conf[k,]) - TP
        if(TP + FP !=0){
          prei = prei + TP / (TP + FP)
        }
      }
      precision = prei / K
    }
    else{
      stop("The option for average method is either micro or macro.")
    }
  }
  else{
    stop("A confusion matrix at least has two rows.")
  }
  return(precision)
}

rec <- function(x.pred, x.true, average="macro"){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K == 2){
    TP <- conf[1,1]
    FN <- conf[2,1]
    recall <- TP / (TP + FN)
  }else if(K > 2){
    if(average == "micro"){
      recall <-  sum(diag(conf)) / sum(conf)
    }else if(average == "macro"){
      reci <- 0
      for(k in 1:K){
        TP = conf[k,k]
        FN = sum(conf[,k]) - TP
        if(TP + FN != 0){
          reci <- reci + TP / (TP + FN)
        }
      }
      recall <- reci / K
    }else{
      stop("The option for average method is either micro or macro.")
    }
  }else{
    stop("A confusion matrix at least has two rows.")
  }
  return(recall)
}

Fscore <- function(x.pred, x.true, average = "macro", beta = 1){
  precision <- pre(x.pred, x.true, average)
  recall <- rec(x.pred, x.true, average)
  Fscore <- (beta^2 + 1) * precision * recall / (beta^2 * precision + recall)
  return(Fscore)
}

errRate <- function(x.pred, x.true){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K > 2){
    erri <- 0
    for(k in 1: K){
      TP <- conf[k,k]
      FP <- sum(conf[k,]) - TP
      FN <- sum(conf[,k]) - TP
      erri <- erri + (FP + FN) / sum(conf)
    }
    errRate <- erri / K
  }else{
    stop("Error Rate is an evaluation measure for multi-class classification.")
  }
  return(errRate)
}

accuracy <- function(x.pred, x.true){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K > 2){
    accui <- 0
    for(k in 1: K){
      TP <- conf[k,k]
      FP <- sum(conf[k,]) - TP
      FN <- sum(conf[,k]) - TP
      TN <- sum(conf) - TP - FP - FN
      accui <- accui + (TP + TN) / sum(conf)
    }
    accuracy <- accui / K
  }else{
    stop("Error Rate is an evaluation measure for multi-class classification.")
  }
  return(accuracy)
}
