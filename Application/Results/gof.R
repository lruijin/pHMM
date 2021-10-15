gof <- function(betaa,pidx,N,M,ni,hs,hs0,X,group){
  niter <- nrow(betaa)
  p <- ncol(X)
  ysimu <- array(NA,dim=c(4,N*ni,niter))
for(r in 1: niter){
  Sigma <- matrix(betaa[r,pidx$Sigma],4,4)
  re <- mvrnorm(N,rep(0,4),Sigma)
  Zmat <- matrix(0,N,ni)
  Zmat[,1] <- sample(1:hs, N, prob=c(betaa[r,pidx$pi0],1-sum(betaa[r,pidx$pi0])),replace = T)
  
  ## P0
  P0 <- betaa[r,pidx$P0]
  P0mat <- matrix(0,hs,hs)
  ctt = 1
  for(i in 1:hs){
    for(j in 1:hs){
      if(i != j){
        P0mat[i,j]=P0[ctt]
        ctt = ctt+1
      }
    }
    P0mat[i,i] = 1 - sum(P0mat[i,])
  }
  ## P1
  P1 <- betaa[r,pidx$P1]
  P1mat <- matrix(0,hs,hs)
  ctt = 1
  for(i in 1:hs){
    for(j in 1:hs){
      if(i != j){
        P1mat[i,j]=P1[ctt]
        ctt = ctt+1
      }
    }
    P1mat[i,i] = 1 - sum(P1mat[i,])
  }
  
  for(t in 2:ni){
    for(i in 1:N){
      if(group[i] == 0){
        Zmat[i,t] = sample(1:hs,1,prob=P0mat[Zmat[i,t-1],])
      }else if(group[i] == 1){
        Zmat[i,t] = sample(1:hs,1,prob=P1mat[Zmat[i,t-1],])
      }else{
        stop("Groups are labeled as 0 and 1")
      }
    }
  }
  
  ## Pm0
  Pm0 <- betaa[r,pidx$Pm0]
  Pm0l <- list()
  ctt <- 0
  for(m in 1:M){
    Pm0mat <- matrix(0,hs,hs0[m])
    hs0tmp <- hs0[m]
    Pm0mat[,1:(hs0tmp-1)] <- matrix(Pm0[ctt +(1:(hs * (hs0tmp-1)))],nr = hs)
    Pm0mat[,hs0tmp] <- 1- rowSums(Pm0mat)
    Pm0l[[m]] <- Pm0mat
    ctt <- ctt + hs * (hs0tmp - 1)
  }
  
  Zm <- matrix(NA, 2, N*ni)
  Zm[1,] <- unlist(lapply(c(Zmat),function(x) sample(1:hs,1, prob=Pm0l[[1]][x,])))
  Zm[2,] <- unlist(lapply(c(Zmat),function(x) sample(1:hs,1, prob=Pm0l[[2]][x,])))
  
  ZZm <- list()
  for(m in 1:2){
    tmpZZm <- matrix(NA,nr=N * ni,nc = hs)
    tmpZZm[,1] <- 1
    for(s in 2:hs){
      tmpZZm[,s] <- as.numeric(1*Zm[m,] > (s-1))
    }
    ZZm[[m]] <- tmpZZm
  }
  
  tauidx <- betaidx <- 0
  for(m in 1:2){
    for(k in 1:2){
      y_m = ZZm[[m]] %*% betaa[r,pidx$tau[tauidx+(1:hs)]] + 
        c(matrix(X %*% betaa[r,pidx$beta[betaidx+(1:p)]],ncol=ni,byrow=T))+
        rep(re[,(m-1)*2+k],ni)
      invsigma <- betaa[r,pidx$invsigma[tauidx + (1:hs)]]
      invsigma <- c(invsigma[1],diff(invsigma))
      y_v = sqrt(1/ (ZZm[[m]] %*%  invsigma))
      ysimu[(m-1)*2 +k,,r] = rnorm(length(y_m),y_m,y_v)
      tauidx <- tauidx + hs
      betaidx <- betaidx + p
    }
  }
}
  return(ysimu)
}
