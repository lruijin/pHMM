rm(list=ls())
setwd("~/Project2_DMHMM/code")
load('../application/var_hs/ori/dhmm1/dhmm_cov_hs3.Rdata')
burnin = 30000
hs = 3
library(ggplot2)
library(ggpubr)
data <- read.csv("C:/Users/lur5.NIH/OneDrive - National Institutes of Health/Project1_LCA/data/Fmod_imputed.csv",header = T)
data.na <- na.omit(data)
A1C <- read.csv("C:/Users/lur5.NIH/OneDrive - National Institutes of Health/Project1_LCA/data/A1C.csv",header = T)
y_c = list(t(matrix(data.na$CxDFR,nr = 5)),t(matrix(log(data.na$CxPCC),nr = 5)))
y_p = list(t(matrix(data.na$PxDFR,nr = 5)),t(matrix(log(data.na$PxPCC),nr = 5)))
X = read.csv("../data/age_sex.csv",header=T)
N <- nrow(y_c[[1]])
ni <- 5
y <- c(y_c,y_p)

y_obs <- matrix(NA,4,N*ni)
for(k in 1:4){
  y_obs[k,] <- c(y[[k]])
}

# plot the data
group <- A1C$GROUP-1
tit <- c("CDFR","PDFR","CPCC","PPCC")
y_tit <- c("CDFR","log(CPCC)","PDFR","log(PPCC)")
for(k in 1:length(y)){
  plot_data <- data.frame(y = c(y[[k]]), Arm = rep(group,ni),
                          x = rep(1:5, each = nrow(y[[k]])))
  plot_data$Arm <- factor(plot_data$Arm)
  levels(plot_data$Arm)= c("Usual-care", "Intervention")
  assign(paste0("p",k), ggboxplot(plot_data, x= "x",y="y",fill ="Arm",palette = "jco") + 
           labs(x = "Visit", y = y_tit[k]) + theme_bw())
}
ggarrange(p1,p3,p2,p4,labels = tit,common.legend = T, vjust = -0.5, legend = "bottom")+
  theme(plot.margin = margin(t=0.8,r = 0, b = 0.2, l =0, unit = "cm"))
betaa <- dhmm$betaa[-(1:burnin),]
library(MASS)
source("gof.R")

y_simu <- gof(betaa,dhmm$pidx,N,2,ni,hs,hs0,as.matrix(X,ncol=2),group) 
save(y_simu, file = "GOF_y_simu.Rdata")
## obtain the simulated distributions:
# obtain the means of the 4 reponse variables
perc <- c(0.025,0.25,0.5,0.75,0.975)
res <- matrix(0,4,length(perc) + 2)

for(i in 1:4){
  ysim <- y_simu[i,,]
  ym_sim <- colMeans(ysim)
  ym_obs <- mean(y_obs[i,])
  res[i,1:(ncol(res)-2)] <- quantile(ym_sim,perc)
  res[i,ncol(res)-1] <- ym_obs
  res[i,ncol(res)] <- sum(ym_sim <= ym_obs) / ncol(ysim)
}
res
# obtain the variances of the 4 response variables
for(i in 1:4){
  ysim <- y_simu[i,,]
  yv_sim <- apply(ysim,2,var)
  yv_obs <- var(y_obs[i,])
  res[i,1:(ncol(res)-2)] <- quantile(yv_sim,perc)
  res[i,ncol(res)-1] <- yv_obs
  res[i,ncol(res)] <- sum(yv_sim <= yv_obs) / ncol(ysim)
}
print(round(res,3))

# obtain the residual plots:
tit <- c("CDFR","CPCC","PDFR","PPCC")
res_plot_data <- NULL
for(i in 1:4){
  ysim <- y_simu[i,,]
  yobs <- y_obs[i,]
  y_merge <- cbind(yobs,ysim)
  x <- apply(y_merge,1,function(x) sum(x[-1] <= x[1])/10000)
  z_obs <- qnorm(unlist(lapply(yobs,function(x) sum(yobs <= x)/1950)))
  z_norm <- qnorm(x)
  res_plot_data = rbind(res_plot_data,data.frame(z_obs,z_norm,construct = tit[i]))
  # plot(z_norm,z_obs,ylab="Observed quantile",xlab="Normal quantile",main=tit[i])
  # abline(0,1,col=2)
}
names(res_plot_data)
ggplot(res_plot_data, aes(x = z_norm, y = z_obs))+
  geom_point()+
  labs(x="Normal quantile", y = "Observed quantile")+
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_wrap(~construct)+
  theme_bw()
source('plot_utils.R')
# lik = dhmm$lik[20001:30000,]
# llik = dhmm$llik[20001:3000,]
# WAIC = -2 * (sum(log(colMeans(lik))) - sum(apply(llik,2,var)))
# WAIC
#load('dhmm_hs3_22.Rdata')
    hs = 3
    hs0 = c(3,3)
    M=2
    burnin = 30000
    dhmm$WAIC
    betaa <- dhmm$betaa[-(1:burnin),]
## pi
    if(hs == 2){
    	pi <- mean(betaa[,dhmm$pidx$pi0])
    }else{
      pi <- colMeans(betaa[,dhmm$pidx$pi0])
    }
    pi <- c(pi,1-sum(pi))
    print(paste(round(pi,2),collapse = " & "))
    
    pi_seq <- cbind(betaa[,dhmm$pidx$pi0,drop=F],1-rowSums(betaa[,dhmm$pidx$pi0,drop=F]))
    pi_lower <- round(apply(pi_seq,2,quantile,0.025),2)
    pi_upper <- round(apply(pi_seq,2,quantile,0.975),2)
    
    print(paste(paste0(round(pi,2)," (",pi_lower,",",pi_upper,")"),collapse = " & "))
## P0
    P0 <- betaa[,dhmm$pidx$P0]
    P0_seq <- array(NA,dim=c(hs,hs,nrow(P0)))
    for(k in 1:nrow(P0)){
      P0mat <- matrix(0,hs,hs)
      ctt = 1
      for(i in 1:hs){
        for(j in 1:hs){
          if(i != j){
            P0mat[i,j]=P0[k,ctt]
            ctt = ctt+1
          }
        }
        P0mat[i,i] = 1 - sum(P0mat[i,])
      }
      P0_seq[,,k] <- P0mat
    }
    
    print(apply(P0_seq,c(1,2),mean))
    print(paste(paste0(t(round(apply(P0_seq,c(1,2),mean),2))," (",
                       t(round(apply(P0_seq,c(1,2),quantile,0.025),2)),",",
                       t(round(apply(P0_seq,c(1,2),quantile,0.975),2)),")"),collapse = " & "))
    
## P1
    P1 <- betaa[,dhmm$pidx$P1]
    P1_seq <- array(NA,dim=c(hs,hs,nrow(P1)))
    for(k in 1:nrow(P1)){
      P1mat <- matrix(0,hs,hs)
      ctt = 1
      for(i in 1:hs){
        for(j in 1:hs){
          if(i != j){
            P1mat[i,j]=P1[k,ctt]
            ctt = ctt+1
          }
        }
        P1mat[i,i] = 1 - sum(P1mat[i,])
      }
      P1_seq[,,k] <- P1mat
    }
    print(paste(paste0(t(round(apply(P1_seq,c(1,2),mean),2))," (",
                       t(round(apply(P1_seq,c(1,2),quantile,0.025),2)),",",
                       t(round(apply(P1_seq,c(1,2),quantile,0.975),2)),")"),collapse = " & "))
## Pm0
    Pm0 <- betaa[,dhmm$pidx$Pm0]
    Pm0l <- list()
    ctt <- 0
    for(m in 1:M){
      hs0tmp <- hs0[m]
      Pm0_seq <- array(NA, dim=c(hs,hs0tmp,nrow(Pm0)))
      for(k in 1:nrow(Pm0)){
        Pm0mat <- matrix(0,hs,hs0tmp)
        Pm0mat[,1:(hs0tmp-1)] <- matrix(Pm0[k,ctt +(1:(hs * (hs0tmp-1)))],nr = hs)
        Pm0mat[,hs0tmp] <- 1- rowSums(Pm0mat)
        Pm0_seq[,,k] <- Pm0mat
      }
	    Pm0l[[m]] <- Pm0_seq
      ctt <- ctt + hs * (hs0tmp - 1)
    }
    for(m in 1:M){
      print(paste(paste0(t(round(apply(Pm0l[[m]],c(1,2),mean),2))," (",
                         t(round(apply(Pm0l[[m]],c(1,2),quantile,0.025),2)),", ",
                         t(round(apply(Pm0l[[m]],c(1,2),quantile,0.975),2)),")"),collapse = " & "))
    }
    

## Pm1
    Pm1 <- round(apply(betaa[,dhmm$pidx$Pm1],2,mean),2)
    Pm1l <- list()
    Pm1mat <- matrix(0,hs,hs)
    ctt <- 0
    for(m in 1:M){
	    Pm1mat[,1:(hs-1)] <- matrix(Pm1[ctt +(1:(hs * (hs-1)))],nr = hs)
      Pm1mat[,hs] <- 1- rowSums(Pm1mat)
      print(apply(Pm1mat,1, function(x) paste(round(x,2),collapse = " & ")))
	    Pm1l[[m]] <- Pm1mat
      ctt <- ctt + hs * (hs - 1)
    	Pm1mat <- matrix(0,hs,hs)
    }
    
## beta
   round(colMeans(dhmm$betaa[-(1:burnin),dhmm$pidx$beta]),3)
   paste0("(",round(apply(dhmm$betaa[-(1:burnin),dhmm$pidx$beta],2,quantile,0.025),3),", ",
          round(apply(dhmm$betaa[-(1:burnin),dhmm$pidx$beta],2,quantile,0.975),3),")")
## tau
    ctt <- 0
    tau <- list()
    for(m in 1:2){
      tau[[m]] <- matrix(NA,hs0[m],2)
      for(i in 1:2){
      print(tau[[m]][,i] <- round(rowMeans(apply(dhmm$betaa[-(1:burnin),dhmm$pidx$tau[ctt + (i-1)*hs0[m]+(1:hs0[m])]],1,cumsum)),2))
      }
      ctt <- ctt + hs0[m]*2
    }
    ctt <- 0
    for(m in 1:2){
      for(i in 1:2){
        print(paste(paste0(tau[[m]][,i]," (",round(apply(apply(dhmm$betaa[-(1:burnin),dhmm$pidx$tau[ctt + (i-1)*hs0[m]+(1:hs0[m])]],1,cumsum),1,quantile,0.025),2),", ",
        round(apply(apply(dhmm$betaa[-(1:burnin),dhmm$pidx$tau[ctt + (i-1)*hs0[m]+(1:hs0[m])]],1,cumsum),1,quantile,0.975),2),")"),collapse=" & "))
      }
      ctt <- ctt + hs0[m]*2
    }
## Sigma
    print(matrix(round(colMeans(dhmm$betaa[-(1:burnin),dhmm$pidx$Sigma]),3),4,4))
## invsigma
    sigma <- 1/sqrt(colMeans(dhmm$betaa[-(1:burnin),dhmm$pidx$invsigma]))
    sigma <- matrix(sigma,nr=hs)
    #print(round(t(sigma),3))
    
    paste(paste0(unlist(tau)," (",round(sigma,3),")"),collapse = " & ")
## use group, Pm and tau to estimate the marginal taus for the families
A1C <- read.csv("C:/Users/lur5.NIH/OneDrive - National Institutes of Health/Project1_LCA/data/A1C.csv",header = T)
group = table(A1C$GROUP)/390
tau_f <- matrix(NA,hs,4)
for(k in 1:2){
	tau_f[,k] = (Pm0l[[1]] * group[1] + Pm1l[[1]] * group[2]) %*% tau[[1]][,k]
}
for(k in 3:4){
	tau_f[,k] = (Pm0l[[2]] * group[1] + Pm1l[[2]] * group[2]) %*% tau[[2]][,k-2]
}


#for(i in 1:2){
#  for(j in 1:3){
#    print(round(mean(y_p[[i]][apply(dhmm$Zmax[30001:40000,],2,function(x) names(which.max(table(x))))==j]),3))
#  }
#}
#mean(y_c[[1]][apply(dhmm$Zmax[30001:40000,],2,function(x) names(which.max(table(x))))==1])

sigma_f <- matrix(NA,hs,4)
for(k in 1:2){
  sigma_f[,k] = (Pm0l[[1]] * group[1] + Pm1l[[1]] * group[2])^2 %*% sigma[,k]^2
}
for(k in 3:4){
  sigma_f[,k] = (Pm0l[[2]] * group[1] + Pm1l[[2]] * group[2])^2 %*% sigma[,k]^2
}
apply(matrix(paste0(round(tau_f,3)," (",round(sqrt(sigma_f),3),")"),nc=hs,byrow = T),1,
      function(x) paste(x,collapse = " & "))


tau_f <- matrix(NA,hs,4)
for(k in 1:2){
	tau_f[,k] = Pm0l[[1]] %*% tau[[1]][,k]
}
for(k in 3:4){
	tau_f[,k] = Pm0l[[2]] %*% tau[[2]][,k-2]
}
apply(t(round(tau_f,2)),1,function(x) paste(x,collapse = " & "))

sigma_f <- matrix(NA,hs,4)
for(k in 1:2){
  sigma_f[,k] = Pm0l[[1]]^2 %*% sigma[(k-1)*hs0[1] +(1:hs0[1])]^2
}
for(k in 3:4){
  sigma_f[,k] = Pm0l[[2]]^2 %*% sigma[2*hs0[1] + (k-2-1)*hs0[2] + (1:hs0[2])]^2
}
#t(round(sqrt(sigma_f),3))

apply(matrix(paste0(round(tau_f,3)," (",round(sqrt(sigma_f),3),")"),nc=hs,byrow = T),1,
      function(x) paste(x,collapse = " & "))

plot_tau_g(dhmm,hs,2)
plot_tau_m(dhmm,hs)


##---------------TPM difference--------------####
library(matrixcalc)
TPM_diff <- array(NA,c(hs,hs,10000))
stat_diff <- matrix(NA,10000,hs)
betaa <- dhmm$betaa[-(1:burnin),]
for(r in 1:10000){
  ## pi
  if(hs == 2){
    pi <- betaa[r,dhmm$pidx$pi0]
  }else{
    pi <- betaa[r,dhmm$pidx$pi0]
  }
  pi <- c(pi,1-sum(pi))
  ## P0
  P0 <- betaa[r,dhmm$pidx$P0]
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
  P1 <-betaa[r,dhmm$pidx$P1]
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
  ## P1- P0
  TPM_diff[,,r] <- P1mat - P0mat
  stat_diff[r,] <- pi%*%(matrix.power(P1mat,5) - matrix.power(P0mat,5))
}
library(plotrix)
#par(mar=c(0.1, 4, 0, 0.8), mfrow=c(2,1),
#    oma = c(4, 4, 0.2, 0.2))
# par(mar = c(0.5,4,2,0.8),mfrow=c(2,1))
# plotCI(x=1:9,y=apply(TPM_diff,c(1,2),mean),li = apply(TPM_diff,c(1,2),quantile,0.025),
#        ui = apply(TPM_diff,c(1,2),quantile,0.975),ylab="Difference between arms",xlab="Transition",
#        xaxt="n", main = "Transition probabilities")
# axis(1, at=1:9, labels=F)
# text(x = 1:9, y = par("usr")[3], labels = expression(paste(rep(c("Discordant","Harmonious","Indifferent"),3), "to", 
#                                                       rep(c("Discordant","Harmonious","Indifferent"),each=3))),
#      xpd = NA, srt = 35, cex = 0.9)
# abline(h=0)
# plotCI(x=c(1.25,2,2.75),y=apply(stat_diff,2,mean),li = apply(stat_diff,2,quantile,0.025),
#        ui = apply(stat_diff,2,quantile,0.975),ylab="Difference between arms",xlab="Hidden States",
#        xaxt="n",main = "Stationary probabilities", ylim=c(-0.3,0.3),xlim=c(1,3))
# axis(1, at=c(1.25,2,2.75), labels=c("Discordant","Harmonious","Indifferent"))
# abline(h=0)

probPlotData_tpm <- data.frame(x = paste(rep(c("Discordant","Harmonious","Indifferent"),3), "to", 
                                         rep(c("Discordant","Harmonious","Indifferent"),each=3)),
                               y = c(apply(TPM_diff,c(1,2),mean)),
                               yl = c(apply(TPM_diff,c(1,2),quantile, 0.025)),
                               yu = c(apply(TPM_diff,c(1,2),quantile, 0.975)),
                               label = "Transition probabilities")
probPlotData_tpm$x <- gsub(" ","\n", probPlotData_tpm$x)
probPlotData_stat <- data.frame(x = c("Discordant","Harmonious","Indifferent"),
                                y = c(apply(stat_diff,2,mean)),
                                yl = c(apply(stat_diff,2,quantile, 0.025)),
                                yu = c(apply(stat_diff,2,quantile, 0.975)),
                                label = "Stationary probabilities")
tmp_plot <- ggplot(probPlotData_tpm, aes(x=x,y = y, ymin = yl, ymax = yu))+
  geom_pointrange()+
  geom_abline(slope = 0, intercept = 0, col="red")+
  labs(y = "Difference between arms", x = " ") + 
  theme_bw()
stat_plot <- ggplot(probPlotData_stat, aes(x=x,y = y, ymin = yl, ymax = yu))+
  geom_pointrange()+
  geom_abline(slope = 0, intercept = 0, col="red")+
  labs(y = "Difference between arms", x = " ") + 
  theme_bw()
#tmp_plot
#stat_plot
ggarrange(tmp_plot,stat_plot,
          labels = c("Transition probabilities","Stationary probabilities"),
          ncol = 1, nrow=2,
          vjust = -0.1)+
  theme(plot.margin = margin(t=1.2,r = 0, b = 0.4, l =0, unit = "cm"))
##---------------Perception differences--------------####
Pm_P_C0 <- array(NA,c(hs,hs0[1],10000))
for(r in 1:10000){
  ## Pm0
  Pm0 <- betaa[r,dhmm$pidx$Pm0]
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
  Pm_P_C0[,,r] <- Pm0l[[2]]-Pm0l[[1]]
}
probPlotData_pm <- data.frame(x = paste(rep(c("Discordant","Harmonious","Indifferent"),3), "to", 
                                         rep(c("Discordant","Harmonious","Indifferent"),each=3)),
                               y = c(apply(Pm_P_C0,c(1,2),mean)),
                               yl = c(apply(Pm_P_C0,c(1,2),quantile, 0.025)),
                               yu = c(apply(Pm_P_C0,c(1,2),quantile, 0.975)))
probPlotData_pm$x <- gsub(" ","\n", probPlotData_pm$x)

pm_plot <- ggplot(probPlotData_pm, aes(x=x,y = y, ymin = yl, ymax = yu))+
  geom_pointrange()+
  geom_abline(slope = 0, intercept = 0, col="red")+
  labs(y = "Difference between parent and child", x = " ") + 
  theme_bw()
pm_plot  

# trace plot fixed effects
## functions are in "plot_utils.R"
plot_tau_m_gg(dhmm,hs,mar=0.5,cex=0.8,ylabs=c("CDFR","log(CPCC)","PDFR","log(PPCC)"),
                          tit=c("CDFR","CPCC","PDFR","PPCC"))
plot_tau_g_gg(dhmm,hs,M,Km=c(2,2),mar=0.5,cex = 0.8,
                          ylabs=c("CDFR","log(CPCC)","PDFR","log(PPCC)"),
                          tit=c("CDFR","CPCC","PDFR","PPCC"))
# trace plot for perception matrix

par(mfrow=c(2,3))
M=2
hs = 3
ctt_pm = 0
for(m in 1:M){
  for(i in 1:hs){
    idx <- ctt_pm + seq(i,((hs-2)*hs+i),by=hs)
    tmpPm <- cbind(dhmm$betaa[,dhmm$pidx$Pm0[idx]],1-rowSums(dhmm$betaa[,dhmm$pidx$Pm0[idx]]))
    plot(tmpPm[,1],type="l",col=1,ylim=c(-0.3,1),ylab=paste0("Percption matrix from Family S",i))
    lines(tmpPm[,2],col=2)
    lines(tmpPm[,3],col=3)
    legend("bottomright",legend = c("S1","S2","S3"), col=1:3, lty=1, border = F, bg="transparent")
    
  }
  ctt_pm = ctt_pm + (hs-1) * hs}
mtext("Child",line=-3,outer = T)
mtext("Parent", line=-34,outer=T)

# initial probability and transition matrices plot
plot_prob(P0mat,P1mat,pi,hs=3)
