rm(list=ls())
load('../application/var_hs/ori/mhmm/mhmm_p_cov_hs3.Rdata')    
hs = 3
K = 4
burnin = 30000
mhmm$WAIC
betaa <- mhmm$betaa[-(1:burnin),]
## pi
if(hs == 2){
  pi <- mean(betaa[,mhmm$pidx$pi0])
}else{
  pi <- colMeans(betaa[,mhmm$pidx$pi0])
}
pi <- c(pi,1-sum(pi))
print(paste(round(pi,2),collapse = " & "))

pi_seq <- cbind(betaa[,mhmm$pidx$pi0,drop=F],1-rowSums(betaa[,mhmm$pidx$pi0,drop=F]))
pi_lower <- round(apply(pi_seq,2,quantile,0.025),2)
pi_upper <- round(apply(pi_seq,2,quantile,0.975),2)

print(paste(paste0(round(colMeans(pi_seq),2)," (",pi_lower,",",pi_upper,")"),collapse = " & "))

## P0
P0 <- betaa[,mhmm$pidx$P0]
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
                   t(round(apply(P0_seq,c(1,2),quantile,0.025),2)),", ",
                   t(round(apply(P0_seq,c(1,2),quantile,0.975),2)),")"),collapse = " & "))

## P1
P1 <- betaa[,mhmm$pidx$P1]
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

## beta
round(colMeans(mhmm$betaa[-(1:burnin),mhmm$pidx$beta]),3)
paste0("(",round(apply(mhmm$betaa[-(1:burnin),mhmm$pidx$beta],2,quantile,0.025),3),", ",
       round(apply(mhmm$betaa[-(1:burnin),mhmm$pidx$beta],2,quantile,0.975),3),")")
## tau
  tau <- matrix(NA,hs,K)
  for(i in 1:K){
    taui <- apply(mhmm$betaa[-(1:burnin),mhmm$pidx$tau[(i-1)*hs+(1:hs)]],1,cumsum)
    print(tau[,i] <- round(rowMeans(taui),3))
    print(round(apply(taui,1,quantile,0.025),3))
    print(round(apply(taui,1,quantile,0.975),3))
  }
## Sigma
print(matrix(round(colMeans(mhmm$betaa[-(1:burnin),mhmm$pidx$Sigma]),3),K,K))
## invsigma
print(round(1/sqrt(colMeans(mhmm$betaa[-(1:burnin),mhmm$pidx$invsigma])),3))

##---------------TPM difference--------------####
library(matrixcalc)
TPM_diff <- array(NA,c(hs,hs,10000))
stat_diff <- matrix(NA,10000,hs)
betaa <- mhmm$betaa[-(1:burnin),]
for(r in 1:10000){
  ## pi
  if(hs == 2){
    pi <- betaa[r,mhmm$pidx$pi0]
  }else{
    pi <- betaa[r,mhmm$pidx$pi0]
  }
  pi <-c(pi,1-sum(pi))
  P0 <- betaa[r,mhmm$pidx$P0]
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
  P1 <- betaa[r,mhmm$pidx$P1]
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


## crosstabulate estiamted hidden classes
z.pred <- apply(mhmm_pc$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
z.c0 <- apply(mhmm_c$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))
z.p0 <- apply(mhmm_p$Zmax[-(1:burnin),],2,function(x) names(which.max(table(x))))

tab.fc <- table(z.pred,z.c0)
tab.fp <- table(z.pred,z.p0)
tab.pc <- table(z.p0,z.c0)
tab.fc;round(tab.fc/rowSums(tab.fc),3)
tab.fp;round(tab.fp/rowSums(tab.fp),3)
tab.pc;round(tab.pc/rowSums(tab.pc),3)

# plot transition matrices
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
trans_theme <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text=element_text(color="black",size=11),
                     axis.title=element_text(size=14))
trans_plot0 <-data.frame(P0mat,From = c("Discordant","Harmonious","Indifferent"))
names(trans_plot0)[1:3] <- c("Discordant","Harmonious","Indifferent")
trans_plot0.melted <- melt(trans_plot0)
hm <- ggplot(data = trans_plot0.melted,aes(x=factor(variable),y=factor(From),fill=value)) + geom_tile() +
  geom_text(aes(label=round(value,2)),color="black",size=rel(3.8))+
  ylab("From") + xlab("To")+ labs(fill="Probability")+
  ggtitle("Usual Care")+
  scale_x_discrete()+scale_y_discrete(limits=c("Indifferent","Harmonious","Discordant"))+
  scale_fill_gradient(low="mistyrose",high="firebrick",space = "Lab",limits=c(0,1),na.value="gray90",guide = "colourbar")+
  trans_theme+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title = element_text(size=10),
        legend.text = element_text(size=8))+
  guides(fill=guide_colourbar(barwidth = unit(6,"cm")))
trans_plot1 <-data.frame(P1mat,From = c("Discordant","Harmonious","Indifferent"))
names(trans_plot1)[1:3] <- c("Discordant","Harmonious","Indifferent")
trans_plot1.melted <- melt(trans_plot1)
hm1 <- ggplot(data = trans_plot1.melted,aes(x=factor(variable),y=factor(From),fill=value)) + geom_tile() +
  geom_text(aes(label=round(value,2)),color="black",size=rel(3.8))+
  ylab("") + xlab("To")+ 
  ggtitle("Intervention")+
  scale_x_discrete()+scale_y_discrete(limits=c("Indifferent","Harmonious","Discordant"))+
  scale_fill_gradient(low="mistyrose",high="firebrick",space = "Lab",limits=c(0,1),na.value="gray90",guide = "colourbar")+
  trans_theme+
  theme(legend.position ="none",plot.margin = margin(0.3,0,1.8,0,"cm"))
bp.y.data <- data.frame(pi=round(pi,2),x=c("Discordant","Harmonious","Indifferent"))
bp.y <- ggplot(data=bp.y.data,aes(x=x,y=pi))+
  geom_bar(stat="identity",aes(fill=pi)) +
  geom_text(aes(x=x,y=pi,label=pi),color="black",size=rel(3.8))+
  coord_flip() + theme_gray()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=14,margin=margin(0,10,0,0),angle = -90),
        legend.position = "none",
        plot.margin = margin(1.6,0,3.1,0,"cm"))+
  scale_fill_gradient(low="mistyrose",high="firebrick",limits=c(0,1))+
  scale_x_discrete(limits=c("Indifferent","Harmonious","Discordant"), position='top')+
  labs(x="Initial Probability")
grid.arrange(hm,hm1,bp.y,nrow=1)

