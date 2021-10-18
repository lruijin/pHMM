library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

# This file includes several functions for plotting
# They are used in the files summarizing posterior samples.

# This function plots the trace plot of fixed effects for the four manifesto variables
plot_tau_m <- function(obj,hs,mar=0.5,cex=0.8,ylabs=c("CDFR","CPCC","PDFR","PPCC")){
  par(mfrow=c(2,2))
  for(i in 1:4){
    tautmp <- obj$betaa[,obj$pidx$tau[(i-1)*hs + (1:hs)]]
    tautmp <- t(apply(tautmp,1,cumsum))
    ylim1 <- min(tautmp) - mar
    ylim2 <- max(tautmp) + mar
    plot(tautmp[,1],type="l",ylim=c(ylim1,ylim2),ylab=ylabs[i],xlab="iteration")
    for(j in 2:hs){
      lines(tautmp[,j],type="l",col=j)
    }
    legend("bottomright",legend = c("Discordant","Harmonious","Indifferent"), col=1:hs, lty=1, cex = cex,border = F, bg="transparent")
  }
}
# This function plots the trace plots of fixed effects for the four manifesto variables
# It is the ggplot2 version, which is used in the paper
plot_tau_m_gg <- function(obj,hs,mar=0.5,cex=0.8,ylabs=c("CDFR","log(CPCC)","PDFR","log(PPCC)"),
                          tit=c("CDFR","CPCC","PDFR","PPCC")){
  for(i in 1:4){
    tautmp <- obj$betaa[,obj$pidx$tau[(i-1)*hs + (1:hs)]]
    tautmp <- t(apply(tautmp,1,cumsum))
    ylim1 <- min(tautmp) - mar
    ylim2 <- max(tautmp) + mar
    plot_data <- data.frame(y = c(tautmp), x = rep(1:nrow(tautmp),3), 
                            Class = rep(c("Discordant","Harmonious","Indifferent"), each = nrow(tautmp)))
    assign(paste0("p",i),ggline(plot_data,x = "x", y = "y", color= "Class",plot_type = "l",palette = "jco")+
      labs(x="Iteration",y = ""))
  }
  ggarrange(p1,p2,p3,p4,labels = tit,common.legend = T, vjust = -0.5, legend = "bottom")+
    theme(plot.margin = margin(t=0.8,r = 0, b = 0.2, l =0, unit = "cm"))
}
# This function plots the trace plots of the fixed effects of the four manifesto variables
# by marginalizing the family members.
plot_tau_g <- function(obj,hs,M,Km=c(2,2),mar=0.5,cex = 0.8,ylabs=c("CDFR","CPCC","PDFR","PPCC")){
  tau_f <- matrix(NA,nrow(obj$betaa),length(obj$pidx$tau))
  ctt_pm = 0
  ctt_tau = 0
  for(m in 1:M){
    for(i in 1:hs){
      idx <- ctt_pm + seq(i,((hs-2)*hs+i),by=hs)
      tmpPm <- cbind(obj$betaa[,obj$pidx$Pm0[idx]],1-rowSums(obj$betaa[,obj$pidx$Pm0[idx]]))
      for(k in 1:Km[m]){
        tau_f[,ctt_tau + (k-1)*hs+i] <- rowSums(tmpPm * t(apply(obj$betaa[,obj$pidx$tau[ctt_tau + (k-1)*hs + (1:hs)]],1,cumsum)))
      }
    }
    ctt_pm = ctt_pm + (hs-1) * hs
    ctt_tau = ctt_tau + Km[m] * hs
  }
  par(mfrow=c(2,2))
  for(i in 1:sum(Km)){
    tautmp <- tau_f[,(i-1)*hs + (1:hs)]
    ylim1 <- min(tautmp) - mar
    ylim2 <- max(tautmp) + mar
    plot(tautmp[,1],type="l",ylim=c(ylim1,ylim2),ylab=ylabs[i],xlab="iteration")
    for(j in 2:hs){
      lines(tautmp[,j],type="l",col=j)
    }
    legend("bottomright",legend = c("Discordant","Harmonious","Indifferent"), col=1:hs, lty=1, border = F,cex=cex, bg="transparent")
  }
}
# This function plots the trace plots of the fixed effects of the four manifesto variables
# by marginalizing the family members.
# This is the ggplot version
plot_tau_g_gg <- function(obj,hs,M,Km=c(2,2),mar=0.5,cex = 0.8,
                          ylabs=c("CDFR","log(CPCC)","PDFR","log(PPCC)"),
                          tit=c("CDFR","CPCC","PDFR","PPCC")){
  tau_f <- matrix(NA,nrow(obj$betaa),length(obj$pidx$tau))
  ctt_pm = 0
  ctt_tau = 0
  for(m in 1:M){
    for(i in 1:hs){
      idx <- ctt_pm + seq(i,((hs-2)*hs+i),by=hs)
      tmpPm <- cbind(obj$betaa[,obj$pidx$Pm0[idx]],1-rowSums(obj$betaa[,obj$pidx$Pm0[idx]]))
      for(k in 1:Km[m]){
        tau_f[,ctt_tau + (k-1)*hs+i] <- rowSums(tmpPm * t(apply(obj$betaa[,obj$pidx$tau[ctt_tau + (k-1)*hs + (1:hs)]],1,cumsum)))
      }
    }
    ctt_pm = ctt_pm + (hs-1) * hs
    ctt_tau = ctt_tau + Km[m] * hs
  }
  for(i in 1:sum(Km)){
    tautmp <- tau_f[,(i-1)*hs + (1:hs)]
    plot_data <- data.frame(y = c(tautmp), x = rep(1:nrow(tautmp),3), 
                            Class = rep(c("Discordant","Harmonious","Indifferent"), each = nrow(tautmp)))
    assign(paste0("p",i),ggline(plot_data,x = "x", y = "y", color= "Class",plot_type = "l",palette = "jco")+
             labs(x="Iteration",y = ""))
  }
  ggarrange(p1,p2,p3,p4,labels = tit,common.legend = T, vjust = -0.5, legend = "bottom")+
    theme(plot.margin = margin(t=0.8,r = 0, b = 0.2, l =0, unit = "cm"))
}

plot_prob <- function(P0mat,P1mat,pi,hs){
trans_theme <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text=element_text(color="black",size=11),
                     axis.title=element_text(size=14))
trans_plot0 <-data.frame(P0mat,From = c("Discordant","Harmonious","Indifferent"))
names(trans_plot0)[1:hs] <- c("Discordant","Harmonious","Indifferent")
trans_plot0.melted <- melt(trans_plot0)
hm <- ggplot(data = trans_plot0.melted,aes(x=factor(variable),y=factor(From),fill=value)) + geom_tile() +
  geom_text(aes(label=value),color="black",size=rel(3.8))+
  ylab("From") + xlab("To")+ labs(fill="Probability")+
  ggtitle("Usual Care")+
  scale_x_discrete()+scale_y_discrete(limits=c("Indifferent","Harmonious","Discordant"))+
  scale_fill_gradient(low="mistyrose",high="firebrick",space = "Lab",limits=c(0,1),na.value="gray90",guide = "colourbar")+
  trans_theme+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title = element_text(size=10),
        legend.text = element_text(size=8))+
  guides(fill=guide_colourbar(barwidth = unit(6,"cm")))
trans_plot1 <-data.frame(P1mat,From = c("Discordant","Harmonious","Indifferent"))
names(trans_plot1)[1:hs] <- c("Discordant","Harmonious","Indifferent")
trans_plot1.melted <- melt(trans_plot1)
hm1 <- ggplot(data = trans_plot1.melted,aes(x=factor(variable),y=factor(From),fill=value)) + geom_tile() +
  geom_text(aes(label=value),color="black",size=rel(3.8))+
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
}

plot_prob_arbi <- function(P0mat,P1mat,pi,hs){
  trans_theme <- theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.text=element_text(color="black",size=11),
                       axis.title=element_text(size=14))
  trans_plot0 <-data.frame(P0mat,From = 1:hs)
  names(trans_plot0)[1:hs] <- 1:hs
  trans_plot0.melted <- melt(trans_plot0)
  hm <- ggplot(data = trans_plot0.melted,aes(x=factor(variable),y=factor(From),fill=value)) + geom_tile() +
    geom_text(aes(label=value),color="black",size=rel(3.8))+
    ylab("From") + xlab("To")+ labs(fill="Probability")+
    ggtitle("Usual Care")+
    scale_x_discrete()+scale_y_discrete(limits=as.character(hs:1))+
    scale_fill_gradient(low="mistyrose",high="firebrick",space = "Lab",limits=c(0,1),na.value="gray90",guide = "colourbar")+
    trans_theme+
    theme(legend.position = "bottom", legend.direction = "horizontal",legend.title = element_text(size=10),
          legend.text = element_text(size=8))+
    guides(fill=guide_colourbar(barwidth = unit(6,"cm")))
  trans_plot1 <-data.frame(P1mat,From = 1:hs)
  names(trans_plot1)[1:hs] <- 1:hs
  trans_plot1.melted <- melt(trans_plot1)
  hm1 <- ggplot(data = trans_plot1.melted,aes(x=factor(variable),y=factor(From),fill=value)) + geom_tile() +
    geom_text(aes(label=value),color="black",size=rel(3.8))+
    ylab("") + xlab("To")+ 
    ggtitle("Intervention")+
    scale_x_discrete()+scale_y_discrete(limits=as.character(hs:1))+
    scale_fill_gradient(low="mistyrose",high="firebrick",space = "Lab",limits=c(0,1),na.value="gray90",guide = "colourbar")+
    trans_theme+
    theme(legend.position ="none",plot.margin = margin(0.3,0,1.8,0,"cm"))
  bp.y.data <- data.frame(pi=round(pi,2),x=1:hs)
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
    scale_x_discrete(limits=as.character(hs:1), position='top')+
    labs(x="Initial Probability")
  grid.arrange(hm,hm1,bp.y,nrow=1)
}
