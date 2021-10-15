####--------Function to obtain the initial values for the constrained model-----######
## sample: the MCMC object from the mhmm model
## burnin: the number of samples to burn-in
## inits: the initial values for sample
## hs: the number of hidden states
## K: the number of variables in the overall model(the sample)
## Km: for which variables are searching for the fixed initial values

get_inits <- function(sample,burnin,inits,hs,K,Km){
  pi <- inits[sample$pidx$pi0]
  P0 <- inits[sample$pidx$P0]
  tau <- colMeans(sample$betaa[-(1:burnin),sample$pidx$tau[((Km[1]-1)*hs+1):(Km[2]*hs)]])
  beta <- rep(0,2*length(Km))
  Sigma <- colMeans(sample$betaa[-(1:burnin),sample$pidx$Sigma])
  Sigma <- matrix(Sigma,K,K)
  Sigma <- as.numeric(Sigma[Km,Km])
  invsigma <- colMeans(sample$betaa[-(1:burnin),sample$pidx$invsigma[((Km[1]-1)*hs+1):(Km[2]*hs)]])
  P1 <- inits[sample$pidx$P1]
  return(c(pi,P0,tau,beta,Sigma,invsigma,P1))
}

# ------ Step 1: read in data ------
data <- read.csv("../Fmod_imputed.csv",header = T)
data.na <- na.omit(data)
# The covariates: age and gender
X <- read.csv("../age_sex.csv",header=T)
# GROUP in A1C.csv indicates the study arm each family belongs to
A1C <- read.csv("../A1C.csv",header = T)
# The four main manifesto variables
y_c = list(t(matrix(data.na$CxDFR,nr = 5)),t(matrix(log(data.na$CxPCC),nr = 5)))
y_p = list(t(matrix(data.na$PxDFR,nr = 5)),t(matrix(log(data.na$PxPCC),nr = 5)))
#y_c <- lapply(y_c,function(x) (x-mean(x))/sd(x))
#y_p <- lapply(y_p,function(x) (x-mean(x))/sd(x))

# load in the functions to implement the method
source('../MHMM_v4.R')

# ------ Step 2: Set up initial values ------
# 2-class model
inits.pc_hs2 = c(0.547,0.353,0.014, 
                34.067,-2.671,3.4925,-0.3385,
                37.197,-3.9,3.274,-0.055,
	  	rep(0,8),
                13.062,0.032,9.663,-0.023,0.032,0.026,0.019,0.01,9.663,0.019,13.259,0.098,-0.023,0.01,0.098,0.029,
                rep(c(0.154,60.897,0.201,52.557),each=2),
                0.459,0.012)
# 3-class model
inits.pc_hs3 = c(0.056,0.815,
		0.223,0.057,0.027,0.39,0.005,0.011, 
                 34.102,-0.008,-2.636,3.791,-0.598,-0.04,
                 37.467,-0.54,-3.63,3.299,-0.05,-0.03,
		 rep(0,8),
 		 13.062,0.032,9.663,-0.023,0.032,0.026,0.019,0.01,9.663,0.019,13.259,0.098,-0.023,0.01,0.098,0.029,
                rep(c(0.154,60.897,0.201,52.557),each=3),
                 0.056,0.077,0.010,0.303,0.005,0.010)
# 4-class model
inits.pc_hs4 = c(0.139,0.408,0.408,
                 0.275,0.275,0.077,0.017,0.14,0.27,0.017,0.14,0.27,0.006,0.009,0.009, 
                 34.102,-0.004,-0.004,-2.636,3.791,-0.299,-0.299,-0.04,
                 37.467,-0.27,-0.27,-3.63,3.299,-0.025,-0.025,-0.03,
		 rep(0,8),
                 13.062,0.032,9.663,-0.023,0.032,0.026,0.019,0.01,9.663,0.019,13.259,0.098,-0.023,0.01,0.098,0.029,
                 rep(c(0.154,60.897,0.201,52.557),each=4),
                 0.356,0.356,0.103,0.014,0.374,0.259,0.014,0.374,0.259,0.006,0.006,0.006)

# 5-class model
inits.pc_hs5 = c(0.139,0.408,0.408,0.023,
                 0.275,0.275,0.039,0.039,0.017,0.14,0.135,0.135,0.017,0.14,0.135,0.135,0.006,0.009,0.009,0.489,0.006,0.009,0.009,0.489,
                 34.102,-0.004,-0.004,-1.318,-1.318,3.791,-0.299,-0.299,-0.02,-0.02,
                 37.467,-0.27,-0.27,-1.82,-1.81,3.299,-0.025,-0.025,-0.015,-0.015,
		 rep(0,8),
                 13.062,0.032,9.663,-0.023,0.032,0.026,0.019,0.01,9.663,0.019,13.259,0.098,-0.023,0.01,0.098,0.029,
                 rep(c(0.154,60.897,0.201,52.557),each=5),
                 0.356,0.356,0.052,0.051,0.014,0.374,0.129,0.130,0.014,0.374,0.13,0.129,0.006,0.006,0.006,0.491,0.006,0.006,0.006,0.491)
# set the number of hidden states as 2, 3, 4, or 5
hs = 3
print(paste("MHMM: The number of hidden states is:",hs))

# ------- Step 3: obtain posterior samples from MHMM model: ------
# The following function fitts MHMM model by setting fitRx as (F,T).
# Thus, transition probabilities are arm-specific.

# Detailed arguments to be found in "MHMM_v4.R"
# Only change "inits." when setting different number of hidden states
# 	e.g: When hs = 3 set inits. = inits.pc_hs3 in the following function
# Note: emission.fix is set to be F to allow different emission models for each family member

mhmm <- simulated_mhmm(y = c(y_c,y_p), X = as.matrix(X,ncol=2),inits. = inits.pc_hs3, report1. = 1000, burnin=30000,emission.fix = F,no.random = F,
                       nsim. = 40000, ksamp. = 1, N. = 390, ni. = rep(5,390), K = 4, hs = hs, rx. = A1C$GROUP-1,
                       fitRx = c(F,T), id = rep(1:390,5), hyperpar= c(5,1,1,1,1,.001, .0002), run. = 1)
save('mhmm','inits.pc_hs3',file = paste0("mhmm_pc_cov_hs",hs,".Rdata"))

Km = 1:2
inits.c <- get_inits(mhmm,30000,inits.pc_hs3,hs,4,Km)
mhmm_c <- simulated_mhmm(y = y_c, X = as.matrix(X,ncol=2),inits. = inits.c, report1. = 1000, burnin=30000,emission.fix = F,no.random = F,
                         nsim. = 40000, ksamp. = 1, N. = 390, ni. = rep(5,390), K = 2, hs = hs, rx. = A1C$GROUP-1,
                         fitRx = c(F,T), id = rep(1:390,5), hyperpar= c(3,1,1,.001, .0002), run. = 1)

save('mhmm_c','inits.c',file = paste0("mhmm_c_hs",hs,".Rdata"))

Km = 3:4
inits.p <- get_inits(mhmm,30000,inits.pc_hs3,hs,4,Km)
mhmm_p <- simulated_mhmm(y = y_p, X = as.matrix(X, ncol=2), inits. = inits.p, report1. = 1000, burnin=30000,
                         emission.fix = F,no.random = F, nsim. = 40000, ksamp. = 1, N. = 390, ni. = rep(5,390),
                         K = 2, hs = hs, rx. = A1C$GROUP-1,
                         fitRx = c(F,T), id = rep(1:390,5), hyperpar= c(3,1,1,.001, .0002), run. = 1)

save('mhmm_p','inits.p',file = paste0("mhmm_p_hs",hs,".Rdata"))
