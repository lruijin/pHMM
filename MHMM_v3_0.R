## Version 0.002
## This code is from Raffa JD and Dubin JA (2014) "Multivariate Longitudinal Data Analysis with Mixed Effects Hidden Markov Models."
## It contains two main functions to simulate multivariate response longitudinal data arises from a hidden Markov Model, 
## and implementation of a MCMC approach to estimate the parameters from such models.  For further details see:
## the original paper, supplementary material and example.R / README files.
## Jesse Raffa (jraffa@uw.edu): http://github.com/jraffa/MHMM
##
#######
#	genSim generates simulated MHMM dataset
#	code written for clarity, maybe not efficient
#	arguments:
# N: Number of study subjects (integer)
# K: Number of response variables
#	ni: Number of observation per subject (fixed for all subjects,, integer)
# hs: Number of hidden states (integer)
#	tau: means for each hidden state for responses (matrix of hs by p)
#	residinvvar: inverse of the residual variance (numeric)
#	reSigma: random effects covariance matrix (K x K matrix)
#	P0: hidden states transition probability matrix (hs x hs matrix) for placebo or baseline group
#	Pi0: initial probability vector (vector, length hs) for placebo or baseline group
#	w: random effects design matrix -- uses default as in paper (separate but correlated REs), requires R list of length 2, with each element containing a 2 x 1 row vector
#	rx: binary vector of treatment indicators of length N; if not specified assumes not treatment effect!
#	P1: hidden state tpm for treatment group (same as P0, requires rx to be specified); if left null, assumes no difference
#	Pi1: hidden state initial probability vector for treatment group (same as Pi0, requires rx to be specified); if left null, assumes no difference
#	fitRx: a logical vector of length two: fitRx[1] is whether to fit intitial probably treatment probabilities separately for Rx and Control; fitRx[2] fits separate tpms
#
# 	Returns as a list everything you need to start computation: 
#	y: generated responses as an array, 1 row per subject, 1 matrix per variable
#	N: Number of subjects, as specified
#	ni: number of observation times per subject
#	Z: generated hidden states (N x ni matrix) -- you don't have these normally
#	randomEffects: subject specific random effects (N x 3 matrix -- you don't have these normally either)
#	rx: treatment dummy variables
#	pars: parameter vector used to simulate data (useful for thinking specifying possible starting values)
#######
mhmmversion <- 0.003;
library(inline)
library(MCMCpack)
genSim_mhmm <- function(N = 354, K=2, ni=6, hs=3, tau=matrix(c(1.2,1.2,1.1,1.2,1.2,1.1),3,2),
                   residinvvar=c(10,10),reSigma=matrix(c(0.1,0.06,0.06,0.1),nr=2),
                   P0=matrix(c(0.85,0.6,0.1,0.1,0.2,0.1,0.05,0.2,0.8),nr=3),Pi0=c(0.6,0.25,0.15),
                   P1=matrix(c(0.85,0.6,0.1,0.1,0.2,0.1,0.05,0.2,0.8),nr=3),Pi1=c(0.6,0.25,0.15),
                   rx.=rep(c(0,1),each = N/2),fitRx=c(FALSE,FALSE),
                   w=list(matrix(c(1,0),nr=1),matrix(c(0,1),nr=1))){
  Res <- matrix(NA,nr=N,nc=dim(reSigma)[1]);
  reSigmachol <- chol(reSigma);
  Z <- matrix(NA,nr=N,nc=ni);
  
  y <- rep(list(matrix(NA, N, ni)), K)
  for(i in 1:N) {
    if(fitRx[1] & !is.null(rx.) & !is.null(Pi1)) {
      if(rx.[i]==1) {
        Z[i,1] <- sample(1:hs,1,replace=TRUE,Pi1); #intial hidden states
      } else {
        Z[i,1] <- sample(1:hs,1,replace=TRUE,Pi0);
      }
    } else {
      Z[i,1] <- sample(1:hs,1,replace=TRUE,Pi0); #intial hidden states
    }
    Res[i,] <- t(reSigmachol)%*%rnorm(dim(reSigma)[1])
    
    for(k in 1:K){
      tmpneta <- sum(tau[1:Z[i,1],k]) + w[[k]]%*%Res[i,];
      y[[k]][i,1]<- rnorm(1,tmpneta,sqrt(1/residinvvar[k]));
    }
  
    for(j in 2:ni) {
      if((fitRx[2] & !is.null(rx.) & !is.null(P1))) {
        if(rx.[i]==1) {
          Z[i,j] <- sample(1:hs,1,replace=TRUE,P1[Z[i,j-1],]);
        } else {
          Z[i,j] <- sample(1:hs,1,replace=TRUE,P0[Z[i,j-1],]);
        }
      } else {
        Z[i,j] <- sample(1:hs,1,replace=TRUE,P0[Z[i,j-1],]);
      }
      
      for(k in 1:K){
        tmpneta <- sum(tau[1:Z[i,j],k]) + w[[k]]%*%Res[i,];
        y[[k]][i,j] <- rnorm(1,tmpneta,sqrt(1/residinvvar[k]));
      }
     
    }
    
  }
  
  if(fitRx[1] & fitRx[2]) { #Both Rx Effects
    mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau,as.numeric(reSigma),residinvvar,as.numeric(t(P1))[-seq(hs^2,1,-(hs+1))],Pi1[1:(hs-1)]);
  }
  if(!fitRx[2] & fitRx[1]) { #Only Initial Prob Rx Effect
    mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau,as.numeric(reSigma),residinvvar,Pi1[1:(hs-1)]);
  }
  if(fitRx[2] & !fitRx[1]) { #Only TPM Rx Effect
    mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau,as.numeric(reSigma),residinvvar,as.numeric(t(P1))[-seq(hs^2,1,-(hs+1))]);
  }
  if(!fitRx[1] & !fitRx[2])  { #No Rx Effect (Default)
    mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau,as.numeric(reSigma),residinvvar);
  }
  return(list(y=y,N=N,ni=rep(ni,N),Z=Z,randomEffects=Res,rx=rx.,hs=hs,pars=mpars,fitRx=fitRx))
}

#######------------ Sample procedure for the hidden states, implemented by Rcpp ----------######
####### This step is to update the prbablity inividual by individual and finally take a log of the result ####
inccode_mhmm <- 'List  MHMMlabl(NumericMatrix y, NumericVector states, NumericVector lmb, NumericMatrix p, NumericVector pi, NumericVector ni) {
//Code adapted from Zucchini and MacDonald (2009) and Raffa and Dubin (2014)

NumericVector nt(ni); NumericVector lambda(lmb); NumericMatrix gamma(p);
NumericVector mt(states); NumericVector delta(pi);
int n = y.ncol(); // number of observations for each individual
int K = y.nrow(); // number of response variables
int m = lambda.size()/K - 1;
NumericMatrix lalpha(m,n);
NumericMatrix v(m,n);
IntegerVector Z(n);

NumericMatrix allprobs(n,m);
arma::mat gma = Rcpp::as<arma::mat>(gamma);
for(int q=0; q<n; q++) {
	for(int r=0; r<m; r++) {
	  double probtemp = ::Rf_dnorm4(y(0,q),lambda[r],lambda[m],false);
	  for(int k=1; k<K; k++){
	  probtemp = probtemp * ::Rf_dnorm4(y(k,q),lambda[k*(m+1) + r], lambda[(k+1)*(m+1)-1],false);
	  }
		allprobs(q,r) = probtemp; // the likelihood of ith sub data if it belongs to the rth class
	}
}
NumericMatrix foo(n,m);
NumericVector lscale(n);
foo.row(0) = delta*allprobs.row(0); // likelihood times the probability of different classes
NumericVector footmp(foo.row(0));
double sumfoo = std::accumulate(footmp.begin(),footmp.end(), 0.0); // scale the probability to make the proportion to the real probability
if(sumfoo==0 && n==1) {
	sumfoo = ::pow(10,-50); //Not usually necessary when using this code in MCMC, but if you have problems, may be worth checking.
	footmp[m-1] = ::pow(10,-50);
}
	
foo.row(0) = footmp/sumfoo;// update the first row with probability
lscale[0] = ::log(sumfoo); // log scale
NumericVector logfootmp(foo.row(0));
std::transform(logfootmp.begin(), logfootmp.end(), logfootmp.begin(), ::log); // log transforamtion
lalpha.column(0) = logfootmp+lscale[0]; 
v.column(0) = logfootmp+lscale[0];

if(n>1) {
  for(int i=1; i<n; i++) {
    NumericVector foa(foo.row(i-1));
    arma::vec lbetatmp(m);
    // iteration for posterior of each hidden state at time i
    for(int j = 0; j < m; j++){
      double maxi = 0;
      for(int k = 0; k < m; k++){
        maxi = std::max(maxi,v(k,i-1)+log(gma(k,j)));
      }
      v(j,i) = log(allprobs(i,j)) + maxi;
    }
    arma::colvec fooa = Rcpp::as<arma::colvec>(foa);
    NumericVector ttt(m);
    ttt = arma::trans(fooa)*gma;
    for(int j=0; j<m; j++) {
      foo(i,j) = ttt[j]*allprobs(i,j);
    }
    NumericVector footmp(foo.row(i));
    double sumfoo = std::accumulate(footmp.begin(),footmp.end(), 0.0);
    lscale[i] = lscale[i-1] + ::log(sumfoo);
    foo.row(i) = footmp/sumfoo;
    NumericVector logfootmp(foo.row(i));
    std::transform(logfootmp.begin(), logfootmp.end(), logfootmp.begin(), ::log);
    lalpha.column(i) = logfootmp+lscale[i];
  }
}
// trace back for the most likely hidden states based on posterior probabilities.
arma::vec tmpv = v.column(n-1);
Z[n-1] = arma::index_max(tmpv) + 1;

for(int i = n-2; i >= 0; i--){
  tmpv = v.column(i);
  Z[i] = arma::index_max(tmpv + log(gma.col(Z[i+1]-1))) + 1;
}


List ret; ret["lalpha"] = lalpha; ret["Z"] = Z;
return ret;
}
';

######-------------Sample procedure for the hidden states, implemented by Rcpp----------------#######
######
code5_mhmm <-'Environment base("package:base"); Function sample = base["sample"];
Rcpp::List y(yl);  // read in the list form of reponses
int K = y.length(); // number of response variables
NumericMatrix y1 = y[1]; // obtain the first matrix for parameter numbers
int N = y1.nrow(); //number of individuals
NumericVector beta(betaa); 
IntegerVector ni(nis); 
NumericMatrix con(ccc);

LogicalVector fitRx(R_fitRx);
NumericVector rx(rx1);
Rcpp::List pidx(paridx);

arma::mat co = Rcpp::as<arma::mat>(con);
IntegerVector states(m);
int me = states[0];
arma::mat pf(arma::zeros(me,me));

int ctt = 0;
IntegerVector ppidx = pidx["P0"];

for(int i=0; i<me; i++) {
	for(int j=0; j<me; j++) {
		if(i!=j) {
			pf(j,i) = beta[ppidx[ctt]-1];
			ctt++;
		}
  }
}

pf = pf.t();
for(int i=0; i<me; i++) {
		pf(i,i) = 1- arma::sum(pf.row(i));
}

NumericMatrix p(as<Rcpp::NumericMatrix>(wrap(pf)));
NumericMatrix pt(me,me);
if(!fitRx[1]) {
	pt = p;
} else {

	arma::mat pft(arma::zeros(me,me));
	int ctt = 0;
	IntegerVector ppidx1 = pidx["P1"];
	for(int i=0; i<me; i++) {
		for(int j=0; j<me; j++) {
			if(i!=j) {
				pft(j,i) = beta[ppidx1[ctt]-1];
				ctt++;
			}
	}
	}
	pft = pft.t();
	for(int i=0; i<me; i++) {
		pft(i,i) = 1- arma::sum(pft.row(i));
	}
	pt = (as<Rcpp::NumericMatrix>(wrap(pft)));
}

NumericVector pint(me);
IntegerVector pintidx = pidx["pi0"];
double pitotal = 0;
for(int i = 0; i<me-1; i++) {
  pint[i] = beta[pintidx[i]-1];
  pitotal = pitotal + pint[i];
}
pint[me-1] = 1-pitotal;
NumericVector pintt(me);
if(!fitRx[0]) {
  pintt = pint;
} else {
  IntegerVector pinttidx = pidx["pi1"];
  double pitotal = 0;
  for(int i = 0; i<me-1; i++) {
    pintt[i] = beta[pinttidx[i]-1];
    pitotal = pitotal + pintt[i];
  }
  pintt[me-1] = 1-pitotal;
}
NumericVector invsigidx = pidx["invsigma"];
NumericVector sdd(K); // use sd instead of inverse variance
for(int k=0;k<K;k++){
  sdd[k] = ::pow(beta[invsigidx[k]-1], -0.5);
}

IntegerVector tauidx = pidx["tau"];

NumericMatrix basenorm(me,K);
for(int k=0; k < K; k++){
  NumericVector normb(me);
  for(int h=0;h<me;h++){
    normb[h]= beta[tauidx[h + k*me] - 1 ];
  }
  arma::colvec nb = Rcpp::as<arma::colvec>(normb);// obtain the parameters for different states
  basenorm.column(k) = as<NumericVector>(wrap(co*nb));
}
NumericVector lmb(K*me+K);
for(int k=0; k<K; k++) {
  for(int i=0; i<me; i++){
    lmb[k*(me+1) + i] = basenorm(i,k);
  }
  lmb[(k+1)*(me+1)-1] = sdd[k];
}

IntegerMatrix Z(N,max(ni));
IntegerMatrix Zmax(N,max(ni));

for(int i = 0; i < N; i++){
  IntegerVector hsseq =	seq_len(me);
  NumericMatrix tmpy(K,ni[i]);
  for(int j=0; j<ni[i]; j++) {
    for(int k=0; k<K; k++){
      NumericMatrix ym = y[k];
      tmpy(k,j) = ym(i,j);
    }
  }
  
  NumericVector lmbtmp(lmb);
  IntegerVector repidx = pidx["re"];
  for(int k=0; k<K; k++) {
    for(int cf=0; cf<me; cf++){
      lmbtmp[cf+k*(me+1)] = basenorm(cf,k) + beta[repidx[k*N + i]-1];
    }
  }  
  
  NumericMatrix pii(me,me);
  NumericVector piii(me);
  if(rx[i]==1 && (fitRx[0] || fitRx[1])) {
    pii = pt;
    piii = pintt;
  } else {
    pii=p;
    piii=pint;
  }
  
  List out = MHMMlabl(tmpy,me,lmbtmp,pii,piii,ni[i]);
  NumericMatrix tmpla = out["lalpha"];
  IntegerVector tmpZmax = out["Z"];
  IntegerVector tmpZZ(ni[i]);
  NumericVector tmpprob(tmpla.column((ni[i]-1)));
  double tmpc = *std::max_element(tmpprob.begin(),tmpprob.end());
  tmpprob = tmpprob-tmpc;
  std::transform(tmpprob.begin(),tmpprob.end(),tmpprob.begin(),::exp);
  double tmpnorm = std::accumulate(tmpprob.begin(),tmpprob.end(), 0.0);
  tmpprob = tmpprob/tmpnorm;
  RNGScope scope;
  IntegerVector t1 = sample(hsseq, 1, false, tmpprob);
  tmpZZ[ni[i]-1] = t1[0];
  for(int qz = ni[i]-2; qz>=0; qz--) {
    NumericVector tmpprobi(tmpla.column(qz));
    double tmpci = *std::max_element(tmpprobi.begin(),tmpprobi.end());
    tmpprobi = tmpprobi-tmpci;
    std::transform(tmpprobi.begin(),tmpprobi.end(),tmpprobi.begin(),::exp);
    double tmpnormi1 = std::accumulate(tmpprobi.begin(),tmpprobi.end(), 0.0);
    tmpprobi = tmpprobi/tmpnormi1;
    NumericVector mtp(pii.column(tmpZZ[qz+1]-1));
    for(int qy = 0; qy<me; qy++) {
      tmpprobi[qy] = tmpprobi[qy]*mtp[qy];
    }
    double tmpnormi = std::accumulate(tmpprobi.begin(),tmpprobi.end(), 0.0);
    tmpprobi = tmpprobi/tmpnormi;
    IntegerVector tx = sample(hsseq, 1, false, tmpprobi);
    tmpZZ[qz] = tx[0];
  }
  Z.row(i) = tmpZZ;
  Zmax.row(i) = tmpZmax;
}
  List ret; ret["Z"] = Z; ret["Zmax"] = Zmax;
return ret;
';

pc.com_mhmm <- cxxfunction(signature( yl = "List", betaa="numeric",m="integer",nis="integer",ccc="matrix",
                                 rx1="integer",paridx="List",R_fitRx="logical"),code5_mhmm,plugin = "RcppArmadillo",
                      includes=c('#include <RcppArmadilloExtensions/sample.h>','#include <cmath>',inccode_mhmm))

#test.out<- pc.com(y,o.betaa,hs,ni.,ccc,rx.,pidx,fitRx)
#test.out


#This code is involved in sampling the random effects b_i

code_mhmm_re <- 'Environment mvtnorm("package:mvtnorm");
	Function rmvnorm = mvtnorm["rmvnorm"];
	Rcpp::List y(yl);  // read in the list form of reponses
	int K = y.length();
	NumericMatrix y1 = y[1];
	NumericMatrix ZZ(ZZZ);
	NumericVector obeta(b1);
	NumericVector nbeta(b2);
	IntegerVector nii(nis);
	arma::mat ZZa = Rcpp::as<arma::mat>(ZZ);
	Rcpp::List pidx(paridx);
	int N = y1.nrow();
	int maxni = y1.ncol();
	
	//obtain the current Sigma: covariance matrix of random effects K * K
	arma::mat sig(K,K);
	IntegerVector reSigidx = pidx["Sigma"];
	for(int i=0; i<K; i++){
	  for(int j=0; j<K; j++){
	    sig(i,j) = nbeta[reSigidx[i*K+j]-1];
	  }
	}
	
	arma::mat invsig = inv(sig);
	
	//obtain the measurement error 1/sigma^2_e
  IntegerVector errorSigidx = pidx["invsigma"];
  arma::mat ninvdiag = maxni * arma::diagmat(as<arma::vec>(obeta[errorSigidx-1]));
	
	// Obetain the new sigma matrix for updated random effect
  arma::mat newSigma = arma::inv(ninvdiag + invsig);
  
  //obtain the coefficients of fixed effects, tau, as a hs * K matrix, each as a column
  NumericVector zeros(K);
  IntegerVector states(m);
  int me = states[0];
  NumericMatrix tmpbetan(me,K);
  IntegerVector tauidx = pidx["tau"];
  for(int k=0; k<K; k++){
    for(int i=0; i<me; i++) {		
      tmpbetan(i,k) = obeta[tauidx[k*me+i]-1];
    }
  }
  arma::mat tmpbetann = Rcpp::as<arma::mat>(tmpbetan);

  arma::mat fixedn(N*maxni,K);
  for(int k=0; k<K; k++){
    fixedn.col(k) = ZZa*tmpbetann.col(k); // obtain the fixed effect for each observation.
  }
 
  NumericMatrix Res(N,K); // claim the matrix for storing the result, each individual takes one row.
  for (int i=0; i<N; i++) {
    arma::colvec tmpresid(K);
    for(int k=0; k<K; k++){
      NumericMatrix tmpy = y[k];
      tmpresid[k] = tmpy(i,0)- fixedn(i,k);
      for(int j=1; j<maxni; j++){
        tmpresid[k] += tmpy(i,j) - fixedn(j*N + i,k);
      }
    }
    arma::vec mean = newSigma * (ninvdiag/maxni) * tmpresid;
    NumericVector sample = rmvnorm(1,as<NumericVector>(wrap(mean)),as<NumericMatrix>(wrap(newSigma)));
    Res(i,_) = sample;
  }
	
	return Res;
';


library(mvtnorm);
re.n_mhmm <- cxxfunction(signature( yl = "List", ZZZ= 'matrix', b1 = "numeric",b2 = "numeric",
                               nis="numeric", paridx="list",m="integer"),
                    code_mhmm_re,
                    plugin = "RcppArmadillo");

#new.res <- re.n(y, ZZZ.r, n.betaa[1:basepars], n.betaa[1:basepars],ni., pidx, hs);
#new.res

code_waic <- '
	Rcpp::List y(yl);  // read in the list form of reponses
	int K = y.length();
	NumericMatrix y1 = y[1];
	NumericMatrix ZZ(ZZZ);
	NumericVector beta(b);
	arma::mat ZZa = Rcpp::as<arma::mat>(ZZ);
	Rcpp::List pidx(paridx);
	int N = y1.nrow();
	int maxni = y1.ncol();
	
	//obtain the measurement error 1/sigma^2_e
	NumericVector invsigidx = pidx["invsigma"];
	NumericVector sdd(K); // use sd instead of inverse variance
  for(int k=0;k<K;k++){
    sdd[k] = ::pow(beta[invsigidx[k]-1], -0.5);
  }
	
  //obtain the coefficients of fixed effects, tau, as a hs * K matrix, each as a column
  IntegerVector states(m);
  int me = states[0];
  NumericMatrix tmpbetan(me,K);
  IntegerVector tauidx = pidx["tau"];
  for(int k=0; k<K; k++){
    for(int i=0; i<me; i++) {		
      tmpbetan(i,k) = beta[tauidx[k*me+i]-1];
    }
  }
  
  arma::mat tmpbetann = Rcpp::as<arma::mat>(tmpbetan);

  arma::mat fixedn(N*maxni,K);
  for(int k=0; k<K; k++){
    fixedn.col(k) = ZZa*tmpbetann.col(k); // obtain the fixed effect for each observation.
  }
 
  IntegerVector repidx = pidx["re"];
  arma::vec Res = arma::ones(N*maxni);
  arma::vec logRes = arma::zeros(N*maxni);
  for(int i=0; i < N; i++){
    for(int k = 0; k < K; k++){
      NumericMatrix tmpy = y[k];
      for(int t = 0; t < maxni; t++){
        Res[t*N+i] = Res[t*N+i] * ::Rf_dnorm4(tmpy(i,t),fixedn(i+N*t,k)+beta[repidx[k*N+i]-1], sdd[k],false);
        logRes[t*N+i] = logRes[t*N+i] + ::Rf_dnorm4(tmpy(i,t),fixedn(i+N*t,k)+beta[repidx[k*N+i]-1], sdd[k],true);
      }
    }
  }
	
	List ret; ret["lik"] = Res; ret["llik"] = logRes;
	return ret;
';


waic <- cxxfunction(signature( yl = "List", ZZZ= 'matrix', b = "numeric",
                               paridx="list",m="integer"),
                    code_waic,
                    plugin = "RcppArmadillo");

#Optional, but may increase speed a little
library(compiler)
enableJIT(3);

###
# simulated: Simulate from the posterior of the MVHMM through Gibbs Sampling
# arguments:
# y : a list of length K (number of responses). Each element is a N * ni matrix (include NAs for missing data after dropout; this software will not currently handle intermittently missing data)
# inits.: initial values to start the simulation.  Can include/exclude the random effects, but the other parameter values must be specified
# nsim: how often to report progress of MCMC
# report1: how often to report progress of MCMC
# ksamp.: thinning parameter: keep every ksamp. samples
# N.: sample size (# of subjects)
# ni.: vector of length N of number of observation times for each subject
# K: number of response variables
# hs:  Number of hidden states in the HMM
# rx.: treatment/covariate vector of length N
# fitRx: logical vector of length two for initial probability vector, and tpm respectively.
# id: vector of length N * ni of ids for all observations; probably just rep(1:N.,max(ni))
# hyperpar: specify hyper-parameters for inverse-Wishart (indices 1-3) and gamma (4-5)
# run.: for keeping tracking chains run in parallel; not supported in this code.
#
# Returns:
# betaa: Matrix of Posterior Samples
# Z: Matrix of hidden state samples	
# pidx: List of Parameter indices
#
###

simulated_mhmm <- function(y,inits.,nsim.,report1.=1000,
                      ksamp.=1, N.,ni., K = 2, emission.fix = F, no.random = F, burnin = 3000,
                      hs=3,rx.=NULL,fitRx=c(FALSE,FALSE),id=rep(1:N.,6),hyperpar=c(3,1,1,.001,.0002),run.=1) {
  
  if(floor(nsim./ksamp.)!=ceiling(nsim./ksamp.)) {
    stop("nsim is not a multiple of ksamp");
  }
  if(sum(!(Reduce('+',is.na(y)) %in% c(0,K))) !=0) {
    stop("responses not matching with respect to missing data");
  }
  
  inc <- !is.na(as.numeric(y[[1]]));
  
  if(sum(fitRx)>0 & is.null(rx.)) {
    stop("No treatment vector passed to fit to");
  }
  
  #Parameter index
  pi0idx <- 1:(hs-1);
  P0idx <- (max(pi0idx)+1):(max(pi0idx)+hs*(hs-1));
  tauidx <- (max(P0idx)+1):(max(P0idx)+hs*K);
  Sigmaidx <- (max(tauidx)+1):(max(tauidx)+K*K);
  invsigmaidx <- (max(Sigmaidx)+1) : (max(Sigmaidx)+K);
  
  if(fitRx[1] & fitRx[2]) {
    P1idx<- (max(invsigmaidx)+1):(max(invsigmaidx)+hs*(hs-1));
    pi1idx <- (max(P1idx)+1):(max(P1idx)+hs-1);
    reidx <- (max(pi1idx)+1):(max(pi1idx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,pi1=pi1idx,re=reidx);# I stopped here.
  } else if(!fitRx[1] & fitRx[2]) {
    P1idx<- (max(invsigmaidx)+1):(max(invsigmaidx)+hs*(hs-1));
    reidx <- (max(P1idx)+1):(max(P1idx)+N.*K)
    pidx <- list(pi0=pi0idx,P0=P0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,re=reidx);
  } else if(fitRx[1] & !fitRx[2]) {
    pi1idx <- (max(invsigmaidx)+1):(max(invsigmaidx)+hs-1);
    reidx <- (max(pi1idx)+1):(max(pi1idx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,pi1=pi1idx,re=reidx);
  } else {
    reidx <- (max(invsigmaidx)+1):(max(invsigmaidx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,re=reidx);
  }
  
  #matrix to hold parameter/REs from MCMC ouput
  betaa <- matrix(NA,nr=nsim./ksamp.,nc=max(pidx$re));
  
  if(!fitRx[1] & !fitRx[2]) {
    basepars <- max(pidx$invsigma);
    rx. <- rep(0,N.);
  } 
  if(fitRx[1] & fitRx[2]) {
    basepars <- max(pidx$pi1)
  }
  if(fitRx[1] & !fitRx[2]) {
    basepars <- max(pidx$pi1)
  }
  if(!fitRx[1] & fitRx[2]) {
    basepars <- max(pidx$P1)
  }
  
  betaa[1,1:basepars] <- inits.[1:basepars];
  if(length(inits.)==(basepars+N.*K)) { #If init values for random effects are passed, use them, otherwise start at zero
    
    betaa[1,c(pidx$re)] <- inits.[pidx$re]
  } else {
    betaa[1,(basepars+1):(basepars+K*N.)] <-0;
  }
  
  o.betaa <- n.betaa_tmp <- betaa[1,]; #Start and save the initial fixed truth 
  n.betaa <- rep(NA,length(betaa[1,]));
  a.param <- hyperpar[2+K];
  b.param <- hyperpar[3+K];

  y.na <- matrix(NA,nr=K,nc=sum(inc))
  for(k in 1:K){
    y.na[k,] <- y[[k]][inc]
  }
  id.na <- id[inc];
  ccc <- matrix(NA,nr=hs,nc=hs);
  for(q in 1:hs) {
    ccc[q,] <- c(rep(1,q),rep(0,hs-q));
  }
  Zmat <- matrix(NA,nr=nsim./ksamp.,nc=max(ni.)*N.)
  Zmaxmat <- matrix(NA,nr=nsim./ksamp.,nc=max(ni.)*N.)
  for(i in 1:nsim.) {
    #Simulate Sigma, following inverse Wishart distribution with hyperpar[1] as degree of freedom
    if(emission.fix){
      n.betaa[Sigmaidx] <- o.betaa[Sigmaidx]
    }else{
      ress <- t(matrix(o.betaa[reidx],nc=K))
      cvvv <- ress%*%t(ress);
      signew <- riwish(N.+hyperpar[1],diag(hyperpar[2:(2+K-1)]) + cvvv);
      n.betaa[Sigmaidx] <- as.numeric(c((signew + t(signew))/2))
    }
    
    #simulated hidden states (Z)
    zz.ret <- pc.com_mhmm(y,o.betaa,hs,ni.,ccc,rx.,pidx,fitRx)
    zz2 <- zz.ret$Z
    zzmax <- zz.ret$Zmax
    zz2[zz2==0] <- NA;
    zz.na <- na.omit(as.numeric(zz2));
    
    #Simulate P and Pi
    if(fitRx[2]) {
      zz2.rx <- zz2[rx.==1,];
      zzt.rx <- table(factor(zz2.rx[,-ncol(zz2.rx)],levels=as.character(c(1:hs))),factor(zz2.rx[,-1],levels=as.character(c(1:hs))));
      zz2.px <- zz2[rx.==0,];
      zzt.px <- table(factor(zz2.px[,-ncol(zz2.px)],levels=as.character(c(1:hs))),factor(zz2.px[,-1],levels=as.character(c(1:hs))));
      ctt <- 1;
      for(ps in 1:hs) {
        n.betaa[P0idx[ctt:(ctt+hs-2)]] <- rdirichlet(1,zzt.px[ps,]+rep(1,hs))[-ps];
        n.betaa[P1idx[ctt:(ctt+hs-2)]] <- rdirichlet(1,zzt.rx[ps,]+rep(1,hs))[-ps];
        ctt <- ctt + hs-1;
      }
    } else {
      zzt <- table(factor(zz2[,-ncol(zz2)],levels=as.character(c(1:hs))),factor(zz2[,-1],levels=as.character(c(1:hs))));
      ctt <- 1;
      for(ps in 1:hs) {
        n.betaa[P0idx[ctt:(ctt+hs-2)]] <- rdirichlet(1,zzt[ps,]+rep(1,hs))[-ps];
        ctt <- ctt + hs-1;
      }
    }
    if(fitRx[1]) {
      zz2.rx <- zz2[rx.==1,];
      zz2.px <- zz2[rx.==0,];
      n.betaa[pi0idx] <- rdirichlet(1,as.numeric(table(factor(zz2.px[,1],levels=as.character(c(1:hs))))+1))[1:(hs-1)];
      n.betaa[pi1idx] <- rdirichlet(1,as.numeric(table(factor(zz2.rx[,1],levels=as.character(c(1:hs))))+1))[1:(hs-1)];
    } else {
      n.betaa[pi0idx] <- rdirichlet(1,as.numeric(table(factor(zz2[,1],levels=as.character(c(1:hs))))+1))[1:(hs-1)];
    }
    ZZZ <- matrix(NA,nr=length(zz.na),nc=hs);
    ZZZ[,1] <- 1;
    for(s in 2:hs) {
      ZZZ[,s] <- as.numeric(1*zz.na>(s-1));
    }
     
    # For testing defining the true states
    # zz <- c(daf$Z)
    # ZZZ <- matrix(NA,nr=length(zz),nc=hs);
    # ZZZ[,1] <- 1;
    # for(s in 2:hs) {
    #   ZZZ[,s] <- as.numeric(1*zz>(s-1));
    # }
      
    #Simulate random effect
    if(no.random){n.betaa[reidx] = 0}else{
      new.res <- re.n_mhmm(y, ZZZ, o.betaa[1:basepars], n.betaa[1:basepars],ni., pidx, hs)
      n.betaa[reidx] = as.vector(new.res)
    }
    
    #Simualte tau
    for(k in 1:K){
      if(emission.fix){
        n.betaa[tauidx[(k-1)*hs + (1:hs)]] <- o.betaa[tauidx[(k-1)*hs + (1:hs)]]
      }else{
        bbcov <- solve(crossprod(ZZZ))
        bbmean <- bbcov%*%t(ZZZ)%*%(y.na[k,]-n.betaa[reidx[(k-1)*N.+id.na]])
        bbcov <- 1/o.betaa[invsigmaidx[k]]*bbcov;
        bbcov <- (bbcov + t(bbcov))/2
        n.betaa[tauidx[(k-1)*hs + (1:hs)]] <- rmvnorm(1,bbmean,bbcov);
      }
    }
    
    
    #Simulate 1/sigma^2_e
    if(emission.fix){
      for(k in 1:K){
        n.betaa[invsigmaidx[k]] <- o.betaa[invsigmaidx[k]]
      }
    }else{
      for(k in 1:K){
        n.betaa[invsigmaidx[k]] <- rgamma(1,(sum(ni.)/2 + a.param),
                                          b.param + sum((y.na[k,]-n.betaa[reidx[(k-1)*N.+id.na]] - ZZZ%*%n.betaa[tauidx[(k-1)*hs + (1:hs)]])^2)/2); #tau
    }}
    #test[i,] <- c(n.betaa[P0idx],n.betaa[pi0idx],n.betaa[Sigmaidx],n.betaa[tauidx],n.betaa[invsigmaidx])
    if(i > burnin){
      llik_s = waic(y,ZZZ,n.betaa,pidx,hs)
      if(i == burnin + 1){
        lik = unlist(llik_s['lik'])
        llik = unlist(llik_s['llik'])
      }else{
        lik = rbind(lik, unlist(llik_s['lik']))
        llik = rbind(llik,unlist(llik_s['llik']))
      }
    }
    
    if(floor(i/ksamp.)==ceiling(i/ksamp.)) {
      betaa[(i-1)/ksamp.+1,] <- n.betaa;
      Zmat[(i-1)/ksamp.+1,] <- as.vector(zz2);
      Zmaxmat[(i-1)/ksamp.+1,] <- as.vector(zzmax);
    }
    o.betaa <- n.betaa;
    rm(n.betaa);
    n.betaa <- rep(NA,length(o.betaa));
    if(ceiling(i/report1.)==floor(i/report1.)) {
      message("iteration: ", i, " ", date(), " Run: ", run.);
    }
  } 
  WAIC = -2 * (sum(log(colMeans(lik))) - sum(apply(llik,2,var)))
  return(list(betaa=betaa,Z=Zmat,Zmax = Zmaxmat,pidx=pidx,lik = lik, llik = llik,WAIC = WAIC));
}




