## Version 0.002
## This code is for the pHMM
##
#######
#	genSim generates simulated pHMM dataset
#	code written for clarity, maybe not efficient
#
#	arguments:
#       N: integer, Number of study groups. e.g.: number of families etc.
# 	M: integer, Number of members in a group, e.g: number of family members is 2.
# 	K: integer vector with M elements, each element is the number of responses from each member.
#	ni: integer, fixed for all groups. Number of observation per group
# 	hs: integer vector with (M+1) elements. Number of hidden states (group-level hidden states, member-level hidden states)
#	tau: a list of M elements, each is a hs[m+1] by K_m matrix. means for each hidden state for responses
#	residinvvar: a list of M elements, each is a each is a hs[m+1] by K_m matrix, inverse of the variance for each hidden state of responses
#	reSigma: \sum K_m x \sum K_m matrix. random effects covariance matrix.
#	P0: hs[1] x hs[1] matrix. hidden states transition probability matrix for placebo or baseline group
#	Pi0: hs[1]-dimensional vector, initial probability vector for placebo or baseline group
#	Pm : a list of M elements, each of which is a hs[1] x hs[m+1] matrix. Perception matrices.
# 	w: \sum K_m by \sum K_m matrix, diagnomal matrix by default, random effects design matrix -- uses default as in paper (separate but correlated REs).
#	rx.: N-d binary vector, treatment indicators of length N; if not specified assumes all come from baseline group
#	P1: hidden state tpm for treatment group (same as P0, requires rx to be specified); if left null, assumes no difference
#	Pi1: hidden state initial probability vector for treatment group (same as Pi0, requires rx to be specified); if left null, assumes no difference
#	fitRx: a logical vector of length two: fitRx[1] is whether to fit intitial probably treatment probabilities separately for Rx and Control; fitRx[2] fits separate tpms
#
# 	Returns as a list everything you need to start computation:
#	y: generated responses as a list of length \sum K_m, each element is an N by ni matrix
#	N: Number of groups, as specified
#	ni: number of observation times per group
#	Z: generated hidden states (N x ni matrix) -- you don't have these normally
# 	Zm : generated hidden states of each member, a list of M and each element is an N x ni matrix
#	randomEffects: subject specific random effects (N x \sum K_m matrix -- you don't have these normally either)
#	rx: treatment dummy variables
#	pars: parameter vector used to simulate data (useful for thinking specifying possible starting values)
#######
dhmmversion <- 0.002;
library(inline)
library(MCMCpack)
genSim <- function(N = 400, M = 2, K = c(2,2), ni = 5, hs = c(3,3,3),
                   tau=list(matrix(c(34.082,0.412,-3.669,3.366,-0.143,-0.095),hs[2],K[1]),matrix(c(38.430,-1.922,-3.227, 3.318,-0.097,-0.063),hs[3],K[2])),
                   residinvvar=list(matrix(c(0.048, 0.436, 0.443, 8.911, 75.614, 145.159),nr=hs[2]),matrix(c(0.139,0.333,0.252,19.069,135.208,62.988),nr=hs[3])),
                   reSigma=matrix(c(13.425,0.051,0.684,-0.027,0.051,0.025,0.033,0.010,0.684,0.033,13.190,0.086,-0.027,0.010,0.086,0.028),nr=sum(K)),
                   P0=matrix(c(0.574,0.029,0.010,0.379,0.465,0.014,0.047,0.506,0.976),nr=3),Pi0=c(0.495,0.437,0.068),
                   Pm0=list(matrix(c(0.531,0.051,0.025,0.433,0.744,0.253,0.036,0.205,0.722),3,3),matrix(c(0.740,0.021,0.015,0.238,0.942,0.023,0.022,0.037,0.962),3,3)),
                   P1=matrix(c(0.487,0.026,0.020,0.465,0.651,0.022,0.048,0.323,0.958),nr=3),Pi1=c(0.495,0.437,0.068),
                   Pm1=list(matrix(c(0.827,0.133,0.000,0.039,0.487,0.009,0.134,0.380,0.991),3,3),matrix(c(0.887,0.001,0.000,0.028,0.970,0.031,0.085,0.092,0.),3,3)),
                   rx.=rep(0,N),fitRx=c(FALSE,FALSE,FALSE), w=diag(sum(K))){
  Res <- matrix(NA,nr=N,nc=dim(reSigma)[1]);
  reSigmachol <- chol(reSigma);
  Z <- matrix(NA,nr=N,nc=ni);
  Zm <- list()
  length(Zm) = M
  for(m in 1:M){
    Zm[[m]] <- matrix(NA,N,ni)
  }

  y <- rep(list(matrix(NA, N, ni)), sum(K))
  for(i in 1:N) {
    if(fitRx[1] & !is.null(rx.) & !is.null(Pi1)) {
      if(rx.[i]==1) {
        Z[i,1] <- sample(1:hs[1],1,replace=TRUE,Pi1); #intial hidden states
      } else {
        Z[i,1] <- sample(1:hs[1],1,replace=TRUE,Pi0);
      }
    } else {
      Z[i,1] <- sample(1:hs[1],1,replace=TRUE,Pi0); #intial hidden states
    }

    Res[i,] <- t(reSigmachol)%*%rnorm(dim(reSigma)[1])
    # Generate hidden states of each member based on the group hidden states
    ctt <- 1
    for(m in 1:M){
      if(fitRx[3] & !is.null(rx.) & !is.null(Pm1)) {
        if(rx.[i]==1) {
          Zm[[m]][i,1] <- sample(1:hs[1+m],1,replace=T,Pm1[[m]][Z[i,1],])
        }else{
          Zm[[m]][i,1] <- sample(1:hs[1+m],1,replace=T,Pm0[[m]][Z[i,1],])
        }
      }else{
        Zm[[m]][i,1] <- sample(1:hs[1+m],1,replace=T,Pm0[[m]][Z[i,1],])
      }
      for(k in 1:K[m]){
        tmpneta <- sum(tau[[m]][1:Zm[[m]][i,1],k]) + w[ctt,]%*%Res[i,];
        y[[ctt]][i,1]<- rnorm(1,tmpneta,sqrt(1/residinvvar[[m]][Zm[[m]][i,1],k]));
        ctt <- ctt + 1
      }
    }

    for(j in 2:ni) {
      if((fitRx[2] & !is.null(rx.) & !is.null(P1))) {
        if(rx.[i]==1) {
          Z[i,j] <- sample(1:hs[1],1,replace=TRUE,P1[Z[i,j-1],]);
        } else {
          Z[i,j] <- sample(1:hs[1],1,replace=TRUE,P0[Z[i,j-1],]);
        }
      } else {
        Z[i,j] <- sample(1:hs[1],1,replace=TRUE,P0[Z[i,j-1],]);
      }
      ctt <- 1
      for(m in 1: M){
        if(fitRx[3] & !is.null(rx.) & !is.null(Pm1)) {
          if(rx.[i]==1) {
            Zm[[m]][i,j] <- sample(1:(hs[1+m]),1,replace=T,Pm1[[m]][Z[i,j],])
          }else{
            Zm[[m]][i,j] <- sample(1:(hs[1+m]),1,replace=T,Pm0[[m]][Z[i,j],])
          }
        }else{
          Zm[[m]][i,j] <- sample(1:(hs[1+m]),1,replace=T,Pm0[[m]][Z[i,j],])
        }
        for(k in 1:K[m]){
          tmpneta <- sum(tau[[m]][1:Zm[[m]][i,j],k]) + w[ctt,]%*%Res[i,];
          y[[ctt]][i,j] <- rnorm(1,tmpneta,sqrt(1/residinvvar[[m]][Zm[[m]][i,j],k]));
          ctt <- ctt + 1
        }
      }

    }

  }

  Pm0.out <- as.numeric(Pm0[[1]][,-hs[2]])
  for(m in 2:M){
    Pm0.out <- c(Pm0.out, as.numeric(Pm0[[m]][,-hs[m+1]]))
  }
  
  Pm1.out <- as.numeric(Pm1[[1]][,-hs[2]])
  for(m in 2:M){
    Pm1.out <- c(Pm1.out, as.numeric(Pm1[[m]][,-hs[m+1]]))
  }

  if(fitRx[1] & fitRx[2] & fitRx[3]) { #All Rx Effects
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar),as.numeric(t(P1))[-seq(hs[1]^2,1,-(hs[1]+1))],Pi1[1:(hs[1]-1)],Pm1.out);
  }
  if(!fitRx[2] & fitRx[1] & !fitRx[3]) { #Only Initial Prob
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar),Pi1[1:(hs[1]-1)]);
  }
  if(fitRx[2] & !fitRx[1] & !fitRx[3]) { #Only TPM Rx Effect
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar),as.numeric(t(P1))[-seq(hs[1]^2,1,-(hs[1]+1))]);
  }
  if(!fitRx[2] & !fitRx[1] & fitRx[3]) { #Only translation Rx Effect
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar),Pm1.out);
  }
  if(fitRx[2] & fitRx[1] & !fitRx[3]) { #Initial Prob and TPM Rx Effect
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar),as.numeric(t(P1))[-seq(hs[1]^2,1,-(hs[1]+1))],Pi1[1:(hs[1]-1)]);
  }
  if(fitRx[2] & !fitRx[1] & fitRx[3]) { #TPM and translation Rx Effect
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar),as.numeric(t(P1))[-seq(hs[1]^2,1,-(hs[1]+1))],Pm1.out);
  }
  if(!fitRx[2] & fitRx[1] & fitRx[3]) { #Initial Prob and translation Rx Effect
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar),Pi1[1:(hs[1]-1)],Pm1.out);
  }
  if(!fitRx[1] & !fitRx[2] & !fitRx[3])  { #No Rx Effect (Default)
    mpars <- c(Pi0[1:(hs[1]-1)],as.numeric(t(P0))[-seq(hs[1]^2,1,-(hs[1]+1))], Pm0.out, unlist(tau), as.numeric(reSigma), unlist(residinvvar));
  }
  return(list(y=y,N=N,ni=rep(ni,N),Z=Z,Zm=Zm,randomEffects=Res,rx=rx.,hs=hs,pars=mpars,fitRx=fitRx))
}

#######------------ Sample procedure for the hidden states of groups, and the members,implemented by Rcpp ----------######
####### This step is to update the prbablity inividual by individual and finally take a log of the result ####
inccode <- 'List  HMMlabl(NumericMatrix y, IntegerVector Kms, int ghs, NumericVector lmb, NumericMatrix p, NumericVector pi, NumericVector ni, List pm) {
//Code adapted from Zucchini and MacDonald (2009) and Raffa and Dubin (2014)
int M = pm.length();
NumericVector nt(ni); NumericVector lambda(lmb); NumericMatrix gamma(p);
NumericVector delta(pi);
NumericVector Km(Kms);
int n = y.ncol(); // number of observations for each individual
NumericMatrix lalpha(ghs,n);
NumericMatrix v(ghs,n);
IntegerVector Z(n);

NumericMatrix allprobs(n,ghs);

arma::mat gma = Rcpp::as<arma::mat>(gamma);


for(int q = 0; q < n; q++){
  for(int r = 0; r < ghs; r++){
    int ctt = 0;
    double probfin = 1.0;
    int lmbs = 0; // lmb start: from which the lambda starts for the current probability
    for(int m=0; m<M;m++){
      NumericMatrix tmppm = pm[m];
      int hsm = tmppm.ncol();
      arma::vec probtmp = arma::ones(hsm);
      for(int i = 0; i < hsm; i++){
        for(int k = 0; k< Km[m]; k++){
          probtmp[i] *= ::Rf_dnorm4(y(ctt+k,q),lambda[lmbs+i+k*hsm*2],lambda[lmbs + 2*hsm*k+hsm+i],false);
        }
        probtmp[i] *= tmppm(r,i);
      }
      ctt = ctt + Km[m];
      probfin *= sum(probtmp);
      lmbs += Km[m] * hsm * 2;
    }
    allprobs(q,r) = probfin;
  }
}

NumericMatrix foo(n,ghs);
NumericVector lscale(n);
foo.row(0) = delta*allprobs.row(0); // likelihood times the probability of different classes
NumericVector footmp(foo.row(0));
double sumfoo = std::accumulate(footmp.begin(),footmp.end(), 0.0); // scale the probability to make the proportion to the real probability
if(sumfoo==0 && n==1) {
  sumfoo = ::pow(10,-50); //Not usually necessary when using this code in MCMC, but if you have problems, may be worth checking.
  footmp[ghs-1] = ::pow(10,-50);
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
    arma::vec lbetatmp(ghs);
    // iteration for posterior of each hidden state at time i
    for(int j = 0; j < ghs; j++){
      double maxi = 0;
      for(int k = 0; k < ghs; k++){
        maxi = std::max(maxi,v(k,i-1)+log(gma(k,j)));
      }
      v(j,i) = log(allprobs(i,j)) + maxi;
    }
    arma::colvec fooa = Rcpp::as<arma::colvec>(foa);
    NumericVector ttt(ghs);
    ttt = arma::trans(fooa)*gma;
    for(int j=0; j<ghs; j++) {
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

//List ret; ret["allprobs"] = allprobs; ret["foo"] = foo; ret["lscale"] = lscale; ret["lalpha"] = lalpha;
List ret; ret["lalpha"] = lalpha; ret["Z"] = Z;
return ret;
}
'

inccode_m <- 'List  mHMMlabl(NumericMatrix y, IntegerVector Kms, IntegerVector ZZ, IntegerVector ZZmax, NumericVector lmb, List pm) {
//Code adapted from Zucchini and MacDonald (2009) and Raffa and Dubin (2014)
Environment base("package:base"); Function sample = base["sample"];
int M = pm.length();
NumericVector lambda(lmb); IntegerVector Z(ZZ); IntegerVector Zmax(ZZmax);
NumericVector Km(Kms);
int n = y.ncol(); // number of observations for each individual

IntegerMatrix tmpZm(M,n);
IntegerMatrix tmpZmmax(M,n);

for(int q=0; q<n; q++){
  int tmpZ = Z[q];
  int tmpZmax = Zmax[q];
  int ctt = 0;
  int lmbs = 0; // lmb start: from which the lambda starts for the current probability
  for(int m=0; m<M; m++){
    NumericMatrix tmppm = pm[m];
    int hsm = tmppm.ncol();
    IntegerVector hsseq =	seq_len(hsm);
    arma::vec probtmp = arma::ones(hsm);
    arma::vec probtmpmax = arma::ones(hsm);
    for(int k=0;k<Km[m];k++){
      for(int i=0;i<hsm;i++){
        probtmp[i] *= ::Rf_dnorm4(y(ctt,q),lambda[lmbs+i],lambda[lmbs + hsm + i],false);
      }
      lmbs += hsm * 2;
      ctt ++;
    }
    for(int i=0;i<hsm;i++){
      probtmp[i] = probtmp[i] * tmppm(tmpZ-1,i);
      probtmpmax[i] = probtmp[i] * tmppm(tmpZmax-1, i);
    }
    probtmp = probtmp / sum(probtmp);
    probtmpmax = probtmpmax / sum(probtmpmax);
    IntegerVector t = sample(hsseq,1,false,probtmp);
    tmpZm(m,q) = t[0];

    tmpZmmax(m,q) = arma::index_max(probtmpmax) + 1;
  }
}
List ret; ret["Zm"] = tmpZm; ret["Zmmax"] = tmpZmmax;

return ret;
}
'

######-------------Sample procedure for the hidden states, implemented by Rcpp----------------#######
######
code5 <-'Environment base("package:base"); Function sample = base["sample"];
Rcpp::List y(yl);  // read in the list form of reponses
NumericMatrix y1 = y[0]; // obtain the first matrix for parameter numbers
int N = y1.nrow(); //number of individuals

NumericVector beta(betaa);
IntegerVector ni(nis);
Rcpp::List con(ccc);

LogicalVector fitRx(R_fitRx);
NumericVector rx(rx1);
Rcpp::List pidx(paridx);

IntegerVector states(m);
IntegerVector Km(Kms);
int me = states[0];
int M = states.length()-1;
int K = sum(Km);
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
  ctt = 0;
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

Rcpp::List pm(M);
IntegerVector pmidx = pidx["Pm0"];
ctt = 0;
for(int m=0; m<M; m++){
  arma::mat tmppm(arma::zeros(me, states[m+1]));
  for(int j=0; j<states[m+1]-1; j++){
    for(int i=0; i<me; i++){
      tmppm(i,j) = beta[pmidx[ctt]-1];
      ctt ++;
    }
  }
  for(int i=0; i<me; i++){
    tmppm(i, states[m+1]-1) = 1-arma::sum(tmppm.row(i));
  }
  pm[m] = as<Rcpp::NumericMatrix>(wrap(tmppm));
}

Rcpp::List pmt(M);
if(!fitRx[2]){
  pmt = pm;
}else{
  IntegerVector pm1idx = pidx["Pm1"];
  ctt = 0;
  for(int m=0; m<M; m++){
    arma::mat tmppm(arma::zeros(me, states[m+1]));
    for(int j=0; j<states[m+1]-1; j++){
      for(int i=0; i<me; i++){
       tmppm(i,j) = beta[pm1idx[ctt]-1];
       ctt ++;
      }
   }
   for(int i=0; i<me; i++){
     tmppm(i, states[m+1]-1) = 1-arma::sum(tmppm.row(i));
   }
    pmt[m] = as<Rcpp::NumericMatrix>(wrap(tmppm));
  }
}

IntegerVector invsigidx = pidx["invsigma"];

IntegerVector tauidx = pidx["tau"];

int lmb_length = 0;
for(int m =0;m < M; m++){
  lmb_length += 2* states[m+1] * Km[m];
}
NumericVector lmb(lmb_length);

ctt = 0;
int lmbidx = 0;
for(int m=0; m<M; m++){
  arma::mat com = Rcpp::as<arma::mat>(con[m]);
  int hsm = states[m+1];
  NumericVector normb(hsm);
  NumericVector sddb(hsm);
  for(int k=0; k<Km[m];k++){
    for(int h=0; h<hsm; h++){
      normb[h] = beta[tauidx[ctt+h]-1];
      sddb[h] =  ::pow(beta[invsigidx[ctt+h]-1], -0.5);
    }
    for(int h=hsm-1; h > 0; h--){
      sddb[h] = sddb[h] - sddb[h-1];
    }
    arma::colvec nb = Rcpp::as<arma::colvec>(normb);// obtain the parameters for different states
    arma::colvec sdd = Rcpp::as<arma::colvec>(sddb);
    NumericVector tmpnb = as<NumericVector>(wrap(com*nb));
    NumericVector tmpsdd = as<NumericVector>(wrap(com*sdd));
    for(int h=0; h<hsm; h++){
      lmb[lmbidx+h] = tmpnb[h];
      lmb[lmbidx+h+hsm] = tmpsdd[h];
    }
    ctt += hsm;
    lmbidx += 2 * hsm;
  }
}

// Initialize the matrix to obtain the labels for each group and each member
IntegerMatrix Z(N,max(ni));
IntegerMatrix Zmax(N,max(ni));
arma::cube Zm(N,max(ni),M);
arma::cube Zmmax(N,max(ni),M);

for(int i = 0;i<N;i++){
  IntegerVector hsseq =	seq_len(me);
  NumericMatrix tmpy(K,ni[i]);
  for(int j=0; j<ni[i]; j++) {
    for(int k=0; k<K; k++){
      NumericMatrix ym = y[k];
      tmpy(k,j) = ym(i,j);
    }
  }
  NumericVector lmbtmp(lmb.length());
  IntegerVector repidx = pidx["re"];
  int ctt = 0;
  int ctk = 0;
  for(int m=0; m<M; m++){
    for(int k=0; k<Km[m]; k++){
      for(int cf =0; cf<states[m+1]; cf++){
        lmbtmp[ctt+cf] = lmb[ctt+cf] + beta[repidx[ctk*N + i]-1];
        lmbtmp[ctt+cf+states[m+1]] = lmb[ctt+cf+states[m+1]];
      }
      ctt = ctt + 2*states[m+1];
      ctk ++;
    }
  }
  NumericMatrix pii(me,me);
  NumericVector piii(me);
  Rcpp::List pmi(M);
  if(rx[i]==1 && (fitRx[0] || fitRx[1] || fitRx[2])) {
    pii = pt;
    piii = pintt;
    pmi = pmt;
  } else {
    pii=p;
    piii=pint;
    pmi = pm;
  }
  
  List out = HMMlabl(tmpy,Km,me,lmbtmp,pii,piii,ni[i],pmi);
  
  NumericMatrix tmpla = out["lalpha"];
  IntegerVector tmpZmax = out["Z"];
  IntegerVector tmpZZ(ni[i]);
  NumericVector tmpprob(tmpla.column((ni[i]-1)));
  double tmpc = *std::max_element(tmpprob.begin(),tmpprob.end());
  tmpprob = tmpprob-tmpc;
  std::transform(tmpprob.begin(),tmpprob.end(),tmpprob.begin(),::exp);
  double tmpnorm = std::accumulate(tmpprob.begin(),tmpprob.end(), 0.0);
  tmpprob = tmpprob/tmpnorm;
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
  List ret = mHMMlabl(tmpy,Km,tmpZZ,tmpZmax,lmbtmp,pmi);
  IntegerMatrix tmpZm = ret["Zm"];
  IntegerMatrix tmpZmmax = ret["Zmmax"];
  for(int m = 0; m< M; m++){
    for(int t=0; t<ni[i];t++){
      Zm(i,t,m) = tmpZm(m,t);
      Zmmax(i,t,m) = tmpZmmax(m,t);
    }
  }
}


//List ret; ret["out"] = out; ret["lmb"]=lmb;ret["lmbtmp"] = lmbtmp; ret["pm"] = pm; ret["tmpy"] = tmpy;//ret["tmpZZ"]=tmpZZ;
List ret; ret["Z"] = Z;ret["Zm"] = Zm; ret["Zmax"] = Zmax; ret["Zmmax"] = Zmmax;
//List rett; rett["out"] = out;rett["ret"] =ret;
//List ret; ret["lmb"] = lmb; ret["ctt"] = ctt;
return ret;
';

pc.com <- cxxfunction(signature( yl = "List", betaa="numeric",m="integer",Kms = "integer", nis="integer",ccc="List",
                                 rx1="integer",paridx="List",R_fitRx="logical"),code5,plugin = "RcppArmadillo",
                      includes=c('#include <cmath>',inccode,inccode_m))

#test.out <- pc.com(y,o.betaa,hs,Km,ni.,ccc,rx.,pidx,fitRx)
#test.out$Pm
#sum(test.out$Z-daf$Z !=0)

code <- 'Environment mvtnorm("package:mvtnorm");
	Function rmvnorm = mvtnorm["rmvnorm"];
	Rcpp::List y(yl);  // read in the list form of reponses
	NumericVector Km(Kms);
	int M = Km.length();
	int K = sum(Km);
	Rcpp::List ZZ(ZZl);
	NumericVector obeta(b1);
	NumericVector nbeta(b2);
	IntegerVector nii(nis);
	Rcpp::List pidx(paridx);

	NumericMatrix y1 = y[1];
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

	
  IntegerVector errorSigidx = pidx["invsigma"];

  //obtain the fixed effects, as an N*maxni by K matrix
  IntegerVector states(m);
  IntegerVector tauidx = pidx["tau"];
  arma::mat fixedn(N*maxni,K);
  arma::mat error(N*maxni,K);
  
  int ctt = 0;
  int ctk = 0;
  for(int m=0; m<M; m++){
    arma::colvec tmpbetan = arma::zeros(states[m+1]);
    arma::colvec tmpsigma = arma::zeros(states[m+1]);
    arma::mat ZZa = as<arma::mat>(ZZ[m]);
    for(int k=0; k<Km[m];k++){
      for(int i=0; i<states[m+1]; i++){
        tmpbetan[i] = obeta[tauidx[ctt+i]-1];
        tmpsigma[i] = obeta[errorSigidx(ctt+i)-1];
      }
      for(int i=states[m+1]-1; i>0; i--){
        tmpsigma[i] = tmpsigma[i] - tmpsigma[i-1];
      }
      ctt += states[m+1];
      fixedn.col(ctk) = ZZa * tmpbetan;
      error.col(ctk) = ZZa* tmpsigma;
      ctk ++;
    }
  }

  NumericMatrix Res(N,K); // claim the matrix for storing the result, each individual takes one row.
  for (int i=0; i<N; i++) {
    //int i = 1;
    arma::colvec tmpresid(K);
    arma::vec tmpmean(K);
    arma::mat tmpSigma = arma::zeros(K,K);
    arma::mat newSigma = arma::zeros(K,K);
    for(int k=0; k<K; k++){
      NumericMatrix tmpy = y[k];
      tmpresid[k] = tmpy(i,0)- fixedn(i,k);
      tmpSigma(k,k) = error(i,k);
    }
    newSigma = tmpSigma;
    tmpmean = tmpSigma * tmpresid;
    
    for(int j=1; j<maxni; j++){
      for(int k = 0; k< K; k++){
        NumericMatrix tmpy = y[k];
        tmpSigma(k,k) = error(i + j*N, k);
        tmpresid[k] = tmpy(i,j) - fixedn(j*N + i,k);
      }
      newSigma = newSigma + tmpSigma;
      tmpmean = tmpmean + tmpSigma * tmpresid;
    }
    
    newSigma = inv(newSigma + invsig);
    arma::vec mean = newSigma * tmpmean;
    NumericVector sample = rmvnorm(1,as<NumericVector>(wrap(mean)),as<NumericMatrix>(wrap(newSigma)));
    Res(i,_) = sample;
  }
    //List ret; ret["newSigma"] = newSigma;ret["invsig"] = invsig;ret["fixedn"] = fixedn; ret["error"] = error;
    //ret["mean"] = mean; ret["sample"] = sample; ret["tmpresid"] = tmpresid;
    //return ret;
  return Res;
';

library(mvtnorm);
re.n <- cxxfunction(signature( yl = "List", Kms = "integer", ZZl= 'List', b1 = "numeric",b2 = "numeric",
                               nis="numeric", paridx="list",m="integer"),
                    code,
                    plugin = "RcppArmadillo")
#new.res <- re.n(y, Km, ZZm, o.betaa[1:basepars], n.betaa[1:basepars],ni., pidx, hs);
code_waic2 <- '
	Rcpp::List y(yl);  // read in the list form of reponses
	NumericVector Km(Kms);
	int M = Km.length();
	int K = sum(Km);
	
	Rcpp::List ZZ(ZZl);
	NumericVector beta(b);
	Rcpp::List pidx(paridx);
	
	NumericMatrix y1 = y[1];
	int N = y1.nrow();
	int maxni = y1.ncol();
	
	//obtain the measurement error 1/sigma^2_e
	NumericVector invsigidx = pidx["invsigma"];
	
  //obtain the coefficients of fixed effects, tau, as a hs * K matrix, each as a column
  IntegerVector states(m);
  IntegerVector tauidx = pidx["tau"];
  
  arma::mat fixedn(N*maxni, K);
  arma::mat error(N*maxni, K);
  int ctt = 0;
  int fixedn_ctt = 0;
  for(int m = 0; m < M; m++){
    int me = states[m+1];
    arma::mat ZZa = as<arma::mat>(ZZ[m]);
    arma::mat tmpbetan(me,Km[m]);
    arma::mat tmpsdd(me,Km[m]);
    for(int k=0; k<Km[m]; k++){
     for(int i=0; i<me; i++) {		
      tmpbetan(i,k) = beta[tauidx[ctt+k*me+i]-1];
      tmpsdd(i,k) = ::pow(beta[invsigidx[ctt + k*me+i]-1], -0.5);
     }
    }
    for(int i=me-1; i > 0; i--){
      tmpsdd.row(i) = tmpsdd.row(i) - tmpsdd.row(i-1);
    }
    for(int k=0; k<Km[m]; k++){
      fixedn.col(fixedn_ctt + k) = ZZa*tmpbetan.col(k); // obtain the fixed effect for each observation.
      error.col(fixedn_ctt + k) = ZZa*tmpsdd.col(k);
    }
    ctt = ctt+Km[m]*me;
    fixedn_ctt = fixedn_ctt+Km[m];
  }
 
  IntegerVector repidx = pidx["re"];
  arma::vec Res = arma::ones(N*maxni);
  arma::vec logRes = arma::zeros(N*maxni);
  for(int i=0; i < N; i++){
    for(int k = 0; k < K; k++){
      NumericMatrix tmpy = y[k];
      for(int t = 0; t < maxni; t++){
        Res[t*N+i] = Res[t*N+i] * ::Rf_dnorm4(tmpy(i,t),fixedn(i+N*t,k)+beta[repidx[k*N+i]-1], error(i+N*t,k),false);
        logRes[t*N+i] = logRes[t*N+i] + ::Rf_dnorm4(tmpy(i,t),fixedn(i+N*t,k)+beta[repidx[k*N+i]-1], error(i+N*t,k),true);
      }
   }
  }
	
	List ret;ret["lik"] = Res; ret["llik"] = logRes;
	return ret;
';


waic2 <- cxxfunction(signature( yl = "List", Kms = "interger", ZZl= 'List', b = "numeric",
                               paridx="list",m="integer"),
                    code_waic2,
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
#
###

simulated <- function(y,inits.,nsim.,report1.=1000,burnin = 3000, no.random=F,
                      ksamp.=1, N.,ni., Km = c(2,2),
                      hs=c(3,3,3), rx.=NULL,fitRx=c(FALSE,FALSE,FALSE),id=rep(1:N.,6),hyperpar=c(5,1,1,1,1,.001,.0002),run.=1) {

  if(floor(nsim./ksamp.)!=ceiling(nsim./ksamp.)) {
    stop("nsim is not a multiple of ksamp");
  }
  K <- sum(Km)
  M <- length(Km)
  if(length(Km) != length(hs)-1){
    stop(paste0("There are ", M, " members in each group. Length of Km and hs should be consistent with this"))
  }
  if(sum(!(Reduce('+',is.na(y)) %in% c(0,K))) !=0) {
    stop("responses not matching with respect to missing data");
  }

  inc <- !is.na(as.numeric(y[[1]]));

  if(sum(fitRx)>0 & is.null(rx.)) {
    stop("No treatment vector passed to fit to");
  }

  #Parameter index
  pi0idx <- 1:(hs[1]-1);
  P0idx <- (max(pi0idx) + 1):(max(pi0idx)+hs[1]*(hs[1]-1));
  Pm0idx <- (max(P0idx) + 1):(max(P0idx) + sum(hs[1] * (hs[-1]-1)))
  tauidx <- (max(Pm0idx)+1):(max(Pm0idx)+sum(Km * hs[-1]))
  Sigmaidx <- (max(tauidx)+1):(max(tauidx)+K*K);
  invsigmaidx <- (max(Sigmaidx)+1) : (max(Sigmaidx)+sum(Km*hs[-1]));

  if(fitRx[1] & fitRx[2] & fitRx[3]) {
    P1idx<- (max(invsigmaidx)+1):(max(invsigmaidx)+hs[1]*(hs[1]-1));
    pi1idx <- (max(P1idx)+1):(max(P1idx)+hs[1]-1);
    Pm1idx <- (max(pi1idx)+1):(max(Pi1dx) + sum(hs[1] * (hs[-1]-1)))
    reidx <- (max(Pm1idx)+1):(max(Pm1idx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,pi1=pi1idx,Pm1 = Pm1idx, re=reidx);
    basepars <- max(pidx$Pm1)
  } else if(!fitRx[1] & fitRx[2] & fitRx[3]) {
    P1idx<- (max(invsigmaidx)+1):(max(invsigmaidx)+hs[1]*(hs[1]-1));
    Pm1idx <- (max(P1idx) + 1) : (max(P1idx) + sum(hs[1] * (hs[-1]-1)))
    reidx <- (max(Pm1idx)+1):(max(Pm1idx)+N.*K)
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx, tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,Pm1=Pm1idx,re=reidx);
    basepars <- max(pidx$Pm1)
  } else if(fitRx[1] & !fitRx[2] & fitRx[3]) {
    pi1idx <- (max(invsigmaidx)+1):(max(invsigmaidx)+hs[1]-1);
    Pm1idx <- (max(pi1idx) + 1) : (max(pi1idx) + sum(hs[1] * (hs[-1]-1)))
    reidx <- (max(Pm1idx)+1):(max(Pm1idx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,pi1=pi1idx,Pm1=Pm1idx,re=reidx);
    basepars <- max(pidx$Pm1)
  } else if(fitRx[1] & fitRx[2] & !fitRx[3]){
    P1idx<- (max(invsigmaidx)+1):(max(invsigmaidx)+hs[1]*(hs[1]-1));
    pi1idx <- (max(P1idx)+1):(max(P1idx)+hs[1]-1);
    reidx <- (max(pi1idx)+1):(max(pi1idx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,pi1=pi1idx,re=reidx);
    basepars <- max(pidx$pi1)
  } else if(fitRx[1] & !fitRx[2] & !fitRx[3]){
    pi1idx <- (max(invsigmaidx)+1):(max(invsigmaidx)+hs[1]-1);
    reidx <- (max(pi1idx)+1):(max(pi1idx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,pi1=pi1idx,re=reidx);
    basepars <- max(pidx$pi1)
  } else if(!fitRx[1] & fitRx[2] & !fitRx[3]){
    P1idx<- (max(invsigmaidx)+1):(max(invsigmaidx)+hs[1]*(hs[1]-1));
    reidx <- (max(P1idx)+1):(max(P1idx)+N.*K)
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx, tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,re=reidx);
    basepars <- max(pidx$P1)
  } else if(!fitRx[1] & !fitRx[2] & fitRx[3]){
    Pm1idx<- (max(invsigmaidx)+1):(max(invsigmaidx)+sum(hs[1] * (hs[-1]-1)));
    reidx <- (max(Pm1idx)+1):(max(Pm1idx)+N.*K)
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx, tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=Pm1idx,re=reidx);
    basepars <- max(pidx$P1)
  } else {
    reidx <- (max(invsigmaidx)+1):(max(invsigmaidx)+N.*K);
    pidx <- list(pi0=pi0idx,P0=P0idx,Pm0 = Pm0idx,tau=tauidx,Sigma=Sigmaidx,invsigma=invsigmaidx,re=reidx);
    rx. <- rep(0,N.)
    basepars <- max(pidx$invsigma)
  }

  #matrix to hold parameter/REs from MCMC ouput
  betaa <- matrix(NA,nr=nsim./ksamp.,nc=max(pidx$re))
  betaa[1,1:basepars] <- inits.[1:basepars];

  if(length(inits.)==(basepars+N.*K)) { #If init values for random effects are passed, use them, otherwise start at zero
    betaa[1,c(pidx$re)] <- inits.[pidx$re]
  } else {
    betaa[1,(basepars+1):(basepars+K*N.)] <-0;
  }

  o.betaa <- betaa[1,]; #Start
  n.betaa <- rep(NA,length(betaa[1,]));
  a.param <- hyperpar[2+K];
  b.param <- hyperpar[3+K];

  y.na <- matrix(NA,nr=K,nc=sum(inc))
  for(k in 1:K){
    y.na[k,] <- y[[k]][inc]
  }
  id.na <- id[inc];
  ccc <- list()
  length(ccc) <- M
  for(m in 1:M){
    ccc[[m]] <- matrix(NA,nr=hs[m+1],nc=hs[m+1])
    for(q in 1:hs[m+1]){
      ccc[[m]][q,] <- c(rep(1,q),rep(0,hs[m+1]-q))
    }
  }

  Zmat <- Zmaxmat <- matrix(NA,nr=nsim./ksamp.,nc=max(ni.)*N.)
  Zmmat <- Zmmaxmat <- list()
  length(Zmmat) <- M
  for(m in 1:M){
    Zmmat[[m]] <- Zmmaxmat[[m]] <- matrix(NA,nr=nsim./ksamp.,nc=max(ni.)*N.)
  }

  for(i in 1:nsim.) {
    #Simulate Sigma, following inverse Wishart distribution with hyperpar[1] as degree of freedom
    if(no.random){
      n.betaa[Sigmaidx] <- 0
    }else{
      ress <- t(matrix(o.betaa[reidx],nc=K)) #daf$randomEffects #new.res
      cvvv <- ress%*%t(ress);
      signew <- riwish(N.+hyperpar[1],diag(hyperpar[2:(2+K-1)]) + cvvv);
      n.betaa[Sigmaidx] <- as.numeric(c((signew + t(signew))/2))
    }
    #simulated family hidden states (Z)
    zz.ret <- pc.com(y,o.betaa,hs,Km,ni.,ccc,rx.,pidx,fitRx)
    zz2 <- zz.ret$Z
    zz2max <- zz.ret$Zmax
    zz2[zz2==0] <- NA;
    zz2max[zz2max==0] <- NA
    zz.na <- na.omit(as.numeric(zz2));

    #Simulate P and Pi
    if(fitRx[2]) {
      zz2.rx <- zz2[rx.==1,];
      zzt.rx <- table(factor(zz2.rx[,-ncol(zz2.rx)],levels=as.character(c(1:hs[1]))),factor(zz2.rx[,-1],levels=as.character(c(1:hs[1]))));
      zz2.px <- zz2[rx.==0,];
      zzt.px <- table(factor(zz2.px[,-ncol(zz2.px)],levels=as.character(c(1:hs[1]))),factor(zz2.px[,-1],levels=as.character(c(1:hs[1]))));
      ctt <- 1;
      for(ps in 1:hs[1]) {
        n.betaa[P0idx[ctt:(ctt+hs[1]-2)]] <- rdirichlet(1,zzt.px[ps,]+rep(1,hs[1]))[-ps];
        n.betaa[P1idx[ctt:(ctt+hs[1]-2)]] <- rdirichlet(1,zzt.rx[ps,]+rep(1,hs[1]))[-ps];
        ctt <- ctt + hs[1]-1;
      }
    } else {
      zzt <- table(factor(zz2[,-ncol(zz2)],levels=as.character(c(1:hs[1]))),factor(zz2[,-1],levels=as.character(c(1:hs[1]))));
      ctt <- 1;
      for(ps in 1:hs[1]) {
        n.betaa[P0idx[ctt:(ctt+hs[1]-2)]] <- rdirichlet(1,zzt[ps,]+rep(1,hs[1]))[-ps];
        ctt <- ctt + hs[1]-1;
      }
    }

    if(fitRx[1]) {
      zz2.rx <- zz2[rx.==1,];
      zz2.px <- zz2[rx.==0,];
      n.betaa[pi0idx] <- rdirichlet(1,as.numeric(table(factor(zz2.px[,1],levels=as.character(c(1:hs[1]))))+1))[1:(hs[1]-1)];
      n.betaa[pi1idx] <- rdirichlet(1,as.numeric(table(factor(zz2.rx[,1],levels=as.character(c(1:hs[1]))))+1))[1:(hs[1]-1)];
    } else {
      n.betaa[pi0idx] <- rdirichlet(1,as.numeric(table(factor(zz2[,1],levels=as.character(c(1:hs[1]))))+1))[1:(hs[1]-1)];
    }

    #simulated member hidden states
    zzm <- zz.ret$Zm
    zzm.na <- matrix(NA,M,length(zz.na))

    for(m in 1:M){
      zzm.na[m,] <- as.vector(zzm[,,m])
    }

    ZZm <- list()
    for(m in 1:M){
      tmpZZm <- matrix(NA,nr=ncol(zzm.na),nc = hs[m+1])
      tmpZZm[,1] <- 1
      for(s in 2:hs[m+1]){
        tmpZZm[,s] <- as.numeric(1*zzm.na[m,] > (s-1))
      }
      ZZm[[m]] <- tmpZZm
    }

    # Simulate Pm
    ctt <- 0
    for(m in 1:M){
      if(fitRx[3]){
        zz2.rx <- zz2[rx. == 1]
        zz2.px <- zz2[rx. == 0]
        zzm.rx <- zzm[rx. == 1,,m]
        zzm.px <- zzm[rx. == 0,,m]
        zzt.rx <- table(factor(zz2.rx,levels=as.character(c(1:hs[1]))),factor(zzm.rx,levels=as.character(c(1:hs[m+1]))))
        zzt.px <- table(factor(zz2.px,levels=as.character(c(1:hs[1]))),factor(zzm.px,levels=as.character(c(1:hs[m+1]))))
        tmppm.rx <-tmppm.px <- matrix(NA,hs[1],hs[m+1])
        for(ps in 1:hs[1]){
          tmppm.rx[ps,] <- rdirichlet(1,zzt.rx[ps,] + rep(1,hs[m+1]))
          tmppm.px[ps,] <- rdirichlet(1,zzt.px[ps,] + rep(1,hs[m+1]))
        }
        n.betaa[Pm0idx[ctt + (1:((hs[m+1]-1)*hs[1]))]]<- as.numeric(tmppm.px[,-hs[m+1]])
        n.betaa[Pm1idx[ctt + (1:((hs[m+1]-1)*hs[1]))]]<- as.numeric(tmppm.rx[,-hs[m+1]])
        ctt <- (hs[m+1]-1)*hs[1]
      }else{
        zzt <- table(factor(zz2,levels=as.character(c(1:hs[1]))),factor(zzm.na[m,],levels=as.character(c(1:hs[m+1]))))
        tmppm <- matrix(NA,hs[1],hs[m+1])
        for(ps in 1:hs[1]){
          tmppm[ps,] <- rdirichlet(1,zzt[ps,]+rep(1,hs[m+1]))
        }
        n.betaa[Pm0idx[ctt + (1:((hs[m+1]-1)*hs[1]))]]<- as.numeric(tmppm[,-hs[m+1]])
        ctt <- (hs[m+1]-1)*hs[1]
      }
    }

    #Simulate random effect
    if(no.random){
      n.betaa[reidx] <- 0
    }else{
      new.res <- re.n(y, Km, ZZm, o.betaa[1:basepars], n.betaa[1:basepars],ni., pidx, hs);
      n.betaa[reidx] = as.vector(new.res)
    }
    #Simualte tau
    ctt <- 1
    ctidx <- 0
    for(m in 1:M){
      ZZZ <- ZZm[[m]]
      me <- hs[m+1]
      for(k in 1:Km[m]){
        invsigma <- o.betaa[invsigmaidx[ctidx + (1:me)]]
        invsigma <- c(invsigma[1],diff(invsigma))
        invsigma <- diag(c(ZZZ%*% matrix(invsigma,nc=1)))
        bbcov <- solve(t(ZZZ) %*% invsigma %*% ZZZ)
        bbmean <- bbcov %*% t(ZZZ) %*% invsigma %*% (y.na[ctt,]-n.betaa[reidx[(ctt-1)*N.+id.na]]);
        bbcov <- (bbcov + t(bbcov))/2
        n.betaa[tauidx[ctidx + (1:me)]] <- rmvnorm(1,bbmean,bbcov);
        ctt <- ctt+1;
        ctidx <- ctidx+me
      }
    }

    #Simulate 1/sigma^2_e
    ctt <- 1
    ctidx <- 0
    for(m in 1:M){
      ZZZ <- ZZm[[m]]
      me <- hs[m+1]
      for(k in 1:Km[m]){
        tau.k <- (y.na[ctt,]-n.betaa[reidx[(ctt-1)*N.+id.na]] - ZZZ%*%n.betaa[tauidx[ctidx + (1:me)]])^2
        for(s in 1:me){
          idx <- which(zzm.na[m,] == s)
          rec <- length(idx)
          n.betaa[invsigmaidx[ctidx + s]] <- rgamma(1,(rec/2 + a.param),
                                                       b.param + sum(tau.k[idx])/2)
        }
        ctidx <- ctidx + me
        ctt <- ctt+1
      }
    }
    
    if(i > burnin){
      llik_s = waic2(y, Km, ZZm,n.betaa,pidx,hs)
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
      Zmaxmat[(i-1)/ksamp.+1,] <- as.vector(zz2max);
      for(m in 1: M){
        Zmmat[[m]][(i-1)/ksamp.+1,] <- as.vector(zzm[,,m])
        Zmmaxmat[[m]][(i-1)/ksamp.+1,] <- as.vector(zz.ret$Zmmax[,,m])
      }
    }
    o.betaa <- n.betaa;
    rm(n.betaa);
    n.betaa <- rep(NA,length(o.betaa));
    if(ceiling(i/report1.)==floor(i/report1.)) {
      message("iteration: ", i, " ", date(), " Run: ", run.);
    }
  }
  WAIC = -2 * (sum(log(colMeans(lik))) - sum(apply(llik,2,var)))
  return(list(betaa=betaa,Z=Zmat,Zmax=Zmaxmat,Zmmat = Zmmat, Zmmax = Zmmaxmat,pidx=pidx,lik = lik, llik = llik,WAIC = WAIC));
}

