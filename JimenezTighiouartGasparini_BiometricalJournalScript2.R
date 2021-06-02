#Article: Cancer phase I trials design using drug combinations when a fraction of dose limiting toxicities is attributable to one or mode agents
#Authors: Jose L. Jimenez, Mourad Tighiouart and Mauro Gasparini.
#E-mail addresses: jose.jimenez@polito.it, Mourad.Tighiouart@cshs.org, gasparini@calvino.polito.it
#File: JimenezTighiouartGasparini_BiometricalJournalScript2.R

########################################################################################
# This program computes and saves posterior probabilities of DLT at MTD in a .csv file.#
########################################################################################

#Code to obtain a random seed (if you want to set a particular seed, comment the lines 11, 12 and 13 and put the seed in line 16)
rm(list=ls())
t <- as.numeric(Sys.time())
seed <- 1e8 * (t - floor(t))

#Now we set the seed (either the random one or the one that we want)
set.seed(seed)


############
#Libraries #
############

library(rjags)

#############
# Functions #
#############

# This function computes the probability of DLT at dose combination (x,y) using Gumbel model

pdlt = function(alpha,beta,gamma,eta,tox,assig,d1,d2,x,y){
  
  #Compute pi10, pi01, pi11, pi00
  
  pi10 = x^alpha*(1-y^beta) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi01 = y^beta*(1-x^alpha) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi11 = x^alpha*y^beta + x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  p = pi10+pi01+pi11
  
  p
}

# Function to calculate MTD set as in equation (10) in the manuscript

MTDy = function(alpha,beta,gamma,x,theta){
  
  xx = x^alpha
  xx2 = x^(2*alpha)
  xxm = x^(-alpha)
  A = (exp(-gamma)-1)/(exp(-gamma)+1)
  root = sqrt((-A*xx2+A*xx+xx-1)^2-4*A*xx*(xx-1)*(theta-xx))
  
  y = 2^(-1/beta)*(xxm/(A*(xx-1))*(A*xx2-A*xx-root-xx+1))^(1/beta)
  
  return(y)
}

# Function to calculate MTD set as in equation (10) in the manuscript but isolating X instead of Y

MTDx = function(alpha,beta,gamma,y,theta){
  
  yy = y^beta
  yy2 = y^(2*beta)
  yym = y^(-beta)
  A = (exp(-gamma)-1)/(exp(-gamma)+1)
  root = sqrt((-A*yy2+A*yy+yy-1)^2-4*A*yy*(yy-1)*(theta-yy))
  
  x = 2^(-1/alpha)*(yym/(A*(yy-1))*(A*yy2-A*yy-root-yy+1))^(1/alpha)
  
  return(x)
}

# Function to calculate posterior probability of DLT at MTD (x,y)

postdlt = function(NN, theta, nx, ny, triali) {
  
  # Prior parameters for alpha, beta, gamma and eta (following Ying and Yuan (2009a))
  a1 = 0.2
  b1 = 2
  a2 = 0.2
  b2 = 2
  a3 = 0.1
  b3 = 0.1
  a4 = 0
  b4 = 1
  
  # Data for a full trial consists of N patients, X doses, Y doses, 
  
  X<-as.numeric(read.table("../dosex.txt")[triali,])
  Y<-as.numeric(read.table("../dosey.txt")[triali,])
  Z<-as.numeric(read.table("../dlt.txt")[triali,])
  A<-as.numeric(read.table("../assign.txt")[triali,])
  D1<-as.numeric(read.table("../assign_d1.txt")[triali,])
  D2<-as.numeric(read.table("../assign_d2.txt")[triali,])
  
  #Compute the posterior distribution of model parameters give data 
  # mcmc parameters for trial conduct
  
  chains<-1
  burn<-10000
  mm<-5000
  
  zeros = rep(0,NN)
  
  j=jags.model('../toxicity_attribution.bug.txt',data=list('zeros'=zeros,'Z'=Z, 'A'=A, 'X'=X,'Y'=Y,'D1'=D1,'D2'=D2,'a1'=a1,'b1'=b1,'a2'=a2,'b2'=b2,'a3'=a3,'b3'=b3,'a4'=a4,'b4'=b4,'N'=NN),n.chains=chains,n.adapt=burn)
  
  s=coda.samples(j,c('alpha','beta','gamma','eta'),mm)
  ss=as.data.frame(s[[1]])
  
  mtdy1 <- MTDy(as.matrix(s[,1]),as.matrix(s[,2]),as.matrix(s[,4]),X[(NN-1)],theta)
  mtdy2 <- MTDy(as.matrix(s[,1]),as.matrix(s[,2]),as.matrix(s[,4]),X[NN],theta)
  
  mtdx1 <- MTDx(as.matrix(s[,1]),as.matrix(s[,2]),as.matrix(s[,4]),Y[(NN-1)],theta)
  mtdx2 <- MTDx(as.matrix(s[,1]),as.matrix(s[,2]),as.matrix(s[,4]),Y[NN],theta)
  
  
  # Calculating Posterior Probability of DLT at MTD (xx,yy)
  temp1 = rep(NA, mm)
  postdlts = data.frame(xx=1:(nx*ny), yy=NA, postdlt=NA, trial=triali)
  
  k=0
  for (xx in seq(0.05,0.3,length.out=nx)) {
    for (yy in seq(0.05,0.3,length.out=ny)) {
      k=k+1
      for (ll in 1:mm) temp1[ll]<-pdlt(ss$alpha[ll],ss$beta[ll],ss$gamma[ll],ss$eta[ll],tox=0,assig=0,d1=0,d2=0,xx,yy)
      postdlts[k,1:3] = c(xx, yy, mean(temp1 > theta+0.1  | temp1 < theta-0.1))
    }
  }
  
  
  return(postdlts)
}

#######################
#Main part of the code#
#######################

# Trial sample size NN
NN<-40

#Target probability of DLT theta
theta<-0.3

M = 1000
nx=4
ny=4

# run for each trial
result = NA
for (jj in 1:M) {
  a= postdlt(NN, theta, nx, ny, jj)
  result = rbind(result, a)
}
result = result[-1,]

###########################################################
#Save posterior probabilities of DLT at MTD in a .csv file#
###########################################################

write.csv(result, file="../postdlts.csv", row.names=F)
