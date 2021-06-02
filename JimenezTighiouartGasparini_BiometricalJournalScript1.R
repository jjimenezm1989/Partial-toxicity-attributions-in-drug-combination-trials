#Article: Cancer phase I trials design using drug combinations when a fraction of dose limiting toxicities is attributable to one or mode agents
#Authors: Jose L. Jimenez, Mourad Tighiouart and Mauro Gasparini.
#E-mail addresses: jose.jimenez@polito.it, Mourad.Tighiouart@cshs.org, gasparini@calvino.polito.it
#File: JimenezTighiouartGasparini_BiometricalJournalScript1.R

##########################################################################################
# This code run the partial toxicity attribution trial with continuous dose combinations # 
##########################################################################################

#Code to obtain a random seed (if you want to set a particular seed, comment the lines 11, 12 and 13 and put the seed in line 16)
rm(list=ls())
t <- as.numeric(Sys.time())
seed <- 1e8 * (t - floor(t))

#Now we set the seed (either the random one or the one that we want)
set.seed(seed)

#############
# Libraries #
#############

library(rjags)

#############
# Functions #
#############

# This function generates all outcomes for 1 patients (DLT, Attribution, Indication of which drug caused the DLT)

outcomegenerator<-function(alpha,beta,gamma,eta,x,y){
  
  pi10 = x^alpha*(1-y^beta) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi01 = y^beta*(1-x^alpha) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi11 = x^alpha*y^beta + x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi00 = (1-x^alpha)*(1-y^beta) + x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  risk = pi10 + pi01 + pi11
  
  tox = rbinom(1,1,risk)
  
  assig = rbinom(1,1,eta)
  
  aux2 = runif(1)
  
  if(aux2<=(1/3)){
    d1=1
    d2=0
  }else if(aux2>(1/3) & aux2<=(2/3)){
    d1=0
    d2=1
  }else{
    d1=1
    d2=1
  }
  
  z = list(tox,assig,d1,d2)
  
}

# This function computes the probability of DLT at dose combination (x,y)

pdlt = function(alpha,beta,gamma,eta,tox,assig,d1,d2,x,y){
  
  #Compute pi10, pi01, pi11, pi00
  
  pi10 = x^alpha*(1-y^beta) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi01 = y^beta*(1-x^alpha) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi11 = x^alpha*y^beta + x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  p = pi01+pi10+pi11
  
  p
}

#######################
#Main part of the code#
#######################

#Maximum "jump" allowed when doing dose escalation
delta1<-0.1

# Number of Trials
M<-1000

# True parameters values.
#For alpha, beta and gamma we put the expected value of the prior distributions just as an example. The value of eta is the one
#We need to change to to increase the fraction of toxicity attribution.
talpha = 0.5*(0.2+2)
tbeta = 0.5*(0.2+2)
tgamma = 0.1/0.1
teta = 0

# Target Probability of DLT 
theta<-0.3

# Number of patients in a trial
NN<-40

# Prior parameters for alpha, beta, gamma, eta following Ying and Yuan (2009a), except for eta.
a1 = 0.2
b1 = 2
a2 = 0.2
b2 = 2
a3 = 0.1
b3 = 0.1
a4 = 0
b4 = 1

#Variables where we are going to store the parameter value at the end of each simulated trial
valpha = numeric()
vbeta = numeric()
vgamma = numeric()
veta = numeric()

# mcmc parameters

chains<-1 #1 chain
burn<-10000 #10000 iteration
mm<-5000 #Burn-in of 5000 

# Declaration of the output statistics
dosex<-mat.or.vec(M,NN)
dosey<-mat.or.vec(M,NN)
dlt<-mat.or.vec(M,NN)
assign = mat.or.vec(M,NN)
assign_d1 = mat.or.vec(M,NN)
assign_d2 = mat.or.vec(M,NN)

for (kk in 1:M) {
  
  
  # Dose combination of the first cohort of two patients.
  X<-c(0.05,0.05)
  Y<-c(0.05,0.05)
  
  # Simulate DLTs, assignment, delta1 and delta2 for patients 1 and 2 from the model
  resp1<-outcomegenerator(talpha,tbeta,tgamma,teta,X[1],Y[1])
  resp2<-outcomegenerator(talpha,tbeta,tgamma,teta,X[2],Y[2])
  Z<-c(resp1[[1]],resp2[[1]])
  A = c(resp1[[2]],resp2[[2]])
  D1 = c(resp1[[3]],resp2[[3]])
  D2 = c(resp1[[4]],resp2[[4]])
  
  n<-2
  
  # Perform MCMC to obtain posterior of the model parameters
  
  zeros = c(0,0) #This is neccesary for the zeros-trick we use. See jags file.
  
  j=jags.model('../toxicity_attribution.bug.txt',data=list('zeros'=zeros,'Z'=Z,'X'=X,'Y'=Y,'A'=A,'D1'=D1,'D2'=D2,'a1'=a1,'b1'=b1,'a2'=a2,'b2'=b2,'a3'=a3,'b3'=b3,'a4'=a4,'b4'=b4,'N'=n),n.chains=chains,n.adapt=burn)
   
  
  s=coda.samples(j,c('alpha','beta','gamma','eta'),mm)
  ss=as.data.frame(s[[1]])
  
  medalpha<-median(ss$alpha)
  medbeta<-median(ss$beta)
  medgamma<-median(ss$gamma)
  medeta<-median(ss$eta)
  
  
  
  #Find dose combination of X give that Y is fixed for second cohort of patients.
  xx1 = optimize(f=function(x){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[1],A[1],D1[1],D2[1],x,Y[1]) - theta)}, c(0.05,0.3),tol=0.0001)$minimum
  
  #Dose escalation is not allowed in the case of a previous DLT in X.
  if(Z[1]==1 & A[1]==1 & D1[1]==1 & (xx1>X[1])){
    
    xx1 = X[1]
    
  }else{
    
    #If escalation is allowed, we limit the it up to delta1
    if ((xx1 - X[1]) > delta1)
      xx1<-X[1]+delta1
    
  }
  
  #Find dose combination of Y give that X is fixed for second cohort of patients.
  yy2 = optimize(f=function(y){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2],A[2],D1[2],D2[2],X[2],y) - theta)}, c(0.05,0.3),tol=0.0001)$minimum
  
  #Dose escalation is not allowed in the case of a previous DLT in Y.
  if(Z[2]==1 & A[2]==1 & D2[2]==1 & (yy2>Y[2])){
    
    yy2 = Y[2]
    
  }else{
    
    #If escalation is allowed, we limit the it up to delta1
    if ((yy2 - Y[2]) > delta1)
      yy2<-Y[2]+delta1
    
  }
  
  # Add dose combinations of second cohort of patients.
  X<-c(X,0.05,xx1)
  Y<-c(Y,yy2,0.05)
  
  # Simulate data for the second cohort of patients.
  resp1<-outcomegenerator(talpha,tbeta,tgamma,teta,X[3],Y[3])
  resp2<-outcomegenerator(talpha,tbeta,tgamma,teta,X[4],Y[4])
  Z<-c(Z,resp1[[1]],resp2[[1]])
  A = c(A,resp1[[2]],resp2[[2]])
  D1 = c(D1,resp1[[3]],resp2[[3]])
  D2 = c(D2,resp1[[4]],resp2[[4]])
  
  n<-4 #Update the number of patients included in the trial so far.
  
  zeros <- c(zeros,0,0) #This is neccesary for the zeros-trick we use in the MCMC (see jags file).
  
  #So far, we have manually dose patients 1,2,3 and 4. Now we enter a loop to do the rest of the patients following the same procedure 
  #we have been following.
  
  for (i in 3:(NN/2)) {
    
    # Perform MCMC to obtain posterior of the model parameters
    
    j=jags.model('../toxicity_attribution.bug.txt',data=list('zeros'=zeros,'Z'=Z,'X'=X,'Y'=Y,'A'=A,'D1'=D1,'D2'=D2,'a1'=a1,'b1'=b1,'a2'=a2,'b2'=b2,'a3'=a3,'b3'=b3,'a4'=a4,'b4'=b4,'N'=n),n.chains=chains,n.adapt=burn)
    
    s=coda.samples(j,c('alpha','beta','gamma','eta'),mm)
    ss=as.data.frame(s[[1]])
    
    medalpha<-median(ss$alpha)
    medbeta<-median(ss$beta)
    medgamma<-median(ss$gamma)
    medeta<-median(ss$eta)
    
    
    #Find dose combination of X give that Y is fixed
    xx1 = optimize(f=function(x){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-3],A[2*i-3],D1[2*i-3],D2[2*i-3],x,Y[2*i-3]) - theta)}, c(0.05,0.3),tol=0.0001)$minimum
    
    #Dose escalation is not allowed in the case of a previous DLT in X.
    if(Z[2*i-3]==1 & A[2*i-3]==1 & D1[2*i-3]==1 & (xx1>X[2*i-3])){
      
      xx1 = X[2*i-3]
      
    }else{
      
      #If escalation is allowed, we limit the it up to delta1
      if ((xx1 - X[2*i-3]) > delta1)
        xx1<-X[2*i-3]+delta1
      
    }
    
    #Find dose combination of X give that Y is fixed
    xx2 = optimize(f=function(x){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-2],A[2*i-2],D1[2*i-2],D2[2*i-2],x,Y[2*i-2]) - theta)}, c(0.05,0.3),tol=0.0001)$minimum
    
    #Dose escalation is not allowed in the case of a previous DLT in X.
    if(Z[2*i-2]==1 & A[2*i-2]==1 & D1[2*i-2]==1 & (xx2>X[2*i-2])){
      
      xx2 = X[2*i-2]
    
    }else{
      
      #If escalation is allowed, we limit the it up to delta1
      if ((xx2 - X[2*i-2]) > delta1)
        xx2<-X[2*i-2]+delta1
    }
    
    #Find dose combination of Y give that X is fixed
    yy1 = optimize(f=function(y){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-3],A[2*i-3],D1[2*i-3],D2[2*i-3],X[2*i-3],y) - theta)}, c(0.05,0.3),tol=0.0001)$minimum
    
    #Dose escalation is not allowed in the case of a previous DLT in Y.
    if(Z[2*i-3]==1 & A[2*i-3]==1 & D2[2*i-3]==1 & (yy1>Y[2*i-3])){
      
      yy1 = Y[2*i-3]
    
    }else{
      
      #If escalation is allowed, we limit the it up to delta1
      if ((yy1 - Y[2*i-3]) > delta1)
        yy1<-Y[2*i-3]+delta1
      
    }
    
    #Find dose combination of Y give that X is fixed
    yy2 = optimize(f=function(y){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-2],A[2*i-2],D1[2*i-2],D2[2*i-2],X[2*i-2],y) - theta)}, c(0.05,0.3),tol=0.0001)$minimum
    
    #Dose escalation is not allowed in the case of a previous DLT in Y.
    if(Z[2*i-2]==1 & A[2*i-2]==1 & D2[2*i-2]==1 & (yy2>Y[2*i-2])){
      
      yy2 = Y[2*i-2]
      
    }else{
      
      #If escalation is allowed, we limit the it up to delta1
      if ((yy2 - Y[2*i-2]) > delta1)
        yy2<-Y[2*i-2]+delta1
      
    }
    
    
    #Add dose combination of the following cohort of patients    
    if (X[2*i-3] == X[2*i-5]) {
      X<-c(X,xx1,X[2*i-2])
      Y<-c(Y,Y[2*i-3],yy2)} else {
        X<-c(X,X[2*i-3],xx2)
        Y<-c(Y,yy1,Y[2*i-2])}
    
    
    n<-n+2 #Update the number of patients included in the trial so far
    
    zeros <- c(zeros,0,0)#This is neccesary for the zeros-trick we use. See jags file.
    
    #Generate data for the next cohort of patients
    resp1<-outcomegenerator(talpha,tbeta,tgamma,teta,X[2*i-1],Y[2*i-1])
    resp2<-outcomegenerator(talpha,tbeta,tgamma,teta,X[2*i],Y[2*i])
    Z<-c(Z,resp1[[1]],resp2[[1]])
    A = c(A,resp1[[2]],resp2[[2]])
    D1 = c(D1,resp1[[3]],resp2[[3]])
    D2 = c(D2,resp1[[4]],resp2[[4]])
    
  }
  
  #Store the data and parameters after each simulated trial.
  dosex[kk,]<-X
  dosey[kk,]<-Y
  dlt[kk,]<-Z
  assign[kk,]=A
  assign_d1[kk,]=D1
  assign_d2[kk,]=D2
  
  valpha[kk]<-median(ss$alpha)
  vbeta[kk]<-median(ss$beta)
  vgamma[kk]<-median(ss$gamma)
  veta[kk]<-median(ss$eta)
  
}


#Save the followig variables. This is neccesary if we want to go to a discrete dose combinations setting

write.table(dosex, "../dosex.txt", sep="\t")
write.table(dosey, "../dosey0.txt", sep="\t")
write.table(assign, "../assign.txt", sep="\t")
write.table(assign_d1, "../assign_d1.txt", sep="\t")
write.table(assign_d2, "../assign_d2.txt", sep="\t")
write.table(dlt, "../dlt.txt", sep="\t")


############################
#We now compute the results#
############################

# we produce results TABLE 2 in the manuscript.

aux = mat.or.vec(1,1000)

for(i in 1:1000){
  
  aux[i] = mean(dlt[i,])
  
}

highToxRate005 = length(which(aux>(theta + 0.05)))/1000 
highToxRate010 = length(which(aux>(theta + 0.1)))/1000 

###################
#Summary of safery#
###################

mean(dlt) #Average % of toxicities
highToxRate005*100 #% of trials with toxicity rate > theta + 0.05
highToxRate010*100 #% of trials with toxicity rate > theta + 0.10