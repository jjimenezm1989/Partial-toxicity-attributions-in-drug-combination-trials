#Article: Cancer phase I trials design using drug combinations when a fraction of dose limiting toxicities is attributable to one or mode agents
#Authors: Jose L. Jimenez, Mourad Tighiouart and Mauro Gasparini.
#E-mail addresses: jose.jimenez@polito.it, Mourad.Tighiouart@cshs.org, gasparini@calvino.polito.it
#File: JimenezTighiouartGasparini_BiometricalJournalScript4.R

####################################################################################################
# This program runs the design and computes operating characteristics under model misspecification #
####################################################################################################

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

# This function simulates a binary DLT response.
outcomegenerator = function(p,eta){
  tox = rbinom(1,1,p)
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
  z
}

# This function computes the probability of DLT at dose combination (x,y)

pdlt = function(alpha,beta,gamma,eta,tox,assig,d1,d2,x,y){
  
  
  pi10 = x^alpha*(1-y^beta) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi01 = y^beta*(1-x^alpha) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi11 = x^alpha*y^beta + x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  p = pi01+pi10+pi11
  
  p
}

# This function rounds a continuous dose to discrete dose downwards
rndown=function(w,x){
  l=length(x)
  if (w==1) {
    y=1} else {
      
      for (i in 1:(l-1)){
        
        if (x[i] <= w & w <x[i+1])
          
          y= x[i] } }
  y}

# This function rounds a continuous dose to the nearest discrete dose
rndnearest=function(w,x){
  l=length(x)
  mid=numeric()
  for(i in 1:(l-1)){
    mid[i]=(x[i]+x[i+1])/2}
  if (w==1) {
    y=1} else 
      if (x[1] <= w & w < mid[1]) {
        y=x[1]} else
          if (w >= mid[l-1]) {
            y=x[l]} else
              
            {
              
              for (i in 2:(l-1)){
                if (mid[i-1] <= w & w <mid[i])
                  y= x[i] } }
  y}

#These 2 functions computes the MTD
mtd_function<-function(alpha,beta,gamma,x,theta){
  y = mat.or.vec(1,length(x))
  for(i in 1:length(x)){
    xx = x[i]^alpha
    xx2 = x[i]^(2*alpha)
    xxm = x[i]^(-alpha)
    A = (exp(-gamma)-1)/(exp(-gamma)+1)
    root = sqrt((-A*xx2+A*xx+xx-1)^2-4*A*xx*(xx-1)*(theta-xx))
    
    y[i] = 2^(-1/beta)*(xxm/(A*(xx-1))*(A*xx2-A*xx-root-xx+1))^(1/beta)
  }
  if(is.nan(y[i])){
    y[i] = 0
  }
  return(y)
}

MTDdoses = function(alpha, beta, gamma, theta, nx, ny) {
  
  # special senarios
  if(mtd_function(alpha,beta,gamma,x=0.075,theta) < 0.075) {
    doses = data.frame(x=0.075, y=0.075)
  } else if(mtd_function(alpha,beta,gamma,x=0.3,theta) > 0.3) {
    doses = data.frame(x=0.3, y=0.3)
  } else {
    
    # set number of descrete levels for drug X and Y
    x1 = seq(0.075,0.3, length.out=nx)
    y1 = seq(0.075,0.3, length.out=ny)
    
    # find intercept
    xmin = 0.075
    xmax = 0.3
    
    if(mtd_function(alpha,beta,gamma,x=0.075,theta) > 0.3) xmin=uniroot(function(x) mtd_function(alpha,beta,gamma,x,theta) - 0.3, c(0.075,0.3),extendInt="downX")$root
    
    #xmin=uniroot(function(x) (log(theta/(1-theta))-log(rho00/(1-rho00))-(log(rho10/(1-rho10)) - log(rho00/(1-rho00)))*x)/(log(rho01/(1-rho01)) - log(rho00/(1-rho00))+eta*x)-1, c(0,1),extendInt="downX")$root
    
    
    
    if(mtd_function(alpha,beta,gamma,x=0.3,theta) < 0.075) xmax=uniroot(function(x) mtd_function(alpha,beta,gamma,x,theta),c(0.075,0.3),extendInt="downX")$root
    
    # Define the function to be minimized
    fdist<-function(x,pt,alpha,beta,gamma,theta){
      # pt is the point
      d<-(pt[1]-x)^2+(pt[2]- mtd_function(alpha,beta,gamma,x,theta))^2
      d<-d^0.5
      d}
    
    ######
    # calculate the mimimum distance for each dose combinations
    ld = matrix(nrow=length(x1), ncol=length(y1))
    
    for (i in 1:length(x1))
      for (j in 1:length(y1)) 
      {
        {
          
          aal<-optimize(f=fdist,c(xmin,xmax),tol=0.0001,c(x1[i],y1[j]),alpha,beta,gamma,theta)
          ld[i,j]<-aal$objective
          
        }
      }
    
    ######
    # find the best combinations
    
    # 1. select x fist
    ld1 = ld
    for (i in 1:length(x1)) ld1[i,which(ld1[i,] != min(ld1[i,],na.rm=T))] = NA
    for (j in 1:length(y1)) {
      if (sum(!is.na(ld1[,j])) > 0)
        ld1[which(ld1[,j] != min(ld1[,j],na.rm=T)),j] = NA
    }
    
    # 2. select y fist
    ld2 = ld
    for (j in 1:length(y1)) ld2[which(ld2[,j] != min(ld2[,j],na.rm=T)),j] = NA
    for (i in 1:length(x1)) {
      if (sum(!is.na(ld2[i,])) > 0)
        ld2[i,which(ld2[i,] != min(ld2[i,],na.rm=T))] = NA
    }
    
    # intersection ld1 and ld2
    temp = data.frame(ld1 = as.numeric(ld1), ld2 = as.numeric(ld2))
    ld = matrix(apply(temp[, 1:2], 1, mean), nrow=nx, ncol=ny) # intersection
    ld[is.nan(ld)]=  NA
    
    # determine the dose combinations for MTD
    dosex = matrix(x1,nrow=length(x1), ncol=length(y1))
    dosey = matrix(y1,byrow=T,nrow=length(x1), ncol=length(y1))
    doses = data.frame(x = dosex[which(!is.na(ld))], y=dosey[which(!is.na(ld))])
    #
    #doses is the dose combinations that are close to the true MTD curve
    return(doses)
  }
}


#######################
#Main part of the code#
#######################


# Threshold for posterior probability of toxicity at minimum dose for stopping trial
delta=0.8

delta1=0.1

# Target Probability of DLT 

theta=0.3
# Number of patients in a trial
NN=40
# Number of Trials
M=100

# Define discrete dose combinations and true probability of DLT

nx=4
ny=4
data1 = data.frame(x=rep(seq(0.075,0.3,length.out=nx),each=ny),y=rep(seq(0.075,0.3,length.out=ny),nx),p=NA)

#Probability of toxicity scenario (see manuscript, section 5 to spot this scenario)
data1$p = c(0.19,0.22,0.25,0.28,0.26,0.30,0.35,0.41,0.34,0.40,0.48,0.55,0.43,0.51,0.60,0.68)

N=nx*ny
x1=seq(0.075,0.3,length.out=nx)
y1=seq(0.075,0.3,length.out=ny)

# Prior parameters for alpha, beta, gamma, eta

a1 = 0.2
b1 = 2
a2 = 0.2
b2 = 2
a3 = 0.1
b3 = 0.1
a4 = 0
b4 = 1
teta=0.40 #Fraction of attribution

valpha = numeric()
vbeta = numeric()
vgamma = numeric()
veta = numeric()

# mcmc parameters for trial conduct

chains=1
burn=8000
mm=4000

# Declaration of intermediate statistics
trcmtdx1=numeric()
trcmtdx2=numeric()
trcmtdy1=numeric()
trcmtdy2=numeric()

# Declaration of the output statistics
dosex=mat.or.vec(M,NN)
dosey=mat.or.vec(M,NN)
dlt=mat.or.vec(M,NN)
assign = mat.or.vec(M,NN)
assign_d1 = mat.or.vec(M,NN)
assign_d2 = mat.or.vec(M,NN)
temp1=numeric()
temp2=numeric()

# Start of main loop for the trials

for (kk in 1:M) {
  
  # data for first cohort of two patients- Z is DLT status
  
  X=c(0.075,0.075)
  Y=c(0.075,0.075)
  
  # Simulate responses for patients 1 and 2 from logistic model
  
  resp1 = outcomegenerator(data1$p[data1$x==X[1] & data1$y==Y[1]],teta)
  resp2 = outcomegenerator(data1$p[data1$x==X[2] & data1$y==Y[2]],teta)
  
  Z=c(resp1[[1]],resp2[[1]])
  A = c(resp1[[2]],resp2[[2]])
  D1 = c(resp1[[3]],resp2[[3]])
  D2 = c(resp1[[4]],resp2[[4]])
  
  n=2
  
  zeros = c(0,0)
  
  # Get dose Y for patient 3 when dose X equals 0 and dose X for patient 4 when dose Y equals 0
  
  j=jags.model('../toxicity_attribution.bug.txt',data=list('zeros'=zeros,'Z'=Z,'X'=X,'Y'=Y,'A'=A,'D1'=D1,'D2'=D2,'a1'=a1,'b1'=b1,'a2'=a2,'b2'=b2,'a3'=a3,'b3'=b3,'a4'=a4,'b4'=b4,'N'=n),n.chains=chains,n.adapt=burn)
  
  s=coda.samples(j,c('alpha','beta','gamma','eta'),mm)
  ss=as.data.frame(s[[1]])
  
  medalpha=median(ss$alpha)
  medbeta=median(ss$beta)
  medgamma=median(ss$gamma)
  medeta=median(ss$eta)
  
  xx1 = optimize(f=function(x){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[1],A[1],D1[1],D2[1],x,Y[1]) - theta)}, c(0.075,0.3),tol=0.0001)$minimum
  
  if(Z[1]==1 & A[1]==1 & D1[1]==1 & (xx1>X[1])){
    
    xx1 = X[1]
    
  }else{
    
    if ((xx1 - X[1]) > delta1)
      xx1=X[1]+delta1
    
  }
  
  xx1=rndnearest(xx1,x1)
  
  yy2 = optimize(f=function(y){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2],A[2],D1[2],D2[2],X[2],y) - theta)}, c(0.075,0.3),tol=0.0001)$minimum
  
  if(Z[2]==1 & A[2]==1 & D2[2]==1 & (yy2>Y[2])){
    
    yy2 = Y[2]
    
  }else{
    
    if ((yy2 - Y[2]) > delta1)
      yy2=Y[2]+delta1
    
  }
  
  yy2=rndnearest(yy2,y1)
  
  # Update data
  
  X=c(X,0.075,xx1)
  Y=c(Y,yy2,0.075)
  
  # Simulate responses for patients 3 and 4 from logistic model
  
  resp1 = outcomegenerator(data1$p[data1$x==X[3] & data1$y==Y[3]],teta)
  resp2 = outcomegenerator(data1$p[data1$x==X[4] & data1$y==Y[4]],teta)
  Z=c(Z,resp1[[1]],resp2[[1]])
  A = c(A,resp1[[2]],resp2[[2]])
  D1 = c(D1,resp1[[3]],resp2[[3]])
  D2 = c(D2,resp1[[4]],resp2[[4]])
  
  n=4
  
  zeros = c(zeros,0,0)
  
  
  for (i in 3:(NN/2)) {
    
    j=jags.model('../toxicity_attribution.bug.txt',data=list('zeros'=zeros,'Z'=Z,'X'=X,'Y'=Y,'A'=A,'D1'=D1,'D2'=D2,'a1'=a1,'b1'=b1,'a2'=a2,'b2'=b2,'a3'=a3,'b3'=b3,'a4'=a4,'b4'=b4,'N'=n),n.chains=chains,n.adapt=burn)
    
    s=coda.samples(j,c('alpha','beta','gamma','eta'),mm)
    ss=as.data.frame(s[[1]])
    
    medalpha=median(ss$alpha)
    medbeta=median(ss$beta)
    medgamma=median(ss$gamma)
    medeta=median(ss$eta)
    
    xx1 = optimize(f=function(x){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-3],A[2*i-3],D1[2*i-3],D2[2*i-3],x,Y[2*i-3]) - theta)}, c(0.075,0.3),tol=0.0001)$minimum
    
    #If there is a DLT, we allow to de-esclate but we won't escalate.
    
    if(Z[2*i-3]==1 & A[2*i-3]==1 & D1[2*i-3]==1 & (xx1>X[2*i-3])){
      
      xx1 = X[2*i-3]
      
      #If there is no DLT, we limit the escalation up to delta1
      
    }else{
      
      if ((xx1 - X[2*i-3]) > delta1)
        xx1=X[2*i-3]+delta1
      
    }
    
    xx1=rndnearest(xx1,x1)
    
    xx2 = optimize(f=function(x){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-2],A[2*i-2],D1[2*i-2],D2[2*i-2],x,Y[2*i-2]) - theta)}, c(0.075,0.3),tol=0.0001)$minimum
    
    if(Z[2*i-2]==1 & A[2*i-2]==1 & D1[2*i-2]==1 & (xx2>X[2*i-2])){
      
      xx2 = X[2*i-2]
      
    }else{
      
      if ((xx2 - X[2*i-2]) > delta1)
        xx2=X[2*i-2]+delta1
    }
    
    xx2=rndnearest(xx2,x1)
    
    yy1 = optimize(f=function(y){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-3],A[2*i-3],D1[2*i-3],D2[2*i-3],X[2*i-3],y) - theta)}, c(0.075,0.3),tol=0.0001)$minimum
    
    if(Z[2*i-3]==1 & A[2*i-3]==1 & D2[2*i-3]==1 & (yy1>Y[2*i-3])){
      
      yy1 = Y[2*i-3]
      
    }else{
      
      if ((yy1 - Y[2*i-3]) > delta1)
        yy1=Y[2*i-3]+delta1
      
    }
    
    yy1=rndnearest(yy1,y1)
    
    yy2 = optimize(f=function(y){abs(pdlt(medalpha,medbeta,medgamma,medeta,Z[2*i-2],A[2*i-2],D1[2*i-2],D2[2*i-2],X[2*i-2],y) - theta)}, c(0.075,0.3),tol=0.0001)$minimum
    
    if(Z[2*i-2]==1 & A[2*i-2]==1 & D2[2*i-2]==1 & (yy2>Y[2*i-2])){
      
      yy2 = Y[2*i-2]
      
    }else{
      
      if ((yy2 - Y[2*i-2]) > delta1)
        yy2=Y[2*i-2]+delta1
      
    }
    
    yy2=rndnearest(yy2,y1)
    
    if (X[2*i-3] == X[2*i-5]) {
      
      X=c(X,xx1,X[2*i-2])
      Y=c(Y,Y[2*i-3],yy2)} else {
        X=c(X,X[2*i-3],xx2)
        Y=c(Y,yy1,Y[2*i-2])}
    
    
    n=n+2
    
    zeros = c(zeros,0,0)
    
    data1$x==X[2*i-1] & data1$y==Y[2*i-1]
    
    resp1 = outcomegenerator(data1$p[data1$x==X[2*i-1] & data1$y==Y[2*i-1]],teta)
    resp2 = outcomegenerator(data1$p[data1$x==X[2*i] & data1$y==Y[2*i]],teta)
    Z=c(Z,resp1[[1]],resp2[[1]])
    A = c(A,resp1[[2]],resp2[[2]])
    D1 = c(D1,resp1[[3]],resp2[[3]])
    D2 = c(D2,resp1[[4]],resp2[[4]])
    
  }
  
  dosex[kk,]=X
  dosey[kk,]=Y
  dlt[kk,]=Z
  assign[kk,]=A
  assign_d1[kk,]=D1
  assign_d2[kk,]=D2
  
  valpha[kk]=median(ss$alpha)
  vbeta[kk]=median(ss$beta)
  vgamma[kk]=median(ss$gamma)
  veta[kk]=median(ss$eta)
  
  cat("Iteration",kk,"\n")
  
  rm(X,Y,Z,A,D1,D2,ss,zeros,n,xx1,yy1,xx2,yy2,resp1,resp2,s)
  
}


#Now we compute the results

result = data.frame(x=NA,y=NA,trial=NA)

for (i in 1:M) {
  
  a = MTDdoses(valpha[i], vbeta[i], vgamma[i], theta, nx, ny)
  
  
  if(nrow(a)>0) a$trial = i
  
  result = rbind(result, a)
  
  print(i)
  
}

############################
#We now compute the results#
############################

# remove NAs (including NA in the first row)
result = result[!is.na(result$x),]


true1 = mat.or.vec(nx*ny,4)
true1[,4]=1:(length(x1)*length(y1))

count=1

for(i in 1:length(x1)){
  for(j in 1:length(y1)){
    true1[count,1] = x1[i]
    true1[count,2] = y1[j]
    count = count + 1
  }
}

true1[,3] = data1$p

true = true1[abs(true1[,3]-theta)<=0.1,]
colnames(true) = c("x","y","prob","id")

#Code to obtain the summary statistcs (percent of times that at leas 25,50,75 and 100% of the recommended doses belong to the true MTD set)

all = mat.or.vec(1,1000)
atleast25 = mat.or.vec(1,1000)
atleast50 = mat.or.vec(1,1000)
atleast75 = mat.or.vec(1,1000)

atleast1 = mat.or.vec(1,1000) #Variable only used to produce Table S3 for the supplementary material

for (i in 1:1000){
  
  num.doses.recommended = nrow(result[result$trial==i,])
  
  num.true.doses.recommended = nrow(merge(result[result$trial==i,],as.data.frame(true[,-3]),by=c("x","y")))
  
  if(num.true.doses.recommended >= 1){
    atleast1[i] = 1
  }
  
  if((num.true.doses.recommended/num.doses.recommended)==1){
    all[i]=1
  }
  
  if((num.true.doses.recommended/num.doses.recommended)>=0.25){
    atleast25[i]=1
  }
  
  if((num.true.doses.recommended/num.doses.recommended)>=0.50){
    atleast50[i]=1
  }
  
  if((num.true.doses.recommended/num.doses.recommended)>=0.75){
    atleast75[i] = 1
  }
  
}

#Summary Safety

aux = mat.or.vec(1,1000)

for(i in 1:1000){
  
  aux[i] = mean(dlt[i,])
  
}

highToxRate005 = length(which(aux>(theta + 0.05)))/1000 
highToxRate010 = length(which(aux>(theta + 0.1)))/1000 

mean(dlt) #Average % of toxicities
highToxRate005*100 #% of trials with toxicity rate > theta + 0.05
highToxRate010*100 #% of trials with toxicity rate > theta + 0.10

####################################
#Summary % of correct MTD selection#
####################################

mean(atleast25)
mean(atleast50)
mean(atleast75)
mean(all)
mean(atleast1)
