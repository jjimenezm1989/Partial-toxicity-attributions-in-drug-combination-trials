#Article: Cancer phase I trials design using drug combinations when a fraction of dose limiting toxicities is attributable to one or mode agents
#Authors: Jose L. Jimenez, Mourad Tighiouart and Mauro Gasparini.
#E-mail addresses: jose.jimenez@polito.it, Mourad.Tighiouart@cshs.org, gasparini@calvino.polito.it
#File: JimenezTighiouartGasparini_BiometricalJournalScript7.R

##################################################################
# Program to obtain Figure 5 (Pointwise % of MTD recommendation) #
##################################################################

#For reference about the scenarios, see Figure 5.

#Load the alpha, beta and gamma estimates for the 3 scenarios. 
#This values were saved in JimenezTighiouartGasparini_BiometricalJournalScript1.R. 
#The names of the parameter files (i.e. the .txt files) are the ones we used and remain here ONLY as an example.
#You should put the .txt files that you obtained.

#Scenario 1

usvalphaeta0 = as.matrix(read.table(paste("../alpha_eta0_menos1sigma.txt", sep="\t")))
usvbetaeta0 = as.matrix(read.table(paste("../beta_eta0_menos1sigma.txt", sep="\t")))
usvgammaeta0 = as.matrix(read.table(paste("../gamma_eta0_menos1sigma.txt", sep="\t")))

usvalphaeta01 = as.matrix(read.table(paste("../alpha_eta010_menos1sigma.txt", sep="\t")))
usvbetaeta01 = as.matrix(read.table(paste("../beta_eta010_menos1sigma.txt", sep="\t")))
usvgammaeta01 = as.matrix(read.table(paste("../gamma_eta010_menos1sigma.txt", sep="\t")))

usvalphaeta025 = as.matrix(read.table(paste("../alpha_eta025_menos1sigma.txt", sep="\t")))
usvbetaeta025 = as.matrix(read.table(paste("../beta_eta025_menos1sigma.txt", sep="\t")))
usvgammaeta025 = as.matrix(read.table(paste("../gamma_eta025_menos1sigma.txt", sep="\t")))

usvalphaeta04 = as.matrix(read.table(paste("../alpha_eta040_menos1sigma.txt", sep="\t")))
usvbetaeta04 = as.matrix(read.table(paste("../beta_eta040_menos1sigma.txt", sep="\t")))
usvgammaeta04 = as.matrix(read.table(paste("../gamma_eta040_menos1sigma.txt", sep="\t")))

#Scenario 2

mvalphaeta0<- as.matrix(read.table(paste("../alpha_eta0_mean.txt", sep="\t")))
mvbetaeta0<- as.matrix(read.table(paste("../beta_eta0_mean.txt", sep="\t")))
mvgammaeta0<- as.matrix(read.table(paste("../gamma_eta0_mean.txt", sep="\t")))

mvalphaeta01<- as.matrix(read.table(paste("../alpha_eta010_mean.txt", sep="\t")))
mvbetaeta01<- as.matrix(read.table(paste("../beta_eta010_mean.txt", sep="\t")))
mvgammaeta01<- as.matrix(read.table(paste("../gamma_eta010_mean.txt", sep="\t")))

mvalphaeta025<- as.matrix(read.table(paste("../alpha_eta025_mean.txt", sep="\t")))
mvbetaeta025<- as.matrix(read.table(paste("../beta_eta025_mean.txt", sep="\t")))
mvgammaeta025<- as.matrix(read.table(paste("../gamma_eta025_mean.txt", sep="\t")))

mvalphaeta04<- as.matrix(read.table(paste("../alpha_eta040_mean.txt", sep="\t")))
mvbetaeta04<- as.matrix(read.table(paste("../beta_eta040_mean.txt", sep="\t")))
mvgammaeta04<- as.matrix(read.table(paste("../gamma_eta040_mean.txt", sep="\t")))

#Scenario 3

dsvalphaeta0 <- as.matrix(read.table(paste("../alpha_eta0_mas1sigma.txt", sep="\t")))
dsvbetaeta0 <- as.matrix(read.table(paste("../beta_eta0_mas1sigma.txt", sep="\t")))
dsvgammaeta0 <- as.matrix(read.table(paste("../gamma_eta0_mas1sigma.txt", sep="\t")))

dsvalphaeta01 <- as.matrix(read.table(paste("../alpha_eta010_mas1sigma.txt", sep="\t")))
dsvbetaeta01 <- as.matrix(read.table(paste("../beta_eta010_mas1sigma.txt", sep="\t")))
dsvgammaeta01 <- as.matrix(read.table(paste("../gamma_eta010_mas1sigma.txt", sep="\t")))

dsvalphaeta025 <- as.matrix(read.table(paste("../alpha_eta025_mas1sigma.txt", sep="\t")))
dsvbetaeta025 <- as.matrix(read.table(paste("../beta_eta025_mas1sigma.txt", sep="\t")))
dsvgammaeta025 <- as.matrix(read.table(paste("../gamma_eta025_mas1sigma.txt", sep="\t")))

dsvalphaeta04 <- as.matrix(read.table(paste("../alpha_eta040_mas1sigma.txt", sep="\t")))
dsvbetaeta04 <- as.matrix(read.table(paste("../beta_eta040_mas1sigma.txt", sep="\t")))
dsvgammaeta04 <- as.matrix(read.table(paste("../gamma_eta040_mas1sigma.txt", sep="\t")))


# This function finds the MTD curve 
mtd_logistic<-function(alpha,beta,gamma,theta,x){
  xx = x^alpha
  xx2 = x^(2*alpha)
  xxm = x^(-alpha)
  A = (exp(-gamma)-1)/(exp(-gamma)+1)
  root = sqrt((-A*xx2+A*xx+xx-1)^2-4*A*xx*(xx-1)*(theta-xx))
  
  y = 2^(-1/beta)*(xxm/(A*(xx-1))*(A*xx2-A*xx-root-xx+1))^(1/beta)
  
  if(is.nan(y)){
    y=0
  }
  
  return(y)
}

# Function to be minimized
fdist<-function(x,pt,alpha,beta,gamma,theta){
  
  y = mtd_logistic(alpha,beta,gamma,theta,x)
  
  d<-(pt[1]-x)^2+(pt[2]-y)^2
  d<-d^0.5
  d
}

M<-1000

tmalpha = 0.5*(0.2+2)
tmbeta = 0.5*(0.2+2)
tmgamma = 0.1/0.1

tusalpha = 0.5*(0.2+2) - sd(0.3^runif(100000,0.2,2))
tusbeta = 0.5*(0.2+2) - sd(0.3^runif(100000,0.2,2))
tusgamma = 0.1/0.1

tdsalpha = 0.5*(0.2+2) + sd(0.3^runif(100000,0.2,2))
tdsbeta = 0.5*(0.2+2) + sd(0.3^runif(100000,0.2,2))
tdsgamma = 0.1/0.1

theta<-0.3

p<-0.2

tmmtdx0 <- mtd_logistic(tmalpha,tmbeta,tmgamma,theta,0.05)
tusmtdx0 <- mtd_logistic(tusalpha,tusbeta,tusgamma,theta,0.05)
tdsmtdx0 <- mtd_logistic(tdsalpha,tdsbeta,tdsgamma,theta,0.05)

x1<-seq(0.05,min(tmmtdx0,0.3),by=0.01)
x2<-seq(0.05,min(tusmtdx0,0.3),by=0.01)
x3<-seq(0.05,min(tdsmtdx0,0.3),by=0.01)

kk1<-length(x1)
kk2<-length(x2)
kk3<-length(x3)

ymmtd<-numeric()
yusmtd<-numeric()
ydsmtd<-numeric()

for(i in 1:kk1){
  ymmtd[i]<-mtd_logistic(tmalpha,tmbeta,tmgamma,theta,x1[i])
}

for(i in 1:kk2){
  yusmtd[i]<-mtd_logistic(tusalpha,tusbeta,tusgamma,theta,x2[i])
}

for(i in 1:kk3){
  ydsmtd[i]<-mtd_logistic(tdsalpha,tdsbeta,tdsgamma,theta,x3[i])
}


x1<-x1[ymmtd >= 0.05 & ymmtd <= 0.3]
ymmtd<-ymmtd[ymmtd >= 0.05 & ymmtd <= 0.3]
kk1<-length(x1)

x2<-x2[yusmtd >= 0.05 & yusmtd <= 0.3]
yusmtd<-yusmtd[yusmtd >= 0.05 & yusmtd <= 0.3]
kk2<-length(x2)

x3<-x3[ydsmtd >= 0.05 & ydsmtd <= 0.3]
ydsmtd<-ydsmtd[ydsmtd >= 0.05 & ydsmtd <= 0.3]
kk3<-length(x3)


lm1<-matrix(nrow=M,ncol=kk1)
lm2=matrix(nrow=M,ncol=kk1)
lm3=matrix(nrow=M,ncol=kk1)
lm4=matrix(nrow=M,ncol=kk1)
lm5=matrix(nrow=M,ncol=kk1)

lus1<-matrix(nrow=M,ncol=kk2)
lus2<-matrix(nrow=M,ncol=kk2)
lus3<-matrix(nrow=M,ncol=kk2)
lus4<-matrix(nrow=M,ncol=kk2)
lus5<-matrix(nrow=M,ncol=kk2)

lds1<-matrix(nrow=M,ncol=kk3)
lds2<-matrix(nrow=M,ncol=kk3)
lds3<-matrix(nrow=M,ncol=kk3)
lds4<-matrix(nrow=M,ncol=kk3)
lds5<-matrix(nrow=M,ncol=kk3)


#Scenario1

for (i in 1:M){
  for (j in 1 :kk1){
    
    #eta=0.00
    
    aal1<- fdist(seq(0.05,0.3,by=0.001),c(x1[j],ymmtd[j]),mvalphaeta0[i,1],mvbetaeta0[i,1],mvgammaeta0[i,1],theta)
    if(!all(is.nan(aal1))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal1)]
      objective = aal1[which.min(aal1)]
      
      lm1[i,j]<- objective
      if (ymmtd[j] > mtd_logistic(mvalphaeta0[i,1],mvbetaeta0[i,1],mvgammaeta0[i,1],theta,minimum)) lm1[i,j] <- -lm1[i,j]
    }
    
    #eta=0.10
    
    aal2<- fdist(seq(0.05,0.3,by=0.001),c(x1[j],ymmtd[j]),mvalphaeta01[i,1],mvbetaeta01[i,1],mvgammaeta01[i,1],theta)
    if(!all(is.nan(aal2))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal2)]
      objective = aal2[which.min(aal2)]
      
      lm2[i,j]<- objective
      if (ymmtd[j] > mtd_logistic(mvalphaeta01[i,1],mvbetaeta01[i,1],mvgammaeta01[i,1],theta,minimum)) lm2[i,j] <- -lm2[i,j]
    }
    
    #eta=0.25
    
    aal3=fdist(seq(0.05,0.3,by=0.001),c(x1[j],ymmtd[j]),mvalphaeta025[i,1],mvbetaeta025[i,1],mvgammaeta025[i,1],theta)
    if(!all(is.nan(aal3))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal3)]
      objective = aal3[which.min(aal3)]
      
      lm3[i,j]<- objective
      if (ymmtd[j] > mtd_logistic(mvalphaeta025[i,1],mvbetaeta025[i,1],mvgammaeta025[i,1],theta,minimum)) lm3[i,j] <- -lm3[i,j]
    }
    
    #eta=0.40
    
    aal4=fdist(seq(0.05,0.3,by=0.001),c(x1[j],ymmtd[j]),mvalphaeta04[i,1],mvbetaeta04[i,1],mvgammaeta04[i,1],theta)
    if(!all(is.nan(aal4))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal4)]
      objective = aal4[which.min(aal4)]
      
      lm4[i,j]<- objective
      if (ymmtd[j] > mtd_logistic(mvalphaeta04[i,1],mvbetaeta04[i,1],mvgammaeta04[i,1],theta,minimum)) lm4[i,j] <- -lm4[i,j]
    }
    
    #eta=1.0
    
    aal5=fdist(seq(0.05,0.3,by=0.001),c(x1[j],ymmtd[j]),mvalphaeta1[i,1],mvbetaeta1[i,1],mvgammaeta1[i,1],theta)
    if(!all(is.nan(aal5))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal5)]
      objective = aal5[which.min(aal5)]
      
      lm5[i,j]<- objective
      if (ymmtd[j] > mtd_logistic(mvalphaeta1[i,1],mvbetaeta1[i,1],mvgammaeta1[i,1],theta,minimum)) lm5[i,j] <- -lm5[i,j]
    }
    
    
  }
}

lm1=lm1[complete.cases(lm1), ]
lm2=lm2[complete.cases(lm2), ]
lm3=lm3[complete.cases(lm3), ]
lm4=lm4[complete.cases(lm4), ]
lm5=lm5[complete.cases(lm4), ]


#Scenario 2

for (i in 1:M){
  for (j in 1 :kk2){
    
    #eta=0.00
    
    aal1<- fdist(seq(0.05,0.3,by=0.001),c(x2[j],yusmtd[j]),usvalphaeta0[i,1],usvbetaeta0[i,1],usvgammaeta0[i,1],theta)
    if(!all(is.nan(aal1))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal1)]
      objective = aal1[which.min(aal1)]
      
      lus1[i,j]<- objective
      if (yusmtd[j] > mtd_logistic(usvalphaeta0[i,1],usvbetaeta0[i,1],usvgammaeta0[i,1],theta,minimum)) lus1[i,j] <- -lus1[i,j]
    }
    
    #eta=0.1
    
    aal2<- fdist(seq(0.05,0.3,by=0.001),c(x2[j],yusmtd[j]),usvalphaeta01[i,1],usvbetaeta01[i,1],usvgammaeta01[i,1],theta)
    if(!all(is.nan(aal2))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal2)]
      objective = aal2[which.min(aal2)]
      
      lus2[i,j]<- objective
      if (yusmtd[j] > mtd_logistic(usvalphaeta01[i,1],usvbetaeta01[i,1],usvgammaeta01[i,1],theta,minimum)) lus2[i,j] <- -lus2[i,j]
    }
    #eta=0.25
    
    aal3=fdist(seq(0.05,0.3,by=0.001),c(x2[j],yusmtd[j]),usvalphaeta025[i,1],usvbetaeta025[i,1],usvgammaeta025[i,1],theta)
    if(!all(is.nan(aal3))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal3)]
      objective = aal3[which.min(aal3)]
      
      lus3[i,j]<- objective
      if (yusmtd[j] > mtd_logistic(usvalphaeta025[i,1],usvbetaeta025[i,1],usvgammaeta025[i,1],theta,minimum)) lus3[i,j] <- -lus3[i,j]
    }
    
    #eta=0.4
    
    aal4=fdist(seq(0.05,0.3,by=0.001),c(x2[j],yusmtd[j]),usvalphaeta04[i,1],usvbetaeta04[i,1],usvgammaeta04[i,1],theta)
    if(!all(is.nan(aal4))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal4)]
      objective = aal4[which.min(aal4)]
      
      lus4[i,j]<- objective
      if (yusmtd[j] > mtd_logistic(usvalphaeta04[i,1],usvbetaeta04[i,1],usvgammaeta04[i,1],theta,minimum)) lus4[i,j] <- -lus4[i,j]
    }
    
    
    #eta=1.0
    
    aal5=fdist(seq(0.05,0.3,by=0.001),c(x2[j],yusmtd[j]),usvalphaeta1[i,1],usvbetaeta1[i,1],usvgammaeta1[i,1],theta)
    if(!all(is.nan(aal5))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal5)]
      objective = aal5[which.min(aal5)]
      
      lus5[i,j]<- objective
      if (yusmtd[j] > mtd_logistic(usvalphaeta1[i,1],usvbetaeta1[i,1],usvgammaeta1[i,1],theta,minimum)) lus5[i,j] <- -lus5[i,j]
    }
    
    
  }
}

lus1=lus1[complete.cases(lus1), ]
lus2=lus2[complete.cases(lus2), ]
lus3=lus3[complete.cases(lus3), ]
lus4=lus4[complete.cases(lus4), ]
lus5=lus5[complete.cases(lus4), ]

#Scenario 3

for (i in 1:M){
  for (j in 1 :kk3){
    
    #eta=0.00
    
    aal1<- fdist(seq(0.05,0.3,by=0.001),c(x3[j],ydsmtd[j]),dsvalphaeta0[i,1],dsvbetaeta0[i,1],dsvgammaeta0[i,1],theta)
    if(!all(is.nan(aal1))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal1)]
      objective = aal1[which.min(aal1)]
      
      lds1[i,j]<- objective
      if (ydsmtd[j] > mtd_logistic(dsvalphaeta0[i,1],dsvbetaeta0[i,1],dsvgammaeta0[i,1],theta,minimum)) lds1[i,j] <- -lds1[i,j]
    }
    
    #eta=0.10
    
    aal2<- fdist(seq(0.05,0.3,by=0.001),c(x3[j],ydsmtd[j]),dsvalphaeta01[i,1],dsvbetaeta01[i,1],dsvgammaeta01[i,1],theta)
    if(!all(is.nan(aal2))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal2)]
      objective = aal2[which.min(aal2)]
      
      lds2[i,j]<- objective
      if (ydsmtd[j] > mtd_logistic(dsvalphaeta01[i,1],dsvbetaeta01[i,1],dsvgammaeta01[i,1],theta,minimum)) lds2[i,j] <- -lds2[i,j]
    }
    
    #eta=0.25
    
    aal3=fdist(seq(0.05,0.3,by=0.001),c(x3[j],ydsmtd[j]),dsvalphaeta025[i,1],dsvbetaeta025[i,1],dsvgammaeta025[i,1],theta)
    if(!all(is.nan(aal3))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal3)]
      objective = aal3[which.min(aal3)]
      
      lds3[i,j]<- objective
      if (ydsmtd[j] > mtd_logistic(dsvalphaeta025[i,1],dsvbetaeta025[i,1],dsvgammaeta025[i,1],theta,minimum)) lds3[i,j] <- -lds3[i,j]
    }
    
    #eta=0.4
    
    aal4=fdist(seq(0.05,0.3,by=0.001),c(x3[j],ydsmtd[j]),dsvalphaeta04[i,1],dsvbetaeta04[i,1],dsvgammaeta04[i,1],theta)
    if(!all(is.nan(aal4))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal4)]
      objective = aal4[which.min(aal4)]
      
      lds4[i,j]<- objective
      if (ydsmtd[j] > mtd_logistic(dsvalphaeta04[i,1],dsvbetaeta04[i,1],dsvgammaeta04[i,1],theta,minimum)) lds4[i,j] <- -lds4[i,j]
    }
    
    #eta=1.0
    
    aal5=fdist(seq(0.05,0.3,by=0.001),c(x3[j],ydsmtd[j]),dsvalphaeta1[i,1],dsvbetaeta1[i,1],dsvgammaeta1[i,1],theta)
    if(!all(is.nan(aal5))){
      minimum = seq(0.05,0.3,by=0.001)[which.min(aal5)]
      objective = aal5[which.min(aal5)]
      
      lds5[i,j]<- objective
      if (ydsmtd[j] > mtd_logistic(dsvalphaeta1[i,1],dsvbetaeta1[i,1],dsvgammaeta1[i,1],theta,minimum)) lds5[i,j] <- -lds5[i,j]
    }
    
  }
}

lds1=lds1[complete.cases(lds1), ]
lds2=lds2[complete.cases(lds2), ]
lds3=lds3[complete.cases(lds3), ]
lds4=lds4[complete.cases(lds4), ]
lds5=lds5[complete.cases(lds4), ]


tol11<-numeric()
tol12<-numeric()
for (j in 1:kk1){
  tol11[j]<-p*((x1[j]^2+ymmtd[j]^2)^0.5)
  tol12[j]<-0.1*((x1[j]^2+ymmtd[j]^2)^0.5)
}

tol21<-numeric()
tol22<-numeric()
for (j in 1:kk2){
  tol21[j]<-p*((x2[j]^2+yusmtd[j]^2)^0.5)
  tol22[j]<-0.1*((x2[j]^2+yusmtd[j]^2)^0.5)
}

tol31<-numeric()
tol32<-numeric()
for (j in 1:kk3){
  tol31[j]<-p*((x3[j]^2+ydsmtd[j]^2)^0.5)
  tol32[j]<-0.1*((x3[j]^2+ydsmtd[j]^2)^0.5)
}


mpercent11<-numeric()
mpercent21<-numeric()
mpercent31<-numeric()
mpercent41<-numeric()
mpercent51<-numeric()
mpercent12<-numeric()
mpercent22<-numeric()
mpercent32<-numeric()
mpercent42<-numeric()
mpercent52<-numeric()

uspercent11<-numeric()
uspercent21<-numeric()
uspercent31<-numeric()
uspercent41<-numeric()
uspercent51<-numeric()
uspercent12<-numeric()
uspercent22<-numeric()
uspercent32<-numeric()
uspercent42<-numeric()
uspercent52<-numeric()

dspercent11<-numeric()
dspercent21<-numeric()
dspercent31<-numeric()
dspercent41<-numeric()
dspercent51<-numeric()
dspercent12<-numeric()
dspercent22<-numeric()
dspercent32<-numeric()
dspercent42<-numeric()
dspercent52<-numeric()


for (i in 1 :kk1){
  mpercent11[i]<-(length(lm1[,i][lm1[,i] <= tol11[i]])/M)*100
  mpercent21[i]<-(length(lm2[,i][lm2[,i] <= tol11[i]])/M)*100
  mpercent31[i]<-(length(lm3[,i][lm3[,i] <= tol11[i]])/M)*100
  mpercent41[i]<-(length(lm4[,i][lm4[,i] <= tol11[i]])/M)*100
  mpercent51[i]<-(length(lm5[,i][lm5[,i] <= tol11[i]])/M)*100
  
  mpercent12[i]<-(length(lm1[,i][lm1[,i] <= tol12[i]])/M)*100
  mpercent22[i]<-(length(lm2[,i][lm2[,i] <= tol12[i]])/M)*100
  mpercent32[i]<-(length(lm3[,i][lm3[,i] <= tol12[i]])/M)*100
  mpercent42[i]<-(length(lm4[,i][lm4[,i] <= tol12[i]])/M)*100
  mpercent52[i]<-(length(lm5[,i][lm5[,i] <= tol12[i]])/M)*100
}

for (i in 1 :kk2){
  uspercent11[i]<-(length(lus1[,i][lus1[,i] <= tol21[i]])/M)*100
  uspercent21[i]<-(length(lus2[,i][lus2[,i] <= tol21[i]])/M)*100
  uspercent31[i]<-(length(lus3[,i][lus3[,i] <= tol21[i]])/M)*100
  uspercent41[i]<-(length(lus4[,i][lus4[,i] <= tol21[i]])/M)*100
  uspercent51[i]<-(length(lus5[,i][lus5[,i] <= tol21[i]])/M)*100
  
  uspercent12[i]<-(length(lus1[,i][lus1[,i] <= tol22[i]])/M)*100
  uspercent22[i]<-(length(lus2[,i][lus2[,i] <= tol22[i]])/M)*100
  uspercent32[i]<-(length(lus3[,i][lus3[,i] <= tol22[i]])/M)*100
  uspercent42[i]<-(length(lus4[,i][lus4[,i] <= tol22[i]])/M)*100
  uspercent52[i]<-(length(lus5[,i][lus5[,i] <= tol22[i]])/M)*100
}


for (i in 1 :kk3){
  dspercent11[i]<-(length(lds1[,i][lds1[,i] <= tol31[i]])/M)*100
  dspercent21[i]<-(length(lds2[,i][lds2[,i] <= tol31[i]])/M)*100
  dspercent31[i]<-(length(lds3[,i][lds3[,i] <= tol31[i]])/M)*100
  dspercent41[i]<-(length(lds4[,i][lds4[,i] <= tol31[i]])/M)*100
  dspercent51[i]<-(length(lds5[,i][lds5[,i] <= tol31[i]])/M)*100
  
  dspercent12[i]<-(length(lds1[,i][lds1[,i] <= tol32[i]])/M)*100
  dspercent22[i]<-(length(lds2[,i][lds2[,i] <= tol32[i]])/M)*100
  dspercent32[i]<-(length(lds3[,i][lds3[,i] <= tol32[i]])/M)*100
  dspercent42[i]<-(length(lds4[,i][lds4[,i] <= tol32[i]])/M)*100
  dspercent52[i]<-(length(lds5[,i][lds5[,i] <= tol32[i]])/M)*100
  
}


# Plot

tiff(file="../figure5.tiff",width=10,height=4,units="in",res=300)  

par(mfrow=c(1,3),mar=c(5.1, 5.1, 4.1, 2.1), mgp=c(3, 1.25, 0))

#Scenario 1

plot(x2,uspercent11,type="l",ylim=c(0,100),xlab="Drug D1", ylab="Percent",cex.lab=2,cex.axis=2, col="black",lwd=2,cex.main=1.5,main="Scenario 1")
points(x2,uspercent21,type="l",lty=1, col="grey40",lwd=4)
points(x2,uspercent31,type="l",lty=1, col="grey70",lwd=4)
points(x2,uspercent41,type="l",lty=1, col="grey90",lwd=4)
points(x2,uspercent12,type="l",lty=2, col="black",lwd=4)
points(x2,uspercent22,type="l",lty=2, col="grey40",lwd=4)
points(x2,uspercent32,type="l",lty=2, col="grey70",lwd=4)
points(x2,uspercent42,type="l",lty=2, col="grey90",lwd=4)

#Scenario 2

plot(x1,mpercent11,type="l",ylim=c(0,100),xlab="Drug D1", ylab="",cex.lab=2,cex.axis=2, col="black",lwd=2,cex.main=1.5,main="Scenario 2")
points(x1,mpercent21,type="l",lty=1, col="grey40",lwd=4)
points(x1,mpercent31,type="l",lty=1, col="grey70",lwd=4)
points(x1,mpercent41,type="l",lty=1, col="grey90",lwd=4)
points(x1,mpercent12,type="l",lty=2, col="black",lwd=4)
points(x1,mpercent22,type="l",lty=2, col="grey40",lwd=4)
points(x1,mpercent32,type="l",lty=2, col="grey70",lwd=4)
points(x1,mpercent42,type="l",lty=2, col="grey90",lwd=4)

#Scenario 3

plot(x3,dspercent11,type="l",ylim=c(0,100),xlab="Drug D1", ylab="",cex.lab=2,cex.axis=2, col="black",lwd=2,cex.main=1.5,main="Scenario 3")
points(x3,dspercent21,type="l",lty=1, col="grey40",lwd=4)
points(x3,dspercent31,type="l",lty=1, col="grey70",lwd=4)
points(x3,dspercent41,type="l",lty=1, col="grey90",lwd=4)
points(x3,dspercent12,type="l",lty=2, col="black",lwd=4)
points(x3,dspercent22,type="l",lty=2, col="grey40",lwd=4)
points(x3,dspercent32,type="l",lty=2, col="grey70",lwd=4)
points(x3,dspercent42,type="l",lty=2, col="grey90",lwd=4)

legend(x = 0.177,y=42,
       legend = c(expression(paste(eta, "=0.00 ")), expression(paste(eta,"=0.10 ")), expression(paste(eta,"=0.25 ")),expression(paste(eta,"=0.40 "))), 
       col=c("black","grey40","grey70","grey90"), lwd=3.5, cex=1.6, horiz = FALSE,box.col = "white")

dev.off()






