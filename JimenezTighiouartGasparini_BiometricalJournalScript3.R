#Article: Cancer phase I trials design using drug combinations when a fraction of dose limiting toxicities is attributable to one or mode agents
#Authors: Jose L. Jimenez, Mourad Tighiouart and Mauro Gasparini.
#E-mail addresses: jose.jimenez@polito.it, Mourad.Tighiouart@cshs.org, gasparini@calvino.polito.it
#File: JimenezTighiouartGasparini_BiometricalJournalScript3.R

################################################################################
#This program computes operating characteristics of the design. See readme.txt #
################################################################################

#############
# functions #
#############

#Function to obtain MTD set as in equation (10) in the manuscript
mtd_set<-function(alpha,beta,gamma,x,theta){
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

#Function to recommend the dose combinations closer to the MTD
MTDdoses = function(alpha, beta, gamma, theta, nx, ny) {
  
  # special senarios
  # do not recommend MTD if rho00 > theta + 0.1
  if(mtd_set(alpha,beta,gamma,x=0.05,theta) < 0.05) {
    doses = data.frame(x=0.05, y=0.05)
  } else if(mtd_set(alpha,beta,gamma,x=0.3,theta) > 0.3) {
    doses = data.frame(x=0.3, y=0.3)
  } else {
    
    # set number of descrete levels for drug X and Y
    x1 = seq(0.05,0.3, length.out=nx)
    y1 = seq(0.05,0.3, length.out=ny)
    
    # find intercept
    xmin = 0.05
    xmax = 0.3
    
    if(mtd_set(alpha,beta,gamma,x=0.05,theta) > 0.3) xmin=uniroot(function(x) mtd_set(alpha,beta,gamma,x,theta) - 0.3, c(0.05,0.3),extendInt="downX")$root
    
    if(mtd_set(alpha,beta,gamma,x=0.3,theta) < 0.05) xmax=uniroot(function(x) mtd_set(alpha,beta,gamma,x,theta),c(0.05,0.3),extendInt="downX")$root
    
    # Define the function to be minimized
    fdist<-function(x,pt,alpha,beta,gamma,theta){
      # pt is the point
      d<-(pt[1]-x)^2+(pt[2]- mtd_set(alpha,beta,gamma,x,theta))^2
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

#Function to compute probability of DLT
probDLT = function(alpha,beta,gamma,x,y){
  
  pi10 = x^alpha*(1-y^beta) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi01 = y^beta*(1-x^alpha) - x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  pi11 = x^alpha*y^beta + x^alpha*(1-x^alpha)*y^beta*(1-y^beta)*((exp(-gamma)-1)/(exp(-gamma)+1))
  
  p = pi10 + pi01 + pi11
  
  return(p)
  
}


#######################
#Main part of the code#
#######################

M = 1000
theta = 0.3
delta = 0.7
nx=4
ny=4

alphas = read.table("../alpha.txt")[,1]
betas = read.table("../beta.txt")[,1]
gammas = read.table("../gamma.txt")[,1]
postdlts = read.csv("../postdlts.csv")
postdlts = postdlts[postdlts[,3]<delta,]

result = data.frame(x=NA,y=NA,trial=NA)


for (i in 1:M) {
  
  a = MTDdoses(alphas[i], betas[i], gammas[i], theta, nx, ny)
  
  if(nrow(a)>0) a$trial = i
  
  result = rbind(result, a)
  
  print(i)
  
}

# remove NAs (including NA in the first row)
result = result[!is.na(result$x),]

# summary
result2 = data.frame(x = 1:(nx*ny), y=NA, percent=NA)
x1 = seq(0.05,0.3, length.out=nx)
y1 = seq(0.05,0.3, length.out=ny)

k=0
for (i in 1:nx) {
  for (j in 1:ny) {
    k=k+1
    result2[k,] = c(x1[i], y1[j],sum(result$x==x1[i] & result$y==y1[j])*100/M)
  }
}


############################
#We now compute the results#
############################

# calculate the probability that a prospective trial will identify at least one MTD

true1 = mat.or.vec(nx*ny,4)
true1[,4]=1:(length(x1)*length(y1))

#True parameter values according to the prior distributions used.
talpha = 0.5*(0.2+2)
tbeta = 0.5*(0.2+2)
tgamma = 0.1/0.1

count=1

for(i in 1:length(x1)){
  for(j in 1:length(y1)){
    true1[count,1] = x1[i]
    true1[count,2] = y1[j]
    true1[count,3] = probDLT(talpha,tbeta,tgamma,x1[i],y1[j])
    count = count + 1
  }
}

true = true1[abs(true1[,3]-theta)<=0.1,]
colnames(true) = c("x","y","prob","id")

#Code to obtain Table 4 in the manuscript

#percent of times that at leas 25,50,75 and 100% of the recommended doses belong to the true MTD set

all = mat.or.vec(1,1000)
atleast25 = mat.or.vec(1,1000)
atleast50 = mat.or.vec(1,1000)
atleast75 = mat.or.vec(1,1000)

for (i in 1:1000){
  
  num.doses.recommended = nrow(result[result$trial==i,])
  
  num.true.doses.recommended = nrow(merge(result[result$trial==i,],as.data.frame(true[,-3]),by=c("x","y")))
  
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

####################################
#Summary % of correct MTD selection#
####################################

mean(atleast25)
mean(atleast50)
mean(atleast75)
mean(all)



