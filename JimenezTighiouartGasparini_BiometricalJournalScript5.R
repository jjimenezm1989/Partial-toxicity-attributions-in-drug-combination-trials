#Article: Cancer phase I trials design using drug combinations when a fraction of dose limiting toxicities is attributable to one or mode agents
#Authors: Jose L. Jimenez, Mourad Tighiouart and Mauro Gasparini.
#E-mail addresses: jose.jimenez@polito.it, Mourad.Tighiouart@cshs.org, gasparini@calvino.polito.it
#File: JimenezTighiouartGasparini_BiometricalJournalScript5.R

####################################################
# Code to obtain Figures 2 and 3 in the manuscript #
####################################################

#This functions computes the MTD curve
MTDcurves = function(alpha,beta,gamma,x,theta,contour){
  y = mat.or.vec(1,length(x))
  for(i in 1:length(x)){
    xx = x[i]^alpha
    xx2 = x[i]^(2*alpha)
    xxm = x[i]^(-alpha)
    A = (exp(-gamma)-1)/(exp(-gamma)+1)
    root = sqrt((-A*xx2+A*xx+xx-1)^2-4*A*xx*(xx-1)*(theta-xx-contour))
    
    y[i] = 2^(-1/beta)*(xxm/(A*(xx-1))*(A*xx2-A*xx-root-xx+1))^(1/beta)
  }
  if(is.nan(y[i])){
    y[i] = 0
  }
  return(y)
}

theta = 0.3
x=seq(0.05,0.3,by=0.0001)
sdev = sd(0.3^runif(10000,0.2,2))

alpha1 = 0.5*(0.2+2)
beta1 = 0.5*(0.2+2)
gamma1 = 0.1/0.1

######
#MEAN#
######

#True MTD curve and contours

y1 = MTDcurves(alpha1,beta1,gamma1,x,theta,contour=0)
y2 = MTDcurves(alpha1,beta1,gamma1,x,theta,contour=0.05)
y3 = MTDcurves(alpha1,beta1,gamma1,x,theta,contour=-0.05)
y4 = MTDcurves(alpha1,beta1,gamma1,x,theta,contour=0.1)
y5 = MTDcurves(alpha1,beta1,gamma1,x,theta,contour=-0.1)

#The values of alpha, beta and gamma are the ones we obtained and remain here ONLY as an example.

#eta=0.00
alpha6 = 1.164443
beta6 = 1.17301
gamma6 = 0.1507451
y6 = MTDcurves(alpha6,beta6,gamma6,x,theta,contour=0)

#eta=0.10
alpha7 = 1.162017
beta7 = 1.156905
gamma7 = 0.0846982
y7 = MTDcurves(alpha7,beta7,gamma7,x,theta,contour=0)

#eta=0.25
alpha8 = 1.134876
beta8 = 1.132421
gamma8 = 0.09548954
y8 = MTDcurves(alpha8,beta8,gamma8,x,theta,contour=0)

#eta=0.40
alpha9 = 1.106895
beta9 = 1.109093
gamma9 = 0.06032829
y9 = MTDcurves(alpha9,beta9,gamma9,x,theta,contour=0)

#eta=1.00
alpha10 = 1.017307
beta10 = 1.013429
gamma10 = 0.01092717
y10 = MTDcurves(alpha10,beta10,gamma10,x,theta,contour=0)

#########
#-1sigma#
#########

alpha11 = 0.5*(0.2+2) - sdev
beta11 = 0.5*(0.2+2) - sdev
gamma11 = 0.1/0.1
y11 = mat.or.vec(1,length(x))

#True MTD curve and contours

y11 = MTDcurves(alpha11,beta11,gamma11,x,theta,contour=0)
y12 = MTDcurves(alpha11,beta11,gamma11,x,theta,contour=0.05)
y13 = MTDcurves(alpha11,beta11,gamma11,x,theta,contour=-0.05)
y14 = MTDcurves(alpha11,beta11,gamma11,x,theta,contour=0.1)
y15 = MTDcurves(alpha11,beta11,gamma11,x,theta,contour=-0.1)

#eta=0.00
alpha16 = 0.9917704
beta16 = 0.9873333
gamma16 = 0.1762382
y16 = MTDcurves(alpha16,beta16,gamma16,x,theta,contour=0)

#eta=0.10
alpha17 = 0.9836879
beta17 = 0.9896146
gamma17 = 0.1259297
y17 = MTDcurves(alpha17,beta17,gamma17,x,theta,contour=0)

#eta=0.25
alpha18 = 0.9540791
beta18 = 0.9514816
gamma18 = 0.08298481
y18 = MTDcurves(alpha18,beta18,gamma18,x,theta,contour=0)

#eta=0.40
alpha19 = 0.9307474
beta19 = 0.9218087
gamma19 = 0.06414504
y19 = MTDcurves(alpha19,beta19,gamma19,x,theta,contour=0)

#eta=1.00
alpha20 = 0.8488039
beta20 = 0.8361285
gamma20 = 0.01482762
y20 = MTDcurves(alpha20,beta20,gamma20,x,theta,contour=0)

#########
#+1sigma#
#########

alpha21 = 0.5*(0.2+2) + sdev
beta21 = 0.5*(0.2+2) + sdev
gamma21 = 0.1/0.1
y21 = mat.or.vec(1,length(x))

y21 = MTDcurves(alpha21,beta21,gamma21,x,theta,contour=0)
y22 = MTDcurves(alpha21,beta21,gamma21,x,theta,contour=0.05)
y23 = MTDcurves(alpha21,beta21,gamma21,x,theta,contour=-0.05)
y24 = MTDcurves(alpha21,beta21,gamma21,x,theta,contour=0.1)
y25 = MTDcurves(alpha21,beta21,gamma21,x,theta,contour=-0.1)

#eta=0.00
alpha26 = 1.330198
beta26 = 1.333697
gamma26 = 0.09265048
y26 = MTDcurves(alpha26,beta26,gamma26,x,theta,contour=0)

#eta=0.10
alpha27 = 1.325436
beta27 = 1.329625
gamma27 = 0.09209786
y27 = MTDcurves(alpha27,beta27,gamma27,x,theta,contour=0)

#eta=0.25
alpha28 = 1.294679
beta28 = 1.305604
gamma28 = 0.09245929
y28 = MTDcurves(alpha28,beta28,gamma28,x,theta,contour=0)

#eta=0.40
alpha29 = 1.277103
beta29 = 1.276137
gamma29 = 0.05051354
y29 = MTDcurves(alpha29,beta29,gamma29,x,theta,contour=0)

#eta=1.00
alpha30 = 1.192004
beta30 = 1.177363
gamma30 = 0.01389495
y30 = MTDcurves(alpha30,beta30,gamma30,x,theta,contour=0)



#######
#PLOTS#
#######

tiff(file="../figure2.tiff",width=10,height=3.75,units="in",res=300)

par(mfrow=c(1,3),mar=c(5.1, 5.1, 4.1, 2.1), mgp=c(3, 1.25, 0))

plot(x[which(y11>=0.05  & y11 <= 0.3)],y11[which(y11>=0.05 & y11 <= 0.3)],type="line",xlim=c(0.05,0.3),ylim=c(0.05,0.3),cex.lab=2,cex.axis=1.5,
     xlab="Drug D1", ylab="Drug D2",lwd=2,cex.main=1.5,xaxs="i",yaxs="i",main="Scenario 1",col="black",lty=2)
lines(x[which(y12>=0.05 & y12 <= 0.3)],y12[which(y12>=0.05 & y12 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y13>=0.05 & y13 <= 0.3)],y13[which(y13>=0.05 & y13 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y14>=0.05 & y14 <= 0.3)],y14[which(y14>=0.05 & y14 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y15>=0.05 & y15 <= 0.3)],y15[which(y15>=0.05 & y15 <= 0.3)],type="line",col="gray",lty=4,lwd=2)


plot(x[which(y1>=0.05  & y1 <= 0.3)],y1[which(y1>=0.05  & y1 <= 0.3)],type="line",xlim=c(0.05,0.3),ylim=c(0.05,0.3),cex.lab=2,cex.axis=1.5,
     xlab="Drug D1", ylab="",lwd=2,cex.main=1.5,xaxs="i",yaxs="i",main="Scenario 2",lty=2)
lines(x[which(y2>=0.05 & y2 <= 0.3)],y2[which(y2>=0.05 & y2 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y3>=0.05  & y3 <= 0.3)],y3[which(y3>=0.05  & y3 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y4>=0.05  & y4 <= 0.3)],y4[which(y4>=0.05  & y4 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y5>=0.05  & y5 <= 0.3)],y5[which(y5>=0.05  & y5 <= 0.3)],type="line",col="gray",lty=4,lwd=2)


plot(x[which(y21>=0.05 & y21 <= 0.3)],y21[which(y21>=0.05 & y21 <= 0.3)],type="line",xlim=c(0.05,0.3),ylim=c(0.05,0.3),cex.lab=2,cex.axis=1.5,
     xlab="Drug D1", ylab="",lwd=2,cex.main=1.5, xaxs="i",yaxs="i",main="Scenario 3",lty=2)
lines(x[which(y22>=0.05  & y22 <= 0.3)],y22[which(y22>=0.05 & y22 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y23>=0.05 & y23 <= 0.3)],y23[which(y23>=0.05 & y23 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y24>=0.05 & y24 <= 0.3)],y24[which(y24>=0.05 & y24 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y25>=0.05 & y25 <= 0.3)],y25[which(y25>=0.05 & y25 <= 0.3)],type="line",col="gray",lty=4,lwd=2)

dev.off()


tiff(file="../figure3.tiff",width=10,height=3.75,units="in",res=300)

par(mfrow=c(1,3),mar=c(5.1, 5.1, 4.1, 2.1), mgp=c(3, 1.25, 0))

plot(x[which(y11>=0.05  & y11 <= 0.3)],y11[which(y11>=0.05 & y11 <= 0.3)],type="line",xlim=c(0.05,0.3),ylim=c(0.05,0.3),cex.lab=2,cex.axis=1.5,
     xlab="Drug D1", ylab="Drug D2",lwd=2,cex.main=1.5,xaxs="i",yaxs="i",main="Scenario 1",col="black",lty=2)
lines(x[which(y12>=0.05 & y12 <= 0.3)],y12[which(y12>=0.05 & y12 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y13>=0.05 & y13 <= 0.3)],y13[which(y13>=0.05 & y13 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y14>=0.05 & y14 <= 0.3)],y14[which(y14>=0.05 & y14 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y15>=0.05 & y15 <= 0.3)],y15[which(y15>=0.05 & y15 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y16>=0.05 & y16 <= 0.3)],y16[which(y16>=0.05 & y16 <= 0.3)],type="line",col="black",lty=1,lwd=4)
lines(x[which(y17>=0.05 & y17 <= 0.3)],y17[which(y17>=0.05 & y17 <= 0.3)],type="line",col="gray40",lty=1,lwd=4)
lines(x[which(y18>=0.05 & y18 <= 0.3)],y18[which(y18>=0.05 & y18 <= 0.3)],type="line",col="gray70",lty=1,lwd=4)
lines(x[which(y19>=0.05 & y19 <= 0.3)],y19[which(y19>=0.05 & y19 <= 0.3)],type="line",col="gray90",lty=1,lwd=4)


plot(x[which(y1>=0.05  & y1 <= 0.3)],y1[which(y1>=0.05  & y1 <= 0.3)],type="line",xlim=c(0.05,0.3),ylim=c(0.05,0.3),cex.lab=2,cex.axis=1.5,
     xlab="Drug D1", ylab="",lwd=2,cex.main=1.5,xaxs="i",yaxs="i",main="Scenario 2",lty=2)
lines(x[which(y2>=0.05 & y2 <= 0.3)],y2[which(y2>=0.05 & y2 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y3>=0.05  & y3 <= 0.3)],y3[which(y3>=0.05  & y3 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y4>=0.05  & y4 <= 0.3)],y4[which(y4>=0.05  & y4 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y5>=0.05  & y5 <= 0.3)],y5[which(y5>=0.05  & y5 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y6>=0.05  & y6 <= 0.3)],y6[which(y6>=0.05  & y6 <= 0.3)],type="line",col="black",lty=1,lwd=4)
lines(x[which(y7>=0.05  & y7 <= 0.3)],y7[which(y7>=0.05  & y7 <= 0.3)],type="line",col="gray40",lty=1,lwd=4)
lines(x[which(y8>=0.05  & y8 <= 0.3)],y8[which(y8>=0.05  & y8 <= 0.3)],type="line",col="gray70",lty=1,lwd=4)
lines(x[which(y9>=0.05  & y9 <= 0.3)],y9[which(y9>=0.05  & y9 <= 0.3)],type="line",col="gray90",lty=1,lwd=4)



plot(x[which(y21>=0.05 & y21 <= 0.3)],y21[which(y21>=0.05 & y21 <= 0.3)],type="line",xlim=c(0.05,0.3),ylim=c(0.05,0.3),cex.lab=2,cex.axis=1.5,
     xlab="Drug D1", ylab="",lwd=2,cex.main=1.5, xaxs="i",yaxs="i",main="Scenario 3",lty=2)
lines(x[which(y22>=0.05  & y22 <= 0.3)],y22[which(y22>=0.05 & y22 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y23>=0.05 & y23 <= 0.3)],y23[which(y23>=0.05 & y23 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y24>=0.05 & y24 <= 0.3)],y24[which(y24>=0.05 & y24 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y25>=0.05 & y25 <= 0.3)],y25[which(y25>=0.05 & y25 <= 0.3)],type="line",col="gray",lty=4,lwd=2)
lines(x[which(y26>=0.05 & y26 <= 0.3)],y26[which(y26>=0.05 & y26 <= 0.3)],type="line",col="black",lty=1,lwd=4)
lines(x[which(y27>=0.05 & y27 <= 0.3)],y27[which(y27>=0.05 & y27 <= 0.3)],type="line",col="gray40",lty=1,lwd=4)
lines(x[which(y28>=0.05 & y28 <= 0.3)],y28[which(y28>=0.05 & y28 <= 0.3)],type="line",col="gray70",lty=1,lwd=4)
lines(x[which(y29>=0.05 & y29 <= 0.3)],y29[which(y29>=0.05 & y29 <= 0.3)],type="line",col="gray90",lty=1,lwd=4)

legend(x = 0.052,y=0.155,
       legend = c(expression(paste(eta, "=0.00 ")), expression(paste(eta,"=0.10 ")), expression(paste(eta,"=0.25 ")),expression(paste(eta,"=0.40 "))), 
       col=c("black","gray40","gray70","gray90"), lwd=3.5, cex=1.6, horiz = FALSE,box.col = "white")


dev.off()