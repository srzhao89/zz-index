######################################################################
## Statistical Inference for Price or Quantity Indices
## Author: Valentin Zelenyuk and Shirong Zhao
## Email:  shironz@163.com
## Dalian, July 13, 2026
######################################################################
rm(list=ls())
require(fdrtool)
#################################
source("./Functions/data.sim.R")
#################################
#### Set Seed ####
if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)

####################################
####################################
index="All"
###################################
####################################
###################################
###################################
q=100
ETA=c(0,0.5,1,2,5,10)
Sigma.varepsilon=c(0.02,0.05,0.1,0.2,0.3)
sigma.p=0.01
diff.res=matrix(0,nrow=length(Sigma.varepsilon)*length(ETA),ncol=7)
###################################
###################################
# first generate the log of the average value for good j (see equation (4.3))
#p1j.log=rnorm(q,0,1) # p1j.log is the log of the average value for good j (see equation (4.3))
p1j.log=log(runif(q,0.5,10))
p1j.log=matrix(p1j.log,ncol=1)


for (sigma.varepsilon.i in 1:length(Sigma.varepsilon)) {
  
  sigma.varepsilon=Sigma.varepsilon[sigma.varepsilon.i]
  
  for (eta.i in 1:length(ETA)) {
    
    eta=ETA[eta.i]
    
###################################
###################################
###################################
###################################
    
# compute the true value using 1e6 observations 
n=1e6
#
data=data.sim(q=q,eta=eta,p1j.log=p1j.log,n=n,sigma.varepsilon=sigma.varepsilon,sigma.p=sigma.p)
#
alpha=data$alpha
p1=data$p1
p2=data$p2
q1=data$q1
q2=data$q2
I1=data$I[,1]
I2=data$I[,2]
### 
L1=apply(p2*q1,1,sum)
L2=apply(p1*q1,1,sum)
L3=apply(p2*q2,1,sum)
L4=apply(p1*q2,1,sum)
#
mu1=mean(L1)
mu2=mean(L2)
mu3=mean(L3)
mu4=mean(L4)

########  Laspeyres
tau.t.log=log(mu1)-log(mu2)
tau.t=exp(tau.t.log)
#
Ti=cbind(L1,L2)
Sigma=cov(Ti)
g1 =  1/mu1
g2 = -1/mu2
g=c(g1,g2)
V.Sigma=t(g)%*%Sigma%*%g
sig.t.log=as.vector(sqrt(V.Sigma))
sig.t=exp(tau.t.log)*sig.t.log
#
L.t=tau.t
L.sig.t=sig.t
#####################


########  Paasche
Paasche.i=L3/L4
L5=Paasche.i*L2
mu5=mean(L5)
#
tau.t.log=log(mu5)-log(mu2)
tau.t=exp(tau.t.log)
#
Ti=cbind(L5,L2)
Sigma=cov(Ti)
g1 =  1/mu5
g2 = -1/mu2
g=c(g1,g2)
V.Sigma=t(g)%*%Sigma%*%g
sig.t.log=as.vector(sqrt(V.Sigma))
sig.t=exp(tau.t.log)*sig.t.log
#
P.t=tau.t
P.sig.t=sig.t
#####################


########  Fisher
#
Fisher.i=((L1/L2)*(L3/L4))^(1/2)
L5=Fisher.i*L2
mu5=mean(L5)
#
tau.t.log=log(mu5)-log(mu2)
tau.t=exp(tau.t.log)
#
Ti=cbind(L5,L2)
Sigma=cov(Ti)
g1 =  1/mu5
g2 = -1/mu2
g=c(g1,g2)
V.Sigma=t(g)%*%Sigma%*%g
sig.t.log=as.vector(sqrt(V.Sigma))
sig.t=exp(tau.t.log)*sig.t.log
#
F.t=tau.t
F.sig.t=sig.t
#####################


########  Tornqvist
#
share2=0.5*(p1*q1/apply(p1*q1,1,sum)+p2*q2/apply(p2*q2,1,sum))
Tornqvist.i.log=apply(share2*log(p2/p1),1,sum)
Tornqvist.i=exp(Tornqvist.i.log)
L6=Tornqvist.i*L2
mu6=mean(L6)
#
tau.t.log=log(mu6)-log(mu2)
tau.t=exp(tau.t.log)
#
Ti=cbind(L6,L2)
Sigma=cov(Ti)
g1 =  1/mu6
g2 = -1/mu2
g=c(g1,g2)
V.Sigma=t(g)%*%Sigma%*%g
sig.t.log=as.vector(sqrt(V.Sigma))
sig.t=exp(tau.t.log)*sig.t.log
#
T.t=tau.t
T.sig.t=sig.t
#####################



########  SV
#
share1.ij=p1*q1/apply(p1*q1,1,sum)
share2.ij=p2*q2/apply(p2*q2,1,sum)
share.ij=(share1.ij-share2.ij)/(log(share1.ij)-log(share2.ij))
####################################################################################
share.ij[abs(share1.ij-share2.ij)<1e-10]=share2.ij[abs(share1.ij-share2.ij)<1e-10]
####################################################################################
share3=share.ij/apply(share.ij, 1, sum)
SV.i.log=apply(share3*log(p2/p1),1,sum)
SV.i=exp(SV.i.log)
L7=SV.i*L2
mu7=mean(L7)
#
tau.t.log=log(mu7)-log(mu2)
tau.t=exp(tau.t.log)
#
Ti=cbind(L7,L2)
Sigma=cov(Ti)
g1 =  1/mu7
g2 = -1/mu2
g=c(g1,g2)
V.Sigma=t(g)%*%Sigma%*%g
sig.t.log=as.vector(sqrt(V.Sigma))
sig.t=exp(tau.t.log)*sig.t.log
#
S.t=tau.t
S.sig.t=sig.t
#####################



########  LM
# calculate cost of living index, i.e., CES index
# share is the expenditure weight on jth good for consumer i.
share=p1*q1/apply(p1*q1,1,sum)
LM.i=(apply(share*(p2/p1)^(1-eta),1,sum))^(1/(1-eta))
# if eta==1, the above equation does not work well, we use the following
if (eta==1){
  LM.i=exp(apply(share*log(p2/p1),1,sum))
}
L5=LM.i*L2
mu5=mean(L5)
#
tau.t.log=log(mu5)-log(mu2)
tau.t=exp(tau.t.log)
#
Ti=cbind(L5,L2)
Sigma=cov(Ti)
g1 =  1/mu5
g2 = -1/mu2
g=c(g1,g2)
V.Sigma=t(g)%*%Sigma%*%g
sig.t.log=as.vector(sqrt(V.Sigma))
sig.t=exp(tau.t.log)*sig.t.log
#
LM.t=tau.t
LM.sig.t=sig.t
#####################



##############################################################
### Save the Data For the Difference
#(L.t-LM.t)/LM.sig.t
LPFT=c(L.t,P.t,F.t,T.t,S.t)
(LPFT-LM.t)/LM.sig.t
row.index=(sigma.varepsilon.i-1)*length(ETA)+eta.i
diff.res[row.index,1]=sigma.varepsilon
diff.res[row.index,2]=eta
diff.res[row.index,3:7]=(LPFT-LM.t)/LM.sig.t

outfile=paste("./Data/diff-true.csv",sep="")
write.csv(diff.res,row.names=FALSE,file=outfile)


tex1=formatC(diff.res[,1],width=4,digits = 2,format = "f")
tex2=formatC(diff.res[,2],width=4,digits = 1,format = "f")
tex=paste(tex1,'&',tex2)
for (k in 3:ncol(diff.res)) {
  tex = paste(tex,"&",formatC(diff.res[,k],width=7,digits = 4,format = "f"))
}
#tex=paste(tex,"&",formatC(LM.t,width=6,digits = 4,format = "f"))
tex = paste(tex,"\\\\")
outfile=paste("./Output-All/difff-true-",index,".tex",sep="")
write(tex,file=outfile)
##############################################################


#####################################
#
N=c(10,20,50,100,200,300,500,1000)
M=1000
#N=c(10)
#M=10
#
L.res=matrix(0,nrow=length(N),ncol=3)
P.res=matrix(0,nrow=length(N),ncol=3)
F.res=matrix(0,nrow=length(N),ncol=3)
T.res=matrix(0,nrow=length(N),ncol=3)
S.res=matrix(0,nrow=length(N),ncol=3)
LM.res=matrix(0,nrow=length(N),ncol=3)
sig=matrix(0,nrow=length(N),ncol=12)
tm=matrix(0,nrow=length(N),ncol=12)
#
####################################
t1=proc.time()
for (i in 1:length(N)) {
  
  for (m in 1:M) {
    
    n=N[i]
    #
    data=data.sim(q=q,eta=eta,p1j.log=p1j.log,n=n,sigma.varepsilon=sigma.varepsilon,sigma.p=sigma.p)
    #
    alpha=data$alpha
    p1=data$p1
    p2=data$p2
    q1=data$q1
    q2=data$q2
    I1=data$I[,1]
    I2=data$I[,2]
    #
    L1=apply(p2*q1,1,sum)
    L2=apply(p1*q1,1,sum)
    L3=apply(p2*q2,1,sum)
    L4=apply(p1*q2,1,sum)
    #
    mu1=mean(L1)
    mu2=mean(L2)
    mu3=mean(L3)
    mu4=mean(L4)
    #
    ########  Laspeyres
    tau.s.log=log(mu1)-log(mu2)
    tau.s=exp(tau.s.log)
    #
    Ti=cbind(L1,L2)
    Sigma=cov(Ti)
    g1 =  1/mu1
    g2 = -1/mu2
    #
    g=c(g1,g2)
    V.Sigma=t(g)%*%Sigma%*%g
    sig.s.log=as.vector(sqrt(V.Sigma))
    sig.s=exp(tau.s.log)*sig.s.log
    #
    L.s=tau.s
    L.sig.s=sig.s
    #####################
    
    ########  Paasche
    #
    Paasche.i=L3/L4
    L5=Paasche.i*L2
    mu5=mean(L5)
    #
    tau.s.log=log(mu5)-log(mu2)
    tau.s=exp(tau.s.log)
    #
    Ti=cbind(L5,L2)
    Sigma=cov(Ti)
    g1 =  1/mu5
    g2 = -1/mu2
    #
    g=c(g1,g2)
    V.Sigma=t(g)%*%Sigma%*%g
    sig.s.log=as.vector(sqrt(V.Sigma))
    sig.s=exp(tau.s.log)*sig.s.log
    #
    P.s=tau.s
    P.sig.s=sig.s
    #####################
    
    ########  Fisher
    #
    Fisher.i=((L1/L2)*(L3/L4))^(1/2)
    L5=Fisher.i*L2
    mu5=mean(L5)
    #
    tau.s.log=log(mu5)-log(mu2)
    tau.s=exp(tau.s.log)
    #
    Ti=cbind(L5,L2)
    Sigma=cov(Ti)
    g1 =  1/mu5
    g2 = -1/mu2
    g=c(g1,g2)
    V.Sigma=t(g)%*%Sigma%*%g
    sig.s.log=as.vector(sqrt(V.Sigma))
    sig.s=exp(tau.s.log)*sig.s.log
    #
    F.s=tau.s
    F.sig.s=sig.s
    #####################
    
    ########   Tornqvist
    share2=0.5*(p1*q1/apply(p1*q1,1,sum)+p2*q2/apply(p2*q2,1,sum))
    Tornqvist.i.log=apply(share2*log(p2/p1),1,sum)
    Tornqvist.i=exp(Tornqvist.i.log)
    L6=Tornqvist.i*L2
    mu6=mean(L6)
    #
    tau.s.log=log(mu6)-log(mu2)
    tau.s=exp(tau.s.log)
    #
    Ti=cbind(L6,L2)
    Sigma=cov(Ti)
    g1 =  1/mu6
    g2 = -1/mu2
    g=c(g1,g2)
    V.Sigma=t(g)%*%Sigma%*%g
    sig.s.log=as.vector(sqrt(V.Sigma))
    sig.s=exp(tau.s.log)*sig.s.log
    #
    T.s=tau.s
    T.sig.s=sig.s
    #####################
    
    ########   SV
    share1.ij=p1*q1/apply(p1*q1,1,sum)
    share2.ij=p2*q2/apply(p2*q2,1,sum)
    share.ij=(share1.ij-share2.ij)/(log(share1.ij)-log(share2.ij))
    ####################################################################################
    share.ij[abs(share1.ij-share2.ij)<1e-10]=share2.ij[abs(share1.ij-share2.ij)<1e-10]
    ####################################################################################
    share3=share.ij/apply(share.ij, 1, sum)
    SV.i.log=apply(share3*log(p2/p1),1,sum)
    SV.i=exp(SV.i.log)
    L7=SV.i*L2
    mu7=mean(L7)
    #
    tau.s.log=log(mu7)-log(mu2)
    tau.s=exp(tau.s.log)
    #
    Ti=cbind(L7,L2)
    Sigma=cov(Ti)
    g1 =  1/mu7
    g2 = -1/mu2
    g=c(g1,g2)
    V.Sigma=t(g)%*%Sigma%*%g
    sig.s.log=as.vector(sqrt(V.Sigma))
    sig.s=exp(tau.s.log)*sig.s.log
    #
    S.s=tau.s
    S.sig.s=sig.s
    #####################
    
    
    ########  LM
    # calculate cost of living index, i.e., CES index
    # share is the expenditure weight on jth good for consumer i.
    share=p1*q1/apply(p1*q1,1,sum)
    LM.i=(apply(share*(p2/p1)^(1-eta),1,sum))^(1/(1-eta))
    # if eta==1, the above equation does not work well, we use the following
    if (eta==1){
      LM.i=exp(apply(share*log(p2/p1),1,sum))
    }
    L5=LM.i*L2
    mu5=mean(L5)
    #
    tau.s.log=log(mu5)-log(mu2)
    tau.s=exp(tau.s.log)
    #
    Ti=cbind(L5,L2)
    Sigma=cov(Ti)
    g1 =  1/mu5
    g2 = -1/mu2
    #
    g=c(g1,g2)
    V.Sigma=t(g)%*%Sigma%*%g
    sig.s.log=as.vector(sqrt(V.Sigma))
    sig.s=exp(tau.s.log)*sig.s.log
    #
    LM.s=tau.s
    LM.sig.s=sig.s
    #####################
    
    
    crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
    
    L.ts=L.sig.s/sqrt(n)
    P.ts=P.sig.s/sqrt(n)
    F.ts=F.sig.s/sqrt(n)
    T.ts=T.sig.s/sqrt(n)
    S.ts=S.sig.s/sqrt(n)
    LM.ts=LM.sig.s/sqrt(n)
    
    L.bounds=matrix((L.s-L.ts*crit),nrow=3,ncol=2)
    P.bounds=matrix((P.s-P.ts*crit),nrow=3,ncol=2)
    F.bounds=matrix((F.s-F.ts*crit),nrow=3,ncol=2)
    T.bounds=matrix((T.s-T.ts*crit),nrow=3,ncol=2)
    S.bounds=matrix((S.s-S.ts*crit),nrow=3,ncol=2)
    LM.bounds=matrix((LM.s-LM.ts*crit),nrow=3,ncol=2)
    
    L.res[i,1:3]=ifelse(L.bounds[1:3,1]<L.t &L.t<L.bounds[1:3,2],
                        L.res[i,1:3]+1,L.res[i,1:3])
    P.res[i,1:3]=ifelse(P.bounds[1:3,1]<P.t & P.t<P.bounds[1:3,2],
                        P.res[i,1:3]+1,P.res[i,1:3])
    F.res[i,1:3]=ifelse(F.bounds[1:3,1]<F.t & F.t<F.bounds[1:3,2],
                        F.res[i,1:3]+1,F.res[i,1:3])
    T.res[i,1:3]=ifelse(T.bounds[1:3,1]<T.t & T.t<T.bounds[1:3,2],
                        T.res[i,1:3]+1,T.res[i,1:3])
    S.res[i,1:3]=ifelse(S.bounds[1:3,1]<S.t & S.t<S.bounds[1:3,2],
                        S.res[i,1:3]+1,S.res[i,1:3])
    LM.res[i,1:3]=ifelse(LM.bounds[1:3,1]<LM.t & LM.t<LM.bounds[1:3,2],
                         LM.res[i,1:3]+1,LM.res[i,1:3])
    
    
    sig[i,]=sig[i,]+c(L.sig.s,L.sig.t,P.sig.s,P.sig.t,F.sig.s,F.sig.t,T.sig.s,T.sig.t,S.sig.s,S.sig.t,LM.sig.s,LM.sig.t)
    tm[i,]=tm[i,]+c(L.s,L.t,P.s,P.t,F.s,F.t,T.s,T.t,S.s,S.t,LM.s,LM.t)
    
    if (m %% 100==0){
      cat("Iteration:", m,"/",M,'of',N[i],'Observations',"\n")
      t2=proc.time()
      cat("Elapsed Time in Hour",(t2[[3]]-t1[[3]])/3600,"\n")
    }
  }
}

#
res=cbind(L.res,P.res,F.res,T.res,S.res,LM.res)
#
res/M
sig/M
tm/M
#


### Save the Data For the Coverages
outfile=paste("./Data/coverage-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".csv",sep="")
write.csv(res/M,row.names=FALSE,file=outfile)

### Save the Data For the Mu
outfile=paste("./Data/mu-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".csv",sep="")
write.csv(tm/M,row.names=FALSE,file=outfile)

### Save the Data For the Sigma
outfile=paste("./Data/sigma-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".csv",sep="")
write.csv(sig/M,row.names=FALSE,file=outfile)






### construct the Table ###
res1=res/M
tex=formatC(N,width=4,digits=0,format="f")
tex=paste(tex,"&",sigma.varepsilon,"&",eta)
for (k in 1:ncol(res1)) {
  tex = paste(tex,"&",formatC(res1[,k],width=5,digits = 3,format = "f"))
}
tex = paste(tex,"\\\\")
outfile=paste("./Output/coverage-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".tex",sep="")
write(tex,file=outfile)


res2=sig/M
tex=formatC(N,width=4,digits=0,format="f")
tex=paste(tex,"&",sigma.varepsilon,"&",eta)
for (k in 1:ncol(res2)) {
  tex = paste(tex,"&",formatC(res2[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")
outfile=paste("./Output/sigma-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".tex",sep="")
write(tex,file=outfile)


res3=tm/M
tex=formatC(N,width=4,digits=0,format="f")
tex=paste(tex,"&",sigma.varepsilon,"&",eta)
for (k in 1:ncol(res3)) {
  tex = paste(tex,"&",formatC(res3[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")
outfile=paste("./Output/mu-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".tex",sep="")
write(tex,file=outfile)

###################################
###################################
###################################
###################################

  }
  
}


