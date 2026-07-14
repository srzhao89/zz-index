######################################################################
## Statistical Inference for Price or Quantity Indices
## Author: Valentin Zelenyuk and Shirong Zhao
## Email:  shironz@163.com
## Dalian, July 13, 2026
######################################################################

rm(list=ls())#

Sigma.varepsilon=c(0.02,0.05,0.1,0.2,0.3)

for (sigma.varepsilon.i in 1:length(Sigma.varepsilon)) {

  sigma.varepsilon=Sigma.varepsilon[sigma.varepsilon.i]



setEPS()
postscript(file = paste("./Figures/plot-All-sigma.varepsilon-",sigma.varepsilon,".eps",sep=""))


par(mfrow = c(3, 2),oma = c(2,1,0,1),cex=0.55,mar=c(5,4,3.5,4))


################
# eta=0
################
eta=0
outfile=paste("./Data/coverage-sigma.varepsilon-",sigma.varepsilon,"-All-eta-",eta,".csv",sep="")
df=read.csv(outfile)
sol1=df$V2
sol2=df$V5
sol3=df$V8
sol4=df$V11
sol5=df$V14
sol6=df$V17


x=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000))
plot(x,sol1, type="l",lty=1,lwd=3,col = 1, xaxt='n',ylim=c(0,1),
     xlab="Sample size", ylab = "Coverage", 
     main = bquote(eta==0))
lines(x,sol2, type="l",lty=2,lwd=3,col = 2)
lines(x,sol3, type="l",lty=3,lwd=3,col = 3)
lines(x,sol4, type="l",lty=4,lwd=3,col = 4)
lines(x,sol5, type="l",lty=5,lwd=3,col = 6)
lines(x,sol6, type="l",lty=6,lwd=3,col = 7)
abline(0.95,0,lty=2,lwd=1)
axis(side=1,at=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000)),
     labels=c("10","20","50","100","200","300","500","1000"),
     cex.axis = 0.82)

axis(side=2,at=c(0.2,0.4,0.6,0.8,1.0),
     labels=c("0.2","0.4","0.6","0.8","1.0"))



################
# eta=0.5
################
eta=0.5
outfile=paste("./Data/coverage-sigma.varepsilon-",sigma.varepsilon,"-All-eta-",eta,".csv",sep="")
df=read.csv(outfile)
sol1=df$V2
sol2=df$V5
sol3=df$V8
sol4=df$V11



x=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000))
plot(x,sol1, type="l",lty=1,lwd=3,col = 1, xaxt='n',ylim=c(0,1),
     xlab="Sample size", ylab = "Coverage", 
     main = bquote(eta==0.5))
lines(x,sol2, type="l",lty=2,lwd=3,col = 2)
lines(x,sol3, type="l",lty=3,lwd=3,col = 3)
lines(x,sol4, type="l",lty=4,lwd=3,col = 4)
lines(x,sol5, type="l",lty=5,lwd=3,col = 6)
lines(x,sol6, type="l",lty=6,lwd=3,col = 7)
abline(0.95,0,lty=2,lwd=1)
axis(side=1,at=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000)),
     labels=c("10","20","50","100","200","300","500","1000"),
     cex.axis = 0.82)

axis(side=2,at=c(0.2,0.4,0.6,0.8,1.0),
     labels=c("0.2","0.4","0.6","0.8","1.0"))




################
# eta=1
################
eta=1
outfile=paste("./Data/coverage-sigma.varepsilon-",sigma.varepsilon,"-All-eta-",eta,".csv",sep="")
df=read.csv(outfile)
sol1=df$V2
sol2=df$V5
sol3=df$V8
sol4=df$V11


x=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000))
plot(x,sol1, type="l",lty=1,lwd=3,col = 1, xaxt='n',ylim=c(0,1),
     xlab="Sample size", ylab = "Coverage", 
     main = bquote(eta==1))
lines(x,sol2, type="l",lty=2,lwd=3,col = 2)
lines(x,sol3, type="l",lty=3,lwd=3,col = 3)
lines(x,sol4, type="l",lty=4,lwd=3,col = 4)
lines(x,sol5, type="l",lty=5,lwd=3,col = 6)
lines(x,sol6, type="l",lty=6,lwd=3,col = 7)
abline(0.95,0,lty=2,lwd=1)
axis(side=1,at=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000)),
     labels=c("10","20","50","100","200","300","500","1000"),
     cex.axis = 0.82)

axis(side=2,at=c(0.2,0.4,0.6,0.8,1.0),
     labels=c("0.2","0.4","0.6","0.8","1.0"))





################
# eta=2
################
eta=2
outfile=paste("./Data/coverage-sigma.varepsilon-",sigma.varepsilon,"-All-eta-",eta,".csv",sep="")
df=read.csv(outfile)
sol1=df$V2
sol2=df$V5
sol3=df$V8
sol4=df$V11


x=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000))
plot(x,sol1, type="l",lty=1,lwd=3,col = 1, xaxt='n',ylim=c(0,1),
     xlab="Sample size", ylab = "Coverage", 
     main = bquote(eta==2))
lines(x,sol2, type="l",lty=2,lwd=3,col = 2)
lines(x,sol3, type="l",lty=3,lwd=3,col = 3)
lines(x,sol4, type="l",lty=4,lwd=3,col = 4)
lines(x,sol5, type="l",lty=5,lwd=3,col = 6)
lines(x,sol6, type="l",lty=6,lwd=3,col = 7)
abline(0.95,0,lty=2,lwd=1)
axis(side=1,at=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000)),
     labels=c("10","20","50","100","200","300","500","1000"),
     cex.axis = 0.82)

axis(side=2,at=c(0.2,0.4,0.6,0.8,1.0),
     labels=c("0.2","0.4","0.6","0.8","1.0"))



################
# eta=5
################
eta=5
outfile=paste("./Data/coverage-sigma.varepsilon-",sigma.varepsilon,"-All-eta-",eta,".csv",sep="")
df=read.csv(outfile)
sol1=df$V2
sol2=df$V5
sol3=df$V8
sol4=df$V11


x=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000))
plot(x,sol1, type="l",lty=1,lwd=3,col = 1, xaxt='n',ylim=c(0,1),
     xlab="Sample size", ylab = "Coverage", 
     main = bquote(eta==5))
lines(x,sol2, type="l",lty=2,lwd=3,col = 2)
lines(x,sol3, type="l",lty=3,lwd=3,col = 3)
lines(x,sol4, type="l",lty=4,lwd=3,col = 4)
lines(x,sol5, type="l",lty=5,lwd=3,col = 6)
lines(x,sol6, type="l",lty=6,lwd=3,col = 7)
abline(0.95,0,lty=2,lwd=1)
axis(side=1,at=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000)),
     labels=c("10","20","50","100","200","300","500","1000"),
     cex.axis = 0.82)

axis(side=2,at=c(0.2,0.4,0.6,0.8,1.0),
     labels=c("0.2","0.4","0.6","0.8","1.0"))




################
# eta=10
################
eta=10
outfile=paste("./Data/coverage-sigma.varepsilon-",sigma.varepsilon,"-All-eta-",eta,".csv",sep="")
df=read.csv(outfile)
sol1=df$V2
sol2=df$V5
sol3=df$V8
sol4=df$V11


x=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000))
plot(x,sol1, type="l",lty=1,lwd=3,col = 1, xaxt='n',ylim=c(0,1),
     xlab="Sample size", ylab = "Coverage", 
     main = bquote(eta==10))
lines(x,sol2, type="l",lty=2,lwd=3,col = 2)
lines(x,sol3, type="l",lty=3,lwd=3,col = 3)
lines(x,sol4, type="l",lty=4,lwd=3,col = 4)
lines(x,sol5, type="l",lty=5,lwd=3,col = 6)
lines(x,sol6, type="l",lty=6,lwd=3,col = 7)
abline(0.95,0,lty=2,lwd=1)
axis(side=1,at=c(log10(10),log10(20),log10(50),log10(100),log10(200),log10(300),log10(500),log10(1000)),
     labels=c("10","20","50","100","200","300","500","1000"),
     cex.axis = 0.82)

axis(side=2,at=c(0.2,0.4,0.6,0.8,1.0),
     labels=c("0.2","0.4","0.6","0.8","1.0"))



######################################################################
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Laspeyres","Paasche","Fisher",expression("Törnqvist"),"SV","LM","0.95"), col = c(1,2,3,4,6,7,1), 
       lty=c(1,2,3,4,5,6,2),lwd =c(3,3,3,3,3,3,1),xpd = TRUE,cex=1.4, ncol=7,bty = 'n')



dev.off()

}
#
