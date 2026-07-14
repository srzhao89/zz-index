######################################################################
## Statistical Inference for Price or Quantity Indices
## Author: Valentin Zelenyuk and Shirong Zhao
## Email:  shironz@163.com
## Dalian, July 13, 2026
######################################################################
rm(list=ls())#
require(fdrtool)
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
index="All"
####################################
## Combining Coverages
ETA=c(0,0.5,1,2,5,10)
Sigma.varepsilon=c(0.02,0.05,0.1,0.2,0.3)
tex=""
for (sigma.varepsilon.i in 1:length(Sigma.varepsilon)) {
  sigma.varepsilon=Sigma.varepsilon[sigma.varepsilon.i]
  for (eta.i in 1:length(ETA)) {
    eta=ETA[eta.i]
    
    outfile=paste("./Output/coverage-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".tex",sep="")
    df=read.delim(outfile, header = FALSE)  
    
    tex=rbind(tex,"\\\\", df)
  }
  
}


tex=unlist(tex)
print(tex)

outfile=paste("./Output-All/coverage-",index,".tex",sep="")
write(tex,file=outfile)


##### Mean Values
tex=""
for (sigma.varepsilon.i in 1:length(Sigma.varepsilon)) {
  sigma.varepsilon=Sigma.varepsilon[sigma.varepsilon.i]
  for (eta.i in 1:length(ETA)) {
    eta=ETA[eta.i]
    
    outfile=paste("./Output/mu-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".tex",sep="")
    df=read.delim(outfile, header = FALSE)  
    
    tex=rbind(tex,"\\\\", df)
  }
  
}


tex=unlist(tex)
print(tex)

outfile=paste("./Output-All/mu-",index,".tex",sep="")
write(tex,file=outfile)



##### Standard Deviations

tex=""
for (sigma.varepsilon.i in 1:length(Sigma.varepsilon)) {
  sigma.varepsilon=Sigma.varepsilon[sigma.varepsilon.i]
  for (eta.i in 1:length(ETA)) {
    eta=ETA[eta.i]
    
    outfile=paste("./Output/sigma-sigma.varepsilon-",sigma.varepsilon,"-",index,"-eta-",eta,".tex",sep="")
    df=read.delim(outfile, header = FALSE)  
    
    tex=rbind(tex,"\\\\", df)
  }
  
}


tex=unlist(tex)
print(tex)

outfile=paste("./Output-All/sigma-",index,".tex",sep="")
write(tex,file=outfile)








