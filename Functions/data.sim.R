### data generation
data.sim <- function(q=2,eta=eta,p1j.log=p1j.log,n=100,sigma.varepsilon=sigma.varepsilon,sigma.p=sigma.p) {
  
  
  #alpha1=matrix(runif(n*(q-1),1/(1.5*q),1/q),nrow=n,ncol=q-1)
  #alpha2=1-apply(alpha1,1,sum)
  #alpha=cbind(alpha1,alpha2)
  
  alpha1=matrix(rexp(n*q),nrow=n,ncol=q)
  alpha=alpha1/apply(alpha1, 1, sum)
  
  p1.log=apply(p1j.log,1,function(x) rnorm(n,x,sigma.p))
  p1=exp(p1.log)
  p2=p1*exp(0.05+rnorm(n*q,0,sigma.varepsilon))
  #epilson.j=rnorm(q,0,0.01)
  #epilson=matrix(rep(epilson.j,n),byrow=TRUE,nrow=n,ncol=q)
  #p2=p1*exp(0.05+epilson)
  
  I1=runif(n,10,100)
  I2=I1*exp(rnorm(n,0,0.02))
  
  I=cbind(I1,I2)
  
  q1=I1*alpha*p1^{-eta}/apply(alpha*p1^{1-eta},1,sum)
  q2=I2*alpha*p2^{-eta}/apply(alpha*p2^{1-eta},1,sum)
  
  #
  res=list(p1=p1,p2=p2,q1=q1,q2=q2,I1=I1,I2=I2,alpha=alpha)
  
  return(res)
  
}