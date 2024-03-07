

#impact curve function

impactcurve<-function(tval,mu=0.5,beta=0.1,M=1,Y=3.5,B=1/(1+exp(mu/beta)),C=(1+exp(-(1-mu)/beta))/(1-B*(1+exp(-(1-mu)/beta)))){
  ifelse(tval>=Y, M, M*C*(1/(1+exp(-(tval/Y-mu)/beta))-B))
}

#function to apply impact curve function to a raster

threatToCond<-function(inval,type=2,mu.val=0.5,beta.val=0.1,M.val=1,Y.val=3.5){
  if(type==1){mu.val=0;beta.val=0.1;M.val=1;Y.val=3.5}
  if(type==3){mu.val=1;beta.val=1;M.val=1;Y.val=3.5}
  if(type==4){mu.val=1;beta.val=0.1;M.val=1;Y.val=3.5}
  impactcurve(tval=inval,mu=mu.val,beta=beta.val,M=M.val,Y=Y.val)
}



#absolute of negative stress values only

absStress<-function(stress,posneg="neg"){   #posneg='neg' to make positive stress zero, posneg='pos' to make negative stress zero
  stressout=abs(stress)
  if(posneg=="neg"){stressout=stressout*(stress<0)}
  if(posneg=="pos"){stressout=stressout*(stress>0)}
  stressout
}

