#Generate Brownian Motion bridge variance-covariance matrix--assumes unit scale (time scale/rate of 1)
#Can be multiplied to work for alternate time scales/rates
.bridge.cov<-function(n){
  pts<-seq(0,1,length.out=n+2)[-c(1,n+2)]
  out<-matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:i){
      out[i,j]<-(1-pts[i])*pts[j]
    }
  }
  tmp<-upper.tri(out)
  out[tmp]<-t(out)[tmp]
  out
}

#Linearly interpolate n equally spaced points between yy[,1,] and yy[,2,]
#Note that it returns the results as a vector
.lin.interp<-function(yy,n){
  y1<-as.vector(yy[,1,])
  y2<-as.vector(yy[,2,])
  xx<-rep(seq(0,1,length.out=n+2)[-c(1,n+2)],each=length(y1))
  c(y1,(y2-y1)*xx+y1,y2)
}