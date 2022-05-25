#make test for getting nums/denoms
vec<-seq(0,1,0.1)[-1]
vec<-sort(c(vec,setNames(runif(2),c('0','1'))))
names(vec)[length(vec)]<-'2'
incl<-nzchar(names(vec))
tmp<-rle(names(vec))
inds<-which(!nzchar(tmp$values))
tmp$values[inds]<-tmp$values[inds+1]
names(vec)<-inverse.rle(tmp)
dt<-c(vec[1],diff(vec))
tmp1<-c(rev(cumsum(rev(dt[-1]))),0)
tmp2<-rle(names(dt))
sub<-tmp1[incl]
tmp2$values<-sub
num<-tmp1-inverse.rle(tmp2)
denom<-num+dt
