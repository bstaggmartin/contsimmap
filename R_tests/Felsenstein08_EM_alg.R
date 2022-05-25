library(phytools)
library(contsimmap)
set.seed(123)
tree<-pbtree(n=100)
Xsig2<-rWishart(1,3,diag(3))[,,1]
Ysig2<-rWishart(1,3,diag(3))[,,1]
X<-sim.conthistory(tree,
                   ntraits=3,Xsig2=Xsig2,
                   nobs=pmax(rpois(100,2),1),Ysig2=Ysig2,
                   res=1,nsims=1)$trait.data[,,1]
plot(X,col=as.factor(rownames(X)))
ntraits<-ncol(X)
CC<-XX<-matrix(ncol=ntraits,nrow=max(tree$edge))
WW<-SS<-VV<-vector('numeric',max(tree$edge))
ntips<-length(tree$tip.label)
#things were getting reordered!
XX[1:ntips,]<-apply(X,2,function(ii) tapply(ii,rownames(X),mean))[tree$tip.label,]
SS[1:ntips]<-1/lengths(split(rownames(X),rownames(X)))[tree$tip.label]
des<-c(list(evorates:::root.edges.phylo(tree)),
       evorates:::des.edges(tree))
#assumes cladewise ordering
prune.seq<-rev(which(lengths(des)>0))
des<-des[prune.seq]
edge<-rbind(c(0,ntips+1),tree$edge)
elen<-c(0,tree$edge.length)
#assumes bifurcating tree...
for(i in seq_along(prune.seq)){
  e<-prune.seq[i]
  des_e<-des[[i]]+1
  n<-edge[e,2]
  des_n<-edge[des_e,2]
  t<-elen[des_e]
  tmp<-t+VV[des_n]
  sum.SS<-sum(SS[des_n])
  CC[n,]<-(XX[des_n[2],]-XX[des_n[1],])/sqrt(sum.SS)
  WW[n]<-sum(tmp)/sum.SS
  wgts<-1/tmp
  sum.wgts<-sum(wgts)
  VV[n]<-1/sum.wgts
  SS[n]<-sum((wgts/sum.wgts)^2*SS[des_n])
  XX[n,]<-VV[n]*.colSums(wgts*XX[des_n,],2,ntraits)
}
foo<-function(x){
  len<-length(x)
  if(len>1){
    inds<-seq_len(len)
    cummeans<-cumsum(x[-len])/inds[-len]
    out<-x[-1]-cummeans
    sqrts<-sqrt(inds)
    sqrts[-len]/sqrts[-1]*out
  }
}
within.CC<-apply(X,2,function(ii) unlist(tapply(ii,rownames(X),foo),use.names=FALSE))
if(is.null(within.CC)){
  within.CC<-matrix(nrow=0,ncol=ntraits)
}

node.seq<-edge[prune.seq,2]
W<-WW[node.seq]
sq.W<-sqrt(W)
C<-lapply(node.seq,function(ii) tcrossprod(CC[ii,]))
A<-diag(3)
P<-diag(3)
tot.len<-length(node.seq)+nrow(within.CC)
AA<-vector('list',tot.len)
PP<-vector('list',tot.len)
tmp.seq<-seq_along(node.seq)
PP[-tmp.seq]<-lapply(seq_len(nrow(within.CC)),function(ii) tcrossprod(within.CC[ii,]))
denom<-lapply(W,function(ii) solve(ii*A+P))
BA<-lapply(tmp.seq,function(ii) sq.W[ii]*A%*%denom[[ii]])
BP<-lapply(denom,function(ii) P%*%ii)
AA[tmp.seq]<-lapply(tmp.seq,function(ii) A+BA[[ii]]%*%(C[[ii]]-(W[ii]*A+P))%*%t(BA[[ii]]))
AA[-tmp.seq]<-list(A)
PP[tmp.seq]<-lapply(tmp.seq,function(ii) P+BP[[ii]]%*%(C[[ii]]-(W[ii]*P+P))%*%t(BP[[ii]]))
prop.A<-Reduce('+',AA)/tot.len
prop.P<-Reduce('+',PP)/tot.len
counter<-1
while(any(abs(A-prop.A)>1e-6)){
  A<-prop.A
  P<-prop.P
  denom<-lapply(W,function(ii) solve(ii*A+P))
  BA<-lapply(tmp.seq,function(ii) sq.W[ii]*A%*%denom[[ii]])
  BP<-lapply(denom,function(ii) P%*%ii)
  AA[tmp.seq]<-lapply(tmp.seq,function(ii) A+BA[[ii]]%*%(C[[ii]]-(W[ii]*A+P))%*%t(BA[[ii]]))
  AA[-tmp.seq]<-list(A)
  PP[tmp.seq]<-lapply(tmp.seq,function(ii) P+BP[[ii]]%*%(C[[ii]]-(W[ii]*P+P))%*%t(BP[[ii]]))
  prop.A<-Reduce('+',AA)/tot.len
  prop.P<-Reduce('+',PP)/tot.len
  counter<-counter+1
}

A;P;counter
#A seems off...it keeps blowing to ridiculous values for some reason...
#I'm probably making some mistake in the contrast algorithm itself...
#FOund a mistake that improved things, but I nonetheless still get A blowing up for some reason...
#A does converge, but to values much higher than the generating ones--hopefully it's not something wrong with simulation!
#Figured it out! Tip labels were just getting shuffled around in some of the beginning steps!
#Can be very slow to converge in some situations, but does work