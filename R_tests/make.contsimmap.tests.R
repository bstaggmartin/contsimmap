rm(list=ls())
set.seed(123)
library(phytools)
library(contsimmap)
tree<-pbtree(n=10,scale=1)
tree$edge.length[8:9]<-0
Q<-matrix(c(-1,2,1,-2),2,2)
colnames(Q)<-rownames(Q)<-c(0,1)
tree<-sim.history(tree,Q,nsim=1)
Xsig2<-asplit(rWishart(2,4,diag(3)),3)
Xsig2[[1]][2,]<-0
Xsig2[[1]][,2]<-0
cont.tree<-sim.conthistory(tree,ntraits=3,nsims=1,res=1,Xsig2=Xsig2[[1]],Ysig2=0.2^2,nobs=rpois(100,5))
trait.data<-matrix(cont.tree$trait.data[,,1],dim(cont.tree$trait.data)[1],dim(cont.tree$trait.data)[2],
                   dimnames=dimnames(cont.tree$trait.data)[-3])
trait.data[sample(length(trait.data),30)]<-NA
Xsig2<-cont.tree$params$Xsig2
Ysig2<-cont.tree$params$Ysig2
trees<-make.simmap(tree,
                   setNames(unlist(lapply(tree$maps,function(ii) names(ii)[length(ii)]))[1:100],tree$tip.label),
                   Q=Q,nsim=50)
tree<-c(tree,trees)
# undebug(make.contsimmap)
test<-make.contsimmap(trait.data,tree,Xsig2=Xsig2,Ysig2=Ysig2,res=1000)
plot(test,sims=20,traits=2,col=rgb(0,0,0,0.1),lty=1)
#need to test where/if make.contsimmap breaks in the cases of...
  #single traits -- check
  #multiple sims per tree -- check
  #single states -- check
  #single trees -- check
  #branches of length 0 (tips in particular) -- not extensively tested, but seems to work

rm(list=ls())
set.seed(123)
library(phytools)
library(contsimmap)
tree<-pbtree(n=100,scale=1)
Q<-matrix(c(-1,2,1,-2),2,2)
colnames(Q)<-rownames(Q)<-c(0,1)
tree<-sim.history(tree,Q,nsim=1)
Xsig2<-asplit(rWishart(2,4,diag(3)),3)
cont.tree<-sim.conthistory(tree,ntraits=3,nsims=1,res=1,Xsig2=Xsig2[[1]],Ysig2=0.2^2,nobs=rpois(100,5))
trait.data<-matrix(cont.tree$trait.data[,,1],dim(cont.tree$trait.data)[1],dim(cont.tree$trait.data)[2],
                   dimnames=dimnames(cont.tree$trait.data)[-3])
trait.data[sample(length(trait.data),300)]<-NA
Xsig2<-cont.tree$params$Xsig2
Ysig2<-cont.tree$params$Ysig2
trees<-make.simmap(tree,
                   setNames(unlist(lapply(tree$maps,function(ii) names(ii)[length(ii)]))[1:100],tree$tip.label),
                   Q=Q,nsim=50)
tree<-c(tree,trees)
test<-make.contsimmap(trait.data,tree,Xsig2=Xsig2,Ysig2=Ysig2,res=1,nsims=5000)
plot(test,col=rgb(0,0,0,0.1),lty=1,traits=2)

library(Rphylopars)
dat<-as.data.frame(trait.data)
dat$species<-rownames(trait.data)
dat<-dat[,c(4,1:3)]
tmp<-Rphylopars::phylopars(dat,tree[[1]],phylocov_fixed=Xsig2[[1]],phenocov_fixed=Ysig2[[1]])
cov(t(test$nodes$'0'))
tmp$anc_cov[[101]]
#it all looks about right...
