rm(list=ls())
set.seed(123)
library(phytools)
library(contsimmap)
tree<-pbtree(n=50,scale=1)
Q<-matrix(c(-1,2,1,-2),2,2)
colnames(Q)<-rownames(Q)<-c(0,1)
tree<-sim.history(tree,Q,nsim=1)
cont.tree<-sim.conthistory(tree,ntraits=3,nsims=1,res=1,Xsig2=asplit(rWishart(2,4,diag(3)),3),Ysig2=0.2^2,nobs=rpois(100,5))
a<-2
cont.tree<-make.traits(cont.tree,r~exp(trait_1-a),c(y,z)~diffusion(y='r',z=1,Xsig2=matrix(c(1,0.9,0.9,1),2,2)))
cont.tree
plot(cont.tree,traits=c('y','z'),Col.by='r',col=hcl.colors(100),sims=1)

trait.data<-matrix(cont.tree$trait.data[,,1],dim(cont.tree$trait.data)[1],dim(cont.tree$trait.data)[2],
                   dimnames=dimnames(cont.tree$trait.data)[-3])
Xsig2<-cont.tree$params$Xsig2
Ysig2<-cont.tree$params$Ysig2
trees<-make.simmap(tree,
                   tree$states,
                   Q=Q,
                   nsim=100)
trait.data[sample(length(trait.data),250)]<-NA
contsimmap2<-make.contsimmap(trait.data,trees,Xsig2=Xsig2,Ysig2=Ysig2,res=100)
test<-get.node.vals(contsimmap2)

X<-fastBM(tree)
sims<-make.contsimmap(as.matrix(X),tree,Xsig2=mean(pic(X,tree)^2),nsims=1000)
test<-get.node.vals(sims)
plot(c(rep(0,50),fastAnc(tree,X,vars=TRUE)$var),apply(test[,1,],1,var))
abline(0,1) #does work

plot(contsimmap,col=rgb(0,0,0,0.1),lty=1)
contsimmap<-make.traits(contsimmap,trait_4~threshold(trait_1,nbreaks=40),trait_5~exp(trait_2))
plot(contsimmap,traits=c(4,5)) #doin' it's thing!
simmap<-as.simmap(contsimmap,nbreaks=39,numeric.names=TRUE)
plot(simmap)
format(contsimmap)