rm(list=ls())
set.seed(123)
library(phytools)
library(contsimmap)
tree<-pbtree(n=200,scale=1)
Q<-matrix(c(-1,2,1,-2),2,2)
colnames(Q)<-rownames(Q)<-c(0,1)
tree<-sim.history(tree,Q,nsim=1)
cont.tree<-sim.conthistory(tree,ntraits=3,nsims=1,res=1,Xsig2=asplit(rWishart(2,4,diag(3)),3),Ysig2=0.2^2,nobs=rpois(100,5))



trait.data<-matrix(cont.tree$trait.data[,,1],dim(cont.tree$trait.data)[1],dim(cont.tree$trait.data)[2],
                   dimnames=dimnames(cont.tree$trait.data)[-3])
Xsig2<-cont.tree$params$Xsig2
Ysig2<-cont.tree$params$Ysig2
trees<-make.simmap(tree,
                   tree$states,
                   Q=Q,
                   nsim=100)
trait.data[sample(length(trait.data),1000)]<-NA
contsimmap<-make.contsimmap(trait.data,trees,Xsig2=Xsig2,Ysig2=Ysig2,res=200)






plot(contsimmap,col=rgb(0,0,0,0.1),lty=1)
contsimmap<-make.traits(contsimmap,trait_4~threshold(trait_1,nbreaks=40),trait_5~exp(trait_2))
plot(contsimmap,traits=c(4,5)) #doin' it's thing!
simmap<-as.simmap(contsimmap,nbreaks=39,numeric.names=TRUE)
plot(simmap)
format(contsimmap)