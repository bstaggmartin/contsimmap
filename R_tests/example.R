tree<-phytools::pbtree(n=1000)
trait.data<-evorates::sim.evorates(tree,Rsig2=0,Xsig2=matrix(c(1,-0.8,0.5,-0.8,1,-0.3,0.5,-0.3,1),3,3),n.obs=rpois(10,5),Ysig2=-0.2^2*matrix(c(-2,-0.8,0.5,-0.8,-2,-0.3,0.5,-0.3,-2),3,3))$trait.data
trait.data[sample(length(trait.data),3000)]<-NA
Xsig2<-matrix(c(1,-0.8,0.5,-0.8,1,-0.3,0.5,-0.3,1),3,3)
Ysig2<- -0.2^2*matrix(c(-2,-0.8,0.5,-0.8,-2,-0.3,0.5,-0.3,-2),3,3)

test<-make.contsimmap(trait.data,tree,Xsig2=Xsig2,Ysig2=Ysig2,nsim=1,res=100)
plot(test,samp=list(1),lty=1)
hist(test$x[[1]][1,1,],breaks=50)

tree<-phytools::pbtree(n=100)
trait.data<-evorates::sim.evorates(tree,Rsig2=0,Xsig2=1,n.obs=rpois(10,5),Ysig2=1)$trait.data
trait.data[sample(length(trait.data),100)]<-NA
test<-make.contsimmap(trait.data,tree,nsim=10,res=100)
