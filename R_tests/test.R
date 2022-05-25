library(phytools)
library(contsimmap)
tree<-pbtree(n=100,scale=1)
Q<-matrix(c(-1,2,1,-2),2,2)
colnames(Q)<-rownames(Q)<-c(0,1)
tree<-sim.history(tree,Q,nsim=1)

sims<-sim.conthistory(tree,ntraits=3,res=1000,
                      Xsig2=list(diag(3),matrix(c(3,-2,-1,
                                                  -2,3,0.5,
                                                  -1,0.5,1),
                                                3,
                                                3)
                      )
) #~2.3 seconds for 100 simmaps with res 100, a little over 3 with res 500
#even 1000 res is pretty fast!
plot(sims,samp=1,trait=c(1,2))
sims<-sim.conthistory(tree,ntraits=2,res=1000,
                      Xsig2=list(matrix(c(2,1,1,2),2,2),
                                 matrix(c(2,-1,-1,2),2,2)),
                      nobs=10,Ysig2=0.1^2,
)
plot(sims,samp=1,trait=c(1,2)) #perfect-looking!
plot(sims$trait.data[,,1],col=rainbow(100)[as.factor(rownames(sims$trait.data))])
plot(sims,samp=1,trait=1)

debug(sim.conthistory)
sims<-sim.conthistory(tree,ntraits=5,res=1000,nobs=0)
debug(contsimmap:::.get.edges)
test<-subset(sims,edges=c(1,100),sims=1)
plot(test,sims=100)
plot(sims,sims=list(100))
#check to see if new method of edge selection/ancestor finding works...
#edge/trait selection seems to work, but SIMS selection is exhibiting some issues
#fixed sims, but now edges is broken again...

#subsetting now working as far as I can tell, but function seems to break with singular inputs...