#Bruce Stagg Martin, 01/18/2022
#Demonstration of 'contsimmap' Package Prototype
rm(list=ls())
set.seed(123)
library(phytools)
library(contsimmap)
library(Rphylopars)

####SIMULATING CONTINUOUS TRAIT DATA####

ntaxa<-100
tree<-pbtree(n=ntaxa,scale=1)
#Xsig2: evolutionary rate matrix
Xsig2<-diag(c(1,3))%*%matrix(c(1,-0.6,-0.6,1),2,2)%*%diag(c(1,3))
#Ysig2: intra-tip variance-covariance matrix
Ysig2<-diag(c(1,1))%*%matrix(c(1,0,0,1),2,2)%*%diag(c(1,1))
#nobs: number of observations per tip
nobs<-pmax(1,rpois(ntaxa,3))
#simulate history of a continuous character
##res is the approximate number of time points to sample trait values along the tree's height
###in this case it's 1, meaning that traits will only be sampled at nodes
##nsims is just the number of simulations to do
sim<-sim.conthistory(tree,ntraits=2,traits=c('x1','x2'),
                     Xsig2=Xsig2,Ysig2=Ysig2,nobs=nobs,
                     res=1,nsims=1)
trait.data<-sim$trait.data[,,1]

####FITTING TRAIT DATA TO MODEL####

#eventually I want to make my own functions for this and potentially do it "under the hood", so to speak
##as in Liam's simmaping functions
##for now, however, phylopars will be sufficient for this example
fit<-phylopars(trait_data=data.frame(species=rownames(trait.data),trait.data),
               tree=tree)
fit.Xsig2<-fit$pars$phylocov
fit.Ysig2<-fit$pars$phenocov

####USING MODEL TO SIMULATE CONTINUOUS STOCHASTIC CHARACTER MAPS (CONTSIMMAPS)####

#pretty simple--only thing to make sure is that trait.data is properly formatted
##trait.data should be a numeric matrix with rownames corresponding to tips and colnames to traits
map<-make.contsimmap(trait.data,tree,
                     Xsig2=fit.Xsig2,Ysig2=fit.Ysig2,
                     res=100,nsims=100)

####VIEWING CONTSIMMAP INFORMATION/PLOTTING####

##PLOTTING
#I made a rudimentary plotting method that could use some improvement in terms of flexibility...
plot(map)
#uses matplot() internally, so different colors/lines correspond to different simulations
plot(map,col=rgb(0,0,0,0.1),lty=1)
#use traits to select different traits
plot(map,col=rgb(0,0,0,0.1),lty=1,traits=2)
plot(map,col=rgb(0,0,0,0.1),lty=1,traits=c(1,2))
#by default, it plots 20 randomly sampled simulations, but we can tell it to sample a different number
plot(map,sims=2)
#or plot specific simulations
plot(map,sims=c(1,2))
#(if you want to select a single, specific simulation, pass the index as a list)
plot(map,sims=list(1))
#you can also select edges
plot(map,edges=c(1,10))

##MORE FINE MANIPULATION/VIEWING
#the plotting function is based on some custom generic methods, namely subset() and format()
#you can use subset to select certain traits, simulations, or edges; compare:
summary(map)
summary(subset(map,sims=list(1),trait=2,edges=c(1,10)))
#the internal structure of contsimmap objects is actually quite complex
##the reasons for this are numerous--preventing redundancy, ease of manipulation, accommodating trees with multiple regimes...
##and also just me generally still figuring out the best way to do all this
##anyways, because of this, I made a format function which makes it easier to pull out trait values and make it more readable
format(map,sims=list(1),trait=2,edges=c(1,10))
#format basically returns lists of 3D arrays of trait values, with each list element corresponding to an edge
##the arrays themselves are arranged with time points along the first dimension, traits second, and simulations third
#you can pull out numeric vectors of the sampled time points along each edge too, if needed
map$ts

####MAKING NEW TRAITS FROM CONTSIMMAPS####

#one of the primary things I'm interested in for later applications is transforming traits in the contsimmap
##as I said, the internal structure of these objects can get kinda complex, so I made a function to do this conveniently
map<-make.traits(map,new_x1~exp(x1),new_x2~plogis(x2),x3~x1+x2)
plot(map,col=rgb(0,0,0,0.1),lty=1,traits=3)
plot(map,col=rgb(0,0,0,0.1),lty=1,traits=4)
#the only warning here is that the formula for the transformation HAS to be "1-to-1"
##if the formula includes anything like mean(), sum(), etc., things will break, badly!
#in the future, I definitely hope to make functions that could summarize across simulations using mean() or quantile()
##but I haven't gotten there yet

####THRESHOLDING AND CONVERTING TO SIMMAPS####

#this brings me to thresholding continuous trait data and converting it to a normal simmap!
#for this, I made a threshold() function which can be plugged into a make.traits() formula
map<-make.traits(map,thresh_x1~threshold(x1,nbreaks=7))
plot(map,traits='thresh_x1')
#info about break points and state names are stored here
map$keys
#nbreaks breaks up the trait into nbreaks+1 equally-spaced bins
##more explicit breaks can be specified with breaks=, as well as different names than the default using state.names=
##(numeric.names=TRUE will make the default names be 0,1,2,... rather than A,B,C,...)
map<-make.traits(map,thresh_x2~threshold(x2,breaks=c(-6,2),state.names=c('lo','med','hi')))
plot(map,traits='thresh_x2')
#the function as.simmap performs thresholding for you, then converts that trait's history to a multiSimmap object based on the thresholds
simmap<-as.simmap(map,trait='x2',breaks=c(-6,2),state.names=c('lo','med','hi'))
plot(simmap[[1]])
#you can still access the break points from the simmaps
simmap[[1]]$breaks

####SOME NOTES####

#while there might be exact ways to simulate the times a BM hits a threshold, it would be complicated to do so
##it also creates conceptual issues with regard to resolution
##in truly continuous time, BM technically crosses a threshold infinitely many times in a single instant (or so I've been told)
##so I simply approximate BM between sampled time points with straight lines and then calculate when these lines cross thresholds
##because of this, simmaps based on contsimmaps might appear very different depending on their resolution

#I made the sim.conthistory and make.contsimmap functions able to handle simmaps/multiSimmaps as well
##in this case, you can pass a labeled list to Xsig2 and simulate trait data under a multi-regime BM model
##you can also pass a labeled list to Ysig2 to allow different tips to have different intra-tip var-cov matrices
##let me know if you're interested in these features and I can demo/explain more
