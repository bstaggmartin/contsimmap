# plot(as.vector(sim_XX[1,,])~rep(ape::node.depth.edgelength(tree),nsim))
# plot.fun<-function(sim_XX,edge,trait=1,nsample=20,...){
#   xx<-edge
#   xx[,]<-ape::node.depth.edgelength(tree)[edge]
#   yy<-matrix(sim_XX[trait,t(edge),sample(dim(sim_XX)[3],nsample)],ncol=2,byrow=TRUE)
#   plot(0,col='white',type='p',xlim=range(xx),ylim=range(yy))
#   segments(x0=xx[,1],x1=xx[,2],y0=yy[,1],y1=yy[,2],...)
# }
# plot.fun(sim_XX,trait=1,edge,col=rgb(0,0,0,0.3))
# plot(t(sim_XX[,5,]))

#need to put some more effort in in argument specification
#abilities to color/change line type by edge, by evolving quantity, and by state all needed
#for now all you get is ability to color/change line type by simulation
#also something

#pass samp as length 1 list if you wish to plot a specific simulation rather than sample simulations randomly

#' @export
plot.contsimmap<-function(contsimmap,traits=1,sims=min(20,NROW(contsimmap$perm.inds)),edges=NULL,
                   xlim=NULL,xlab=NULL,ylim=NULL,ylab=NULL,...){
  #formatting
  if(length(traits)>1){
    traits<-traits[c(1,2)]
  }
  contsimmap<-subset(contsimmap,edges,sims,traits)
  contsimmap[['x']]<-format(contsimmap)
  dims<-dim(contsimmap$x[[1]])
  dimnms<-dimnames(contsimmap$x[[1]])
  nsims<-dims[3]
  k<-dims[2]
  trait.names<-dimnms[[2]]
  
  #plotting parameters
  if(is.null(xlim)&k==1) xlim<-range(unlist(contsimmap$ts))
  if(is.null(xlim)) xlim<-range(unlist(lapply(contsimmap$x,function(ii) ii[,1,])))
  if(is.null(xlab)&k==1) xlab<-'time'
  if(is.null(xlab)){
    xlab<-trait.names[1]
  }
  if(is.null(ylim)&k==1) ylim<-range(unlist(contsimmap$x))
  if(is.null(ylim)) ylim<-range(unlist(lapply(contsimmap$x,function(ii) ii[,2,])))
  if(is.null(ylab)&k==1){
    ylab<-trait.names[1]
  }
  if(is.null(ylab)&k==2){
    ylab<-trait.names[2]
  }
  
  #put it all together
  ne<-length(contsimmap$x)
  if(k==1){
    for(e in seq_len(ne)){
      matplot(x=contsimmap$ts[[e]],y=contsimmap$x[[e]][,1,],
              xlim=xlim,xlab=xlab,ylim=ylim,ylab=ylab,
              type='l',add=if(e==1) FALSE else TRUE,...)
    }
  }else{
    for(e in seq_len(ne)){
      matplot(x=contsimmap$x[[e]][,1,],y=contsimmap$x[[e]][,2,],
              xlim=xlim,xlab=xlab,ylim=ylim,ylab=ylab,
              type='l',add=if(e==1) FALSE else TRUE,...)
    }
  }
}
