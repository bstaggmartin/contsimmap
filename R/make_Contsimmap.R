#2/13/23 Parameter formatting:
# - Generally expects list arrays with each column corresponding to a state (mu, Xsig2), tip (Ysig2, nobs), or trait (X0)
#    - Also accepts matrices in the case of X0/nobs
#    - Also accepts lists of lists, where each sublist represents a column
# - Unnamed elements readily recycled to fill missing portions of inputs, but if no unnamed elements are found, will be
#   filled in with defaults (identity for Xsig2, 0s for mu/Ysig2/nobs/X0)
#    - Entries of nobs corresponding to internal nodes are never derived from recycling and always default to 0 if unspecified

#' @export
sim.conthistory<-function(tree,ntraits=1,traits=paste0('trait_',seq_len(ntraits)),nobs=NULL,
                          nsims=100,res=100,
                          X0=NULL,Xsig2=NULL,Ysig2=NULL,mu=NULL,
                          verbose=FALSE){
  #takes care of all the necessary formatting work...
  list2env(.prep.contsimmap(tree,nsims,res,X0,Xsig2,Ysig2,mu,conditional=FALSE,ntraits,traits,nobs),envir=environment())
  
  #at some point want to generalize to have multiple nobs per call to simulate differently-structured trait.data...
  #a little confusing at this point since the trait.data lookup element refers to X0 here, but trait.data for make.contsimmap
  #but that can be resolved manually for now
  #^2/14/23: need to look up if the above still holds
  #Can now handle multiple nobs, but things are still a bit confusing to me with how lookups are handled
  #Should eventually become more concrete/organized as code matures
  tmp<-.uncond.traversals(prune.seq,anc,tree[[1]][['edge']],Ntip(tree[[1]]),
                          maps,
                          X0,nobs,Xsig2,Ysig2,mu,lookup,
                          nts,seed,
                          verbose=verbose)
  
  #format output
  x<-tmp[[1]]
  trait.data<-tmp[[2]]
  attr(x,'ts')<-ts
  attr(x,'tree')<-tree
  attr(x,'maps')<-maps
  attr(x,'treeID')<-treeID
  traits<-colnames(Xsig2[[1]])
  attr(x,'traits')<-setNames(rep(1,length(traits)),traits)
  #need to first reconfigure lookup to refer to include newly simulated trait.data...
  counter<-1
  for(i in seq_along(lookup)){
    tmp.n<-length(lookup[[i]][['matches']])
    lookup[[i]][['table']]<-lookup[[i]][['table']][lookup[[i]][['matches']],,drop=FALSE]
    lookup[[i]][['matches']]<-seq_len(tmp.n)
    lookup[[i]][['table']][,1]<-counter:(counter+tmp.n-1)
    counter<-counter+tmp.n
  }
  attr(x,'params')<-matrix(list(trait.data,Xsig2,Ysig2,mu,lookup,list(fxn="sim.conthistory")),
                           6,1,
                           dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  .uncompress(x)
}

#may want to make some function for auto-estimation of Xsig2/Ysig2
#' @export
make.contsimmap<-function(tree,trait.data,
                          nsims=100,res=100,
                          Xsig2=NULL,Ysig2=NULL,mu=NULL,
                          verbose=FALSE){
  #takes care of all the necessary formatting work...
  list2env(.prep.contsimmap(tree,nsims,res,trait.data,Xsig2,Ysig2,mu,TRUE),envir=environment())
  
  #does the post/preorder traversals, returns final simulated trait data
  x<-.cond.traversals(prune.seq,anc,des,ndes,
                      maps,
                      parsed.obs,parsed.mis,nobs,Xsig2,Ysig2,mu,lookup,
                      nts,NTS,t1s,seed,x,v,dx,dv,
                      verbose=verbose)
  
  #format output
  attr(x,'ts')<-ts
  attr(x,'tree')<-tree
  attr(x,'maps')<-maps
  attr(x,'treeID')<-treeID
  traits<-colnames(Xsig2[[1]])
  attr(x,'traits')<-setNames(rep(1,length(traits)),traits)
  attr(x,'params')<-matrix(list(trait.data,Xsig2,Ysig2,mu,lookup,list(fxn="make.contsimmap")),
                           6,1,
                           dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  .uncompress(x)
}