quick.sim<-function(contsimmap,ntraits=1,traits=paste0('trait_',seq_len(ntraits)),
                    Xvar=rep(list(1),ntraits),Xsig2=diag(ntraits),
                    X0=0,
                    nobs=0,Ysig2=0){
  
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(contsimmap[['nodes']])==1){
      stop('The quick.sim() function will support single subtrees in the future, but not yet')
    }else{
      stop('The quick.sim() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  
  states<-.get.states(contsimmap[['tree']])
  nstates<-length(states)
  Xsig2<-.fix.mat.list(Xsig2,ntraits,traits,nstates,states)
  traits<-rownames(Xsig2[[1]])
  tips<-tree[[1]]$tip.label
  ntips<-length(tips)
  out.trait.data<-if(!sum(nobs)) FALSE else TRUE
  if(out.trait.data){
    Ysig2<-.fix.mat.list(Ysig2,ntraits,traits,ntips,tips)
    nobs<-.fix.nobs(nobs,ntips,tips)
  }
  X0<-matrix(rep(X0,length.out=ntraits*nsims),ntraits,nsims)
}