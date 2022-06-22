#function to put in make.traits formulae to threshold traits!
#' @export
threshold<-function(x,breaks=NULL,nbreaks=10,
                    state.names=NULL,numeric.names=FALSE,
                    ...){
  #nbreaks processing
  if(is.null(breaks)){
    if(is.null(nbreaks)){
      stop('PLACEHOLDER ABOUT NO BREAKS SPECIFIED')
    }
    ran<-range(x,na.rm=TRUE)
    breaks<-seq(ran[1],ran[2],length.out=nbreaks+2)[-c(1,nbreaks+2)]
  }
  #naming
  term.nm<-NA
  if(!is.null(state.names)){
    nnms<-length(state.names)
    nbreaks<-length(breaks)
    if(nnms==nbreaks){
      names(breaks)<-state.names
    }else{
      if(nnms==nbreaks+1){
        names(breaks)<-state.names[-nnms]
        term.nm<-state.names[nnms]
      }else{
        stop('PLACEHOLDER ABOUT STATE.NAMES BEING OF WRONG LENGTH')
      }
    }
  }
  nms<-names(breaks)
  if(is.null(nms)){
    nms<-rep(NA,length(breaks))
  }
  ord<-order(breaks)
  breaks<-breaks[ord]
  nms<-nms[ord]
  nms<-c(nms,term.nm) #for last, always unnamed, group
  prob.nms<-!nzchar(nms)|is.na(nms)
  len<-length(nms)
  if(any(prob.nms)){
    if(numeric.names){
      nn<-ceiling(log(len,10))
      tmp<-do.call(expand.grid,rep(list(as.character(c(0,seq_len(9)))),nn))
    }else{
      nn<-ceiling(log(len,26))
      tmp<-do.call(expand.grid,rep(list(LETTERS),nn))
    }
    tmp<-as.matrix(tmp)[seq_len(len),,drop=FALSE]
    nms[prob.nms]<-do.call(paste,c(rev(asplit(tmp,2)),sep=''))[prob.nms]
  }
  #output
  out<-findInterval(x,breaks,...)+1
  attr(breaks,'state_names')<-nms
  attr(out,'keys')<-unname(breaks)
  out
}

#function for fixing lists of (potentially variable) rates, used in diffusion() and quick.sim()
.fix.var.list<-function(Xvar,traits,avail.traits){
  var.rates<-unlist(lapply(Xvar,is.character),use.names=FALSE)
  con.rates<-unlist(lapply(Xvar,is.numeric),use.names=FALSE)
  #check variable rates match up to a trait...needs a warning as well
  tmp<-match(unlist(Xvar[var.rates],use.names=FALSE),
             avail.traits)
  var.rates[var.rates]<-!is.na(tmp)
  #should probably give a warning about defaulting to 1...
  Xvar[!(var.rates|con.rates)]<-list(1)
  #should do length check and warning before doing this, probably...
  Xvar<-lapply(Xvar,'[',1)
  nms<-names(Xvar)
  if(is.null(nms)|all(!nzchar(nms))){
    names(Xvar)<-traits
  }
  Xvar[traits]
}

#more efficient alternative to traits2rates in the case that some other trait is only affecting the rate at which
#a new trait evolves
#may be nice to make this eventually store the extra parameters created
#as it stands right now, I think this might be messed up by reordering simulations?
#I think this is the way to go potentially--can do anything traits2rates does if Xcor doesn't vary as well
#ALMOST works now--just need to figure out a way to make .extract.traits work nicely with multiple outputs...
#works fully, but I need to remove the "double formula" effect to be more consistent and less prone to bugs
#Instead, just make it so that rates can be tied to already made traits through some kind of list indicating vectors to tie to other traits
#Eventually, extend this to be able to condition on trait values to make contsimmaps under things like evorates models!
#6/3: In a future update, make it such that ... can only be numeric or character arguments
#   Parse out list(...) [assigning new trait names as necessary?], default unmatched arguments to rate 1 like elsewhere
#   Fill out Xsig2 accordingly; numeric elements of list(...) can just multiply the appropriate entries of Xsig2
#   Should unmatched traits in Xsig2 be added? Would be difficult, probs not worth it...
#   Then you can just go through the contsimmap like before without the extra make.traits call
#   And you won't have the weird double formula effect
#   Hmmm...only issue is how you recognize the new traits to construct
#   Right now, you rely on the formula format to return multiple traits in a single formula
#   You'd have to figure out some way to grab the trait names in this function
#   One method might be keeping the formula format, but doing checks to make sure they only include numeric/character variables...
#   Another, better method might be using multiple response formulae of the form c(x,y,z) and using this to inform ntraits and traits...
#' @export
diffusion<-function(...,Xsig2=NULL,X0=0){
  #grab contsimmap and response variables from most recent make.traits call
  #this might be REALLY fragile...but it seems to get the job done for now
  tmp<-unlist(lapply(sys.calls(),function(ii) as.character(ii[[1]])),use.names=FALSE)
  tmp<-sys.frame(which.max(tmp=='make.traits'))
  contsimmap<-get('contsimmap',envir=tmp)
  traits<-get('tmp.resp',envir=tmp)
  ntraits<-length(traits)
  
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(contsimmap[['nodes']])==1){
      stop('The diffusion() function will support single subtrees in the future, but not yet')
    }else{
      stop('The diffusion() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  
  #parse out list(...)
  Xvar<-.fix.var.list(list(...),traits,contsimmap[['traits']])
  var.rates<-unlist(lapply(Xvar,is.character),use.names=FALSE)
  con.rates<-unlist(lapply(Xvar,is.numeric),use.names=FALSE)
  
  #fix up Xsig2
  #still need to remember to make better X0 handling procedures...
  states<-.get.states(contsimmap[['tree']])
  nstates<-length(states)
  if(is.null(Xsig2)){
    Xsig2<-diag(ntraits)
  }
  Xsig2<-.fix.mat.list(Xsig2,ntraits,traits,nstates,states)
  if(any(con.rates)){
    scalars<-sqrt(unlist(Xvar[con.rates],use.names=FALSE))
    foo<-function(x){
      x[con.rates,]<-sweep(x[con.rates,,drop=FALSE],1,scalars,'*')
      x[,con.rates]<-sweep(x[,con.rates,drop=FALSE],2,scalars,'*')
      x
    }
    Xsig2<-lapply(Xsig2,foo)
  }
  #should probably do a check that the final Xsig2 is valid (no negative variances)
  nsims<-nsim(contsimmap)
  X0<-matrix(rep(X0,length.out=ntraits*nsims),ntraits,nsims)
  
  #simulate new diffusion traits
  tmp.ntraits<-ntrait(contsimmap)
  treeID<-contsimmap[['perm.inds']][,1,drop=FALSE]
  treeID<-as.integer(factor(treeID,unique(treeID)))
  sims.per.tree<-tabulate(treeID)
  nts<-do.call(cbind,lapply(contsimmap[['x']],function(ii) lengths(ii)/tmp.ntraits/sims.per.tree))
  seeds.per.tree.edge<-nts*sims.per.tree
  seeds.per.edge<-colSums(seeds.per.tree.edge)
  seed<-.get.seed(seeds.per.edge,ntraits)
  anc<-anc.edges(contsimmap[['tree']][[1]])
  anc[lengths(anc)==0]<-0
  anc<-as.character(unlist(anc))
  names(anc)<-seq_along(anc)
  cladewise.ord<-.get.cladewise.ord(anc)
  rate.inds<-rep(NA,ntraits)
  rate.inds[var.rates]<-match(unlist(Xvar[var.rates],use.names=FALSE),contsimmap[['traits']])
  new.trait.inds<-match(traits,contsimmap[['traits']])
  contsimmap[['x']]<-.accumulate.seed(seed,
                                      anc,cladewise.ord,contsimmap$maps,
                                      Xsig2,X0,
                                      ntraits,traits,nstates,states,
                                      treeID,nts,sims.per.tree,seeds.per.tree.edge,
                                      scalar.trait=TRUE,original=contsimmap[['x']],rate.inds=rate.inds,scalar.inds=new.trait.inds)
  contsimmap[['nodes']][['0']][new.trait.inds,]<-X0
  
  #unlist the new diffusion traits for re-integration
  #could probably be made more efficient
  nodes<-contsimmap[['nodes']]
  nnodes<-length(nodes)
  node.array<-array(unlist(nodes,use.names=FALSE),
                    c(dim(nodes[[1]]),nnodes),
                    list(NULL,NULL,names(nodes)))
  tmp.nodes<-asplit(node.array[new.trait.inds,,,drop=FALSE],1)
  edge.seq<-seq_along(contsimmap[['x']])
  tree.seq<-seq_along(sims.per.tree)
  out<-setNames(rep(contsimmap['x'],ntraits),traits)
  for(i in seq_len(ntraits)){
    for(j in edge.seq){
      for(k in tree.seq){
        out[[i]][[j]][[k]]<-out[[i]][[j]][[k]][,new.trait.inds[i],,drop=FALSE]
      }
    }
    out[[i]]<-c(tmp.nodes[[i]],unlist(out[[i]],use.names=FALSE))
  }
  
  out
}