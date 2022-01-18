#I think some of this might need some cleaning up...just the organization of indices alone is getting crazy!
#Perhaps not the most efficient functions...but I think it all works reasonably fast

.get.perm.inds<-function(treeID){
  treeID<-as.integer(factor(treeID,unique(treeID)))
  nsims<-length(treeID)
  out<-matrix(treeID,nrow=nsims,ncol=3)
  per.tree<-tabulate(treeID)
  ntrees<-length(per.tree)
  coarse.inds<-c(0,cumsum(per.tree[-ntrees]))[treeID]
  fine.inds<-lapply(per.tree,seq_len)
  for(i in seq_along(fine.inds)){
    inds<-treeID==i
    out[inds,2]<-fine.inds[[i]]
    out[inds,3]<-coarse.inds[inds]+out[inds,2]
  }
  out
}

.get.cladewise.ord<-function(anc){
  anc<-as.numeric(anc)
  ord<-order(anc)
  #this was the fastest I could come up with
  for(i in seq_along(ord)){
    inds<-which(anc==ord[i])
    if(length(inds)){
      pos<-match(inds,ord)
      after<-i-sum(pos<i)
      ord<-append(ord[-pos],inds,after)
    }
  }
  ord
}

#assume numbers correspond to edge names, rather than indices-->I think this will be less confusing
.get.edges<-function(x,contsimmap){
  nedges<-length(contsimmap[['ts']])
  edges<-names(contsimmap[['ts']])
  if(is.numeric(x)){
    x<-as.character(x)
  }
  if(is.character(x)){
    out.of.range<-is.na(match(x,edges))
    if(all(out.of.range)){
      stop("all specified edges don't exist; edge indices must match those in contsimmap exactly")
    }
    if(any(out.of.range)){
      warning("some specified edges don't exist; edge indices must match those in contsimmap exactly")
      x<-x[!out.of.range]
    }
  }else{
    stop("invalid edges specified; please provide integers or names specifying edges to select")
  }
  contsimmap[['anc']]<-contsimmap[['anc']][x]
  nodes.to.store<-contsimmap[['anc']][!(contsimmap[['anc']]%in%x)&!(contsimmap[['anc']]%in%names(contsimmap[['nodes']]))]
  if(length(nodes.to.store)){
    out<-matrix(nrow=length(contsimmap[['traits']]),ncol=nrow(contsimmap[['perm.inds']]))
    foo<-function(ii){
      out[]<-unlist(lapply(ii,function(jj) jj[dim(jj)[1],,,drop=FALSE]))
      out[,contsimmap[['perm.inds']][,3],drop=FALSE]
    }
    contsimmap[['nodes']][nodes.to.store]<-setNames(lapply(contsimmap[['x']][nodes.to.store],foo),nodes.to.store)
    contsimmap[['nodes']]<-contsimmap[['nodes']][names(contsimmap[['nodes']])%in%contsimmap[['anc']]]
  }
  for(i in c('x','maps','ts')){
    contsimmap[[i]]<-contsimmap[[i]][x]
  }
  contsimmap
}

.get.sims<-function(x,contsimmap){
  nsims<-nrow(contsimmap[['perm.inds']])
  if(is.numeric(x)&length(x)>1){
    x<-list(x)
  }
  if(is.numeric(x)){
    if(!(x>0&x<=nsims)){
      stop("can only sample between 1 and ",nsims," simulations")
    }
    x<-sample(nsims,x)
  }else if(is.list(x)){
    x<-unlist(x)
    if(!is.numeric(x)){
      stop("invalid simulations specified; please provide a single integer specifying a number of simulations to randomly sample or a list of indices")
    }
    out.of.range<-!(x>0&x<=nsims)
    if(all(out.of.range)){
      stop("all specified simulations don't exist; only indices between 1 and ",nsims," are allowed")
    }
    if(any(out.of.range)){
      warning("some specified simulations don't exist; only indices between 1 and ",nsims," are allowed")
      x<-x[!out.of.range]
    }
  }else{
    stop("invalid simulations specified; please provide a single integer specifying a number of simulations to randomly sample or a list of indices")
  }
  if(!is.null(contsimmap[['trait.data']])){
    if(dim(contsimmap[['trait.data']])[3]>1){
      contsimmap[['trait.data']]<-contsimmap[['trait.data']][,,x,drop=FALSE]
    }
  }
  inds<-contsimmap[['perm.inds']]
  inds<-inds[x,,drop=FALSE]
  treeID<-unique(inds[,1])
  tmp<-split(inds[,2],inds[,1])
  tmp<-tmp[as.character(treeID)]
  tree.seq<-seq_along(tmp)
  foo<-function(ii){
    for(j in tree.seq){
      ii[[j]]<-ii[[j]][,,tmp[[j]],drop=FALSE]
    }
    ii
  }
  contsimmap[['x']]<-lapply(contsimmap[['x']],function(ii) foo(ii[treeID]))
  contsimmap[['nodes']]<-lapply(contsimmap[['nodes']],function(ii) ii[,x,drop=FALSE])
  contsimmap[['maps']]<-lapply(contsimmap[['maps']],'[',treeID)
  contsimmap[['tree']]<-contsimmap[['tree']][treeID]
  contsimmap[['perm.inds']]<-.get.perm.inds(inds[,1])
  contsimmap
}

.get.traits<-function(x,contsimmap){
  traits<-contsimmap[['traits']]
  ntraits<-length(traits)
  nas<-is.na(x)
  if(is.numeric(x)){
    out.of.range<-!(x>0&x<=ntraits)
    out.of.range[nas]<-FALSE
    if(all(out.of.range)){
      stop("all specified traits don't exist; only trait indices between 1 and ",ntraits," are allowed")
    }
    if(any(out.of.range)){
      warning("some specified traits don't exist; only trait indices between 1 and ",ntraits," are allowed")
      x<-x[!out.of.range]
    }
  }else if(is.character(x)){
    x<-pmatch(x,traits)
    out.of.range<-is.na(x)
    out.of.range[nas]<-FALSE
    if(all(out.of.range)){
      stop("all specified traits don't exist")
    }
    if(any(out.of.range)){
      warning("some specified traits don't exist")
      x<-x[!out.of.range]
    }
  }else{
    stop("invalid traits specified; please provide integers or names specifying traits to select")
  }
  contsimmap[['x']]<-lapply(contsimmap[['x']],function(ii) lapply(ii,function(jj) jj[,x,,drop=FALSE]))
  contsimmap[['nodes']]<-lapply(contsimmap[['nodes']],function(ii) ii[x,,drop=FALSE])
  contsimmap[['traits']]<-contsimmap[['traits']][x]
  if(!is.null(contsimmap[['trait.data']])){
    contsimmap[['trait.data']]<-contsimmap[['trait.data']][,x,,drop=FALSE]
  }
  contsimmap
}

#' @export
subset.contsimmap<-function(contsimmap,edges=NULL,sims=NULL,traits=NULL){
  if(!is.null(edges)){
    contsimmap<-.get.edges(edges,contsimmap)
  }
  if(!is.null(sims)){
    contsimmap<-.get.sims(sims,contsimmap)
  }
  if(!is.null(traits)){
    contsimmap<-.get.traits(traits,contsimmap)
  }
  contsimmap
}