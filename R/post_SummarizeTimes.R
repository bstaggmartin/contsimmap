
#convert to a "traits through time" plot, averaging not across simulations (like summarize.traits)
#but across edges!
#' @export
summarize.times<-function(contsimmap,traits=NULL,res=100,FUN="mean",...,
                          simplify=TRUE,linear.interpolation=TRUE){
  if(!is.null(traits)){
    contsimmap<-contsimmap[,traits,]
  }
  
  #reshaping--basically the opposite of .uncompress while also adding in the first time point for each edge
  #dealing with tree stuff
  treeID<-contsimmap:::.get.treeID(contsimmap)
  ntrees<-length(attr(contsimmap,'tree'))
  tree.seq<-seq_len(ntrees)
  tree.inds<-lapply(tree.seq,'==',treeID)
  sims.per.tree<-unlist(lapply(tree.inds,sum),use.names=FALSE)
  #dealing with edge stuff
  tmp<-unclass(contsimmap)
  edge.inds<-dimnames(tmp)[[1]]
  not.dups<-!duplicated(edge.inds)
  tmp<-tmp[not.dups,,,drop=FALSE]
  edge.inds<-edge.inds[not.dups]
  nodes<-grepl("^N",edge.inds)
  dimnames(tmp)[[1]][nodes]<-substr(dimnames(tmp)[[1]][nodes],2,nchar(dimnames(tmp)[[1]][nodes]))
  edge.inds[nodes]<-NA
  edge.inds<-as.numeric(edge.inds)
  nts<-contsimmap:::.get.ns(contsimmap)[edge.inds,,drop=FALSE]
  nts[nodes,]<-1
  rownames(nts)<-dimnames(tmp)[[1]]
  edges<-attr(contsimmap,'tree')[[1]][["edge"]]
  edge.inds[nodes]<-0
  anc<-as.character(match(edges[,1,drop=FALSE],edges[,2,drop=FALSE],nomatch=0L)[edge.inds])
  edge.inds<-as.character(edge.inds[!nodes])
  #trait stuff-simple
  ntraits<-dim(tmp)[2]
  #functions to actually reshape everything
  outer.foo<-function(i){
    a<-anc[i]
    e<-edge.inds[i]
    lapply(tree.seq,inner.foo,a=a,e=e)
  }
  inner.foo<-function(j,a,e){
    x0<-matrix(unlist(lapply(tmp[a,,tree.inds[[j]],drop=FALSE],'[',nts[a,j]),use.names=FALSE),ntraits,sims.per.tree[j])
    x1<-aperm(array(unlist(tmp[e,,tree.inds[[j]],drop=FALSE],use.names=FALSE),c(nts[e,j],ntraits,sims.per.tree[j])),c(2,3,1))
    aperm(array(c(x0,x1),c(ntraits,sims.per.tree[j],nts[e,j]+1)),c(3,1,2))
  }
  tmp.mat<-matrix(unlist(lapply(seq_along(edge.inds),outer.foo),recursive=FALSE,use.names=FALSE),
                         length(edge.inds),ntrees,
                  dimnames=list(edge.inds,NULL))
  
  
  erans<-edge.ranges(contsimmap)[as.numeric(edge.inds),,drop=FALSE]
  tpts<-seq(min(erans),max(erans),length.out=max(1,res+1))
  tips<-as.numeric(edge.inds)%in%tip.edges(contsimmap)
  # tips<-as.numeric(edge.inds)%in%match(1:Ntip(contsimmap),attr(contsimmap,'tree')[[1]][['edge']][,2])
  erans[tips,2]<-erans[tips,2,drop=FALSE]+1e-8 #generally seems to makes sure tips aren't truncated too early
  which.edges<-outer(erans[,1],tpts,FUN="<")&outer(erans[,2],tpts,FUN=">=")
  which.edges[,1]<-abs(tpts[1]-erans[,1])<1e-8
  nedges<-length(edge.inds)
  no.tpts<-.rowSums(which.edges,nedges,res+1)
  nts<-nts[edge.inds,,drop=FALSE]
  ts<-contsimmap:::.get.maps(contsimmap,'ts')[as.numeric(edge.inds),,drop=FALSE]
  ts[]<-lapply(seq_along(ts),function(ii) c(erans[(ii-1)%%nedges+1,1],ts[[ii]]))
  nsims<-dim(contsimmap)[3]
  inner.foo<-function(i,j,tmp.tpts){
    tmp.ts<-ts[[i,j]]
    tmp.xs<-tmp.mat[[i,j]]
    tmp.inds<-findInterval(tmp.tpts,tmp.ts)
    tmp.inds[tmp.inds>=nts[i,j]]<-nts[i,j]-1
    tmp.inds[tmp.inds<1]<-1
    tmp.inds.p1<-tmp.inds+1
    if(linear.interpolation){
      tmp.x0<-tmp.xs[tmp.inds,,,drop=FALSE]
      tmp.t0<-tmp.ts[tmp.inds]
      (tmp.xs[tmp.inds.p1,,,drop=FALSE]-tmp.x0)/(tmp.ts[tmp.inds.p1]-tmp.t0)*(tmp.tpts-tmp.t0)+tmp.x0
    }else{
      tmp.inds.inds<-abs(tmp.ts[tmp.inds]-tmp.tpts)>abs(tmp.ts[tmp.inds.p1]-tmp.tpts)
      tmp.inds[tmp.inds.inds]<-tmp.inds.p1[tmp.inds.inds]
      tmp.xs[tmp.inds,,,drop=FALSE]
    }
  }
  outer.foo<-function(i){
    out<-array(dim=c(no.tpts[i],ntraits,nsims))
    tmp<-lapply(tree.seq,inner.foo,i=i,tmp.tpts=tpts[which.edges[i,]])
    for(j in tree.seq){
      out[,,tree.inds[[j]]]<-tmp[[j]]
    }
    out
  }
  tmp.res<-lapply(seq_len(nedges),outer.foo)
  
  which.edges<-apply(which.edges,2,which)
  no.edges<-lengths(which.edges)
  counters<-rep(1,nedges)
  foo<-function(j,tmp){
    tmp.res[[j]][counters[j],,,drop=FALSE]
  }
  for(i in seq_len(res+1)){
    if(i==1){
      tmp<-apply(array(unlist(lapply(which.edges[[i]],foo,tmp=tmp),use.names=FALSE),c(ntraits,nsims,no.edges[i])),
                 c(1,2),FUN,...)
      if(length(dim(tmp))==3){
        dim4<-dim(tmp)[1]
        dimnms4<-dimnames(tmp)[[1]]
        tmp.reshape<-function(x){
          aperm(x,c(2,3,1))
        }
      }else{
        dim4<-1
        dimnms4<-NULL
        tmp.reshape<-function(x){
          array(x,c(ntraits,nsims,1))
        }
      }
      if(is.null(dimnms4)){
        dimnms4<-paste0(if(is.character(FUN)) FUN else deparse(substitute(FUN)),if(dim4>1) seq_len(dim4) else "")
      }
      out<-array(dim=c(res+1,ntraits,nsims,dim4),
                 dimnames=c(list(NULL),dimnames(contsimmap)[-1],FUN=list(dimnms4)))
      out[i,,,]<-tmp.reshape(tmp)
    }else if(no.edges[i]>0){
      out[i,,,]<-tmp.reshape(apply(array(unlist(lapply(which.edges[[i]],foo,tmp=tmp),use.names=FALSE),c(ntraits,nsims,no.edges[i])),
                                   c(1,2),FUN,...))
    }
    counters[which.edges[[i]]]<-counters[which.edges[[i]]]+1
  }
  if(simplify){
    dims<-dim(out)
    dimnms<-dimnames(out)
    out<-drop(out)
    for(i in which(dims==1)){
      attr(out,names(dimnms)[i])<-dimnms[[i]]
    }
  }
  attr(out,'ts')<-tpts
  attr(out,'FUN')<-FUN
  attr(out,'args')<-list(...)
  out
}