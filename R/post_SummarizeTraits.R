#still have some work to do--need to update indexing/plotting/make.traits function to account for
#the weird ways states get processed, as well as the different style of simulation naming!
#3/10/23: the above is mostly done as far as I know, but I may have missed something

#' @export
summarize.traits<-function(contsimmap,traits=NULL,FUN="mean",...){
  if(!is.null(traits)){
    contsimmap<-contsimmap[,traits,]
  }
  #perhaps not the most efficient approach (especially if there's only 1 tree for example)
  #but works for now
  #this part creates the basic array
  dims<-dim(contsimmap)
  incl<-.get.maps(contsimmap,'incl',uncompress=TRUE)
  tmp.rep<-unlist(lapply(incl[,1,drop=FALSE],sum))
  new.dim1<-sum(tmp.rep)
  incl<-unlist(incl[,rep(seq_len(dims[3]),each=dims[2]),drop=FALSE],use.names=FALSE)
  res<-array(unlist(contsimmap,use.names=FALSE)[incl],c(new.dim1,dims[-1]))
  res<-apply(res,c(1,2),FUN,...)
  if(length(dim(res))==3){
    res<-aperm(res,c(2,3,1))
  }else{
    res<-array(res,c(dim(res),1))
  }
  new.dim3<-dim(res)[3]
  if(is.null(dimnames(res)[[3]])){
    tmp<-if(is.character(FUN)) FUN else deparse(substitute(FUN))
    dimnames(res)[[3]]<-paste0(tmp,if(new.dim3>1) seq_len(new.dim3) else "")
  }
  out<-unclass(contsimmap)[,,seq_len(new.dim3),drop=FALSE]
  tmp.seq<-seq_len(dims[1]*dims[2]*new.dim3)
  splits<-rep(factor(tmp.seq,levels=tmp.seq),rep(tmp.rep,dims[2]*new.dim3))
  out[]<-split(res,splits)
  dimnames(out)[[3]]<-paste0('tree1_',dimnames(res)[[3]])
  #now you need to reconfigure attributes...
  #ts
  ts<-attr(contsimmap,'ts')
  #tree
  tree<-attr(contsimmap,'tree')[[1]]
  tree[['maps']]<-tree[['mapped.edge']]<-tree[['node.states']]<-tree[['states']]<-tree[['Q']]<-tree[['logL']]<-NULL
  #maps
  tmp.seq<-seq_along(ts)
  lens<-lengths(ts)
  dts<-lapply(tmp.seq,function(ii) rep(ts[[ii]][2]-ts[[ii]][1],lens[ii]-1))
  tmp<-1/(1:max(lens))
  bb.dts<-lapply(tmp.seq,function(ii) tmp[lens[ii]:1])
  rev.ts<-lapply(ts,function(ii) rev(ii)-ii[1])
  bb.sds<-lapply(tmp.seq,function(ii) sqrt(dts[[ii]]*rev.ts[[ii]][-1]/rev.ts[[ii]][-lens[ii]]))
  incl<-.get.maps(contsimmap,'incl')
  state<-.get.maps(contsimmap,'state')
  tmp.state<-sort(unique(unlist(state,use.names=FALSE)))
  nstate<-length(tmp.state)
  ntree<-length(attr(contsimmap,'tree'))
  state[]<-lapply(seq_along(state),function(ii) state[[ii]][incl[[ii]]])
  foo<-function(x){
    tmp<-do.call(cbind,x)
    tmp[]<-match(tmp,tmp.state)
    mode(tmp)<-"numeric"
    tmp<-as.matrix(apply(tmp,1,tabulate,nbins=nstate)/ntree)
    if(nstate>1){
      tmp<-t(tmp)
    }
    colnames(tmp)<-tmp.state
    tmp
  }
  state<-apply(state,1,foo)
  incl<-lapply(tmp.seq,function(ii) rep(TRUE,lens[ii]))
  inds<-as.list(lens)
  tmp.dims<-dim(attr(contsimmap,'maps'))
  tmp.dims[2]<-1
  maps<-array(c(as.list(tree[['edge.length']]),lapply(ts,'[',-1),dts,bb.dts,bb.sds,state,incl,inds),
              tmp.dims,dimnames(attr(contsimmap,'maps')))
  #traits
  traits<-attr(contsimmap,'traits')
  traits[]<-1
  #params
  #should you deparse? Leaning against it for now...
  if(!is.character(FUN)) FUN<-deparse(substitute(FUN))
  params<-matrix(list(NULL,NULL,NULL,NULL,NULL,list(fxn="summarize.traits",FUN=FUN,args=list(...))),
                 6,1,
                 dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  attr(out,'ts')<-ts
  attr(out,'tree')<-list(tree)
  attr(out,'maps')<-maps
  attr(out,'traits')<-traits
  attr(out,'params')<-params
  class(out)<-c("contsimmap","summarized_contsimmap")
  out
}
