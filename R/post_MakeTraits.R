.fix.formulae<-function(formulae,traits){
  lens<-lengths(formulae)
  probs<-lens<3
  if(any(probs)){
    all.traits<-c(unique(unlist(lapply(formulae,all.vars),use.names=FALSE)),
                  traits)
    def.nms<-grepl('^trait_\\d+$',all.traits)
    if(any(def.nms)){
      def.nms.offset<-max(as.numeric(gsub('^trait_(\\d+)','\\1',all.traits[def.nms])))
    }else{
      def.nms.offset<-0
    }
    new.def.nms<-paste0('trait_',seq_len(sum(probs))+def.nms.offset)
    deparsed.formulae<-unlist(lapply(formulae[probs],deparse),use.names=FALSE)
    formulae[probs]<-lapply(paste0(new.def.nms,deparsed.formulae),str2lang)
  }
  formulae
}

.get.call.envir<-function(foc.traits,lhs,contsimmap,summarized.flag){
  flags<-setNames(rep(TRUE,4),paste0(c('n','t','dt','s'),'_RES'))
  inds<-rep(TRUE,length(foc.traits))
  for(i in names(flags)){
    tmp<-foc.traits==i
    inds<-inds&!tmp
    if(i!='n_RES'){
      flags[i]<-any(tmp)
    }
  }
  foc.traits<-foc.traits[inds]
  lhs<-lhs[inds]
  traits<-dimnames(contsimmap)[[2]]
  traits.ind<-match(foc.traits,traits)
  #use &lhs to exclude non-existent traits that aren't being created in call
  #(likely represent variables stored in global environment!)
  prob.traits<-is.na(traits.ind)&lhs
  if(any(prob.traits)){
    new.traits<-unique(foc.traits[prob.traits])
    n.new.traits<-length(new.traits)
    contsimmap<-contsimmap[,c(seq_along(traits),rep(NA,n.new.traits)),]
    traits<-dimnames(contsimmap)[[2]]<-c(traits,new.traits)
    traits.ind<-match(foc.traits,traits)
  }
  inds<-!is.na(traits.ind)
  foc.traits<-foc.traits[inds]
  traits.ind<-traits.ind[inds]
  lhs<-lhs[inds]
  nts<-.get.ns(contsimmap,uncompress=TRUE)
  tmp.seq<-seq_along(nts)
  splits<-factor(rep(tmp.seq,nts),levels=tmp.seq)
  n_RES<-length(splits)
  out<-setNames(vector("list",length(foc.traits)),foc.traits)
  tmp<-unclass(contsimmap)
  for(i in seq_along(foc.traits)){
    out[i]<-list(unlist(tmp[,traits.ind[i],,drop=FALSE],use.names=FALSE))
  }
  probs<-lengths(out)!=n_RES
  out[probs]<-list(rep(NA,n_RES))
  
  if(flags["n_RES"]){
    out["n_RES"]<-n_RES
  }
  focs<-names(flags[-1])[flags[-1]]
  for(i in focs){
    #could technically be made more efficient by uncompressing maps elements in one go...but probably makes only a small difference
    tmp<-.get.maps(contsimmap,switch(i,t_RES='ts',dt_RES='dts',s_RES='state'),uncompress=TRUE)
    if(summarized.flag&i=="s_RES"){
      out[[i]]<-do.call(rbind,tmp)
    }else{
      out[[i]]<-unlist(tmp,use.names=FALSE)
    }
  }
  c(out,contsimmap_RES=list(contsimmap),splits_RES=list(splits))
}

#for making new traits given arbitrary formula of existing traits
#' @export
make.traits<-function(contsimmap,...){
  summarized.flag<-inherits(contsimmap,"summarized_contsimmap")
  formulae<-.fix.formulae(list(...),contsimmap[['traits']])
  resps<-lapply(formulae,function(ii) all.vars(ii[[2]]))
  calls<-lapply(formulae,'[[',3)
  foc.traits<-unique(unlist(lapply(formulae,all.vars),use.names=FALSE))
  lhs<-foc.traits%in%unlist(resps,use.names=FALSE)
  call.envir<-.get.call.envir(foc.traits,lhs,contsimmap,summarized.flag)
  check.length<-function(x){
    len<-length(x)
    if(len==1){
      x<-rep(x,call.envir[['n_RES']])
    }else if(len!=call.envir[['n_RES']]){
      stop("New trait is of wrong length: double-check provided formulae!")
    }
    x
  }
  #depending on how things go, may want re-integrate traits into contsimmap within this loop...
  for(i in seq_along(formulae)){
    tmp.resp<-resps[[i]]
    call.envir[["tmp.resp_RES"]]<-tmp.resp
    tmp.res<-eval(calls[[i]],envir=call.envir,enclos=parent.frame())
    if(length(tmp.resp)==1){
      if(is.list(tmp.res)){
        tmp.res<-tmp.res[[1]]
      }
      call.envir[[tmp.resp]]<-check.length(tmp.res)
    }else{
      if(!is.list(tmp.res)){
        tmp.res<-list(tmp.res)
      }
      tmp.res<-lapply(tmp.res,check.length)
      #maybe should make this check stricter--tmp.res is recycled to match length of tmp.resp
      call.envir[tmp.resp]<-tmp.res
    }
  }
  
  #re-forming contsimmap
  #deciding to drop keys and just store the appropriate traits as character vectors for now...may revisit later
  #changed my mind--threshold returns numerics again
  #BUT instead of having the key system, breaks and state_names are just stored in call_info in the params attribute
  #having threshold return character vectors just made indexing and plotting functions more annoying
  #for now, though, not doing any checks to enforce numeric data in contsimmaps
  #this provides more flexibility, but users just need to be aware that various problems will be caused by non-numeric trait data
  contsimmap<-call.envir[['contsimmap_RES']]
  old.traits<-attr(contsimmap,'traits')
  new.traits<-unique(dimnames(contsimmap)[[2]])
  new.traits<-setNames(old.traits[new.traits],new.traits)
  dups<-split(!duplicated(unlist(resps,use.names=TRUE),fromLast=TRUE),rep(seq_along(resps),lengths(resps)))
  any.dups<-unlist(lapply(dups,any),use.names=FALSE)
  old.params<-attr(contsimmap,'params')
  new.params<-matrix(list(NULL),
                     6,ncol(old.params)+sum(any.dups),
                     dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  new.params[,seq_len(ncol(old.params))]<-old.params
  contsimmap<-unclass(contsimmap)
  for(i in seq_along(resps)){
    if(any.dups[i]){
      ii<-resps[[i]][dups[[i]]]
      tmp.ind<-max(new.traits,na.rm=TRUE)+1
      new.traits[ii]<-tmp.ind
      #check if threshold call
      if(!is.null(attr(call.envir[[ii[1]]],'breaks'))){
        new.params[["call_info",tmp.ind]]<-list(fxn="threshold",
                                                state_names=attr(call.envir[[ii[1]]],'state_names'),
                                                breaks=attr(call.envir[[ii[1]]],'breaks'),
                                                formula=formulae[[i]])
        
        
      }else if(!is.null(attr(call.envir[[ii[1]]],'params'))){
        new.params[,tmp.ind]<-attr(call.envir[[ii[1]]],'params')
        new.params[["call_info",tmp.ind]]<-c(new.params[["call_info",tmp.ind]],
                                             formula=formulae[i])
      }else{
        new.params[["call_info",tmp.ind]]<-list(fxn="make.traits",
                                                formula=formulae[[i]])
      }
      for(jj in ii){
        contsimmap[,jj,]<-split(call.envir[[jj]],call.envir[["splits_RES"]])
      }
    }
  }
  traits.levs<-sort(unique(new.traits))
  new.params<-new.params[,traits.levs,drop=FALSE]
  new.traits<-setNames(match(new.traits,traits.levs),names(new.traits))
  attr(contsimmap,'traits')<-new.traits
  attr(contsimmap,'params')<-new.params
  if(summarized.flag){
    class(contsimmap)<-c("contsimmap","summarized_contsimmap")
  }else{
    class(contsimmap)<-"contsimmap"
  }
  contsimmap
}
