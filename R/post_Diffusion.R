#function for fixing lists of (potentially variable) rates/mu, used in diffusion()
.fix.mods.list<-function(mods,traits,avail.traits){
  #form default output
  out<-setNames(rep(list(1,0),length(traits)),paste0(rep(traits,each=2),c("_Xsig2","_mu")))
  if(!length(mods)){
    mods<-list(NULL)
  }else{
    #start checking input
    var.mods<-unlist(lapply(mods,is.character),use.names=FALSE)
    con.mods<-unlist(lapply(mods,is.numeric),use.names=FALSE)
    #check variable rates match up to a trait...probably should come with a warning
    tmp<-match(unlist(mods[var.mods],use.names=FALSE),
               avail.traits)
    var.mods[var.mods]<-!is.na(tmp)
    #NEED names for now--maybe I'll think of a better system later
    inds<-var.mods|con.mods
    tmp<-match(names(mods)[inds],names(out))
    tmp<-tmp[!is.na(tmp)]
    #probably should do better length check/warning
    out[tmp]<-lapply(mods[inds],'[',1)
  }
  out
}

#one helpful thing might be not reordering edges (traversal functions technically equipped to deal with out-of-order edges)
#the only issue is if any edges are duplicated, the traversal functions will break...
#could be fixed with a simple match() though...
#actually, matching with names ONLY grabs the first matching name--it might all work...
#except for the assignment bit, where it will only assign to the first matching name too
#might be worth for ease of coding to fix this somewhere down the line...but not too important for now, I think
#' @export
diffusion<-function(...,Xsig2=NULL,Ysig2=NULL,mu=NULL,
                    trait.data=NULL,nobs=NULL,X0=NULL,
                    verbose=FALSE){
  #now using "reserved" (RES) variables in call.envir rather than messing with scoping too much
  #ugh, for some reason you need to manually add call.envir (should always be parent.frame() I think to)
  #why do I have to do this here but not for other functions like threshold()???
  #works for now, but there's something I don't understand here!
  #oohhh, because the variables aren't directly called in the function I think...
  #I think dynGet is the "safest" thing to use here...
  contsimmap<-dynGet("contsimmap_RES")
  splits<-dynGet("splits_RES")
  traits<-dynGet("tmp.resp_RES")
  ntraits<-length(traits)
  
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(.stored.nodes(contsimmap))==1){
      stop('The diffusion() function will support single subtrees in the future, but not yet')
    }else{
      stop('The diffusion() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  
  #parse out list(...)
  mods<-.fix.mods.list(list(...),traits,dimnames(contsimmap)[[2]])
  
  #some manual prep necessary since input could either specify conditional or unconditional contsimmap...
  if(!is.null(trait.data)){
    tips<-attr(contsimmap,'tree')[[1]][['tip.label']]
    nodes<-attr(contsimmap,'tree')[[1]][['node.label']]
    trait.data<-.fix.trait.data(trait.data,tips,nodes,traits)
    #add X0 to trait.data if trait.data is included (implies conditionality)
    if(!is.null(X0)){
      root<-as.character(Ntip(contsimmap)+1)
      foo<-function(x){
        probs<-rownames(x)==root
        if(any(probs)){
          warning("PLACEHOLDER warning about X0 overwriting any stuff in trait.data")
        }
        x[!probs,,drop=FALSE]
      }
      trait.data<-lapply(trait.data,foo)
      X0<-.fix.nobs.X0(X0,traits)
      n.trait.data<-length(trait.data)
      n.X0<-ncol(X0)
      colnames(X0)<-rep(root,n.X0)
      nn<-max(n.trait.data,n.X0)
      trait.data<-trait.data[rep(seq_len(n.trait.data),length.out=nn)]
      X0<-t(X0[,rep(seq_len(n.X0),length.out=nn),drop=FALSE])
      trait.data<-lapply(seq_len(nn),function(ii) rbind(trait.data[[ii]],X0[ii,,drop=FALSE]))
      #also ensure associated Ysig2 is 0!
      Ysig2<-.fix.param.list(Ysig2,traits,c(tips,nodes))
      Ysig2[,root]<-list(matrix(0,length(traits),length(traits),dimnames=list(traits,traits)))
    }
    #now set X0/nobs
    #I think .prep.contsimmap is fine with non-NULL ntraits and traits now, but I could be wrong...
    X0<-nobs<-NULL
    conditional<-TRUE
  }else{
    #swap trait.data for X0
    trait.data<-X0
    X0<-NULL
    conditional<-FALSE
  }
  
  #now I think we're good for this with a few edits...
  list2env(.prep.contsimmap(contsimmap,nsims=NULL,res=NULL,trait.data,Xsig2,Ysig2,mu,conditional=conditional,ntraits,traits,nobs,diffusion=TRUE),envir=environment())
  
  #eventually might want to generalize so that mods only affect particular states, but not a big deal for now I think...
  #find constant mods and incorporate into param list...
  con.mods<-unlist(lapply(mods,is.numeric),use.names=FALSE)
  Xsig2.con.mods<-which(con.mods[c(TRUE,FALSE)])
  if(length(Xsig2.con.mods)){
    tmp<-sqrt(unlist(mods[Xsig2.con.mods*2-1],use.names=FALSE))
    foo<-function(x){
      x[Xsig2.con.mods,]<-sweep(x[Xsig2.con.mods,,drop=FALSE],1,tmp,'*')
      x[,Xsig2.con.mods]<-sweep(x[,Xsig2.con.mods,drop=FALSE],2,tmp,'*')
      x
    }
    Xsig2[]<-lapply(Xsig2,foo)
  }
  mu.con.mods<-which(con.mods[c(FALSE,TRUE)])
  if(length(mu.con.mods)){
    tmp<-unlist(mods[mu.con.mods*2],use.names=FALSE)
    foo<-function(x){
      x[mu.con.mods]<-x[mu.con.mods]+tmp
      x
    }
    mu[]<-lapply(mu,foo)
  }
  #form lists variable mods
  #some initialization
  holder1<-array(unclass(contsimmap)[,1,,drop=FALSE],dim(contsimmap)[-2],dimnames(contsimmap)[-2])
  holder2<-array(list(),dim(seed),dimnames(seed))
  tree.seq<-dimnames(seed)[[2]]
  edge.seq<-seq_len(dim(seed)[1])
  tree.inds<-lapply(tree.seq,'==',treeID)
  edge.inds<-match(dimnames(seed)[[1]],dimnames(holder1)[[1]])
  foo<-function(x){
    x<-x[edge.inds,,drop=FALSE]
    for(i in seq_along(tree.seq)){
      holder2[,i]<-lapply(edge.seq,function(ii) do.call(cbind,x[ii,tree.inds[[i]]]))
    }
    holder2
  }
  var.mods<-unlist(lapply(mods,is.character),use.names=FALSE)
  #Xsig2 variable mods
  Xsig2.var.mods<-which(var.mods[c(TRUE,FALSE)])
  if(length(Xsig2.var.mods)){
    Xsig2.mods<-setNames(rep(list(NULL),ntraits),traits)
    nms<-unlist(mods[Xsig2.var.mods*2-1])
    for(i in seq_along(Xsig2.var.mods)){
      tmp<-dynGet(nms[i],ifnotfound=NULL)
      if(is.null(tmp)){
        holder1[]<-lapply(unclass(contsimmap)[,nms[i],,drop=FALSE],sqrt)
      }else{
        holder1[]<-split(sqrt(tmp),splits)
      }
      Xsig2.mods[[Xsig2.var.mods[i]]]<-foo(holder1)
    }
  }else{
    Xsig2.mods<-NULL
  }
  #mu variable mods
  mu.var.mods<-which(var.mods[c(FALSE,TRUE)])
  if(length(mu.var.mods)){
    mu.mods<-setNames(rep(list(NULL),ntraits),traits)
    nms<-unlist(mods[mu.var.mods*2])
    for(i in seq_along(mu.var.mods)){
      tmp<-dynGet(nms[i],ifnotfound=NULL)
      if(is.null(tmp)){
        holder1[]<-unclass(contsimmap)[,nms[i],,drop=FALSE]
      }else{
        holder1[]<-split(tmp,splits)
      }
      mu.mods[[mu.var.mods[i]]]<-foo(holder1)
    }
  }else{
    mu.mods<-NULL
  }
  
  #do traversal
  if(conditional){
    #yeah, definitely might be beneficial to have some univariate shortcut and/or a verbose setting to at least see progress...
    tmp<-.cond.traversals(prune.seq,anc,des,ndes,
                          maps,
                          parsed.obs,parsed.mis,nobs,Xsig2,Ysig2,mu,lookup,
                          nts,NTS,t1s,seed,x,v,dx,dv,
                          Xsig2.mods,mu.mods,
                          verbose=verbose) #rather slow, but whatcha gonna do?
  }else{
    tmp<-.uncond.traversals(prune.seq,anc,tree[[1]][['edge']],Ntip(tree[[1]]),
                            maps,
                            X0,nobs,Xsig2,Ysig2,mu,lookup,
                            nts,seed,
                            Xsig2.mods,mu.mods,
                            verbose=verbose)
    trait.data<-tmp[[2]]
    tmp<-tmp[[1]]
    counter<-1
    for(i in seq_along(lookup)){
      tmp.n<-length(lookup[[i]][['matches']])
      lookup[[i]][['table']]<-lookup[[i]][['table']][lookup[[i]][['matches']],,drop=FALSE]
      lookup[[i]][['matches']]<-seq_len(tmp.n)
      lookup[[i]][['table']][,1]<-counter:(counter+tmp.n-1)
      counter<-counter+tmp.n
    }
  }
  
  #forming final unlisted versions to return to make.traits...
  out<-tmp<-aperm(array(unlist(apply(tmp,1,function(ii) do.call(cbind,lapply(ii,asplit,c(2,3)))),recursive=FALSE),
                        c(ntraits,dim(contsimmap)[3],dim(tmp)[1]),list(NULL,NULL,dimnames(tmp)[[1]])),
                  c(3,1,2))
  base<-rep(FALSE,length(tree))
  foo<-function(x){
    base[x]<-TRUE
    rep(base,sims.per.tree)
  }
  out.lookup<-lookup<-do.call(rbind,lapply(lookup,function(ii) ii[['table']][ii[['matches']],seq_len(4),drop=FALSE]))
  tree.inds2<-lapply(seq_along(tree),foo)
  for(i in seq_along(tree)){
    out[,,tree.inds[[i]]]<-tmp[,,tree.inds2[[i]],drop=FALSE]
    out.lookup[tree.inds[[i]],]<-lookup[tree.inds2[[i]],,drop=FALSE]
  }
  out<-out[dimnames(holder1)[[1]],,,drop=FALSE]
  mods[con.mods]<-as.character(NA)
  mods<-unlist(mods,use.names=FALSE)
  Xsig2.mods<-setNames(mods[seq_len(3)*2-1],traits)
  mu.mods<-setNames(mods[seq_len(3)*2],traits)
  params<-matrix(list(trait.data,Xsig2,Ysig2,mu,out.lookup,list(fxn="diffusion",Xsig2.mods=Xsig2.mods,mu.mods=mu.mods)),
                 6,1,
                 dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  foo<-function(x){
    out<-unlist(out[,x,,drop=FALSE],use.names=FALSE)
    attr(out,'params')<-params
    out
  }
  lapply(seq_len(dim(out)[2]),foo)
}