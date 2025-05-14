#' @export
accumulate<-function(x,FUN="+",...,forget.prev=FALSE,root.val=NULL){
  #now using "reserved" (RES) variables in call.envir rather than messing with scoping too much
  #ugh, for some reason you need to manually add call.envir (should always be parent.frame() I think to)
  #why do I have to do this here but not for other functions like threshold()???
  #works for now, but there's something I don't understand here!
  #oohhh, because the variables aren't directly called in the function I think...
  #I think dynGet is the "safest" thing to use here...
  contsimmap<-dynGet("contsimmap_RES")
  splits<-dynGet("splits_RES")
  
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(.stored.nodes(contsimmap))==1){
      stop('The accumulate() function will support single subtrees in the future, but not yet')
    }else{
      stop('The accumulate() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  
  ntraits<-length(x)
  
  #try to implement better argument handling going forward...
  args.ls<-c(rep(list(NA),ntraits*2),list(...))
  
  #some initialization to get things in the right shape...
  #lots of cannibalized code from diffusion()...
  holder1<-array(unclass(contsimmap)[,1,,drop=FALSE],dim(contsimmap)[-2],dimnames(contsimmap)[-2])
  
  treeID<-.get.treeID(contsimmap)
  tree.seq<-sort(unique(treeID))
  
  maps<-attr(contsimmap,"maps")
  dims<-dim(maps)[-3]
  dims[1]<-dims[1]+1
  dimnms<-dimnames(maps)[-3]
  dimnms[[1]]<-c("N0",dimnms[[1]])
  holder2<-array(list(),dims,dimnms)
  
  edge.seq<-seq_len(dims[1])
  
  tree.inds<-lapply(tree.seq,'==',treeID)
  edge.inds<-match(dimnms[[1]],dimnames(holder1)[[1]])
  foo<-function(x){
    x<-x[edge.inds,,drop=FALSE]
    for(i in seq_along(tree.seq)){
      holder2[,i]<-lapply(edge.seq,function(ii) do.call(cbind,x[ii,tree.inds[[i]]]))
    }
    holder2
  }
  
  tmp.xs<-vector("list",ntraits)
  for(i in seq_along(tmp.xs)){
    tmp<-dynGet(x[i],ifnotfound=NULL)
    if(is.null(tmp)){
      holder1[]<-unclass(contsimmap)[,x[i],,drop=FALSE]
    }else{
      holder1[]<-split(tmp,splits)
    }
    tmp.xs[[i]]<-foo(holder1)
  }
  
  tree<-attr(contsimmap,"tree")[[1]]
  attr(tree,"order")<-NULL
  prune.seq<-c(reorder(tree,order="postorder",index.only=TRUE)+1,1)
  anc<-anc.edges(contsimmap)
  anc[lengths(anc)==0]<-list(NA)
  anc<-unlist(anc,use.names=FALSE)
  anc<-anc+1
  anc[is.na(anc)]<-1
  anc<-c(0,anc)
  nts<-rbind(1,.get.ns(contsimmap))
  
  #hack to more easily specify root/initial values...
  if(forget.prev){
    #now do this BEFORE specifying root values to ensure root value is properly "forgotten"
    #could be more efficient here for sure
    #because all tmp.xs elements other than the first have forget.prev=TRUE implicitly...
    tmp.xs.old<-tmp.xs
  }
  if(!is.null(root.val)){
    root.val<-rep(root.val,length.out=length(treeID))
    for(t in tree.seq){
      tmp.xs[[1]][[1,t]][]<-root.val[tree.inds[[t]]]
    }
  }
  
  if(forget.prev){
    for(e in rev(prune.seq)[-1]){
      for(t in tree.seq){
        for(l in seq_len(ntraits)){
          args.ls[[l*2-1]]<-tmp.xs.old[[l]][[anc[[e]],t]][nts[anc[[e]],t],,drop=FALSE]
          args.ls[[l*2]]<-tmp.xs.old[[l]][[e,t]][1,,drop=FALSE]
        }
        tmp.xs[[1]][[e,t]][1,]<-do.call(FUN,args.ls)
        for(k in seq_len(nts[e,t])[-1]){
          for(l in seq_len(ntraits)){
            args.ls[[l*2-1]]<-tmp.xs.old[[l]][[e,t]][k-1,,drop=FALSE]
            args.ls[[l*2]]<-tmp.xs.old[[l]][[e,t]][k,,drop=FALSE]
          }
          tmp.xs[[1]][[e,t]][k,]<-do.call(FUN,args.ls)
        }
      }
    }
  }else{
    for(e in rev(prune.seq)[-1]){
      for(t in tree.seq){
        for(l in seq_len(ntraits)){
          args.ls[[l*2-1]]<-tmp.xs[[l]][[anc[[e]],t]][nts[anc[[e]],t],,drop=FALSE]
          args.ls[[l*2]]<-tmp.xs[[l]][[e,t]][1,,drop=FALSE]
        }
        tmp.xs[[1]][[e,t]][1,]<-do.call(FUN,args.ls)
        for(k in seq_len(nts[e,t])[-1]){
          for(l in seq_len(ntraits)){
            args.ls[[l*2-1]]<-tmp.xs[[l]][[e,t]][k-1,,drop=FALSE]
            args.ls[[l*2]]<-tmp.xs[[l]][[e,t]][k,,drop=FALSE]
          }
          tmp.xs[[1]][[e,t]][k,]<-do.call(FUN,args.ls)
        }
      }
    }
  }
  
  out<-matrix(list(),
              dims[1],length(treeID),
              dimnames=list(dimnames(tmp.xs[[1]])[[1]],NULL))
  sims.per.tree<-unlist(lapply(tree.inds,sum),use.names=FALSE)
  #reshape one tree index at a time...will be easiest this way
  for(i in tree.seq){
    out[,tree.inds[[i]]]<-matrix(unlist(lapply(tmp.xs[[1]][,i],asplit,2),use.names=FALSE,recursive=FALSE),
                                 dims[1],sims.per.tree[i],byrow=TRUE)
  }
  out<-out[dimnames(holder1)[[1]],,drop=FALSE]
  
  #need to improve param info based on additional functionality (2/14/25)
  params<-matrix(list(NULL,NULL,NULL,NULL,NULL,list(fxn="accumulate",FUN=FUN,args=list(...))),
                 6,1,
                 dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  out<-unlist(out,use.names=FALSE)
  attr(out,'params')<-params
  out
}

# accumulate<-function(x,FUN="+",...,forget.prev=FALSE){
#   #now using "reserved" (RES) variables in call.envir rather than messing with scoping too much
#   #ugh, for some reason you need to manually add call.envir (should always be parent.frame() I think to)
#   #why do I have to do this here but not for other functions like threshold()???
#   #works for now, but there's something I don't understand here!
#   #oohhh, because the variables aren't directly called in the function I think...
#   #I think dynGet is the "safest" thing to use here...
#   contsimmap<-dynGet("contsimmap_RES")
#   splits<-dynGet("splits_RES")
#   
#   if(nedge(contsimmap)<Nedge(contsimmap)){
#     if(length(.stored.nodes(contsimmap))==1){
#       stop('The accumulate() function will support single subtrees in the future, but not yet')
#     }else{
#       stop('The accumulate() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
#     }
#   }
#   
#   args.ls<-c(NA,NA,list(...))
#   
#   #some initialization to get things in the right shape...
#   #lots of cannibalized code from diffusion()...
#   holder1<-array(unclass(contsimmap)[,1,,drop=FALSE],dim(contsimmap)[-2],dimnames(contsimmap)[-2])
#   maps<-attr(contsimmap,"maps")
#   dims<-dim(maps)[-3]
#   dims[1]<-dims[1]+1
#   dimnms<-dimnames(maps)[-3]
#   dimnms[[1]]<-c("N0",dimnms[[1]])
#   holder2<-array(list(),dims,dimnms)
#   tree.seq<-seq_len(dims[2])
#   edge.seq<-seq_len(dims[1])
#   treeID<-.get.treeID(contsimmap)
#   tree.inds<-lapply(tree.seq,'==',treeID)
#   edge.inds<-match(dimnms[[1]],dimnames(holder1)[[1]])
#   foo<-function(x){
#     x<-x[edge.inds,,drop=FALSE]
#     for(i in seq_along(tree.seq)){
#       holder2[,i]<-lapply(edge.seq,function(ii) do.call(cbind,x[ii,tree.inds[[i]]]))
#     }
#     holder2
#   }
#   
#   tmp<-dynGet(x,ifnotfound=NULL)
#   if(is.null(tmp)){
#     holder1[]<-unclass(contsimmap)[,x,,drop=FALSE]
#   }else{
#     holder1[]<-split(tmp,splits)
#   }
#   tmp.x<-foo(holder1)
#   
#   tree<-attr(contsimmap,"tree")[[1]]
#   attr(tree,"order")<-NULL
#   prune.seq<-c(reorder(tree,order="postorder",index.only=TRUE)+1,1)
#   anc<-anc.edges(contsimmap)
#   anc[lengths(anc)==0]<-list(NA)
#   anc<-unlist(anc,use.names=FALSE)
#   anc<-anc+1
#   anc[is.na(anc)]<-1
#   anc<-c(0,anc)
#   nts<-rbind(1,.get.ns(contsimmap))
#   if(forget.prev){
#     tmp.x.old<-tmp.x
#     for(e in rev(prune.seq)[-1]){
#       for(t in tree.seq){
#         args.ls[[1]]<-tmp.x.old[[anc[[e]],t]][nts[anc[[e]],t],,drop=FALSE]
#         args.ls[[2]]<-tmp.x.old[[e,t]][1,,drop=FALSE]
#         tmp.x[[e,t]][1,]<-do.call(FUN,args.ls)
#         for(k in seq_len(nts[e,t])[-1]){
#           args.ls[[1]]<-tmp.x.old[[e,t]][k-1,,drop=FALSE]
#           args.ls[[2]]<-tmp.x.old[[e,t]][k,,drop=FALSE]
#           tmp.x[[e,t]][k,]<-do.call(FUN,args.ls)
#         }
#       }
#     }
#   }else{
#     for(e in rev(prune.seq)[-1]){
#       for(t in tree.seq){
#         args.ls[[1]]<-tmp.x[[anc[[e]],t]][nts[anc[[e]],t],,drop=FALSE]
#         args.ls[[2]]<-tmp.x[[e,t]][1,,drop=FALSE]
#         tmp.x[[e,t]][1,]<-do.call(FUN,args.ls)
#         for(k in seq_len(nts[e,t])[-1]){
#           args.ls[[1]]<-tmp.x[[e,t]][k-1,,drop=FALSE]
#           args.ls[[2]]<-tmp.x[[e,t]][k,,drop=FALSE]
#           tmp.x[[e,t]][k,]<-do.call(FUN,args.ls)
#         }
#       }
#     }
#   }
#   
#   out<-tmp.x<-matrix(unlist(lapply(tmp.x,asplit,2),use.names=FALSE,recursive=FALSE),
#                      dims[1],length(treeID),byrow=TRUE,
#                      dimnames=list(dimnames(tmp.x)[[1]],NULL))
#   sims.per.tree<-unlist(lapply(tree.inds,sum),use.names=FALSE)
#   base<-rep(FALSE,dims[2])
#   foo<-function(x){
#     base[x]<-TRUE
#     rep(base,sims.per.tree)
#   }
#   tree.inds2<-lapply(tree.seq,foo)
#   for(i in tree.seq){
#     out[,tree.inds[[i]]]<-tmp.x[,tree.inds2[[i]],drop=FALSE]
#   }
#   out<-out[dimnames(holder1)[[1]],,drop=FALSE]
#   params<-matrix(list(NULL,NULL,NULL,NULL,NULL,list(fxn="accumulate",FUN=FUN,args=list(...))),
#                  6,1,
#                  dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
#   out<-unlist(out,use.names=FALSE)
#   attr(out,'params')<-params
#   out
# }