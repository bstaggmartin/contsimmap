#much faster than old function, but doesn't do as much attempting to coerce trees to be identical
#makes sure node, edge, and tips information is all equivalent for a given set of trees
.check.matching.tops<-function(tree,nodes,nnodes,edges,nedges,tips,ntips,elen){
  # #old version
  # all(unlist(lapply(tree[-1],function(ii) all.equal(ii,tree[[1]]))))
  out<-TRUE
  for(i in seq_along(tree)[-1]){
    foc<-tree[[i]]
    tmp.nodes<-foc[['node.label']]
    tmp.nnodes<-foc[['Nnode']]
    if(tmp.nnodes!=nnodes){
      out<-FALSE
      break
    }
    if(!is.null(tmp.nodes)){
      #current has node labels, but first tree doesn't
      if(is.null(nodes)){
        out<-FALSE
        break
        #both trees have node labels, but they don't match
      }else if(any(tmp.nodes!=nodes)){
        out<-FALSE
        break
      }
      #otherwise, both have matching node labels
      #current doesn't have node labels, but first tree does
    }else if(!is.null(nodes)){
      out<-FALSE
      break
    }
    #otherwise both don't have node labels
    tmp.edges<-foc[['edge']]
    tmp.nedges<-nrow(tmp.edges)
    if(tmp.nedges!=nedges){
      out<-FALSE
      break
    }
    if(any(tmp.edges!=edges)){
      out<-FALSE
      break
    }
    tmp.elen<-foc[['edge.length']]
    if(any(tmp.elen!=elen)){
      out<-FALSE
      break
    }
    tmp.tips<-foc[['tip.label']]
    tmp.ntips<-length(tmp.tips)
    if(tmp.ntips!=ntips){
      out<-FALSE
      break
    }
    if(any(tmp.tips!=tips)){
      out<-FALSE
      break
    }
  }
  out
}

.get.states<-function(tree){
  sort(unique(unlist(lapply(tree,function(ii) colnames(ii$mapped.edge)))))
}

.get.tree.info<-function(tree,nsim){
  #checking to see if it's at least a phylo object
  if(!inherits(tree,'multiPhylo')&!inherits(tree,'phylo')){
    stop('tree must be of class phylo or simmap')
  }
  #checking to see if it includes multiple phylogenies, coercing it if it does
  if(!inherits(tree,'multiPhylo')){
    tree<-as.multiPhylo(tree)
    if(!is.null(tree[[1]]$maps)){
      class(tree)<-c(class(tree),'multiSimmap')
    }
  }
  if(any(!unlist(lapply(tree,is.rooted),use.names=FALSE))){
    stop('All phylogenies must be rooted')
  }
  if(any(unlist(lapply(tree,function(ii) is.null(ii[['edge.length']]))))){
    stop('All phylogenies must have edge lengths (ideally in units of time)')
  }
  nodes<-tree[[1]][['node.label']]
  if(!is.null(nodes)){
    if(any(duplicated(nodes,incomparables=c(NA,'')))){
      stop('Phylogenies cannot have duplicate node labels')
    }
  }
  nnodes<-tree[[1]][['Nnode']]
  edges<-tree[[1]][['edge']]
  nedges<-nrow(edges)
  tips<-tree[[1]][['tip.label']]
  if(any(duplicated(tips))){
    stop('Phylogenies cannot have duplicate tip labels')
  }
  ntips<-length(tips)
  elen<-tree[[1]][['edge.length']]
  if(length(tree)>1){
    if(!.check.matching.tops(tree,nodes,nnodes,edges,nedges,tips,ntips,elen)){
      stop('This function only accepts multiple phylogenies with identical topologies')
    }
  }
  #checking to see if it's a simmap or not
  if(is.null(tree[[1]]$maps)){
    states<-'0'
    maps<-as.list(elen)
    for(i in seq_along(maps)){
      maps[[i]]<-setNames(maps[[i]],'0')
    }
    mapped.edge<-as.matrix(unlist(maps,use.names=FALSE))
    colnames(mapped.edge)<-'0'
    for(i in seq_along(tree)){
      tree[[i]]$maps<-maps
      tree[[i]]$mapped.edge<-mapped.edge
    }
  }else{
    states<-.get.states(tree)
  }
  if(is.null(nodes)){
    nodes<-rep(NA,nnodes)
  }
  probs<-is.na(nodes)|!nzchar(nodes)
  nodes[probs]<-as.character(((ntips+1):(ntips+nnodes))[probs])
  for(i in seq_along(tree)){
    tree[[i]][['node.label']]<-nodes
  }
  ntrees<-length(tree)
  # treeID<-paste0('tree',rep(seq_len(ntrees),length.out=nsim),'_sim',rep(seq_len(ceiling(nsim/ntrees)),each=ntrees,length.out=nsim))
  treeID<-rep(seq_len(ntrees),length.out=nsim)
  list(tree=tree,treeID=treeID,ntrees=ntrees,
       nodes=nodes,
       edges=edges,
       tips=tips,
       states=states)
}

#"special function" that gets anc, des, sis, tips, edge lengths, and pruning sequence in one fell swoop
#also assigns root edge an index of 1 and offsets things automatically!
.get.topo.info<-function(phy){
  anc2<-match(phy[['edge']][,1],phy[['edge']][,2])+1
  anc2[is.na(anc2)]<-1
  anc<-c(list(NULL),as.list(anc2))
  anc2<-c(NA,anc2)
  len<-length(anc2)
  sis<-des<-rep(list(integer(0)),len)
  inds<-!is.na(anc)
  anc2<-anc2[inds]
  tmp<-seq_len(len)
  des2<-split(tmp[inds],anc2)
  des[as.numeric(names(des2))]<-des2
  ndes<-lengths(des2)
  di<-ndes==2
  di.des<-unlist(des2[di])
  if(length(di.des)>0){
    odds<-seq.int(1,length(di.des),2)
    evens<-odds+1
    odds<-di.des[odds]
    evens<-di.des[evens]
    sis[odds]<-evens
    sis[evens]<-odds
  }
  poly.des<-des2[!di]
  if(length(poly.des)>0){
    ndes<-ndes[!di]
    unlist.poly.des<-unlist(poly.des,use.names=FALSE)
    sis[unlist.poly.des]<-rep(poly.des,ndes)
    foo<-function(x){
      tmp<-sis[[x]]
      tmp[tmp!=x]
    }
    sis[unlist.poly.des]<-lapply(unlist.poly.des,foo)
  }
  tip<-which(!lengths(des))
  len<-c(0,phy[['edge.length']])
  attr(phy,'order')<-NULL
  prune.seq<-c(reorder(phy,'postorder',index.only=TRUE)+1,1)
  # prune.seq<-prune.seq[lengths(des[prune.seq])>0]
  list(anc=anc,des=des,sis=sis,tip=tip,len=len,prune.seq=prune.seq)
}

.get.edge.info<-function(tree,treeID,res){
  ntrees<-length(tree)
  nsims<-length(treeID)
  len<-tree[[1]][['edge.length']]
  nedges<-length(len)
  edgerans<-edge.ranges(tree[[1]])
  dt<-max(edgerans)/max(1,res)
  nts<-ceiling(len/dt)
  foo<-function(x){
    t1<-edgerans[x,1]
    t2<-edgerans[x,2]
    if(t2==t1){
      rep(t1,2)
    }else{
      seq(t1,t2,length.out=nts[x]+1)
    }
  }
  tmp<-seq_len(nedges)
  ts<-setNames(lapply(tmp,foo),tmp)
  dts<-len/nts
  dts[is.nan(dts)]<-0
  coarse.maps<-do.call(cbind,lapply(tree,'[[','maps'))
  single.state<-lengths(coarse.maps)==1
  foo<-function(x){
    ts<-ts[[x]]
    nt<-nts[x]
    dt<-dts[x]
    maps<-lapply(coarse.maps[x,],function(ii) cumsum(ii)+ts[1])
    single.state<-single.state[x,]
    out<-vector('list',ntrees)
    if(any(single.state)){
      tmp<-rev(ts-ts[1])[-1]
      bb.dts<-dt/(dt+tmp)
      out[single.state]<-lapply(maps[single.state],
                                function(ii) list(ts=ts[-1], #nearly 10-fold increase from storing as lists rather than data.frames
                                                  dts=rep(dt,nt),
                                                  bb.dts=bb.dts,
                                                  bb.sds=sqrt(tmp*bb.dts),
                                                  state=rep(names(ii),nt),
                                                  incl=rep(TRUE,nt),
                                                  inds=nt))
    }
    if(any(!single.state)){
      maps<-maps[!single.state]
      maps.seq<-seq_along(maps)
      nms<-rep(maps.seq,lengths(maps))
      ts<-sort(c(ts,unlist(maps)),index.return=TRUE)
      tmp.inds<-ts[['ix']]<(nt+2)
      inds<-lapply(maps.seq,function(ii) ts[['ix']]%in%(which(nms==ii)+nt+1)|tmp.inds)
      ts<-ts[['x']]
      for(j in maps.seq){
        map<-maps[[j]]
        t<-ts[inds[[j]]]
        tlen<-length(t)
        #switch in case rounding error makes last point not a state switch point
        if(!nzchar(names(t)[tlen])){
          names(t)[tlen-c(0,1)]<-names(t)[tlen-c(1,0)]
        }
        state<-rle(names(t))
        incl<-state
        incl[['values']]<-!nzchar(incl[['values']])
        state[['values']][incl[['values']]]<-state[['values']][-1][incl[['values']]]
        dt<-t[-1]-t[-length(t)]
        dups<-dt<1e-13 #tolerance seems to break down sometimes
        #may have to come up with better solution
        tmp.len<-sum(!dups)+1
        excl.inds<- -(which(dups)+1)
        dt<-t[excl.inds][-1]-t[excl.inds][-tmp.len]
        excl.inds<-c(-1,excl.inds)
        t<-t[excl.inds]
        state<-inverse.rle(state)[excl.inds]
        incl<-inverse.rle(incl)[excl.inds]
        tmp<-t[tmp.len-1]-t
        bb.dts<-rle(state)
        bb.dts[['values']]<-c(tmp[!incl],0)
        tmp<-tmp-inverse.rle(bb.dts)
        bb.dts<-dt/(dt+tmp)
        out[!single.state][[j]]<-list(ts=t,
                                      dts=dt,
                                      bb.dts=bb.dts,
                                      bb.sds=sqrt(tmp*bb.dts),
                                      state=state,
                                      incl=incl,
                                      inds=c(which(!incl),tmp.len-1))
      }
    }
    rbind(coarse.maps[x,,drop=FALSE],do.call(cbind,out))
  }
  maps<-aperm(array(unlist(lapply(tmp,foo),recursive=FALSE),
                    c(8,ntrees,nedges),
                    list(c('coarse','ts','dts','bb.dts','bb.sds','state','incl','inds'),NULL,tmp)),
              c(3,2,1))
  list(ts=ts,
       maps=maps)
}

.fix.trait.data<-function(trait.data,tips,nodes,traits=NULL){
  ntips<-length(tips)
  nnodes<-length(nodes)
  labs<-c(tips,nodes)
  if(is.data.frame(trait.data)|!is.list(trait.data)){
    trait.data<-list(trait.data)
  }
  for(i in seq_along(trait.data)){
    #initial processing into matrix
    if(length(dim(trait.data[[i]]))<2&is.numeric(trait.data[[i]])){
      input.code<-1
      tmp<-as.matrix(trait.data[[i]])
    }else if(is.data.frame(trait.data[[i]])){
      input.code<-2
      tmp<-as.matrix(trait.data[[i]][unlist(lapply(trait.data[[i]],function(ii) is.numeric(ii)),use.names=FALSE)])
      labels.cols<-unlist(lapply(trait.data[[i]],function(ii) is.character(ii)|is.factor(ii)),use.names=FALSE)
      labels<-trait.data[[i]][labels.cols]
      if(is.null(names(labels))){
        names(labels)<-which(labels.cols)
      }
      labels<-c(list(rownames=rownames(trait.data[[i]])),
                labels)
      labels<-labels[lengths(labels)>0]
      if(length(labels)>1){
        matches<-unlist(lapply(labels,function(ii) sum(!is.na(match(ii,labs)))),use.names=FALSE)
        if(max(matches)==0){
          stop('')
        }
        labels<-labels[which.max(matches)]
        warning('Multiple columns of label-like data in trait.data[[',i,
                ']]; column ',names(labels),' assumed to correspond to tip/node labels based on number of matches')
      }
      if(length(labels)>0){
        rownames(tmp)<-as.character(labels[[1]])
      }
    }else if(length(dim(trait.data[[i]]))==2&is.numeric(trait.data[[i]])){
      input.code<-3
      tmp<-trait.data[[i]]
    }else{
      stop('Format of trait.data[[',i,
           ']] not recognized')
    }
    #row checks--tip/node labels, completely missing observations...
    if(is.null(rownames(tmp))){
      nn<-nrow(tmp)
      if(nn==ntips+nnodes){
        rownames(tmp)<-c(tips,nodes)
        warning('No labels found for trait.data[[',
                i,
                ']]; assumed entries correspond to each tip and node in tree')
      }else if(nn==ntips){
        rownames(tmp)<-tips
        warning('No labels found for trait.data[[',
                i,
                ']]; assumed entries correspond to each tip in tree')
      }else{
        stop('No labels found for trait.data[[',
             i,
             ']]; please add a ',
             c('names attribute','column','rownames attribute')[input.code],
             ' providing labels to be matched with tips/nodes in tree')
      }
    }
    probs<-apply(tmp,1,function(ii) all(is.na(ii)))
    if(any(probs)){
      if(all(probs)){
        stop('All observations appear to be missing in trait.data[[',i,
             "]]; formatting error perhaps? If you're trying to simulate data unconditional on observations use sim.conthistory() instead")
      }
      warning(.report.names(c('Observation','Observations'),
                            which(probs),
                            paste0(c('is','are'),
                                   ' completely missing in trait.data[[',i,']] and ',
                                   c('was','were'),' excluded')))
      tmp<-tmp[!probs,,drop=FALSE]
    }
    probs<-is.na(match(rownames(tmp),labs))
    if(any(probs)){
      if(all(probs)){
        stop('No labels for trait.data[[',i,
             ']] matched with any tip/node labels in tree')
      }
      warning(.report.names(c('Label','Labels'),unique(rownames(tmp)[probs]),c("doesn't","didn't")),
              'match any tip/node labels in tree; associated observations excluded')
      tmp<-tmp[!probs,,drop=FALSE]
    }
    #column checks...
    colnames(tmp)<-.fix.trait.names(if(is.null(traits)) ncol(tmp) else traits,
                                    colnames(tmp))
    trait.data[[i]]<-tmp
  }
  if(is.null(traits)){
    traits<-unique(unlist(lapply(trait.data,colnames),use.names=FALSE))
  }
  trait.data<-lapply(trait.data,function(ii) ii[,traits,drop=FALSE])
  # #last check for invalid columns
  # #wait...does it work with "invalid" columns?
  # probs<-apply(do.call(rbind,trait.data),1,function(ii) all(is.na(ii)))
  # if(any(probs)){
  #   warning(.report.names(c('Trait','Traits'),
  #                         traits[probs],
  #                         paste0('completely ',
  #                                c('consists','consist'),
  #                                ' of missing observations and ',
  #                                c('was','were'),' excluded')))
  #   trait.data<-lapply(trait.data,function(ii) ii[,!probs,drop=FALSE])
  # }
  trait.data
}

.fix.trait.names<-function(k,nms){
  if(is.character(k)){
    #probably should have some warning messages in here...
    tmp<-match(k,nms)
    probs<-is.na(nms)|!nzchar(nms)
    probs<-probs&cumsum(probs)<=sum(is.na(tmp))
    nms[probs]<-k[is.na(tmp)]
  }else{ #assume numeric
    nms<-c(nms,rep(NA,k-length(nms)))
    nas<-which(is.na(nms)|!nzchar(nms))
    nms[nas]<-paste0('trait_',nas)
  }
  nms
}

.fix.nobs.X0<-function(l,states,ntips=NULL){
  nobs<-!is.null(ntips)
  nice.name<-if(nobs) "nobs" else "X0"
  nstates<-length(states)
  if(nobs){
    tmp.seq<-(ntips+1):nstates
  }
  if(is.null(l)){
    #new default behavior to be 1 if representing observations--otherwise no trait data generated!
    l<-if(nobs) 1 else 0
  }
  if(is.matrix(l)){
    l<-asplit(l,1)
  }else if(length(dim(l))>2|is.data.frame(l)){
    stop(nice.name," should either be a list of vectors of ",
         if(nobs) "numbers of observations" else "trait values at the root",
         " or a matrix with columns corresponding to ",
         if(nobs) "tips/nodes" else "traits")
  }
  if(!is.list(l)){
    l<-list(l)
  }
  #recycled from .fix.mu mainly
  base.vec<-setNames(rep(NA,nstates),states)
  out<-vector('list',length(l))
  for(i in seq_along(l)){
    if(!is.null(l[[i]])){
      out[[i]]<-base.vec
      if(length(dim(l[[i]]))>1){
        stop(nice.name,'[[',i,']] is multidimensional but should be a vector')
      }
      nms<-names(l[[i]])
      n<-length(l[[i]])
      if(!length(nms)){
        names(l[[i]])<-rep('',n)
        nms<-names(l[[i]])
        prob.nms<-rep(TRUE,n)
        has.prob.nms<-TRUE
        quiet.flag<-TRUE
      }else{
        prob.nms<-!(nms%in%states)
        has.prob.nms<-any(prob.nms)
        quiet.flag<-FALSE
      }
      if(any(!prob.nms)){
        out[[i]][nms[!prob.nms]]<-l[[i]][!prob.nms]
      }
      #set unspecified number of observations for internal nodes to 0
      #no warning for now--may want to change in the future
      if(nobs){
        out[[i]][tmp.seq][is.na(out[[i]][tmp.seq])]<-0
      }
      probs<-is.na(out[[i]])
      if(any(probs)){
        if(has.prob.nms){
          out[[i]][probs]<-l[[i]][prob.nms]
          message<-paste0("; missing entries of ",nice.name,"[[",i,"]] recycled from other entries of ",nice.name,"[[",i,"]]")
        }else{
          out[[i]][prob.traits]<-0
          message<-paste0("; missing entries of ",nice.name,"[[",i,"]] filled with 0s")
        }
        if(!quiet.flag){
          warning(.report.names(if(nobs) c('Tip','Tips') else c('Trait','Traits'),states[probs],c('has','have')),
                  ' no associated entries in ',nice.name,'[[',i,']]',
                  message,
                  immediate.=TRUE)
        }
      }
    }
  }
  probs<-!lengths(out)
  if(any(probs)){
    if(all(probs)){
      tmp<-list(setNames(rep(0,nstates),states))
      message<-paste0("; missing elements of ",nice.name," filled in with 0 vectors")
    }else{
      tmp<-out[!probs]
      message<-paste0("; missing elements of ",nice.name," recycled from other entries of ",nice.name)
    }
    out[probs]<-rep(tmp,length.out=sum(probs))
    warning("Some elements of ",nice.name," are NULL",message,
            immediate.=TRUE)
  }
  do.call(cbind,out)
}

.get.trait.info<-function(trait.data,tips,nodes,edges,
                          ntraits,traits,nobs,
                          conditional){
  if(conditional){
    trait.data<-.fix.trait.data(trait.data,tips,nodes,traits)
    ntraits<-ncol(trait.data[[1]])
    nms<-c(nodes[1],c(tips,nodes)[edges[,2]])
    mis<-lapply(trait.data,is.na)
    parsed.mis<-do.call(cbind,lapply(mis,function(ii) split(ii,rownames(ii))[nms]))
    foo<-function(x){
      trait.data[[x]][mis[[x]]]<-0
      split(trait.data[[x]],rownames(trait.data[[x]]))[nms]
    }
    parsed.obs<-do.call(cbind,lapply(seq_along(trait.data),foo))
    nobs<-lengths(parsed.obs)/ntraits
    tmp<-which(nobs>0)
    parsed.obs[tmp]<-lapply(tmp,function(ii) matrix(parsed.obs[[ii]],nobs[ii],ntraits))
    parsed.mis[tmp]<-lapply(tmp,function(ii) matrix(parsed.mis[[ii]],ntraits,nobs[ii],byrow=TRUE))
    list(traits=colnames(trait.data[[1]]),
         ntrait.data=length(trait.data),
         trait.data=trait.data,
         nobs=nobs,
         parsed.obs=parsed.obs,
         parsed.mis=parsed.mis)
  }else{
    nobs<-.fix.nobs.X0(nobs,c(tips,nodes),length(tips))
    traits<-.fix.trait.names(ntraits,traits)
    X0<-.fix.nobs.X0(trait.data,traits)
    list(traits=rownames(X0),
         ntrait.data=1,
         nobs=nobs,
         X0=X0)
  }
}


# recycling explanation: Generally speaking, named entries of input will be matched up to specified trait/state names, then any entries with either no
# or unmatched names are recycled in their given order to form an output of appropriate dimensions. When there are unmatched or unnamed entries, missing
# entries are filled in with default values: 1 for variances and 0 for covariances (i.e., the identity matrix).
# 
# In the case of matrices specifically, the function first attempts to "symmetrize" the input. For vector inputs, this creates a diagonal matrix. For
# matrix inputs, rows/columns with missing names are labelled either with their corresponding column/row names (if they exist) or a numerical index.
# Extra rows/columns of NAs are tacked on to the matrix so each row has a corresponding column and vice versa, the matrix is rearranged such that its row
# and column names are identical, and all NA entries i,j are replaced with non-NA entries j,i. From here, the output is formed from the symmetrized input
# as described above. Unmatched rows/columns of the input matrix are recycled into a block-diagonal matrix, and any unspecified covariances remaining are
# set to 0.

.symmetrize.mat<-function(mat,nice.name,Ysig2.flag){
  ndims<-length(dim(mat))
  if(ndims>2){
    stop(nice.name,' has more than two dimensions; multiple covariance matrices should be stored in one or two-dimensional lists, not arrays')
  }else if(ndims<2){
    if(length(mat)>1){
      nms<-names(mat)
      mat<-diag(mat)
      rownames(mat)<-colnames(mat)<-nms
    }else{
      mat<-matrix(mat,1,1,dimnames=rep(list(names(mat))))
    }
  }
  dims<-dim(mat)
  #get names
  nms<-dimnames(mat)
  if(!length(nms)) nms<-rep(list(NULL),2)
  for(i in c(1,2)){
    if(is.null(nms[[i]])) nms[[i]]<-rep('',dims[i])
  }
  minseq<-seq_len(min(dims))
  offset<-0
  for(i in c(1,2)){
    inds<-!nzchar(nms[[i]])
    other<-nms[[c(2,1)[i]]]
    other.inds<-nzchar(other)
    #fill in with names from other dimension, if they exist
    nms[[i]][minseq][inds[minseq]&other.inds[minseq]]<-other[other.inds[minseq]]
    inds<-!nzchar(nms[[i]])
    nms[[i]][inds]<-paste0('TEMPORARY_',seq_len(sum(inds))+offset)
    offset<-sum(inds)
  }
  dimnames(mat)<-nms
  #append rows corresponding to cols
  missing.rows<-!(nms[[2]]%in%nms[[1]])
  mat<-do.call(rbind,c(list(mat),rep(NA,sum(missing.rows))))
  rownames(mat)<-c(nms[[1]],nms[[2]][missing.rows])
  nms<-dimnames(mat)
  #append cols corresponding to rows
  missing.cols<-!(nms[[1]]%in%nms[[2]])
  mat<-do.call(cbind,c(list(mat),rep(NA,sum(missing.cols))))
  colnames(mat)<-c(nms[[2]],nms[[1]][missing.cols])
  #sort
  nms<-sort(nms[[1]])
  mat<-mat[nms,nms,drop=FALSE]
  #reflect
  lo<-lower.tri(mat)
  up<-upper.tri(mat)
  tmat<-t(mat)
  inds<-is.na(mat)
  inds<-inds&t(!inds)
  mat[lo&inds]<-tmat[lo&inds]
  mat[up&inds]<-tmat[up&inds]
  #fill in the rest with defaults
  diagonal<-diag(mat)
  probs<-is.na(diagonal)
  if(any(probs)){
    warning('PLACEHOLDER: diagonal fill-ins')
  }
  diagonal[probs]<-if(Ysig2.flag) 0 else 1
  diag(mat)<-diagonal
  probs<-is.na(mat)
  if(any(probs)){
    warning('PLACEHOLDER: off-diagonal fill-ins')
  }
  mat[probs]<-0
  mat
}

.fix.mat<-function(mat,nice.name,k,traits,Ysig2.flag){
  #form output
  base.mat<-matrix(NA,k,k,dimnames=rep(list(traits),2))
  out<-vector('list',length(mat))
  for(i in seq_along(mat)){
    if(!is.null(mat[[i]])){
      out[[i]]<-base.mat
      tmp.mat<-.symmetrize.mat(mat[[i]],nice.name(i),Ysig2.flag)
      nms<-dimnames(tmp.mat)[[1]]
      prob.nms<-grepl('^TEMPORARY_\\d+$',nms)
      if(all(prob.nms)){
        has.prob.nms<-TRUE
        quiet.flag<-TRUE
      }else{
        prob.nms<-prob.nms|!(nms%in%traits)
        has.prob.nms<-any(prob.nms)
        quiet.flag<-FALSE
      }
      if(any(!prob.nms)){
        out[[i]][nms[!prob.nms],nms[!prob.nms]]<-tmp.mat[!prob.nms,!prob.nms]
      }
      diagonal<-diag(out[[i]])
      prob.traits<-is.na(diagonal)
      if(any(prob.traits)){
        if(has.prob.nms){
          #take elements with no matching names and fill in
          nrem<-sum(prob.traits)
          tmp<-kronecker(diag(ceiling(nrem/sum(prob.nms))),tmp.mat[prob.nms,prob.nms])
          n<-nrow(tmp)
          incl.inds<-seq_len(n)<=nrem
          tmp<-tmp[incl.inds,incl.inds]
          out[[i]][prob.traits,prob.traits]<-tmp
          message<-paste0("; missing entries of ",nice.name(i)," recycled from other entries of ",nice.name(i))
        }else{
          #insert diagonal of 1s
          diagonal[prob.traits]<-if(Ysig2.flag) 0 else 1
          diag(out[[i]])<-diagonal
          message<-paste0("; missing entries of ",nice.name(i)," filled in with ",
                          if(Ysig2.flag) "0s" else "entries of identity matrix")
        }
        out[[i]][is.na(out[[i]])]<-0
        if(!quiet.flag){
          warning(.report.names(c('Trait','Traits'),traits[prob.traits],c('has','have')),
                  ' no associated entries in ',
                  nice.name(i),message,
                  immediate.=TRUE)
        }
      }
      if(!isSymmetric(out[[i]])|any(eigen(out[[i]])$values<0)){
        stop('Failed to create proper covariance matrix from ',nice.name(i))
      }
    }
  }
  probs<-!lengths(out)
  if(any(probs)){
    if(all(probs)){
      tmp<-diag(k)
      dimnames(tmp)<-list(traits,traits)
      if(Ysig2.flag) tmp<-0*tmp
      tmp<-list(tmp)
      message<-paste0("; missing elements of ",nice.name('')," filled in with ",
                      if(Ysig2.flag) "0" else "identity"," matrices")
    }else{
      tmp<-out[!probs]
      message<-paste0("; missing elements of ",nice.name('')," recycled from other entries of ",nice.name(''))
    }
    out[probs]<-rep(tmp,length.out=sum(probs))
    warning("Some elements of ",nice.name('')," are NULL",message,
            immediate.=TRUE)
  }
  out
}

.fix.mu<-function(mu,nice.name,k,traits){
  base.vec<-setNames(rep(NA,k),traits)
  out<-vector('list',length(mu))
  for(i in seq_along(mu)){
    if(!is.null(mu[[i]])){
      out[[i]]<-base.vec
      if(length(dim(mu[[i]]))>1){
        stop(nice.name(i),' is multidimensional; multiple mu vectors should be stored in one or two-dimensional lists, not matrices/arrays')
      }
      nms<-names(mu[[i]])
      n<-length(mu[[i]])
      if(!length(nms)){
        names(mu[[i]])<-rep('',n)
        nms<-names(mu[[i]])
        prob.nms<-rep(TRUE,n)
        has.prob.nms<-TRUE
        quiet.flag<-TRUE
      }else{
        prob.nms<-!(nms%in%traits)
        has.prob.nms<-any(prob.nms)
        quiet.flag<-FALSE
      }
      if(any(!prob.nms)){
        out[[i]][nms[!prob.nms]]<-mu[[i]][!prob.nms]
      }
      prob.traits<-is.na(out[[i]])
      if(any(prob.traits)){
        if(has.prob.nms){
          out[[i]][prob.traits]<-mu[[i]][prob.nms]
          message<-paste0("; missing entries of ",nice.name(i)," recycled from other entries of ",nice.name(i))
        }else{
          out[[i]][prob.traits]<-0
          message<-paste0("; missing entries of ",nice.name(i)," filled in with 0s")
        }
        if(!quiet.flag){
          warning(.report.names(c('Trait','Traits'),traits[prob.traits],c('has','have')),
                  ' no associated entries in ',
                  nice.name(i),message,
                  immediate.=TRUE)
        }
      }
    }
  }
  probs<-!lengths(out)
  if(any(probs)){
    if(all(probs)){
      tmp<-list(setNames(rep(0,k),traits))
      message<-paste0("; missing elements of ",nice.name('')," filled in with 0 vectors")
    }else{
      tmp<-out[!probs]
      message<-paste0("; missing elements of ",nice.name('')," recycled from other elements of ",nice.name(''))
    }
    out[probs]<-rep(tmp,length.out=sum(probs))
    warning("Some elements of ",nice.name('')," are NULL",message,
            immediate.=TRUE)
  }
  out
}

#may want a warning for empty state/trait names and trait names of pattern 'TEMPORARY_\\d+'
#now takes things in either a list-of-list format (not unlike data.frames) or a matrix list
#diff parameters are assumed to be columns, diff values for same parameters are assumed to be rows
#if given a vector list, assume it is for diff parameters unless there is only 1 parameter!
.fix.param.list<-function(l,traits,states){
  ntraits<-length(traits)
  nstates<-length(states)
  type<-deparse(substitute(l))
  if(is.null(l)){
    l<-switch(type,
              Xsig2=diag(1,ntraits),
              Ysig2=0,
              mu=0)
  }
  if(!is.list(l)){
    l<-list(l)
  }
  if(length(dim(l))<2){
    #check to see if it's in list of lists format!
    l.checks<-unlist(lapply(l,is.list),use.names=FALSE)
    if(any(l.checks)){
      l[!l.checks]<-lapply(l[!l.checks],function(ii) list(ii))
      max.len<-max(lengths(l))
      l<-lapply(l,function(ii) c(ii,rep(list(NULL),max.len-length(ii))))
      l<-do.call(cbind,l)
    }else{
      l<-matrix(l,1,length(l),
                dimnames=list(NULL,names(l)))
      #automatically assume l is supposed to be a column vector when there's only 1 state!
      if(nstates==1){
        l<-t(l)
      }
    }
  }
  nms<-colnames(l)
  n<-ncol(l)
  if(!length(nms)){
    colnames(l)<-rep('',n)
    nms<-colnames(l)
    prob.nms<-rep(TRUE,n)
    has.prob.nms<-TRUE
    quiet.flag<-TRUE
  }else{
    prob.nms<-!(nms%in%states)
    has.prob.nms<-any(prob.nms)
    quiet.flag<-FALSE
  }
  for(i in seq_len(n)){
    nm<-nms[i]
    if(nzchar(nm)){
      report.i<-paste0("'",nm,"'")
    }else{
      report.i<-i
    }
    if(type=='Xsig2'|type=='Ysig2'){
      l[,i]<-.fix.mat(l[,i],function(j) paste0(type,'[[',j,',',report.i,']]'),
                      ntraits,traits,
                      type=='Ysig2')
    }else{
      l[,i]<-.fix.mu(l[,i],function(j) paste0(type,'[[',j,',',report.i,']]'),
                     ntraits,traits)
    }
  }
  out<-matrix(list(),nrow(l),nstates,
              dimnames=list(NULL,states))
  if(any(!prob.nms)){
    out[,nms[!prob.nms]]<-l[,!prob.nms]
  }
  prob.states<-!lengths(out[1,])
  if(any(prob.states)){
    if(has.prob.nms){
      #take elements with no matching names and fill in
      tmp<-l[,prob.nms,drop=FALSE]
      message<-paste0("; missing ",type,"'s recycled from elements of ",type," list")
    }else{
      #take default Xsig2/Ysig2/mu element and fill in
      tmp<-matrix(list(switch(type,
                              Xsig2=matrix(diag(ntraits),ntraits,ntraits,
                                           dimnames=list(traits,traits)),
                              Ysig2=matrix(0,ntraits,ntraits,
                                           dimnames=list(traits,traits)),
                              mu=setNames(rep(0,ntraits),traits))),
                  nrow(l),1)
      message<-paste0("; missing ",type,"'s filled in with ",
                      if(type=='Ysig2'|type=='mu') "0 " else "identity ",
                      if(type=='mu') "vectors" else "matrices")
    }
    out[,prob.states]<-tmp[,rep(seq_len(ncol(tmp)),length.out=sum(prob.states))]
    if(!quiet.flag){
      prefix<-if(type=='Ysig2') c('Tip/Node','Tips/Nodes') else c('State','States')
      warning(.report.names(prefix,states[prob.states],c('has','have')),
              ' no associated elements in ',
              type,message,
              immediate.=TRUE)
    }
  }
  out
}

#one slight issue is that things will get pointlessly divided up for differing values of X0 and nobs, despite them not affecting the main
#preorder traversal calculations...
.get.lookup<-function(treeID,ntrait.data,nXsig2,nYsig2,nmu,nnobs=NULL,nX0=NULL){
  tree.inds<-unique(treeID)
  sub.ntrees<-length(tree.inds)
  lookup<-vector('list',sub.ntrees)
  lens<-c(ntrait.data,nXsig2,nYsig2,nmu,nnobs,nX0)
  len<-length(lens)
  mods<-lapply(lens,function(ii) (seq_len(ii)-1)%%sub.ntrees+1)
  for(i in seq_len(sub.ntrees)){
    inds<-treeID==tree.inds[i]
    tmp.nrow<-sum(inds)
    lookup[[i]][['table']]<-matrix(nrow=tmp.nrow,ncol=len)
    for(j in seq_len(len)){
      lookup[[i]][['table']][,j]<-rep(which(mods[[j]]%%lens[j]==i%%lens[j]),length.out=tmp.nrow)
    }
    if(len>4){
      lookup[[i]][['nobs.X0']]<-lookup[[i]][['table']][,c(5,6),drop=FALSE]
    }
    nms<-do.call(paste,asplit(lookup[[i]][['table']][,seq_len(4),drop=FALSE],2))
    lookup[[i]][['table']]<-lookup[[i]][['table']][!duplicated(nms),,drop=FALSE]
    lookup[[i]][['matches']]<-match(nms,nms)
    lookup[[i]][['inds']]<-do.call(rbind,lapply(seq_len(max(lookup[[i]][['matches']])),function(ii) lookup[[i]][['matches']]==ii))
  }
  lookup
}

.init.arrays<-function(treeID,ntraits,ts,inds,conditional=FALSE){
  nedges<-nrow(ts)
  ntrees<-ncol(ts)
  tree.inds<-unique(treeID)
  sub.ntrees<-length(tree.inds)
  sims.per.tree<-unlist(lapply(tree.inds,function(ii) sum(treeID==ii)),use.names=FALSE)
  nts<-matrix(1,nrow=nedges+1,ncol=sub.ntrees)
  if(conditional){
    t1s<-NTS<-nts
    NTS[-1,]<-lengths(inds[,tree.inds,drop=FALSE])
    t1s[-1,]<-unlist(lapply(inds[,tree.inds,drop=FALSE],'[[',1),use.names=FALSE)
  }
  nts[-1,]<-lengths(ts[,tree.inds,drop=FALSE])
  seeds.per.tree.edge<-ntraits*sweep(nts,2,sims.per.tree,'*')
  tmp.seed<-rnorm(sum(seeds.per.tree.edge))
  seed<-matrix(list(),nedges+1,sub.ntrees,
               dimnames=list('edge'=c('N0',seq_len(nedges)),'tree'=tree.inds))
  if(conditional){
    dv<-dx<-v<-x<-seed
    inf.I<-diag(Inf,nrow=ntraits)
  }
  counter<-0
  for(i in seq_len(nedges+1)){
    for(j in seq_len(sub.ntrees)){
      seed[[i,j]]<-array(tmp.seed[(counter+1):(counter+seeds.per.tree.edge[i,j])],
                         dim=c(nts[i,j],ntraits,sims.per.tree[j]))
      if(conditional){
        x[[i,j]]<-array(0,
                        dim=c(nts[i,j],ntraits,sims.per.tree[j]))
        dx[[i,j]]<-array(0,
                         dim=c(NTS[i,j],ntraits,sims.per.tree[j]))
        dv[[i,j]]<-v[[i,j]]<-array(inf.I,
                                   dim=c(ntraits,ntraits,NTS[i,j],sims.per.tree[j]))
      }
      counter<-counter+seeds.per.tree.edge[i,j]
    }
  }
  out<-list(nts=nts,
            sims.per.tree=sims.per.tree,
            seed=seed)
  if(conditional){
    out[c('NTS','t1s','x','v','dx','dv')]<-list(NTS=NTS,
                                                t1s=t1s,
                                                x=x,
                                                v=v,
                                                dx=dx,
                                                dv=dv)
  }
  out
}

.contsimmap2treeinfo<-function(contsimmap){
  tree<-attr(contsimmap,'tree')
  tree1<-tree[[1]]
  list(tree=tree,
       treeID=.get.treeID(contsimmap),
       ntrees=length(tree),
       nodes=tree1[['node.label']],
       nnodes=tree1[['Nnode']],
       edges=tree1[['edge']],
       nedges=nrow(tree1[['edge']]),
       tips=tree1[['tip.label']],
       ntips=length(tree1[['tip.label']]),
       states=colnames(tree1[['mapped.edge']]),
       nstates=ncol(tree1[['mapped.edge']]))
}

#can I reuse this for prep conthistory?
#looks that way
.prep.contsimmap<-function(tree,nsims,res,
                           trait.data,Xsig2,Ysig2,mu,
                           conditional,
                           ntraits=NULL,traits=NULL,nobs=NULL,
                           diffusion=FALSE){
  
  if(diffusion){
    tree.info<-.contsimmap2treeinfo(tree)
  }else{
    ####INITIAL INFO EXTRACTION####
    #outputs...
    #tree = list of trees in multiSimmap format
    #treeID = which tree each simulation will correspond to
    #ntrees = number of trees
    #nodes = node labels
    #nnodes = number of internal nodes
    #edges = edge matrix
    #nedges = number of edges
    #tips = tip labels
    #ntips = number of tips
    #states = state labels
    #nstates = number of states
    tree.info<-.get.tree.info(tree,nsims)
  }
  #outputs...
  #
  topo.info<-.get.topo.info(tree.info[['tree']][[1]])
  if(diffusion){
    edge.info<-list(ts=attr(tree,'ts'),maps=attr(tree,'maps'))
  }else{
    #outputs...
    #ts = list of final timepoints along each edge (excluding timepoints added to account for discrete state shifts)
    #maps = list including pertinent information on increments along each edge in each tree
    edge.info<-.get.edge.info(tree.info[['tree']],tree.info[['treeID']],res)
  }
  
  ####FORMATTING PARAMETERS####
  #trait.data
  #outputs...
  #
  trait.info<-.get.trait.info(trait.data,tree.info[['tips']],tree.info[['nodes']],tree.info[['edges']],
                              ntraits,traits,nobs,conditional)
  traits<-trait.info[[1]]
  ntraits<-length(traits)
  ntrait.data<-trait.info[[2]]
  #Xsig2
  #takes absurdly long when Xsig2 is just a list of scalars...
  Xsig2<-.fix.param.list(Xsig2,traits,tree.info[['states']])
  #Ysig2
  Ysig2<-.fix.param.list(Ysig2,traits,c(tree.info[['tips']],tree.info[['nodes']]))
  #mu
  mu<-.fix.param.list(mu,traits,tree.info[['states']])
  #make lookup table for parameters
  lookup<-.get.lookup(tree.info[['treeID']],ntrait.data,nrow(Xsig2),nrow(Ysig2),nrow(mu),
                      if(conditional) NULL else ncol(trait.info[['nobs']]),
                      if(conditional) NULL else ncol(trait.info[['X0']]))
  
  ####INITIALIZING OUTPUT ARRAYS####
  #outputs...
  #seed = 
  #nts = 
  #sims.per.tree = 
  #and potentially also...
  #NTS = 
  #t1s = 
  #x = 
  #v = 
  #dx = 
  #dv = 
  arrs<-.init.arrays(tree.info[['treeID']],ntraits,
                     matrix(edge.info[['maps']][,,'ts',drop=FALSE],ncol=tree.info[['ntrees']]),
                     matrix(edge.info[['maps']][,,'inds',drop=FALSE],ncol=tree.info[['ntrees']]),
                     conditional)
  
  c(tree.info[c('tree','treeID')],
    Yperm=list(c(length(tree.info[["tips"]])+1,tree.info[["edges"]][,2])), #so Ysig2 can be rearranged into edgewise format
    topo.info[c('anc','des','prune.seq')], #technically don't need des/ndes for unconditional version...
    ndes=list(lengths(topo.info[['des']])),
    edge.info,
    trait.info,
    Xsig2=list(Xsig2),
    Ysig2=list(Ysig2),
    mu=list(mu),
    lookup=list(lookup),
    arrs)
}

.uncompress<-function(x){
  traits<-names(attr(x,'traits'))
  edges<-rownames(x)
  nsims<-length(attr(x,'treeID'))
  trees<-colnames(x)
  sims.per.tree<-unlist(lapply(x[1,],function(ii) dim(ii)[3]),use.names=FALSE)
  in.treeID<-rep(trees,sims.per.tree)
  out.treeID<-as.character(attr(x,'treeID'))
  for(i in seq_along(trees)){
    inds<-in.treeID==i
    in.treeID[inds]<-paste0('tree',in.treeID[inds],'_sim',seq_len(sum(inds)))
    inds<-out.treeID==i
    out.treeID[inds]<-paste0('tree',out.treeID[inds],'_sim',seq_len(sum(inds)))
  }
  matches<-match(out.treeID,in.treeID)
  out<-aperm(array(unlist(apply(x,1,function(ii) do.call(cbind,lapply(ii,asplit,c(2,3)))),recursive=FALSE),
                   c(length(traits),nsims,length(edges)),list(trait=traits,sim=in.treeID,edge=edges)),
             c(3,1,2))[,,matches,drop=FALSE]
  
  #have to reconfigure lookup
  #have to figure out better way of indicating WHICH lookups need to be merged!
  attr(x,'params')[['lookup',1]]<-do.call(rbind,
                                          lapply(attr(x,'params')[['lookup',1]],
                                                 function(ii) ii[['table']][ii[['matches']],,drop=FALSE])
  )[matches,seq_len(4),drop=FALSE]

  for(i in c('ts','tree','maps','traits','params')){
    attr(out,i)<-attr(x,i)
  }
  class(out)<-'contsimmap'
  out
}