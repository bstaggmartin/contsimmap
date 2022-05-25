.extract.traits<-function(contsimmap,foc.traits){
  ntraits<-ntrait(contsimmap)
  nts<-lapply(contsimmap[['x']],function(ii) lengths(ii)/ntraits)
  edge.seq<-seq_along(nts)
  edge.split.inds<-rep(edge.seq,unlist(lapply(nts,sum),use.names=FALSE))
  tree.seq<-seq_along(contsimmap[['tree']])
  tree.split.inds<-lapply(nts,function(ii) rep(tree.seq,ii))
  trait.ind<-match(foc.traits,contsimmap[['traits']])
  trait.seq<-seq_along(foc.traits)
  prob.traits<-is.na(trait.ind)
  if(any(prob.traits)){
    n.prob.traits<-sum(prob.traits)
    #use subset to extend data to include new traits
    contsimmap<-subset(contsimmap,traits=c(seq_len(ntraits),rep(1,n.prob.traits)))
    tmp.seq<-seq_len(n.prob.traits)+ntraits
    contsimmap[['traits']][tmp.seq]<-foc.traits[prob.traits]
    trait.ind[prob.traits]<-tmp.seq
  }
  nodes<-contsimmap[['nodes']]
  nnodes<-length(nodes)
  node.array<-array(unlist(nodes,use.names=FALSE),
                    c(dim(nodes[[1]]),nnodes),
                    list(NULL,NULL,names(nodes)))
  tmp.nodes<-asplit(node.array[trait.ind,,,drop=FALSE],1)
  out<-setNames(rep(contsimmap['x'],length(foc.traits)),foc.traits)
  for(i in seq_along(foc.traits)){
    for(j in edge.seq){
      for(k in tree.seq){
        out[[i]][[j]][[k]]<-out[[i]][[j]][[k]][,trait.ind[i],,drop=FALSE]
      }
    }
    out[[i]]<-c(tmp.nodes[[i]],unlist(out[[i]],use.names=FALSE))
  }
  init.split<-rep.int(2,length(out[[1]]))
  init.split[seq_along(tmp.nodes[[1]])]<-1
  resplit.and.integrate<-function(contsimmap,foc.traits){
    for(i in trait.seq){
      init<-split(foc.traits[[i]],init.split)
      edges<-split(init[[2]],edge.split.inds)
      nodes<-init[[1]]
      for(j in edge.seq){
        tmp<-split(edges[[j]],tree.split.inds[[j]])
        for(k in tree.seq){
          contsimmap[['x']][[j]][[k]][,trait.ind[i],]<-tmp[[k]]
        }
      }
      node.array[trait.ind[i],,]<-nodes
    }
    contsimmap[['nodes']]<-asplit(node.array,3)
    contsimmap
  }
  out<-c(out,fxn=resplit.and.integrate)
  out
}

#start of function for naming unnamed categories of breaks
#simply fill in unnamed categories with respective entry from here
.get.names<-function(breaks){
  if(is.character(breaks)){
    len<-as.numeric(breaks)+1
  }else{
    len<-length(breaks)+1
  }
  out<-names(breaks)
  if(is.null(out)){
    out<-vector('character',len)
  }else{
    out<-c(out,'')
  }
  prob.nms<-!nzchar(out)|is.na(out)
  if(any(prob.nms)){
    nn<-ceiling(log(len,26))
    tmp<-do.call(expand.grid,rep(list(LETTERS),nn))
    tmp<-as.matrix(tmp)[seq_len(len),,drop=FALSE]
    out[prob.nms]<-do.call(paste,c(rev(asplit(tmp,2)),sep=''))[prob.nms]
  }
  out
}

#need some convention for naming the breaks themselves...see .get.names()
.fix.breaks<-function(breaks,traits){
  ntraits<-length(traits)
  if(!is.list(breaks)){
    breaks<-list(breaks)
  }
  if(is.null(names(breaks))&length(breaks)==ntraits){
    names(breaks)<-traits
  }
  nms<-names(breaks)
  prob.nms<-!(nms%in%traits)
  if(all(prob.nms)){
    stop(report.names(c('Name','Names'),nms[prob.nms]),' did not match to any trait names in contsimmap.')
  }
  if(any(prob.nms)){
    breaks<-breaks[!prob.nms]
    warning(report.names(c('Name','Names'),nms[prob.nms]),' did not match to any trait names in contsimmap: these breaks were ignored.')
  }
  for(i in seq_along(breaks)){
    
    attr(breaks[[i]],'nms')<-.get.names(breaks[[i]])
  }
  breaks
}

#enter breaks as a list of break points for a given trait
#this will only threshold the trait
#if no breaks are specified, defaults to using nbreaks (input as named vector)
#if an nbreaks is specified for a name in breaks, nbreaks will replace corresponding entry in breaks
#' @export
thresh.traits<-function(contsimmap,breaks=NULL,nbreaks=10,...){
  traits<-contsimmap[['traits']]
  if(is.null(breaks)){
    nbreaks<-setNames(rep(nbreaks,length.out=length(traits)),traits)
  }
  if(!is.null(names(nbreaks))){
    prob.nms<-!(names(nbreaks)%in%traits)
    if(any(prob.nms)){
      warning('PLACEHOLDER WARNING ABOUT NBREAKS WITH UNRECOGNIZED NAMES')
    }
    nbreaks<-nbreaks[!prob.nms]
    overwrite.inds<-match(names(nbreaks),names(breaks))
    overwrites<-!is.na(overwrite.inds)
    overwrite.inds<-overwrite.inds[overwrites]
    nbreaks<-setNames(as.list(as.character(nbreaks)),names(nbreaks))
    if(sum(overwrites)){
      breaks[overwrite.inds]<-nbreaks[overwrites]
      warning('PLACEHOLDER WARNING ABOUT OVERWRITING BREAKS WITH NBREAKS')
    }
    breaks<-c(breaks,nbreaks[!overwrites])
  }
  breaks<-.fix.breaks(breaks,traits)
  foc.traits<-.extract.traits(contsimmap,names(breaks))
  resplit.and.integrate<-foc.traits[['fxn']]
  breaks.seq<-seq_along(breaks)
  foc.traits<-foc.traits[breaks.seq]
  #may need to do something about breaks that imply intervals of 0 length--could really screw things up
  for(i in breaks.seq){
    if(is.character(breaks[[i]])){
      nms<-attr(breaks[[i]],'nms')
      ran<-range(c(foc.traits[[i]]),na.rm=TRUE)
      nn<-as.numeric(breaks[[i]])
      breaks[[i]]<-seq(ran[1],ran[2],length.out=nn+2)[-c(1,nn+2)]
      attr(breaks[[i]],'nms')<-nms
    }
    foc.traits[[i]]<-findInterval(foc.traits[[i]],breaks[[i]])+1
  }
  contsimmap<-resplit.and.integrate(contsimmap,foc.traits)
  #add cut-points/names as a separate attribute...
  contsimmap[['keys']]<-breaks
  contsimmap
}

#!!!needs to be fixed to account for cases when multiple thresholds have been crossed in a single time-step!!!!
#this function will wrap around thresh.traits and convert it to a phytools simmap object
#' @export
as.simmap<-function(contsimmap,breaks=NULL,nbreaks=10,...){
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(contsimmap[['nodes']])==1){
      #convert tree object to single subtree
    }else{
      stop('The as.simmap() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  thresh<-thresh.traits(contsimmap,breaks,nbreaks)
  breaks<-thresh[['keys']]
  keys<-lapply(breaks,function(ii) attr(ii,'nms'))
  thresh.nodes<-thresh[['nodes']]
  thresh<-thresh[['x']]
  out<-contsimmap[['tree']][[1]]
  ehgt<-evorates:::edge.ranges(out)
  elen<-out[['edge.length']]
  out[['maps']]<-vector('list',nrow(out[['edge']]))
  out[['node.states']]<-out[['states']]<-out[['mapped.edge']]<-NULL
  out<-matrix(list(out),length(breaks),nsim(contsimmap),
              dimnames=list('trait'=names(breaks),'sim'=NULL))
  edge.inds<-as.numeric(names(contsimmap[['x']]))
  tree.seq<-seq_along(contsimmap[['tree']])
  out.inds<-contsimmap[['perm.inds']][,1]
  out.inds<-split(seq_along(out.inds),out.inds)
  traits<-contsimmap[['traits']]
  trait.inds<-match(names(breaks),traits)
  trait.seq<-seq_along(trait.inds)
  for(i in seq_along(contsimmap[['x']])){
    edge<-contsimmap[['x']][[i]]
    maps<-contsimmap[['maps']][[i]]
    t0<-ehgt[i,1]
    tend<-ehgt[i,2]
    tlen<-elen[i]
    edge.thresh<-thresh[[i]]
    anc<-contsimmap[['anc']][i]
    anc.edge<-contsimmap[['x']][[anc]]
    if(is.null(anc.edge)){
      anc.node<-TRUE
      anc.edge<-contsimmap[['nodes']][[anc]]
      anc.edge.thresh<-thresh.nodes[[anc]]
    }else{
      anc.node<-FALSE
      anc.edge.thresh<-thresh[[anc]]
    }
    for(j in tree.seq){
      sub.edge<-edge[[j]]
      dims<-dim(sub.edge)
      len<-dims[1]+1
      sub.seq<-seq_len(dims[3])
      tmp.out.inds<-out.inds[[j]]
      ts<-c(t0,maps[[j]][['ts']])
      sub.edge.thresh<-edge.thresh[[j]]
      if(anc.node){
        sub.anc<-anc.edge[,tmp.out.inds,drop=FALSE]
        sub.anc.thresh<-anc.edge.thresh[,tmp.out.inds,drop=FALSE]
      }else{
        sub.anc<-anc.edge[[j]]
        tmp.dims<-dim(sub.anc)
        sub.anc<-matrix(sub.anc[tmp.dims[1],,],
                        tmp.dims[2],tmp.dims[3])
        sub.anc.thresh<-matrix(anc.edge.thresh[[j]][tmp.dims[1],,],
                               tmp.dims[2],tmp.dims[3])
      }
      for(k in trait.seq){
        trait.ind<-trait.inds[k]
        tmp<-matrix(sub.edge[,trait.ind,],dims[1],dims[3])
        tmp<-rbind(sub.anc[trait.ind,],tmp)
        tmp.thresh<-matrix(sub.edge.thresh[,trait.ind,],dims[1],dims[3])
        tmp.thresh<-rbind(sub.anc.thresh[trait.ind,],tmp.thresh)
        starts<-keys[[k]][tmp.thresh[1,]]
        shifts<-lapply(asplit(tmp.thresh,2),function(ii) which(ii[-len]!=ii[-1]))
        nshifts<-lengths(shifts)
        no.shifts<-!nshifts
        if(any(!no.shifts)){
          tmp.inds<-cbind(unlist(shifts,use.names=FALSE),rep(sub.seq,nshifts))
          x1<-tmp[tmp.inds]
          thresh1<-tmp.thresh[tmp.inds]
          t1<-ts[tmp.inds[,1]]
          tmp.inds[,1]<-tmp.inds[,1]+1
          x2<-tmp[tmp.inds]
          thresh2<-tmp.thresh[tmp.inds]
          t2<-ts[tmp.inds[,1]]
          polarity<-(thresh2>thresh1)-1 #0 for positive crosses, -1 for negative crosses
          thresholds<-breaks[[k]][(thresh2-polarity)-1] #I think this works...
          shiftpts<-split(setNames((thresholds-x1)/(x2-x1)*(t2-t1)+t1,keys[[k]][thresh2]),tmp.inds[,2])
          tmp.starts<-starts[!no.shifts]
          for(l in seq_along(shiftpts)){
            tmp.pts<-shiftpts[[l]]
            shiftpts[[l]]<- -c(setNames(t0,tmp.starts[l]),tmp.pts)+c(tmp.pts,tend)
          }
          shifts[!no.shifts]<-shiftpts
        }
        shifts[no.shifts]<-lapply(starts[no.shifts],function(ii) setNames(tlen,ii))
        for(l in sub.seq){
          out[[k,tmp.out.inds[l]]][['maps']][[i]]<-shifts[[l]]
        }
      }
    }
  }
  rownms<-out[[1]][['edge']]
  rownms<-paste(rownms[,1],rownms[,2],sep=',')
  nodes<-t(out[[1]][['edge']])
  tips<-out[[1]][['tip.label']]
  ntips<-length(tips)
  tips<-setNames(vector('character',ntips),tips)
  tip.inds<-which(nodes[2,]<=ntips)
  tip.inds<-tip.inds[order(nodes[2,tip.inds])]
  tmp.seq<-seq_len(ncol(out))
  for(i in trait.seq){
    tmp.key<-keys[[i]]
    holder<-matrix(nrow=length(tmp.key),ncol=length(elen),
                   dimnames=list('state'=tmp.key,rownms))
    foo<-function(simmap){
      holder[]<-unlist(lapply(simmap[['maps']],function(ii) 
        unlist(lapply(tmp.key,function(jj) 
          sum(ii[names(ii)==jj])),use.names=FALSE)),
        use.names=FALSE)
      t(holder)
    }
    for(j in tmp.seq){
      nodes[]<-unlist(lapply(out[[i,j]][['maps']],
                             function(ii) names(ii)[c(1,length(ii))]))
      tips[]<-nodes[2,tip.inds]
      out[[i,j]][c('mapped.edge','node.states','states')]<-list(foo(out[[i,j]]),t(nodes),tips)
    }
  }
  out<-setNames(asplit(out,1),names(breaks))
  for(i in trait.seq){
    class(out[[i]])<-c('multiSimmap','multiPhylo')
  }
  out
}

#it may be best to turn thresh.trait into a wrapper for findInterval that can be plugged into a formula for the function below...
#yeah, I need a thresh() function, but it comes with a couple difficulties
  #needs to precompute and cache info regarding break naming
  #needs to compute breaks for nbreaks argument from BOTH foc.traits and nodes...
  #I think I need to combine nodes with foc.traits after all...
#for making new traits given arbitrary formula of existing traits
#' @export
make.traits<-function(contsimmap,...){
  formulae<-list(...)
  out.traits<-lapply(formulae,function(ii) as.character(ii[[2]]))
  calls<-lapply(formulae,function(ii) ii[[3]])
  #maybe don't check for variable names, just let calls run their course and throw an error if they don't find the right variable
  #make sure to have some warning about how formulae are applied in sequence regarding new traits created in same function, replacements, etc...
  #actually not sure what the priority is regarding eval() with multiple variables named the same thing...
  # lapply(calls,all.vars)
  foc.traits<-.extract.traits(contsimmap,out.traits)
  resplit.and.integrate<-foc.traits[['fxn']]
  form.seq<-seq_along(formulae)
  foc.traits<-foc.traits[form.seq]
  for(i in form.seq){
    foc.traits[[i]]<-eval(calls[[i]],envir=as.environment(foc.traits))
  }
  contsimmap<-resplit.and.integrate(contsimmap,foc.traits)
  #check if any foc.traits have a 'keys' attribute and add to contsimmap
  keys<-lapply(foc.traits,function(ii) attr(ii,'keys'))
  keys<-keys[lengths(keys)>0]
  contsimmap[['keys']]<-keys
  contsimmap
}

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
  attr(out,'breaks')<-unname(breaks)
  out
}
