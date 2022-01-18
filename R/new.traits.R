#helper function to merge lists in contsimmap by trait
#then automatically makes a function that resplits and integrates those merged lists back into contsimmap object
.extract.traits<-function(contsimmap,foc.traits,lhs){
  ntraits<-ntrait(contsimmap)
  nts<-lapply(contsimmap[['x']],function(ii) lengths(ii)/ntraits)
  edge.seq<-seq_along(nts)
  edge.split.inds<-rep(edge.seq,unlist(lapply(nts,sum),use.names=FALSE))
  tree.seq<-seq_along(contsimmap[['tree']])
  tree.split.inds<-lapply(nts,function(ii) rep(tree.seq,ii))
  trait.seq<-seq_along(foc.traits)
  trait.ind<-match(foc.traits,contsimmap[['traits']])
  prob.traits<-is.na(trait.ind)
  if(any(prob.traits)){
    new.traits<-unique(foc.traits[prob.traits])
    n.new.traits<-length(new.traits)
    #use subset to extend data to include new traits
    contsimmap<-subset(contsimmap,traits=c(seq_len(ntraits),rep(NA,n.new.traits)))
    contsimmap[['traits']][seq_len(n.new.traits)+ntraits]<-new.traits
    trait.ind<-match(foc.traits,contsimmap[['traits']])
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
  resplit.and.integrate<-function(foc.traits){
    for(i in trait.seq[lhs]){
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

#for making new traits given arbitrary formula of existing traits
#' @export
make.traits<-function(contsimmap,...){
  formulae<-list(...)
  out.traits<-unlist(lapply(formulae,function(ii) as.character(ii[[2]])),use.names=FALSE)
  calls<-lapply(formulae,function(ii) ii[[3]])
  other.traits<-unlist(lapply(calls,all.vars),use.names=FALSE)
  foc.traits<-c(out.traits,other.traits)
  lhs<-rep(c(TRUE,FALSE),c(length(out.traits),length(other.traits)))
  #maybe don't check for variable names, just let calls run their course and throw an error if they don't find the right variable
  #make sure to have some warning about how formulae are applied in sequence regarding new traits created in same function, replacements, etc...
  #actually not sure what the priority is regarding eval() with multiple variables named the same thing...
  # lapply(calls,all.vars)
  foc.traits<-.extract.traits(contsimmap,foc.traits,lhs)
  resplit.and.integrate<-foc.traits[['fxn']]
  foc.traits<-foc.traits[-length(foc.traits)]
  for(i in seq_along(formulae)){
    foc.traits[[i]]<-with(foc.traits,eval(calls[[i]]))
  }
  contsimmap<-resplit.and.integrate(foc.traits)
  #check if any foc.traits have a 'keys' attribute and add to contsimmap
  keys<-lapply(foc.traits,function(ii) attr(ii,'keys'))
  keys<-keys[lengths(keys)>0]
  if(length(keys)){
    matches<-match(names(keys),names(contsimmap[['keys']]))
    new<-is.na(matches)
    contsimmap[['keys']][matches[!new]]<-keys[!new]
    contsimmap[['keys']]<-c(contsimmap[['keys']],keys[new])
  }
  contsimmap
}

#convert contSimmap to simmap via thresholding
#for simplicity, let's make it only do 1 trait at a time for now...
#... is extra arguments to pass to threshold()
#now should work, even when multiple thresholds are crossed in a single time interval
#' @export
as.simmap<-function(contsimmap,trait=1,...){
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(contsimmap[['nodes']])==1){
      stop('The as.simmap() function will support single subtrees in the future, but not yet')
    }else{
      stop('The as.simmap() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  #make appropriate formula
  trait<-trait[1]
  in.trait<-trait
  if(is.character(trait)){
    trait<-pmatch(trait,contsimmap[['traits']])
  }
  trait.ind<-trait
  trait<-contsimmap[['traits']][trait]
  if(is.na(trait)){
    stop('PLACEHOLDER:',in.trait,'did not match up with any available traits in contsimmap')
  }
  #threshold trait
  ex.args<-paste(paste(names(list(...)),list(...),sep='='),collapse=',')
  form<-as.formula(paste0('thresH~threshold(',trait,',',ex.args,')'))
  contsimmap<-make.traits(contsimmap,form)
  trait.ind<-c(trait.ind,which(contsimmap[['traits']]=='thresH'))
  breaks<-contsimmap[['keys']][['thresH']]
  keys<-attr(breaks,'state_names')
  #lots of initialization for main loop
  out<-contsimmap[['tree']][[1]]
  class(out)<-c('simmap','phylo')
  ehgt<-edge.ranges(out)
  elen<-out[['edge.length']]
  out[['maps']]<-vector('list',nrow(out[['edge']]))
  out[['node.states']]<-out[['states']]<-out[['mapped.edge']]<-NULL
  nsims<-nsim(contsimmap)
  out<-rep(list(out),nsims)
  edge.inds<-as.numeric(names(contsimmap[['x']]))
  tree.seq<-seq_along(contsimmap[['tree']])
  perm.inds<-contsimmap[['perm.inds']]
  out.inds<-perm.inds[,1]
  out.inds<-split(seq_along(out.inds),out.inds)
  sims.per.tree<-lengths(out.inds)
  sim.seqs<-lapply(sims.per.tree,seq_len)
  perm.inds<-perm.inds[,3]
  holder<-matrix(nrow=2,ncol=nsims)
  #main loop: preparing 'maps' elements
  for(i in seq_along(contsimmap[['x']])){
    edge<-contsimmap[['x']][[i]]
    maps<-contsimmap[['maps']][[i]]
    t0<-ehgt[i,1]
    tend<-ehgt[i,2]
    tlen<-elen[i]
    anc.ind<-contsimmap[['anc']][i]
    anc<-contsimmap[['x']][[anc.ind]]
    if(is.null(anc)){
      anc<-contsimmap[['nodes']][[anc.ind]][trait.ind,,drop=FALSE]
    }else{
      holder[,perm.inds]<-unlist(lapply(anc,function(ii) ii[dim(ii)[1],trait.ind,,drop=FALSE]),use.names=FALSE)
      anc<-holder
    }
    for(j in tree.seq){
      ts<-c(t0,maps[[j]][['ts']])
      nts<-length(ts)
      tmp.out.inds<-out.inds[[j]]
      sims<-sims.per.tree[j]
      sim.seq<-sim.seqs[[j]]
      tx<-xx<-matrix(nrow=nts,ncol=sims)
      tmp.edge<-edge[[j]][,trait.ind,,drop=FALSE]
      xx[1,]<-anc[1,tmp.out.inds,drop=FALSE]
      xx[-1,]<-tmp.edge[,1,,drop=FALSE]
      tx[1,]<-anc[2,tmp.out.inds,drop=FALSE]
      tx[-1,]<-tmp.edge[,2,,drop=FALSE]
      starts<-keys[tx[1,]]
      shifts<-lapply(asplit(tx,2),function(ii) which(ii[-nts]!=ii[-1]))
      nshifts<-lengths(shifts)
      no.shifts<-!nshifts
      if(any(!no.shifts)){
        tmp.inds<-cbind(unlist(shifts,use.names=FALSE),rep(sim.seq,nshifts))
        xx1<-xx[tmp.inds]
        tx1<-tx[tmp.inds]
        t1<-ts[tmp.inds[,1]]
        tmp.inds[,1]<-tmp.inds[,1]+1
        xx2<-xx[tmp.inds]
        tx2<-tx[tmp.inds]
        t2<-ts[tmp.inds[,1]]
        cross<-tx2-tx1
        polarity<-sign(cross)
        ncross<-abs(cross)
        nncross<-ncross>1
        #handle multiple crosses in single timestep...
        if(any(nncross)){
          xx1<-rep(xx1,ncross)
          t1<-rep(t1,ncross)
          xx2<-rep(xx2,ncross)
          t2<-rep(t2,ncross)
          tx1<-rep(tx1,ncross)
          ii<-which(nncross)
          tmp.seqs<-unlist(lapply(ii,function(ii) polarity[ii]*seq_len(ncross[ii]-1)),use.names=FALSE)
          ii<-rep(ii,ncross[nncross]-1)
          ii<-ii+seq_along(ii)
          tx1[ii]<-tx1[ii]+tmp.seqs
          polarity<-rep(polarity,ncross)
          tmp.inds<-rep(tmp.inds[,2],ncross)
          #rep'ing tx2 is unnecessary all needed info is in tx1 at this point
        }else{
          tmp.inds<-tmp.inds[,2]
        }
        thresholds<-breaks[tx1-(polarity<0)]
        shiftpts<-split(setNames((thresholds-xx1)/(xx2-xx1)*(t2-t1)+t1,keys[tx1+polarity]),tmp.inds)
        tmp.starts<-starts[!no.shifts]
        for(k in seq_along(shiftpts)){
          tmp.pts<-shiftpts[[k]]
          shiftpts[[k]]<- -c(setNames(t0,tmp.starts[k]),tmp.pts)+c(tmp.pts,tend)
        }
        shifts[!no.shifts]<-shiftpts
      }
      shifts[no.shifts]<-lapply(starts[no.shifts],function(ii) setNames(tlen,ii))
      for(k in sim.seq){
        out[[tmp.out.inds[k]]][['maps']][[i]]<-shifts[[k]]
      }
    }
  }
  #get extra simmap elements --> mapped.edge, node.states, states
  #also added breaks so folks could keep track of break points...
  edge<-out[[1]][['edge']]
  nodes<-t(edge)
  tips<-out[[1]][['tip.label']]
  ntips<-length(tips)
  tips<-setNames(vector('character',ntips),tips)
  tip.inds<-which(nodes[2,]<=ntips)
  tip.inds<-tip.inds[order(nodes[2,tip.inds])]
  holder<-matrix(nrow=length(keys),ncol=length(elen),
                 dimnames=list('state'=keys,
                               paste(edge[,1],edge[,2],sep=',')))
  foo<-function(simmap){
    holder[]<-unlist(lapply(simmap[['maps']],function(ii) 
      unlist(lapply(keys,function(jj) 
        sum(ii[names(ii)==jj])),use.names=FALSE)),
      use.names=FALSE)
    t(holder)
  }
  for(i in seq_len(nsims)){
    nodes[]<-unlist(lapply(out[[i]][['maps']],
                           function(ii) names(ii)[c(1,length(ii))]))
    tips[]<-nodes[2,tip.inds]
    out[[i]][c('mapped.edge','node.states','states','breaks')]<-list(foo(out[[i]]),t(nodes),tips,breaks)
  }
  class(out)<-c('multiSimmap','multiPhylo')
  out
}