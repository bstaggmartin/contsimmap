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
  out[['node.states']]<-out[['states']]<-out[['mapped.edge']]<-tree[['Q']]<-tree[['logL']]<-NULL
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
                           function(ii) names(ii)[c(1,length(ii))]),use.names=FALSE)
    tips[]<-nodes[2,tip.inds]
    out[[i]][c('mapped.edge','node.states','states','breaks')]<-list(foo(out[[i]]),t(nodes),tips,breaks)
  }
  class(out)<-c('multiSimmap','multiPhylo')
  out
}