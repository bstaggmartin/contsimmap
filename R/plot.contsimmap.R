.parse.bys<-function(choices,
                     Col.by,Alpha.by,Mix.by,Wgt.by,Lty.by,Lwd.by){
  #determine bys
  bys<-list(col=Col.by,alpha=Alpha.by,mix=Mix.by,wgt=Wgt.by,lty=Lty.by,lwd=Lwd.by)
  trait.flag<-setNames(rep(FALSE,6),names(bys))
  for(i in seq_len(6)){
    if(is.null(bys[[i]])){
      if(i==1){
        bys[[i]]<-'simulation'
      }else if(i==3){
        bys[[i]]<-bys[[i-2]]
      }else if(i==5){
        bys[[i]]<-bys[[i-4]]
      }else{
        bys[[i]]<-bys[[i-1]]
      }
    }else{
      bys[[i]]<-pmatch(bys[[i]][1],choices)
      if(is.na(bys[[i]])){
        if(i==1){
          bys[[i]]<-'simulation'
        }else if(i==3){
          bys[[i]]<-bys[[i-2]]
        }else if(i==5){
          bys[[i]]<-bys[[i-4]]
        }else{
          bys[[i]]<-bys[[i-1]]
        }
        warning(paste0(names(bys)[i],'.by'),' should be one of the following:',.report.names(nms=choices,combine='or'),
                ': it was set to ',choices[bys[[i]]],' by default')
      }
      if(bys[[i]]>4){
        trait.flag[i]<-TRUE
      }
    }
  }
  bys<-unlist(bys)
  bys<-setNames(choices[bys],names(bys))
  attr(bys,'trait.flag')<-trait.flag
  bys
}

.parse.contsimmap<-function(contsimmap,edges,sims,traits,bys,join.lines){
  #put focal traits into z1 (and optionally z0) caches for easy retrieval
  #grabs time and states as necessary too
  #only puts things into z0 caches if necessary because the quantity is being plotted along x or y
  #might want to generalize to >2 dimensions, eventually
  #some redundancy if a trait is specified for a plotting dimension AND in a by argument, but overall pretty good...
  if(length(traits)>1){
    doub.trait<-TRUE
    traits<-traits[c(1,2)]
  }else{
    doub.trait<-FALSE
    traits<-c(NA,traits)
  }
  if(is.numeric(traits)){
    #may want to put some more checks here...
    traits<-contsimmap[['traits']][traits]
  }
  extra.traits<-unique(bys[attr(bys,'trait.flag')])
  contsimmap<-subset(contsimmap,edges=edges,sims=sims,traits=c(traits,extra.traits))
  if(!doub.trait){
    hgts<-edge.ranges(contsimmap[['tree']][[1]])
    for(i in names(contsimmap[['nodes']])){
      contsimmap[['nodes']][[i]][1,]<-if(i=='0') 0 else hgts[as.numeric(i),2]
    }
  }
  trait.nms<-c('@xXXx@','@yYYy@',contsimmap[['traits']][-c(1,2)])
  ntraits<-length(trait.nms)
  contsimmap[['z1']]<-setNames(rep(contsimmap['x'],ntraits),trait.nms)
  if(!join.lines){
    contsimmap[['z0']]<-contsimmap[['z1']][c(1,2)]
  }
  ####THIS CHUNK OF CODE IS PERHAPS GETTING REPETITIVE####
  ntrees<-length(contsimmap[['tree']])
  nedges<-length(contsimmap[['x']])
  nts<-matrix(lengths(unlist(contsimmap[['x']],recursive=FALSE,use.names=FALSE)),
              ntrees,nedges)
  nts<-nts/ntrait(contsimmap)
  denoms<-unlist(lapply(contsimmap[['x']][[1]],function(ii) dim(ii)[3]),use.names=FALSE)
  nts<-sweep(nts,1,denoms,'/')
  if(join.lines){
    nts<-nts+2
  }
  contsimmap[['nts']]<-nts
  ########
  edge.seq<-seq_len(nedges)
  tree.seq<-seq_len(ntrees)
  out.inds<-contsimmap[['perm.inds']][,1]
  out.inds<-split(seq_along(out.inds),out.inds)
  sims.per.tree<-lengths(out.inds)
  trait.seq<-seq_len(ntraits)
  #might be easier to have separate join.lines and not loops...
  for(i in edge.seq){
    edge<-contsimmap[['x']][[i]]
    if(!doub.trait){
      maps<-contsimmap[['maps']][[i]]
    }
    node.flag<-FALSE
    anc.ind<-contsimmap[['anc']][i]
    anc<-contsimmap[['x']][[anc.ind]]
    if(is.null(anc)){
      node.flag<-TRUE
      anc<-contsimmap[['nodes']][[anc.ind]]
    }
    anc.ind<-as.numeric(anc.ind)
    for(j in tree.seq){
      foc.edge<-edge[[j]]
      tmp.out<-out.inds[[j]]
      holder<-matrix(nrow=nts[j,i],ncol=sims.per.tree[j])
      if(join.lines){
        tmp.inds<-c(1,seq_len(nts[j,i]-2),nts[j,i]-2)
      }
      if(node.flag){
        foc.anc<-anc[,tmp.out,drop=FALSE]
      }else{
        if(join.lines){
          foc.anc<-matrix(anc[[j]][nts[j,anc.ind]-2,,,drop=FALSE],k,sims.per.tree[j])
        }else{
          foc.anc<-matrix(anc[[j]][nts[j,anc.ind],,,drop=FALSE],k,sims.per.tree[j])
        }
      }
      for(k in trait.seq){
        if(k==1&!doub.trait){
          if(join.lines){
            holder[]<-maps[[j]][['ts']][tmp.inds]
          }else{
            holder[]<-maps[[j]][['ts']]
          }
        }else{
          if(join.lines){
            holder[]<-foc.edge[tmp.inds,k,,drop=FALSE]
          }else{
            holder[]<-foc.edge[,k,,drop=FALSE]
          }
        }
        contsimmap[['z1']][[k]][[i]][[j]]<-holder
        if(k<3){
          if(k==1&!doub.trait){
            tmp.foc.anc<-hgts[i,1,drop=FALSE]
          }else{
            tmp.foc.anc<-foc.anc[k,,drop=FALSE]
          }
          if(join.lines){
            contsimmap[['z1']][[k]][[i]][[j]][1,]<-tmp.foc.anc
            contsimmap[['z1']][[k]][[i]][[j]][nts[j,i],]<-NA
          }else{
            holder[-1,]<-holder[-nts[j,i],,drop=FALSE]
            holder[1,]<-tmp.foc.anc
            contsimmap[['z0']][[k]][[i]][[j]]<-holder
          }
        }
      }
    }
  }
  #grab any time/state stuff...
  time.flag<-FALSE
  if(any(bys=='time')){
    contsimmap[['z1']][['@tTTt@']]<-contsimmap[['z1']][['@xXXx@']]
    if(doub.trait|join.lines){
      time.flag<-TRUE
    }
  }
  state.flag<-FALSE
  if(any(bys=='state')){
    contsimmap[['z1']][['@sSSs@']]<-contsimmap[['z1']][['@xXXx@']]
    state.flag<-TRUE
  }
  if(time.flag|state.flag){
    for(i in edge.seq){
      maps<-contsimmap[['maps']][[i]]
      for(j in tree.seq){
        foc<-maps[[j]]
        if(join.lines){
          tmp.inds<-c(1,seq_len(nts[j,i]-2),nts[j,i]-2)
        }
        if(time.flag){
          if(join.lines){
            contsimmap[['z1']][['@tTTt@']][[i]][[j]][]<-foc[['ts']][tmp.inds]
          }else{
            contsimmap[['z1']][['@tTTt@']][[i]][[j]][]<-foc[['ts']]
          }
        }
        if(state.flag){
          if(join.lines){
            contsimmap[['z1']][['@sSSs@']][[i]][[j]][]<-foc[['state']][tmp.inds]
          }else{
            contsimmap[['z1']][['@sSSs@']][[i]][[j]][]<-foc[['state']]
          }
        }
      }
    }
  }
  #use nodes to finalize ranges as it comes up in .parse.args...
  contsimmap[['x']]<-NULL
  if(!doub.trait){
    contsimmap[['traits']][1]<-'time'
  }
  contsimmap
}

.parse.args<-function(contsimmap,args.ls,bys,scale.res,join.lines){
  if(join.lines){
    args.ls[c('x','y')]<-
      lapply(contsimmap[['z1']][c(1,2)],unlist,use.names=FALSE)
  }else{
    args.ls[c('x0','y0','x1','y1')]<-
      lapply(c(contsimmap[['z0']][c(1,2)],contsimmap[['z1']][c(1,2)]),unlist,use.names=FALSE)
  }
  #set defaults
  defaults<-list('col'=palette(),
                 'alpha'=NA,
                 'mix'=NA,
                 'wgt'=0.5,
                 'lty'=par('lty'),
                 'lwd'=par('lwd'),
                 'xlab'=contsimmap[['traits']][1],
                 'ylab'=contsimmap[['traits']][2])
  for(i in names(defaults)){
    if(is.null(args.ls[[i]])){
      args.ls[[i]]<-defaults[[i]]
    }
  }
  #should I put nodes together in a single array first?
  if(is.null(args.ls[['xlim']])){
    tmp<-unlist(lapply(contsimmap[['nodes']],function(ii) ii[1,,drop=FALSE]),use.names=FALSE)
    if(join.lines){
      args.ls[['xlim']]<-range(tmp,args.ls[['x']],na.rm=TRUE)
    }else{
      args.ls[['xlim']]<-range(tmp,args.ls[['x1']],na.rm=TRUE)
    }
  }
  if(is.null(args.ls[['ylim']])){
    tmp<-unlist(lapply(contsimmap[['nodes']],function(ii) ii[2,,drop=FALSE]),use.names=FALSE)
    if(join.lines){
      args.ls[['ylim']]<-range(tmp,args.ls[['y']],na.rm=TRUE)
    }else{
      args.ls[['ylim']]<-range(tmp,args.ls[['y1']],na.rm=TRUE)
    }
  }
  #form final argument list/recycling
  #could code in some "shortcuts" if no alpha or mix are specified...
  #could also save time/space if any arg is of length 1!
  inds<-vector('list')
  nts<-contsimmap[['nts']]
  nedges<-length(contsimmap[['z1']][[1]])
  ntrees<-length(contsimmap[['z1']][[1]][[1]])
  sims.per.tree<-unlist(lapply(contsimmap[['z1']][[1]][[1]],ncol),use.names=FALSE)
  nsims<-sum(sims.per.tree)
  has.alpha<-TRUE
  has.mix<-TRUE
  for(i in names(bys)){
    if(i=='alpha'){
      if(all(is.na(args.ls[[i]]))){
        args.ls[[i]]<-NA
        has.alpha<-FALSE
      }
    }
    if(i=='mix'){
      if(all(is.na(args.ls[[i]]))){
        args.ls[[i]]<-NA
        has.mix<-FALSE
      }
    }
    if(i=='wgt'){
      if(all(!args.ls[[i]])){
        args.ls[[i]]<-0
        has.mix<-FALSE
      }
    }
    if(length(args.ls[[i]])>1){
      foc<-bys[i]
      if(foc=='simulation'){
        if(is.null(inds[['simulation']])){
          tmp.nts<-nts[rep.int(seq_len(ntrees),sims.per.tree),]
          inds[['simulation']]<-rep.int(rep.int(contsimmap[['perm.inds']][,3],nedges),tmp.nts)
        }
        args.ls[[i]]<-rep(args.ls[[i]],length.out=nsims)[inds[['simulation']]]
      }else if(foc=='state'){
        if(is.null(inds[['state']])){
          tmp<-unlist(contsimmap[['z1']][['@sSSs@']],use.names=FALSE)
          attr(tmp,'states')<-sort(unique(tmp))
          inds[['state']]<-tmp
        }
        states<-attr(inds[['state']],'states')
        nms<-names(args.ls[[i]])
        if(is.null(nms)){
          args.ls[[i]]<-rep(args.ls[[i]],length.out=length(states))
          names(args.ls[[i]])<-states
        }else{
          tmp<-args.ls[[i]][states]
          prob.nms<-!(nms%in%states)
          nas<-is.na(args.ls[[i]])
          tmp[nas]<-if(any(prob.nms)) args.ls[[i]][prob.nms] else NA
          args.ls[[i]]<-tmp
        }
        args.ls[[i]]<-args.ls[[i]][inds[['state']]]
      }else if(foc=='edge'){
        if(is.null(inds[['edge']])){
          tmp.nts<-.colSums(nts*sims.per.tree,ntrees,nedges)
          inds[['edge']]<-rep.int(seq_len(nedges),tmp.nts)
        }
        args.ls[[i]]<-rep(args.ls[[i]],length.out=nedges)[inds[['edge']]]
      }else{
        if(is.null(inds[[foc]])){
          if(foc=='time'){
            tmp<-unlist(contsimmap[['z1']][['@tTTt@']],use.names=FALSE)
            lim<-c(0,max(tmp,na.rm=TRUE))
            tmp<-tmp/lim[2]*(scale.res-1)+1
            attr(tmp,'lim')<-lim
            inds[['time']]<-round(tmp)
          }else{
            tmp<-unlist(contsimmap[['z1']][[foc]],use.names=FALSE)
            tmp.ind<-match(foc,contsimmap[['traits']])
            tmp.node<-unlist(lapply(contsimmap[['nodes']],function(ii) ii[tmp.ind,,drop=FALSE]),use.names=FALSE)
            lim<-range(tmp.node,tmp,na.rm=TRUE)
            tmp<-(tmp-lim[1])/(lim[2]-lim[1])*(scale.res-1)+1
            attr(tmp,'lim')<-lim
            inds[[foc]]<-round(tmp)
          }
        }
        if(i=='col'){
          args.ls[[i]]<-colorRampPalette(args.ls[[i]],alpha=TRUE)(scale.res)
        }else{
          if(is.numeric(args.ls[[i]])){
            args.ls[[i]]<-approx(args.ls[[i]],n=scale.res)$y
          }else{
            args.ls[[i]]<-args.ls[[i]][round(approx(seq_along(args.ls[[i]]),n=scale.res)$y)]
          }
        }
        args.ls[[i]]<-args.ls[[i]][inds[[foc]]]
      }
    }
  }
  if(has.alpha|has.mix){
    args.ls[['col']]<-alter.cols(args.ls[['col']],args.ls[['alpha']],args.ls[['mix']],args.ls[['wgt']])
  }
  args.ls[c('alpha','mix','wgt')]<-NULL
  args.ls
}

#instead of using segments, it might instead be better to break up things into passes to lines() to prevent
##overlapping segments when partial transparency is used...
##using a mix of split() and NAs to break up lines by arguments...
##would just require a simple alteration to .parse.contsimmap where there is no z0--instead, z1's just have starting
###points "grafted on"
#a bit more tricky than anticipated, but done--join.lines works (and is honestly faster and probs should be the default...)
#okay, join.lines is faster, but leads to some weird artifacts in how lines are "prioritized" over each other...
#essentially, with the segments approach, things go in terms of edges, then sims, then time points
#so with the segments approach, you have control over how lines overlap by rearranging edges/sims
#I can't think of an intuitive way to do this in the case of join.lines, however...
#So I'm making segments the default approach again, but keeping join.lines as an interesting option
#there might be a compromise here where you could parse out the joined lines in some intelligent way based on simulations/edges/times...
#using the machinery in parse.args to get the edge, sim, and time indices...
#yeah, because the segments approach with 20, 200-tip, 200-res maps is taking FOREVER...
#' @export
plot.contsimmap<-function(contsimmap,traits=1,sims=min(20,nrow(contsimmap$perm.inds)),edges=NULL,
                          Col.by=c('simulation','state','edge','time',contsimmap[['traits']]),
                          Alpha.by=NULL,Mix.by=NULL,Wgt.by=NULL,Lty.by=NULL,Lwd.by=NULL,
                          scale.res=100,
                          add=FALSE,
                          join.lines=FALSE,
                          ...){
  #some stuff going very wrong here
  #things not matching up correctly AND running into errors when sims=1
  bys<-.parse.bys(c('simulation','state','edge','time',contsimmap[['traits']]),
                  Col.by,Alpha.by,Mix.by,Wgt.by,Lty.by,Lwd.by)
  contsimmap<-.parse.contsimmap(contsimmap,edges,sims,traits,bys,join.lines)
  args.ls<-.parse.args(contsimmap,list(...),bys,scale.res,join.lines)
  if(join.lines){
    #might be a way to make this more efficient...but it works
    par.fac<-do.call(paste,c(args.ls[c('col','lty','lwd')],sep='@'))
    tmp<-rle(par.fac)
    rep1s<-cumsum(tmp[['lengths']])
    rep0s<-c(1,rep1s[-length(rep1s)])
    times<-rep.int(1,length(par.fac))
    times[rep1s]<-2
    times[rep0s]<-times[rep0s]+1
    tmp[['lengths']]<-tmp[['lengths']]+2
    na.inds<-cumsum(tmp[['lengths']])
    na.inds<-c(1,na.inds[-length(na.inds)])
    tmp<-inverse.rle(tmp)
    args.ls[['x']]<-rep.int(args.ls[['x']],times)
    args.ls[['x']][na.inds]<-NA
    args.ls[['y']]<-rep.int(args.ls[['y']],times)
    args.ls[['y']][na.inds]<-NA
    args.ls[['x']]<-split(args.ls[['x']],tmp)
    args.ls[['y']]<-split(args.ls[['y']],tmp)
    new.args<-strsplit(names(args.ls[['x']]),'@')
    args.ls[['col']]<-lapply(new.args,'[[',1)
    args.ls[['lty']]<-lapply(new.args,'[[',2)
    args.ls[['lwd']]<-lapply(new.args,'[[',3)
    pruned.args<-args.ls[!(names(args.ls)%in%c('x','y','col','lty','lwd'))]
    if(!add){
      do.call(plot,c('x'=0,pch=NA,pruned.args))
    }
    for(i in seq_along(args.ls[['x']])){
      do.call(lines,
              c('x'=list(args.ls[['x']][[i]]),
                'y'=list(args.ls[['y']][[i]]),
                'col'=args.ls[['col']][i],
                'lty'=args.ls[['lty']][i],
                'lwd'=args.ls[['lwd']][i],
                pruned.args))
    }
  }else{
    if(!add){
      do.call(plot,c('x'=0,pch=NA,args.ls[!(names(args.ls)%in%c('x0','y0','x1','y1'))]))
    }
    do.call(segments,args.ls)
  }
}
