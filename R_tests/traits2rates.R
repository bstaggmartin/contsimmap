#make simmaps with Xsig2 regimes based off contsimmap traits
#should make this more convenient down the line
#while I would love to let correlations continuously change, it's too easy to get impossible correlation structures...
#for now, set to numeric constant, but eventually allow for folks to call discrete states and new trait construction functions
#and let correlations vary between states
#' @export
traits2rates<-function(contsimmap,Xvar=NA,Xcor=diag(length(Xvar)),collapse=TRUE){
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(contsimmap[['nodes']])==1){
      stop('The traits2rates() function will support single subtrees in the future, but not yet')
    }else{
      stop('The traits2rates() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  if(!is.list(Xvar)){
    Xvar<-list(Xvar)
  }
  k<-length(Xvar)
  out.traits<-trimws(unlist(lapply(strsplit(as.character(Xvar),'~'),'[[',1),use.names=FALSE))
  out.traits<-paste0('TEMPORARY_',.fix.trait.names(k,out.traits)) #to avoid accidental replacements...
  calls<-lapply(Xvar,function(ii) ii[[length(ii)]])
  other.traits<-unlist(lapply(calls,all.vars),use.names=FALSE)
  foc.traits<-c(out.traits,other.traits)
  lhs<-rep(c(TRUE,FALSE),c(k,length(other.traits)))
  #maybe don't check for variable names, just let calls run their course and throw an error if they don't find the right variable
  #make sure to have some warning about how formulae are applied in sequence regarding new traits created in same function, replacements, etc...
  #actually not sure what the priority is regarding eval() with multiple variables named the same thing...
  foc.traits<-.extract.traits(contsimmap,foc.traits,lhs,traits2rates=TRUE)
  regime.map<-foc.traits[['fxn']]
  foc.traits<-foc.traits[-length(foc.traits)]
  for(i in seq_along(Xvar)){
    foc.traits[[i]]<-with(foc.traits,eval(calls[[i]]))
  }
  Xsds<-sqrt(do.call(rbind,foc.traits[-length(foc.traits)][lhs]))
  nms<-gsub('^TEMPORARY_','',rownames(Xsds))
  Xsig2<-array(Xcor,c(k,k,foc.traits[['placeholder_n']]),
               list(nms,nms,seq_len(foc.traits[['placeholder_n']])))
  Xsig2<-sweep(Xsig2,c(1,3),Xsds,'*')
  Xsig2<-sweep(Xsig2,c(2,3),Xsds,'*')
  Xsig2<-asplit(Xsig2,3)
  
  #lots of initialization for main loop
  out<-contsimmap[['tree']][[1]]
  class(out)<-c('simmap','phylo')
  out[['maps']]<-vector('list',nrow(out[['edge']]))
  out[['node.states']]<-out[['states']]<-out[['mapped.edge']]<-out[['Q']]<-out[['logL']]<-NULL
  nsims<-nsim(contsimmap)
  out<-rep(list(out),nsims)
  tree.seq<-seq_along(contsimmap[['tree']])
  perm.inds<-contsimmap[['perm.inds']]
  out.inds<-perm.inds[,1]
  out.inds<-split(seq_along(out.inds),out.inds)
  sims.per.tree<-lengths(out.inds)
  sim.seqs<-lapply(sims.per.tree,seq_len)
  perm.inds<-perm.inds[,3]
  mapped.edge.nms<-vector('list',nsims)
  nedges<-nedge(contsimmap)
  edge.seq<-seq_len(nedges)
  #main loop: preparing 'maps' elements
  for(i in edge.seq){
    edge<-contsimmap[['x']][[i]]
    maps<-contsimmap[['maps']][[i]]
    for(j in tree.seq){
      tmp.out.inds<-out.inds[[j]]
      tmp.map<-maps[[j]][['dts']]
      tmp.nms<-regime.map[[i]][[j]]
      tmp.nms[]<-as.character(tmp.nms)
      sim.seq<-sim.seqs[[j]]
      for(k in sim.seq){
        tmp.tmp.nms<-tmp.nms[,,k,drop=FALSE]
        out[[tmp.out.inds[k]]][['maps']][[i]]<-setNames(tmp.map,tmp.tmp.nms)
        #this could be made more efficient by pre-allocating elements using numbers computed in .extract.traits...but fine for now
        mapped.edge.nms[[tmp.out.inds[k]]]<-c(mapped.edge.nms[[tmp.out.inds[k]]],tmp.tmp.nms)
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
  ncols<-lengths(mapped.edge.nms)
  edge.nms<-paste(edge[,1],edge[,2],sep=',')
  foo<-function(i){
    holder<-matrix(0,nedges,ncols[i],
                   dimnames=list(edge.nms,'state'=mapped.edge.nms[[i]]))
    counter<-0
    nstates<-lengths(out[[i]][['maps']])
    for(j in edge.seq){
      holder[j,counter+seq_len(nstates[j])]<-out[[i]][['maps']][[j]]
      counter<-counter+nstates[j]
    }
    holder
  }
  for(i in seq_len(nsims)){
    nodes[]<-unlist(lapply(out[[i]][['maps']],
                           function(ii) names(ii)[c(1,length(ii))]),use.names=FALSE)
    tips[]<-nodes[2,tip.inds]
    out[[i]][c('mapped.edge','node.states','states')]<-list(foo(i),t(nodes),tips)
  }
  class(out)<-c('multiSimmap','multiPhylo')
  
  if(collapse){
    tmp<-collapse.traits2rates(out,Xsig2)
    out<-tmp$tree
    Xsig2<-tmp$Xsig2
  }
  list(tree=out,Xsig2=Xsig2)
}

#function to simplify above maps to only have 1 state per branch for quicker simulation!
#' @export
collapse.traits2rates<-function(tree,Xsig2){
  #need to add some better checks here
  #nice thing is this function should work with differing topologies, unlike other functions!
  nedges<-unlist(lapply(tree,Nedge),use.names=FALSE)
  sum.nedges<-sum(nedges)
  out.Xsig2<-setNames(rep(Xsig2[1],sum.nedges),seq_len(sum.nedges))
  dims<-dim(Xsig2[[1]])
  counter<-1
  for(i in seq_along(tree)){
    tree[[i]][['Q']]<-tree[[i]][['logL']]<-NULL
    elen<-tree[[i]][['edge.length']]
    tmp<-diag(elen)
    rownames(tmp)<-rownames(tree[[i]]$mapped.edges)
    edge.seq<-seq_along(elen)
    states<-counter+edge.seq-1
    colnames(tmp)<-states
    tree[[i]]$mapped.edge<-tmp
    tree[[i]]$node.states[]<-states
    tree[[i]]$states[]<-counter+which(tree[[i]]$edge[,2]<=Ntip(tree[[i]]))-1
    #need to handle edge-cases involving 0 branch lengths still!
    maps<-tree[[i]][['maps']]
    nts<-lengths(maps)
    dts<-unlist(maps)
    dts<-dts/rep(elen,nts)
    tmp<-array(unlist(Xsig2[names(dts)],use.names=FALSE),c(dims,sum(nts)))
    tmp<-sweep(tmp,3,dts,'*')
    j.counter<-0
    for(j in seq_along(nts)){
      out.Xsig2[[counter]][]<-eval(str2lang(paste0('tmp[,,',j.counter+seq_len(nts[j]),',drop=FALSE]',collapse='+')))
      counter<-counter+1
      j.counter<-j.counter+nts[j]
    }
    tree[[i]]$maps<-unname(split(setNames(elen,states),edge.seq))
  }
  list(tree=tree,Xsig2=out.Xsig2)
}