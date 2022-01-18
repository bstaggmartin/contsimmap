#old version
# .check.matching.tops<-function(tree){
#   all(unlist(lapply(tree[-1],function(ii) all.equal(ii,tree[[1]]))))
# }

#much faster than old function, but doesn't do as much attempting to coerce trees to be identical
#makes sure node, edge, and tips information is all equivalent for a given set of trees
.check.matching.tops<-function(tree,nodes,nnodes,edges,nedges,tips,ntips,elen){
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

.proc.tree<-function(tree,nsim){
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
  nnodes<-tree[[1]][['Nnode']]
  edges<-tree[[1]][['edge']]
  nedges<-nrow(edges)
  tips<-tree[[1]][['tip.label']]
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
    nstates<-1
    maps<-as.list(tree[[1]]$edge.length)
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
    nstates<-length(states)
  }
  treeID<-rep(seq_len(length(tree)),length.out=nsim)
  list(tree=tree,treeID=treeID,
       nodes=nodes,nnodes=nnodes,
       edges=edges,nedges=nedges,
       tips=tips,ntips=ntips,
       elen=elen,
       states=states,nstates=nstates)
}
