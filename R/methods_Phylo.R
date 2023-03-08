
####EDGE RANGES####

.edge.ranges<-function(phy){
  tmp<-node.depth.edgelength(phy)
  matrix(tmp[phy[['edge']]],ncol=2)
}

#' @export edge.ranges
edge.ranges<-function(phy){
  UseMethod('edge.ranges')
}

#' @export
#' @method edge.ranges phylo
edge.ranges.phylo<-function(phy){
  .edge.ranges(phy)
}

#' @export
#' @method edge.ranges multiPhylo
edge.ranges.multiPhylo<-function(phy){
  lapply(phy,.edge.ranges)
}

####ANCESTRAL EDGES####

.anc.edges<-function(phy,commonformat=TRUE){
  out<-match(phy$edge[,1],phy$edge[,2])
  if(commonformat){
    nulls<-is.na(out)
    out<-as.list(out)
    out[nulls]<-list(integer(0))
  }
  out
}

#' @export anc.edges
anc.edges<-function(phy){
  UseMethod('anc.edges')
}

#' @export
#' @method anc.edges phylo
anc.edges.phylo<-function(phy){
  .anc.edges(phy)
}

#' @export
#' @method anc.edges multiPhylo
anc.edges.multiPhylo<-function(phy){
  lapply(phy,.anc.edges)
}

####DESCENDANT EDGES####

.des.edges<-function(phy){
  anc<-.anc.edges(phy,commonformat=FALSE)
  len<-length(anc)
  out<-rep(list(integer(0)),len)
  inds<-!is.na(anc)
  anc<-anc[inds]
  tmp<-split(seq_len(len)[inds],anc)
  out[as.numeric(names(tmp))]<-tmp
  out
}

#' @export des.edges
des.edges<-function(phy){
  UseMethod('des.edges')
}

#' @export
#' @method des.edges phylo
des.edges.phylo<-function(phy){
  .des.edges(phy)
}

#' @export
#' @method des.edges multiPhylo
des.edges.multiPhylo<-function(phy){
  lapply(phy,.des.edges)
}

####SISTER EDGES####

.sis.edges<-function(phy){
  anc<-.anc.edges(phy,commonformat=FALSE)
  len<-length(anc)
  out<-rep(list(integer(0)),len)
  inds<-!is.na(anc)
  anc2<-anc[inds]
  tmp<-seq_len(len)
  des<-c(list(which(!inds)),split(tmp[inds],anc2))
  ndes<-lengths(des)
  di<-ndes==2
  di.des<-unlist(des[di])
  if(length(di.des)>0){
    odds<-seq.int(1,length(di.des),2)
    evens<-odds+1
    odds<-di.des[odds]
    evens<-di.des[evens]
    out[odds]<-evens
    out[evens]<-odds
  }
  poly.des<-des[!di]
  if(length(poly.des)>0){
    ndes<-ndes[!di]
    unlist.poly.des<-unlist(poly.des,use.names=FALSE)
    out[unlist.poly.des]<-rep(poly.des,ndes)
    foo<-function(x){
      tmp<-out[[x]]
      tmp[tmp!=x]
    }
    out[unlist.poly.des]<-lapply(unlist.poly.des,foo)
  }
  out
}

#' @export sis.edges
sis.edges<-function(phy){
  UseMethod('sis.edges')
}

#' @export
#' @method sis.edges phylo
sis.edges.phylo<-function(phy){
  .sis.edges(phy)
}

#' @export
#' @method sis.edges multiPhylo
sis.edges.multiPhylo<-function(phy){
  lapply(phy,.sis.edges)
}

####ROOT EDGES####

.root.edges<-function(phy){
  which(phy$edge[,1]==length(phy$tip.label)+1)
}

#' @export root.edges
root.edges<-function(phy){
  UseMethod('root.edges')
}

#' @export
#' @method root.edges phylo
root.edges.phylo<-function(phy){
  .root.edges(phy)
}

#' @export
#' @method root.edges multiPhylo
root.edges.multiPhylo<-function(phy){
  lapply(phy,.root.edges)
}

####TIP EDGES####

.tip.edges<-function(phy,include.names){
  tips<-phy$tip.label
  out<-match(seq_along(tips),phy$edge[,2])
  if(include.names){
    names(out)<-tips
  }
  out
}

#' @export tip.edges
tip.edges<-function(phy,...){
  UseMethod('tip.edges')
}

#' @export
#' @method tip.edges phylo
tip.edges.phylo<-function(phy,include.names=TRUE){
  .tip.edges(phy,include.names)
}

#' @export
#' @method tip.edges multiPhylo
tip.edges.multiPhylo<-function(phy,include.names=TRUE){
  lapply(phy,.tip.edges,include.names=include.names)
}

