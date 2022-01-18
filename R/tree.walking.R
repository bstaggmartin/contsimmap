#' @export
edge.ranges<-function(phy){
  tmp<-ape::node.depth.edgelength(phy)
  matrix(tmp[phy$edge],ncol=2)
}

#' @export
root.edges<-function(phy){
  which(phy$edge[,1]==length(phy$tip.label)+1)
}

#not using this one as much as I should
#' @export
tip.edges<-function(phy,include.names=TRUE){
  tips<-phy$tip.label
  out<-match(seq_along(tips),phy$edge[,2])
  if(include.names){
    names(out)<-tips
  }
  out
}

#' @export
anc.edges<-function(phy,commonformat=TRUE){
  out<-match(phy$edge[,1],phy$edge[,2])
  if(commonformat){
    nulls<-is.na(out)
    out<-as.list(out)
    out[nulls]<-list(integer(0))
  }
  out
}

#' @export
des.edges<-function(phy){
  anc<-anc.edges(phy,commonformat=FALSE)
  len<-length(anc)
  out<-rep(list(integer(0)),len)
  inds<-!is.na(anc)
  anc<-anc[inds]
  tmp<-split(seq_len(len)[inds],anc)
  out[as.numeric(names(tmp))]<-tmp
  out
}

#' @export
sis.edges<-function(phy){
  anc<-anc.edges(phy,commonformat=FALSE)
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