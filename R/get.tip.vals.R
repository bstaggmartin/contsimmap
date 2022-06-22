#generalize to arbitrary node values (done)
#add ability to select specific traits/simulations
#make get.tip.vals a wrapper for get.node.vals

#' @export
get.tip.vals<-function(contsimmap){
  ntips<-length(contsimmap[['tree']][[1]][['tip.label']])
  edge.mat<-contsimmap[['tree']][[1]][['edge']][as.numeric(names(contsimmap[['x']])),,drop=FALSE]
  tips<-which(edge.mat[,2]<(ntips+1))
  names(tips)<-contsimmap[['tree']][[1]][['tip.label']][edge.mat[tips,2]]
  ntips<-length(tips)
  ntraits<-ntrait(contsimmap)
  nsims<-nsim(contsimmap)
  nms<-list('tip'=names(tips),'trait'=contsimmap[['traits']],'sim'=NULL)
  out<-array(dim=c(ntips,ntraits,nsims),dimnames=nms)
  for(i in seq_len(ntips)){
    out[i,,]<-unlist(lapply(contsimmap[['x']][[tips[i]]],function(ii) ii[dim(ii)[1],,]))
    out[i,,]<-out[i,,contsimmap[['perm.inds']][,3],drop=FALSE]
  }
  out
}

#' @export
get.node.vals<-function(contsimmap,nodes=NULL){
  ntips<-Ntip(contsimmap)
  edge.mat<-contsimmap[['tree']][[1]][['edge']]
  avail.edges<-avail.nodes<-vector('list',2)
  edges<-as.numeric(names(contsimmap[['nodes']]))
  tmp<-as.numeric(names(contsimmap[['nodes']]))
  tmp[tmp==0]<-NA
  avail.edges[[1]]<-tmp
  avail.edges[[2]]<-as.numeric(names(contsimmap[['x']]))
  for(i in seq_len(2)){
    tmp<-edge.mat[avail.edges[[i]],2]
    tmp[is.na(tmp)]<-ntips+1
    avail.nodes[[i]]<-tmp
  }
  lens<-lengths(avail.nodes)
  offset<-lens[1]
  avail.nodes<-unlist(avail.nodes,use.names=FALSE)
  if(is.null(nodes)){
    nodes<-sort(avail.nodes)
  }else{
    #probably return a warning about ignoring unavailable nodes...
    nodes<-nodes[nodes%in%avail.nodes]
  }
  nnodes<-length(nodes)
  inds<-match(nodes,avail.nodes)
  ntraits<-ntrait(contsimmap)
  nsims<-nsim(contsimmap)
  nms<-nodes
  nms[nms<=ntips]<-contsimmap[['tree']][[1]][['tip.label']]
  nms<-list('node'=nms,'trait'=contsimmap[['traits']],'sim'=NULL)
  out<-array(dim=c(nnodes,ntraits,nsims),dimnames=nms)
  for(i in seq_len(nnodes)){
    ii<-inds[i]
    if(ii<=offset){
      out[i,,]<-contsimmap[['nodes']][[ii]]
    }else{
      out[i,,]<-unlist(lapply(contsimmap[['x']][[ii-offset]],function(ii) ii[dim(ii)[1],,]),use.names=FALSE)
      out[i,,]<-out[i,,contsimmap[['perm.inds']][,3],drop=FALSE]
    }
  }
  out
}

# #extract "exact" values of contsimmap across all edges for arbitrary time point
# #make exact by simulating according to Brownian Bridge rules?
# #maybe not, since it assumes that assumes all traits are untransformed...
# #it would really be helpful to have some sort of means of integration/drawing "geodesic" lines,
# ##but I have no idea how I'd do that for arbitrary user functions without serious overhead...
# get.val.at.t<-function(contsimmap,t){
#   edgerans<-edge.ranges(contsimmap[['tree']][[1]])
#   foc.edges<-t>=edgerans[,1]&t<=edgerans[,2]
#   
# }