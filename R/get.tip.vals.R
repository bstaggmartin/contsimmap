#generalize to arbitrary node values and ability to select specific traits/simulations
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