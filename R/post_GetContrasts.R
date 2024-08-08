
#needs better handling of subsetted contsimmaps
#as well as contsimmaps with polytomies
#' @export
get.contrasts<-function(contsimmap,traits=NULL){
  
  x<-get.nodes(contsimmap,traits=traits,internal.only=FALSE)
  tree<-attr(contsimmap,"tree")[[1]]
  
  e2n<-tree[["edge"]][,2]
  dd<-des.edges(tree)
  inds<-which(lengths(dd)>0)
  tmp<-c(Ntip(tree)+1,tree$edge[sapply(dd[inds],"[[",1),1])
  tmp<-tree[["node.label"]][tmp-Ntip(tree)]
  out<-array(dim=c(tree$Nnode,dim(x)[-1]),
             dimnames=c(list(tmp),dimnames(x)[-1]))
  d<-root.edges(tree)
  out[1,,]<-(x[e2n[d[2]],,]-x[e2n[d[1]],,])/
    sqrt(sum(tree$edge.length[d]))
  for(i in seq_len(tree$Nnode-1)){
    d<-dd[[inds[i]]]
    out[i+1,,]<-(x[e2n[d[2]],,]-x[e2n[d[1]],,])/
      sqrt(sum(tree$edge.length[d]))
  }
  
  out[tree[["node.label"]],,,drop=FALSE]
}