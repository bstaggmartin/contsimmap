#add name matching for X0
#' @export
sim.conthistory<-function(tree,ntraits=1,traits=paste0('trait_',seq_len(ntraits)),
                          nsims=100,res=100,
                          Xsig2=diag(ntraits),
                          X0=0,
                          nobs=0,Ysig2=0){
  #formatting trees and parameters
  list2env(.proc.tree(tree,nsims),envir=environment())
  #OUTPUT:
  #tree = list of trees in multiSimmap format; treeID = which tree each simulation will correspond to;
  #nodes = node labels/NULL if they don't exist; nnodes = number of internal nodes;
  #edges = edge matrix; nedges = number of edges;
  #tips = tip labels; ntips = number of tips;
  #elen = edge lengths;
  #states = state labels; nstates = number of states
  Xsig2<-.fix.mat.list(Xsig2,ntraits,traits,nstates,states)
  traits<-rownames(Xsig2[[1]])
  tips<-tree[[1]]$tip.label
  ntips<-length(tips)
  out.trait.data<-if(!sum(nobs)) FALSE else TRUE
  if(out.trait.data){
    Ysig2<-.fix.mat.list(Ysig2,ntraits,traits,ntips,tips)
    nobs<-.fix.nobs(nobs,ntips,tips)
  }
  X0<-matrix(rep(X0,length.out=ntraits*nsims),ntraits,nsims)
  #initialize edge stuff
  list2env(.get.edge.maps(tree,treeID,nsims,res),envir=environment()) #outputs maps, ts, nts, sims.per.tree, seeds.per.tree.edge, and seeds.per.edge
  #nts is a matrix of the number of time points for each tree and edge, with rows corresponding to tree and cols to edge
  #sims.per.tree is a vector of the number of simulations for each tree
  #seeds.per.tree.edge is a matrix of the total number of seeds (i.e., random multinormal draws) for each tree and edge
  #dimensions of seeds.per.tree.edge match that of nts
  #seeds.per.edge is a vector of the total number of seeds per edge across all trees
  seed<-.get.seed(seeds.per.edge,ntraits)
  anc<-anc.edges(tree[[1]])
  anc[lengths(anc)==0]<-0
  anc<-as.character(unlist(anc))
  names(anc)<-seq_along(anc)
  cladewise.ord<-.get.cladewise.ord(anc)
  out<-.accumulate.seed(seed,
                        anc,cladewise.ord,maps,
                        Xsig2,X0,
                        ntraits,traits,nstates,states,
                        treeID,nts,sims.per.tree,seeds.per.tree.edge)
  out<-list(x=out,
            maps=maps,
            ts=ts,
            tree=tree,
            traits=traits,
            perm.inds=.get.perm.inds(treeID),
            anc=anc,
            nodes=list('0'=X0),
            params=list(Xsig2=Xsig2))
  if(out.trait.data){
    tip.vals<-aperm(get.tip.vals(out),c(2,3,1))
    dimnms<-dimnames(tip.vals)
    dimnms[[3]]<-rep(tips,nobs)
    seeds.per.tip<-nsims*nobs
    trait.matrix<-matrix(rnorm(ntraits*sum(seeds.per.tip)),ntraits)
    trait.data<-array(dim=c(ntraits,nsims,sum(nobs)),dimnames=dimnms)
    cholY<-lapply(Ysig2,function(ii) t(.pseudo.chol(ii)))
    j.counter<-0
    k.counter<-0
    for(i in tips){
      j<-nobs[i]
      if(j){
        k<-seeds.per.tip[i]
        array.inds<-seq_len(j)+j.counter
        matrix.inds<-seq_len(k)+k.counter
        trait.data[,,array.inds]<-as.vector(tip.vals[,,i])+cholY[[i]]%*%trait.matrix[,matrix.inds]
        j.counter<-j.counter+j
        k.counter<-k.counter+k
      }
    }
    out[['trait.data']]<-aperm(trait.data,c(3,1,2))
    out[['params']][['Ysig2']]<-Ysig2
  }
  class(out)<-'contsimmap'
  out
}
