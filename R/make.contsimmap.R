#may want to make some function for auto-estimation of Xsig2/Ysig2
#also should add some means of conditionining on specific X0s
#need better warning/feedback about averaging behavior in cases when observations with infinite precision conflict with one another
  #this produces nonsensical results in many cases, I would think...
  #probably need a disclaimer that absolute certainty in observations, though not necessarily 0 branch lengths, should be avoided
  #in general, the parts of this function that deal with infinite precision could probably be improved in terms of efficiency
#' @export
make.contsimmap<-function(trait.data,tree,nsims=100,res=100,
                          Xsig2=diag(NCOL(trait.data)),
                          Ysig2=0){
  #formatting tree, data, and parameters
  list2env(.proc.tree(tree,nsims),envir=environment())
  #OUTPUT:
  #tree = list of trees in multiSimmap format; treeID = which tree each simulation will correspond to;
  #nodes = node labels/NULL if they don't exist; nnodes = number of internal nodes;
  #edges = edge matrix; nedges = number of edges;
  #tips = tip labels; ntips = number of tips;
  #elen = edge lengths;
  #states = state labels; nstates = number of states
  #temporary, bare bones processing of trait.data--eventually want to have functions to check/format this as well...
  ntraits<-ncol(trait.data)
  traits<-colnames(trait.data)
  Xsig2<-.fix.mat.list(Xsig2,ntraits,traits,nstates,states)
  traits<-rownames(Xsig2[[1]])
  Ysig2<-.fix.mat.list(Ysig2,ntraits,traits,ntips,tips)
  nobs<-tapply(rownames(trait.data),rownames(trait.data),length)
  nobs<-.fix.nobs(nobs,ntips,tips)
  
  #initialize arrays
  nn<-ntips+nnodes
  ntrees<-length(tree)
  PP<-array(0,
            dim=c(ntraits,ntraits,nn,ntrees),
            dimnames=list(traits,traits,c(tips,seq.int(ntips+1,nn)),NULL))
  VV<-PP
  XX<-array(dim=c(ntraits,nn,ntrees),dimnames=dimnames(PP)[-1])
  
  #species means + variance
  infs<-matrix(0,ntraits,ntraits)
  diag(infs)<-Inf
  for(i in tips){
    des_X<-t(trait.data[rownames(trait.data)==i,,drop=FALSE])
    des_n<-ncol(des_X)
    if(des_n){
      #could be sped up by avoiding calculations for traits with no observations/no intraspecific variance
      missing.inds<-which(is.na(des_X),arr.ind=TRUE)
      des_X[missing.inds]<-0
      des_V<-array(Ysig2[[i]],c(ntraits,ntraits,des_n))
      des_V[missing.inds[,c(1,1,2),drop=FALSE]]<-Inf
      des_P<-lapply(asplit(des_V,3),.pseudo.solve)
      PP[,,i,]<-Reduce('+',des_P)
      VV[,,i,]<-.pseudo.solve(PP[,,i,1])
      tmp.XX<-VV[,,i,1]%*%Reduce('+',lapply(seq_len(des_n),function(ii) des_P[[ii]]%*%des_X[,ii]))
      tmp.XX[is.nan(tmp.XX)]<-0
      XX[,i,]<-tmp.XX
      inf.prec<-!diag(Ysig2[[i]])
      if(any(inf.prec)&des_n>1){
        warning('PLACEHOLDER ABOUT MULTIPLE OBSERVATIONS WITH 0 INTRASPECIFIC VARIANCE')
        XX[inf.prec,i,]<-.rowMeans(des_X[inf.prec,,drop=FALSE],sum(inf.prec),des_n)
      }
    }else{
      VV[,,i,]<-infs
      XX[,i,]<-0
    }
  }
  
  #tree info
  foo<-function(tree){
    maps<-tree[['mapped.edge']]
    probs<-!(states%in%colnames(maps))
    sum.probs<-sum(probs)
    if(sum.probs){
      maps<-do.call(cbind,c(list(probs),setNames(rep(list(0),sum.probs),states[probs])))
    }
    maps[,states,drop=FALSE]
  }
  maps<-array(unlist(lapply(tree,foo),use.names=FALSE),
              c(nedges,
                nstates,
                ntrees))
  maps<-aperm(maps,c(2,3,1))
  Xsig2.array<-array(unlist(Xsig2,use.names=FALSE),c(ntraits,ntraits,nstates,ntrees))
  des_e<-des.edges(tree[[1]])
  root_e<-root.edges(tree[[1]])
  anc_n<-edges[,2]
  attr(tree[[1]],'order')<-NULL
  prune.seq<-c(reorder(tree[[1]],'pruningwise',index.only=TRUE),0)
  tree.seq<-seq_len(ntrees)
  #check that branch lengths of 0 propagate correctly...
  #I think I need a manual check for any infinite precisions that might crop up!
  for(i in prune.seq){
    if(!i){
      des<-root_e
      anc<-ntips+1
    }else{
      des<-des_e[[i]]
      anc<-anc_n[i]
    }
    des_n<-length(des)
    if(des_n){
      des_N<-anc_n[des]
      des_V<-VV[,,des_N,,drop=FALSE]
      des_P<-des_V
      tmp.XX<-des_X<-XX[,des_N,,drop=FALSE]
      des_Xsig2<-array(Xsig2.array,c(dim(Xsig2.array),des_n))
      des_maps<-maps[,,des,drop=FALSE]
      des_Xsig2<-sweep(des_Xsig2,c(3,4,5),des_maps,'*')
      des_Xsig2<-Reduce('+',asplit(des_Xsig2,3))
      des_P[]<-unlist(lapply(asplit(des_V+aperm(des_Xsig2,c(1,2,4,3)),c(3,4)),.pseudo.solve),use.names=FALSE)
      #need to figure out some way to ensure 0s in off-diagonal entries with 1 infinity along corresponding diagonal entries STAY 0
      #I think below works
      PP[,,anc,]<-Reduce('+',asplit(des_P,3))
      inf.prec<-which(is.infinite(PP[,,anc,,drop=FALSE]),arr.ind=TRUE,useNames=FALSE)
      inf.prec[,3]<-anc
      PP[inf.prec[,1],,anc,inf.prec[,4]]<-0
      PP[,inf.prec[,1],anc,inf.prec[,4]]<-0
      PP[inf.prec]<-Inf
      VV[,,anc,]<-unlist(lapply(asplit(PP[,,anc,,drop=FALSE],c(3,4)),.pseudo.solve),use.names=FALSE)
      for(j in seq_len(des_n)){
        tmp.XX[,j,]<-unlist(lapply(tree.seq,function(ii) des_P[,,j,ii]%*%des_X[,j,ii]),use.names=FALSE)
      }
      tmp.XX<-Reduce('+',asplit(tmp.XX,2))
      #need this for 0-length branches
      tmp.XX[is.infinite(tmp.XX)]<-0
      tmp.XX[]<-unlist(lapply(tree.seq,function(ii) VV[,,anc,ii]%*%tmp.XX[,ii]),use.names=FALSE)
      tmp.XX[is.nan(tmp.XX)]<-0
      #0-length branches...
      inf.prec<-which(is.infinite(des_P),arr.ind=TRUE,useNames=FALSE)
      tmp<-split(des_X[inf.prec[,-1,drop=FALSE]],list(inf.prec[,1],inf.prec[,4]))
      out.inds<-do.call(rbind,strsplit(names(tmp),'\\.'))
      mode(out.inds)<-'numeric'
      tmp.XX[out.inds]<-unlist(lapply(tmp,sum),use.names=FALSE)/lengths(tmp)
      XX[,anc,]<-tmp.XX
    }
  }
  
  #I think it perhaps makes sense to do all the preorder stuff and Brownian Bridge sampling in "one step", so to speak...
  list2env(.get.edge.maps(tree,treeID,nsims,res),envir=environment()) #outputs maps, ts, nts, sims.per.tree, seeds.per.tree.edge, and seeds.per.edge
  seed<-.get.seed(seeds.per.edge,ntraits)
  anc<-anc.edges(tree[[1]])
  anc[lengths(anc)==0]<-0
  anc<-as.character(unlist(anc))
  names(anc)<-seq_along(anc)
  cladewise.ord<-.get.cladewise.ord(anc)
  #new for conditional contsimmap
  coarse.maps<-lapply(tree,'[[','maps')
  #perhaps could use some speeding up...but sufficient for now
  out<-.accumulate.seed.conditional(seed,
                                    anc,cladewise.ord,maps,
                                    Xsig2,XX,PP,VV,edges,ntips,coarse.maps,
                                    ntraits,traits,nstates,states,
                                    treeID,nts,sims.per.tree,seeds.per.tree.edge)
  perm.inds<-.get.perm.inds(treeID)
  X0<-out[['X0']][,perm.inds[,3],drop=FALSE]
  out<-out[['out']]
  out<-list(x=out,
            maps=maps,
            ts=ts,
            tree=tree,
            traits=traits,
            perm.inds=perm.inds,
            anc=anc,
            nodes=list('0'=X0),
            params=list(Xsig2=Xsig2),
            trait.data=array(trait.data,c(dim(trait.data),1),c(dimnames(trait.data),'sim'=NULL)))
  class(out)<-'contsimmap'
  out
}
