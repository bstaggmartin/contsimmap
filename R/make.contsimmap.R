# #new approach-->work with codes!
# #0 --> missing data, 1 --> observed data, 2 --> observed data with infinite precision/0 variance
# #when combining, higher numbers take priority (use max())
# #indices for pseudo.solve operations valid whenever code is 1...
# 
# #need to make this more efficient
# #0s and 2s will only ever pop up for traits in which an observation has a 0 or 2...
# 
# .gen.codes<-function(codes){
#   codes<-cbind(1,codes)
#   out<-as.matrix(do.call(expand.grid,apply(codes,1,unique)))
#   code.lab<-do.call(paste0,asplit(out,2))
#   list(code.num=out,code.ind=out==1,code.inf=out==2,code.lab=code.lab)
# }



#may want to make some function for auto-estimation of Xsig2/Ysig2
#also should add some means of conditionining on specific X0s
#need better warning/feedback about averaging behavior in cases when observations with infinite precision conflict with one another
  #this produces nonsensical results in many cases, I would think...
  #probably need a disclaimer that absolute certainty in observations, though not necessarily 0 branch lengths, should be avoided
  #in general, the parts of this function that deal with infinite precision could probably be improved in terms of efficiency
  #new approach-->keep track of perfect and missing observations with codes, index to avoid directly dealing with them
  #missing should be constant across sims, but perfect could change (0 entries in Xsig2...)
  #the key is to keep track of infinite precisions in separate variance/precision matrices, with entries only being a function of branch lengths
  #(essentially layers a new "unit-BM" process on top that takes "priority" over existing one)
  #above approach may seem desirable, but has odd, recursive feature where observations connected by 0 branch lengths will still pose
  #serious issues...it might be better to instead to introduce some "fudge factor" for 0s, like adding 1e-6 to X/Ysig2s with 0 entries along
  #diagonal!
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
  
  #initialize arrays
  nn<-ntips+nnodes
  ntrees<-length(tree)
  VV<-PP<-array(0,
                dim=c(ntraits,ntraits,nn,ntrees),
                dimnames=list(traits,traits,c(tips,seq.int(ntips+1,nn)),NULL))
  inf.I<-matrix(0,ntraits,ntraits)
  diag(inf.I)<-Inf
  VV[]<-inf.I
  XX<-array(0,dim=c(ntraits,nn,ntrees),dimnames=dimnames(PP)[-1])
  
  ##CAN MAKE THIS ALL MORE EFFICIENT THROUGH IMPLICIT Infs/0s VIA KEEPING TRACK OF CODES##
  ##ALSO PAVES THE WAY FOR MORE ELEGANT HANDLING OF CONFLICTING OBSERVATIONS##
  ##WILL REQUIRE A LARGE, TIME-CONSUMING OVERHAUL##
  # inf.PP<-PP #I'll start working on this later...
  # inf.VV<-VV
  # CC<-matrix(0L,nn,ntrees) #keeps track of codes...
  # #species means + variance
  # #I like where this is going, but need to deal with 0 variance, infinite precision cases...
  # codes<-Yvar<-Y<-t(trait.data)
  # tip.matches<-match(colnames(Y),tips)
  # dims<-c(ntraits,ntraits)
  # Yvar[]<-unlist(Ysig2[tip.matches],use.names=FALSE)[.row(dims)==.col(dims)]
  # codes[]<-0
  # codes[is.na(Y)]<-1
  # codes[Yvar==0]<-2
  # 
  # #initialize pseudo-solve operations
  # list2env(.gen.codes(codes),envir=environment())
  # #OUTPUT:
  # #code.num = matrix of numeric integers corresponding to whether trait value is missing (0), known with error (1), or known without error (2)
  # #code.ind = matrix of TRUE/FALSE indicating if a trait value should be included in pseudo matrix operations (only if trait value is known with error)
  # #code.inf = matrix of TRUE/FALSE indicating if a trait value is known without error
  # #code.lab = vector of labels for all codes (just concatenated code.num)
  # 
  # #this function combines the codes for multiple observation
  # com.codes<-function(codes){
  #   match(paste0(apply(codes.num[codes,,drop=FALSE],2,max),collapse=''),code.lab)
  # }
  # #this functions pseudo-inverts a single matrix given corresponding code.ind
  # holder<-matrix(0,ntraits,ntraits)
  # single.pseudo.solve<-function(mat,code.ind){
  #   tmp<-mat[code.ind,code.ind,drop=FALSE]
  #   if(length(tmp)){
  #     holder[code.ind,code.ind]<-solve(mat[code.ind,code.ind,drop=FALSE])
  #   }
  #   holder
  # }
  # #this function pseudo-inverts many matrices given a list of matrices and vector of associated code.lab's
  # pseudo.solve<-function(mats,code.lab){
  #   code.ind<-code.ind[code.lab,,drop=FALSE]
  #   lapply(seq_along(mats),function(ii) single.pseudo.solve(mats[[ii]],code.ind[ii,]))
  # }
  # 
  # Y<-split(Y,rep(tip.matches,each=ntraits))
  # Y<-lapply(Y,function(ii) matrix(ii,nrow=ntraits))
  # codes<-split(codes,tip.matches)
  # 
  # for(i in as.numeric(names(codes))){
  #   tmp.tip<-tips[i]
  #   tmp.Y<-Y[[i]]
  #   tmp.codes<-codes[[i]]
  #   unique.codes<-unique(tmp.codes)
  #   unique.matches<-match(tmp.codes,unique.codes)
  #   unique.n<-length(unique.codes)
  #   tmp.des_P<-pseudo.solve(rep(Ysig2[tmp.tip],unique.n),unique.codes)[unique.matches]
  #   tip.code<-com.codes(unique.codes)
  #   CC[tmp.tip,]<-tip.code
  #   tip.ind<-code.ind[tip.code,]
  #   PP[tip.ind,tip.ind,tmp.tip,]<-Reduce('+',tmp.des_P)[tip.ind,tip.ind,drop=FALSE]
  #   VV[tip.ind,tip.ind,]<-solve(PP[tip.ind,tip.ind,tmp.tip,1])
  #   
  #   #infinite precision observations
  #   tip.inf<-code.inf[tip.code,]
  #   PP[tip.inf,tip.inf,tmp.tip,]<-inf.I[tip.inf,tip.inf]
  #   VV[tip.inf,tip.inf,tmp.tip,]<-0
  # }
  
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
      }else{
        XX[inf.prec,i,]<-des_X
      }
    }else{
      VV[,,i,]<-inf.I
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
