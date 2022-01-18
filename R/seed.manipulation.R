.get.seed<-function(seeds.per.edge,ntraits){
  e<-length(seeds.per.edge)
  tmp.seed<-matrix(rnorm(ntraits*sum(seeds.per.edge)),nrow=ntraits)
  seed<-vector('list',e)
  counter<-0
  for(i in seq_len(e)){
    seed[[i]]<-tmp.seed[,counter+seq_len(seeds.per.edge[i]),drop=FALSE]
    counter<-counter+seeds.per.edge[i]
  }
  seed
}

.accumulate.seed<-function(seed,
                           ancs,cladewise.ord,maps,
                           Xsig2,X0,
                           ntraits,traits,nstates,states,
                           treeID,nts,sims.per.tree,seeds.per.tree.edge){
  ancs<-as.numeric(ancs)
  tree.seq<-seq_along(maps[[1]])
  cholX<-lapply(Xsig2,function(ii) t(.pseudo.chol(ii)))
  X0<-split(X0,rep(treeID,each=ntraits))
  X0<-lapply(tree.seq,function(ii) matrix(X0[[ii]],ntraits,sims.per.tree))
  e<-length(seed)
  out<-setNames(vector('list',e),seq_len(e))
  for(i in cladewise.ord){
    map<-maps[[i]]
    state<-unlist(rep(lapply(map,'[[','state'),sims.per.tree))
    dx<-seed[[i]]
    #covariance structure
    for(j in states){
      inds<-state==j
      dx[,inds]<-cholX[[j]]%*%dx[,inds]
    }
    #scale
    dts<-unlist(rep(lapply(map,function(ii) sqrt(ii[['dts']])),sims.per.tree))
    dx<-sweep(dx,2,dts,'*')
    #accumulate
    anc<-ancs[i]
    is.root<-anc==0
    x<-lapply(tree.seq,function(j) array(dim=c(ntraits,nts[j,i],sims.per.tree[j])))
    counter<-0
    for(j in tree.seq){
      x[[j]][]<-dx[,counter+seq_len(seeds.per.tree.edge[j,i])]
      counter<-counter+seeds.per.tree.edge[j,i]
      if(is.root){
        x[[j]][,1,]<-x[[j]][,1,]+X0[[j]]
      }else{
        x[[j]][,1,]<-x[[j]][,1,]+out[[anc]][[j]][nts[j,anc],,]
      }
      if(dim(x[[j]])[2]>1){ #prevents dimension dropping
        x[[j]]<-apply(x[[j]],c(1,3),cumsum) #faster than Reduce()
        #must use c(1,3) over -2, the latter of which can cause further dimension dropping
        # tmp<-aperm(x,c(1,3,2))
        # tmp<-unlist(Reduce(foo,asplit(x,2),accumulate=TRUE))
        # x[[j]][]<-aperm(tmp,c(1,3,2))
        #where foo is a bridge construction function!
      }else{
        x[[j]]<-aperm(x[[j]],c(2,1,3))
      }
    }
    out[[i]]<-x
  }
  out
}

#some parts of this function should probably be part of edge mapping functions instead
  #will improve code cleanliness and efficiency
.accumulate.seed.conditional<-function(seed,
                                       ancs,cladewise.ord,maps,
                                       Xsig2,XX,PP,VV,edge,ntips,coarse.maps,
                                       ntraits,traits,nstates,states,
                                       treeID,nts,sims.per.tree,seeds.per.tree.edge){
  ancs<-as.numeric(ancs)
  ntrees<-length(maps[[1]])
  tree.seq<-seq_len(ntrees)
  cholX<-lapply(Xsig2,function(ii) t(.pseudo.chol(ii)))
  XX<-asplit(XX,3)
  PP<-asplit(PP,4)
  VV<-asplit(VV,4)
  X0<-lapply(XX,function(ii) ii[,ntips+1])
  cholX0<-lapply(VV,function(ii) t(.pseudo.chol(ii[,,ntips+1])))
  X0.seed<-lapply(sims.per.tree,function(ii) matrix(rnorm(ii*ntraits),ntraits,ii))
  X0<-lapply(tree.seq,function(ii) cholX0[[ii]]%*%X0.seed[[ii]]+X0[[ii]])
  e<-length(seed)
  edge.seq<-seq_len(e)
  coarse.maps<-lapply(edge.seq,function(ii) lapply(coarse.maps,'[[',ii))
  out<-setNames(vector('list',e),edge.seq)
  I<-diag(ntraits)
  for(i in cladewise.ord){
    #basic  stuff...
    map<-maps[[i]]
    coarse.map<-coarse.maps[[i]]
    foo<-function(ii){
      out<-!ii[['incl']]
      out[length(out)]<-TRUE
      out
    }
    incl.ls<-lapply(map,foo)
    incl<-unlist(rep(incl.ls,sims.per.tree))
    x<-lapply(tree.seq,function(j) array(dim=c(ntraits,nts[j,i],sims.per.tree[j])))
    #get dx's...
    dx<-seed[[i]]
    nodes.dx<-dx[,incl,drop=FALSE]
    incl<-!incl
    dx<-dx[,incl,drop=FALSE]
    #get nodes...
    anc<-ancs[i]
    is.root<-anc==0
    new.seeds.per.tree.edge<-seeds.per.tree.edge
    counter<-0
    for(j in tree.seq){
      tmp_map<-c(0,coarse.map[[j]])
      tmp_n<-length(tmp_map)
      tmp_nn<-(tmp_n-1)*sims.per.tree[j]
      new.seeds.per.tree.edge[j,i]<-seeds.per.tree.edge[j,i]-tmp_nn
      tmp_dx<-tmp_XX<-array(XX[[j]][,edge[i,2]],
                            c(ntraits,tmp_n,sims.per.tree[j]))
      if(is.root){
        tmp_XX[,1,]<-X0[[j]]
      }else{
        tmp_XX[,1,]<-aperm(out[[anc]][[j]][nts[j,anc],,,drop=FALSE],c(2,1,3))
      }
      tmp_VV<-tmp_PP<-array(PP[[j]][,,edge[i,2]],c(ntraits,ntraits,tmp_n))
      tmp_VV[]<-VV[[j]][,,edge[i,2]]
      tmp_PP[,,1]<-PP[[j]][,,edge[i,1]]
      tmp_VV[,,1]<-VV[[j]][,,edge[i,1]]
      tmp_dx<-tmp_dx[,-1,,drop=FALSE]
      tmp_dx[]<-nodes.dx[,counter+seq_len(tmp_nn)]
      counter<-counter+tmp_nn
      tmp_seq<-seq_len(tmp_n)
      anc_XX<-matrix(nrow=ntraits,ncol=sims.per.tree[j])
      for(k in rev(tmp_seq)[-c(1,tmp_n)]){
        tmp_VV[,,k]<-tmp_VV[,,k+1]+tmp_map[k+1]*Xsig2[[names(tmp_map)[k+1]]]
        tmp_PP[,,k]<-.pseudo.solve(tmp_VV[,,k])
        #XX is just descendant XX
      }
      for(k in tmp_seq[-1]){
        tmp1<-tmp_map[k]*Xsig2[[names(tmp_map)[k]]]
        tmp2<-tmp1%*%tmp_PP[,,k]
        #if any tmp_PP is infinite, then current node has infinite precision and must remain the same
        #if any of tmp1 is 0, then ancestral node is equivalent to current node and current node should be set to ancestral node
        anc.inf.prec<-!(diag(tmp1))
        cur.inf.prec<-is.infinite(diag(tmp_PP[,,k]))
        #correct for any infinite precisions to prevent NaN propagation
        tmp2[anc.inf.prec,]<-0
        tmp2[,anc.inf.prec]<-0
        tmp2[cur.inf.prec,]<-0
        tmp2[,cur.inf.prec]<-0
        diag(tmp2)[cur.inf.prec]<-Inf
        cur_PP<-tmp2%*%.pseudo.solve(I+tmp2)
        #correct for any infinite precisions to prevent NaN propagation
        cur_PP[anc.inf.prec,]<-0
        cur_PP[,anc.inf.prec]<-0
        cur_PP[cur.inf.prec,]<-0
        cur_PP[,cur.inf.prec]<-0
        diag(cur_PP)[cur.inf.prec]<-Inf
        anc_XX[]<-tmp_XX[,k-1,]
        cur_XX<-tmp_XX[,k,1]
        tmp_XX[,k,]<-rep(cur_PP%*%cur_XX,sims.per.tree[j])+anc_XX-cur_PP%*%anc_XX
        #correct for any infinite precisions to prevent NaN propagation
        tmp_XX[anc.inf.prec,k,]<-anc_XX[anc.inf.prec,]
        tmp_XX[cur.inf.prec,k,]<-cur_XX[cur.inf.prec]
        #average in case of both? not sure...this produces very odd results in some cases
        #intuitively, branch lengths should still weight observations here, but doing so would be complicated and of little benefit...
        tmp_inds<-anc.inf.prec&cur.inf.prec
        tmp_XX[tmp_inds,k,]<-(anc_XX[tmp_inds,]+cur_XX[tmp_inds])/2
        cur_PP<-tmp_PP[,,k]+.pseudo.solve(tmp1)
        #off-diagonals automatically converted to 0 via .pseudo.solve
        tmp_XX[,k,]<-t(.pseudo.chol(.pseudo.solve(cur_PP)))%*%tmp_dx[,k-1,]+tmp_XX[,k,]
      }
      x[[j]][,incl.ls[[j]],]<-tmp_XX[,-1,]
    }
    state.ls<-lapply(map,'[[','state')
    dt.ls<-lapply(map,'[[','dts')
    #might consider making this part of map.increments function...
    foo<-function(j){
      dt<-dt.ls[[j]]
      num<-c(rev(cumsum(rev(dt[-1]))),0)
      sub<-num[incl.ls[[j]]]
      tmp<-rle(state.ls[[j]])
      tmp$values<-sub
      num<-num-inverse.rle(tmp)
      denom<-num+dt
      out<-dt/denom
      out<-list(out,
                sqrt(out*num))
      out
    }
    dt.ls<-lapply(tree.seq,foo)
    state<-unlist(rep(state.ls,sims.per.tree))[incl]
    dts<-unlist(rep(lapply(dt.ls,'[[',2),sims.per.tree))[incl]
    #covariance structure
    for(j in states){
      inds<-state==j
      dx[,inds]<-cholX[[j]]%*%dx[,inds]
    }
    #scale
    dx<-sweep(dx,2,dts,'*')
    #accumulate
    j.counter<-0
    #I think I need two counters...
    for(j in tree.seq){
      if(new.seeds.per.tree.edge[j,i]){
        tmp_incl<-incl.ls[[j]]
        x[[j]][,!tmp_incl,]<-dx[,j.counter+seq_len(new.seeds.per.tree.edge[j,i])]
        tmp_dx<-x[[j]][,!tmp_incl,,drop=FALSE]
        tmp_dt<-dt.ls[[j]][[1]][!tmp_incl]
        j.counter<-j.counter+new.seeds.per.tree.edge[j,i]
        runs<-rle(state.ls[[j]])
        if(is.root){
          x0<-X0[[j]]
        }else{
          x0<-out[[anc]][[j]][nts[j,anc],,]
        }
        tmp_XX<-c(list(x0),asplit(x[[j]][,tmp_incl,,drop=FALSE],2))
        k.counter<-0
        for(k in seq_along(tmp_XX)[-1]){
          tmp_n<-runs$lengths[k-1]-1
          if(tmp_n){
            cur_XX<-tmp_XX[[k]]
            tmp_seq<-seq_len(tmp_n)
            inds<-k.counter+tmp_seq
            k.counter<-k.counter+tmp_n
            cur_dx<-tmp_dx[,inds,,drop=FALSE]
            cur_dt<-tmp_dt[inds]
            #can't use Reduce because it requires 3 arguments...
            cur_dx[,1,]<-tmp_XX[[k-1]]+(cur_XX-tmp_XX[[k-1]])*cur_dt[1]+cur_dx[,1,]
            for(tt in tmp_seq[-1]){
              cur_dx[,tt,]<-cur_dx[,tt-1,]+(cur_XX-cur_dx[,tt-1,])*cur_dt[tt]+cur_dx[,tt,]
            }
            tmp_dx[,inds,]<-cur_dx
          }
        }
        x[[j]][,!tmp_incl,]<-tmp_dx
      }
      x[[j]]<-aperm(x[[j]],c(2,1,3))
    }
    out[[i]]<-x
  }
  list(out=out,X0=matrix(unlist(X0,use.names=FALSE),nrow=ntraits,ncol=sum(sims.per.tree)))
}
