.subdivide.edges<-function(tree,res){
  edgerans<-evorates::edge.ranges(tree)
  dt<-max(edgerans)/res
  nts<-ceiling(tree$edge.length/dt)
  ts<-lapply(seq_len(nrow(edgerans)),function(ii) seq(edgerans[ii,1],edgerans[ii,2],length.out=nts[ii]+1))
  ts
}

.map.increments<-function(tree,ts,nts,dts,nstates,states){
  ntrees<-length(tree)
  maps<-do.call(cbind,lapply(tree,'[[','maps'))
  single.state<-lengths(maps)==1
  foo2<-function(j){
    if(single.state[j]){
      out<-data.frame(ts=ts,dts=dt,state=names(map[j,][[1]]),incl=TRUE)
      out[['dt']][1]<-0
    }else{
      map<-maps[[j]]
      ts<-sort(c(ts,map))
      tmp<-rle(nzchar(names(ts)))
      tmp$values<-!tmp$values
      incl<-inverse.rle(tmp)
      inds<-cumsum(tmp$lengths)
      nn<-length(inds)+1
      inds<-matrix(c(1,inds+c(0,1))[-nn],nrow=2)
      for(k in seq_len((nn-1)/2)){
        names(ts)[seq.int(inds[1,k],inds[2,k])]<-names(ts)[inds[2,k]+1]
      }
      dups<-duplicated(ts)
      incl<-incl[!dups]
      ts<-ts[!dups]
      data.frame(ts=ts,dts=c(0,diff(ts)),state=names(ts),incl=incl)
    }
  }
  foo1<-function(i){
    ts<-ts[[i]]
    nt<-nts[i]
    dt<-dts[i]
    maps<-lapply(maps[i,],function(ii) cumsum(ii)+ts[1])
    single.state<-single.state[i,]
    out<-vector('list',ntrees)
    
    
    
    out<-matrix(list(NULL),nrow=ntrees,ncol=nt)
    if(any(single.state)){
      nn<-sum(single.state)
      which.states<-unlist(lapply(maps[single.state],names))
      out[single.state,]<-as.list(setNames(rep(dt,nn),which.states))
    }
    
    
    
    ts<-ts-ts[1]
    ts<-ts[-(nt+1)]
    tmp<-sort(c(maps,ts),method='quick')
    out<-matrix(0,nrow=nt,ncol=nstates)
    colnames(out)<-states
    nn<-length(tmp)
    nms<-names(tmp)
    cur.state<-nms[nn]
    cur.ind<-nn-1
    for(t in rev(seq_len(nt))){
      nm<-nms[cur.ind]
      if(nzchar(nm)){
        while(nzchar(nm)){
          out[t,cur.state]<-out[t,cur.state]+tmp[cur.ind+1]-tmp[cur.ind]
          cur.state<-nm
          cur.ind<-cur.ind-1
          nm<-nms[cur.ind]
        }
        out[t,cur.state]<-out[t,cur.state]+tmp[cur.ind+1]-tmp[cur.ind]
        cur.ind<-cur.ind-1
      }else{
        out[t,cur.state]<-dt
        cur.ind<-cur.ind-1
      }
    }
    out
  }
  lapply(seq_along(ts),foo) #a tad slow--but not enough to worry about yet
  #I think this is actually pretty good
}

.init.edges<-function(tree,ntraits,traits,
                      res,nsims,
                      nstates,states){
  ts<-.subdivide.edges(tree[[1]],res)
  nts<-lengths(ts)-1
  elen<-tree[[1]]$edge.length
  dts<-elen/nts
  maps<-.map.increments(tree,ts,nts,dts,nstates,states)
  
  
  dims<-c(ntraits,sum(nts),nsims)
  seed<-array(rnorm(prod(dims)),dims,
              list(trait=trait.names,time=NULL,sim=NULL))
  inds<-asplit(matrix(c(1,rep(cumsum(nts[-length(elen)]),each=2)-c(1,0),dims[1]),nrow=2),2)
  out<-seed<-lapply(inds,function(ii) seed[,seq.int(ii[1],ii[2]),,drop=FALSE])
}
