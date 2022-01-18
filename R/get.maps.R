.subdivide.edges<-function(tree,res){
  edgerans<-edge.ranges(tree)
  dt<-max(edgerans)/max(1,res)
  nts<-ceiling(tree$edge.length/dt)
  e.seq<-seq_len(nrow(edgerans))
  foo<-function(x){
    t1<-edgerans[x,1]
    t2<-edgerans[x,2]
    if(t2==t1){
      rep(t1,2)
    }else{
      seq(t1,t2,length.out=nts[x]+1)
    }
  }
  setNames(lapply(e.seq,foo),e.seq)
}

#might consider run length encoding to save on space!
.map.increments<-function(tree,ts,nts,dts){
  ntrees<-length(tree)
  maps<-do.call(cbind,lapply(tree,'[[','maps'))
  single.state<-lengths(maps)==1
  out<-vector('list',ntrees)
  foo<-function(i){
    ts<-ts[[i]]
    nt<-nts[i]
    dt<-dts[i]
    maps<-lapply(maps[i,],function(ii) cumsum(ii)+ts[1])
    single.state<-single.state[i,]
    if(any(single.state)){
      out[single.state]<-lapply(maps[single.state],
                                function(ii) list(ts=ts[-1], #nearly 10-fold increase from storing as lists rather than data.frames
                                                  dts=rep(dt,nt),
                                                  state=rep(names(ii),nt),
                                                  incl=rep(TRUE,nt)))
    }
    if(any(!single.state)){
      maps<-maps[!single.state]
      maps.seq<-seq_along(maps)
      nms<-rep(maps.seq,lengths(maps))
      ts<-sort(c(ts,unlist(maps)),index.return=TRUE)
      tmp.inds<-ts$ix<(nt+2)
      inds<-lapply(maps.seq,function(ii) ts$ix%in%(which(nms==ii)+nt+1)|tmp.inds)
      ts<-ts$x
      for(j in maps.seq){
        map<-maps[[j]]
        t<-ts[inds[[j]]]
        tlen<-length(t)
        #switch in case rounding error makes last point not a state switch point
        if(!nzchar(names(t)[tlen])){
          names(t)[tlen-c(0,1)]<-names(t)[tlen-c(1,0)]
        }
        state<-rle(names(t))
        incl<-state
        incl$values<-!nzchar(incl$values)
        state$values[incl$values]<-state$values[-1][incl$values]
        ####SOMETING WRONG HERE####
        #sometimes incl gets "clipped" to shorter than it should be???
        #more numerically stable than duplicated() which can thrown off by rounding errors...
        # dt<-t[-1]-t[-length(t)]
        # dups<-c(FALSE,dt<1e-14) #maybe allow users to specify tolerance...
        # dt<-dt[!dups[-1]]
        # t<-t[!dups]
        #less errors (and probs negligible increase in later computations) if you just get rid of the last point
        #which you know is equivalent to preceding--just keep the rest
        #maybe make a special flag for when regimes switch at every time point (e.g., contsimmap for evorate models)
        #turned out the issue was a rounding error causing the positions of the last two entries to be switched!
        #see if statement to fix this above
        #still think it's best to probably just not search for duplicates
        excl.inds<- -c(1,tlen)
        dt<-t[excl.inds]-t[excl.inds[2]+c(0,1)]
        t<-t[excl.inds]
        out[!single.state][[j]]<-list(ts=t,
                                      dts=dt,
                                      state=inverse.rle(state)[excl.inds],
                                      incl=inverse.rle(incl)[excl.inds])
      }
    }
    list(dim=unlist(lapply(out,function(ii) length(ii[['ts']]))),info=out)
  }
  lapply(seq_along(ts),foo)
}

.get.edge.maps<-function(tree,treeID,nsims,res){
  ts<-.subdivide.edges(tree[[1]],res)
  nts<-lengths(ts)-1
  elen<-tree[[1]]$edge.length
  dts<-elen/nts
  maps<-.map.increments(tree,ts,nts,dts)
  nts<-do.call(cbind,lapply(maps,'[[','dim')) #rows are trees, columns are edges
  sims.per.tree<-tabulate(treeID)
  seeds.per.tree.edge<-nts*sims.per.tree
  seeds.per.edge<-.colSums(seeds.per.tree.edge,length(tree),length(elen))
  maps<-lapply(maps,'[[','info')
  names(maps)<-names(ts)
  list(maps=maps,
       ts=ts,
       nts=nts,
       sims.per.tree=sims.per.tree,
       seeds.per.tree.edge=seeds.per.tree.edge,
       seeds.per.edge=seeds.per.edge)
}

#quick viz test
# tmp<-lapply(out,'[[',1)
# plot(0,ylim=range(unlist(tmp)),xlim=range(unlist(ts)))
# tmp.ts<-lapply(maps,function(ii) ii[[1]][['ts']])
# for(i in seq_along(tmp)){
#   if(length(tmp.ts[[i]])>1){
#     matplot(tmp[[i]][1,,c(1,2)],x=tmp.ts[[i]],type='l',lty=1,add=TRUE)
#   }
# }
#looks right

