#Invert matrices including infinite and zero entries on diagonal
#Basically, 0s become Infs and Infs become 0s along diagonals, both have no covariance with other dimensions
.pseudo.solve<-function(mat){
  mat<-as.matrix(mat)
  diagonal<-diag(mat)
  infs<-which(is.infinite(diagonal))
  zeros<-which(diagonal==0)
  if(sum(length(infs),length(zeros))==0){
    solve(mat)
  }else{
    valids<-(1:nrow(mat))[-c(infs,zeros)]
    out<-mat
    out[,]<-0
    out[cbind(zeros,zeros)]<-Inf
    if(length(valids)>0){
      out[valids,valids]<-solve(mat[valids,valids])
    }
    out
  }
}

#Take cholesky composition of matrix including zero entries on diagonal (and infinites, I guess?)
#Returns expected result: t(out)%*%out returns original matrix (EXCEPT with infinites, but this shouldn't be
##issue in practice since they should never occur during Simmap construction)
#Keeping infinites cause proper propagation of Infs along diagonal, but also propagates NaNs during matrix
##multiplication, so should be avoided here if it can be helped
#The nice thing with the current construction is that both 0s and Infs cause diagonal entries to be 0 with 0
##correlation with other variables--this means that, in the context of multiplying uncorrelated variables to be
##correlated, dimensions with 0 and infinite variance will stay the same and not affect other dimensions; this
##also works in the case of infinite variance since we want these means to stay at 0
#Note that there are multiple decompostions for a matrices with 0 on diagonal in general (which is why normal
##Cholesky decompostion doesn't work); I avoid this problem here by making sure dimensions with 0 on the diagonal
##are filled up with 0s, ensuring a unique solution (well, unique according to what I read on Wikipedia)
.pseudo.chol<-function(mat){
  mat<-as.matrix(mat)
  diagonal<-diag(mat)
  infs<-which(is.infinite(diagonal))
  zeros<-which(diagonal==0)
  if(sum(length(infs),length(zeros))==0){
    chol(mat)
  }else{
    valids<-(1:nrow(mat))[-c(infs,zeros)]
    out<-mat
    out[,]<-0
    if(length(valids)>0){
      out[valids,valids]<-chol(mat[valids,valids])
    }
    out
  }
}

#Generate Brownian Motion bridge variance-covariance matrix--assumes unit scale (time scale/rate of 1)
#Can be multiplied to work for alternate time scales/rates
.bridge.cov<-function(n){
  pts<-seq(0,1,length.out=n+2)[-c(1,n+2)]
  out<-matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:i){
      out[i,j]<-(1-pts[i])*pts[j]
    }
  }
  tmp<-upper.tri(out)
  out[tmp]<-t(out)[tmp]
  out
}

#Linearly interpolate n equally spaced points between yy[,1,] and yy[,2,]
#Note that it returns the results as a vector
.lin.interp<-function(yy,n){
  y1<-as.vector(yy[,1,])
  y2<-as.vector(yy[,2,])
  xx<-rep(seq(0,1,length.out=n+2)[-c(1,n+2)],each=length(y1))
  c(y1,(y2-y1)*xx+y1,y2)
}

.get.samp<-function(samp,nsim){
  if(is.numeric(samp)&length(samp)>1){
    samp<-list(samp)
  }
  if(is.numeric(samp)){
    if(!(samp>0&samp<=nsim)){
      stop("can only sample between 1 and ",nsim," contsimmaps")
    }
    samp<-sample(nsim,samp)
  }else if(is.list(samp)){
    samp<-unlist(samp)
    if(!is.numeric(samp)){
      stop("invalid sample specified; please provide a single integer of contsimmaps to randomly sample or a list of indices")
    }
    out.of.range<-!(samp>0&samp<=nsim)
    if(all(out.of.range)){
      stop("all specified contsimmaps don't exist; only indices between 1 and ",nsim," are allowed")
    }
    if(any(out.of.range)){
      warning("some specified contsimmaps don't exist; only indices between 1 and ",nsim," are allowed")
      samp<-samp[!out.of.range]
    }
  }else{
    stop("invalid sample specified; please provide a single integer of contsimmaps to randomly sample or a list of indices")
  }
  samp
}

.get.trait<-function(trait,k,trait.names){
  if(length(trait)>1){
    trait<-trait[1:2]
  }
  if(is.numeric(trait)){
    out.of.range<-!(trait>0&trait<=k)
    if(all(out.of.range)){
      stop("all specified traits don't exist; only trait indices between 1 and ",k," are allowed")
    }
    if(any(out.of.range)){
      warning("some specified traits don't exist; only traits between 1 and ",k," are allowed")
      trait<-trait[!out.of.range]
    }
  }else if(is.character(trait)){
    out.of.range<-!grepl(trait,trait.names)
    if(all(out.of.range)){
      stop("all specified traits don't exist; trait names must match those in contsimmap exactly")
    }
    if(any(out.of.range)){
      warning("some specified traits don't exist; trait names must match those in contsimmap exactly")
      trait<-trait[!out.of.range]
    }
  }
  trait
}

.fix.mat<-function(mat,k,trait.names){
  #recycle from coercion functions in evorates
}

.subdivide.edges<-function(tree,res){
  edgerans<-evorates::edge.ranges(tree)
  dt<-max(edgerans)/res
  nts<-ceiling(tree$edge.length/dt)
  ts<-lapply(seq_len(nrow(edgerans)),function(ii) seq(edgerans[ii,1],edgerans[ii,2],length.out=nts[ii]+1))
  ts
}

.map.increments<-function(tree,ts,nts,dts,nstates,states){
  maps<-tree$maps
  foo<-function(i){
    maps<-cumsum(maps[[i]])
    ts<-ts[[i]]
    nt<-nts[i]
    dt<-dts[i]
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

.check.matching.tops<-function(tree){
  all(unlist(lapply(tree[-1],function(ii) all.equal(ii,tree[[1]]))))
}

.get.states<-function(tree){
  
}

.proc.tree<-function(tree){
  #checking to see if it's at least a phylo object
  if(!inherits(tree,'multiPhylo')&!inherits(tree,'phylo')){
    stop('tree must be of class phylo or simmap')
  }
  #checking to see if it includes multiple phylogenies, coercing it if it does
  if(!inherits(tree,'multiPhylo')){
    tree<-as.multiPhylo(tree)
    if(inherits(tree[[1]],'simmap')){
      class(tree)<-c(class(tree),'multiSimmap')
    }
  }else{
    if(!.check.matching.tops(tree)){
      stop('This function only accepts multiple phylogenies with identical topologies')
    }
  }
  #checking to see if it's a simmap or not
  if(!inherits(tree,'multiSimmap')){
    states<-'0'
    nstates<-1
    maps<-as.list(setNames(tree$edge.length,'0'))
    for(i in seq_along(tree)){
      tree$maps<-maps
    }
  }else{
    states<-unlist(lapply(tree,function(ii) unique(c(getStates(tree,'nodes'),getStates(tree,'tips')))))
    states<-sort(unique(states))
    nstates<-length(states)
  }
  treeID<-rep(seq_len(length(tree)),length.out=nsim)
}