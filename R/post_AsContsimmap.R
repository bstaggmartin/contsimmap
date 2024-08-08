#based on these initial tests, coming up with some kind of built-in "shortcut" for processing inputs might be good
#either based on things be univariate and/or some argument like check.inputs=FALSE

#' @export
as.contsimmap<-function(fit,res=100,tol=0.001,max.iter=1e5,verbose=FALSE){
  chains<-fit[['chains']]
  dims<-dim(chains)
  if(is.null(dims)){
    chains<-array(chains,c(length(chains),1,1),list(iterations=NULL,parameters=attr(chains,"parameters"),chains="chain 1"))
  }else if(length(dims)<3){
    chains<-array(chains,c(dims,1),c(dimnames(chains),chains="chain 1"))
  }
  R0<-as.list(chains[,'R_0',,drop=FALSE])
  nsims<-length(R0)
  Rsig2<-if(any(dimnames(chains)[[2]]=='R_sig2')) as.list(chains[,'R_sig2',,drop=FALSE]) else 0
  Rmu<-if(any(dimnames(chains)[[2]]=='R_mu')) as.list(chains[,'R_mu',,drop=FALSE]) else 0
  #set up "dummy" contsimmap
  list2env(.prep.contsimmap(tree=fit[['call']][['tree']],
                                         nsims=nsims,
                                         res=res,
                                         trait.data=R0,
                                         Xsig2=Rsig2,
                                         Ysig2=0,
                                         mu=Rmu,
                                         conditional=FALSE,
                                         ntraits=1,
                                         traits='R',
                                         nobs=0),envir=environment())
  #some initialization...
  nhgt<-c(0,edge.ranges(tree[[1]])[,2])
  if(length(Rsig2)==1&Rsig2[1]==0){
    seed[[1,1]][1,,]<-X0[,lookup[[1]][['nobs.X0']][,2],drop=FALSE]
    for(e in seq_len(dim(seed)[1])[-1]){
      tmp.dim1<-dim(seed[[e,1]])[1]
      tmp.ts<-ts[[e-1]]
      seed[[e,1]][]<-unlist(lapply(R0,rep,each=tmp.dim1),use.names=FALSE)
      if(!(length(Rmu)==1&Rmu[1]==0)){
        tmp.rmu<-unlist(lapply(Rmu,rep,each=tmp.dim1),use.names=FALSE)
        seed[[e,1]]<-seed[[e,1]]+log(abs(exp(tmp.rmu*tmp.ts[-1])-exp(tmp.rmu*tmp.ts[-(tmp.dim1+1)])))-log(abs(tmp.rmu))-log(diff(tmp.ts))
      }
    }
    scores<-NA
  }else{
    elen<-c(0,tree[[1]]$edge.length)
    Rsig2<-unlist(Rsig2,use.names=FALSE)
    Rsig<-sqrt(Rsig2)
    Rp<-1/Rsig2
    Rmu<-unlist(Rmu,use.names=FALSE)
    dt.Rs<-sweep(cbind(0,matrix(aperm(remove.trend(fit,simplify=FALSE),c(1,3,2)),nsims,Nedge(fit))),2,elen,'/',check.margin=FALSE)
    t.Rs<-sweep(cbind(0,matrix(aperm(get.R(fit,simplify=FALSE),c(1,3,2)),nsims,Nedge(fit))),2,log(elen),'+',check.margin=FALSE)
    
    #helper functions...
    get.ndes<-function(){
      ndes[e]
    }
    calc.v<-function(){
      1/(2/elen[e]+3*sum(1/elen[des[[e]]]))
    }
    calc.r<-function(){
      tmp.v*(2*dt.Rs[,e]+3*.rowSums(dt.Rs[,des[[e]],drop=FALSE],nsims,tmp.ndes))+Rmu*nhgt[e]
    }
    get.nt<-function(){
      nts[e,1]
    }
    get.dt<-function(){
      maps[[e-1,1,'dts']][1]
    }
    get.chol<-function(){
      tmp<-.col(c(tmp.nt,tmp.nt))
      ttmp<-t(tmp)
      tmp.inds<-tmp>ttmp
      tmp[tmp.inds]<-ttmp[tmp.inds]
      tmp.seq<-tmp.nt:1
      ttmp<-tmp[tmp.seq,tmp.seq]
      tmp.seq<-seq(tmp.dt,length.out=tmp.nt,by=tmp.dt)/sqrt(elen[e])
      chol(matrix(tmp.seq[tmp]*tmp.seq[ttmp],tmp.nt,tmp.nt)+tcrossprod(seq(0,sqrt(tmp.v),length.out=tmp.nt+1)[-1]))
    }
    get.mu<-function(){
      anc.r<-seed[[anc[[e]],1]][nts[anc[[e]],1],,]
      tcrossprod(tmp.r-anc.r,seq(0,1,length.out=tmp.nt+2)[-c(1,tmp.nt+2)])+anc.r
    }
    
    scores<-matrix(0,nsims,dim(seed)[1])
    
    #main traversal...
    #first do R0!
    seed[[1,1]][1,,]<-X0[,lookup[[1]][['nobs.X0']][,2],drop=FALSE]
    if(verbose){
      counter<-cur.prog<-0
      tot<-length(prune.seq)-1
      cat("Sampling rate maps...\n")
      cat(.prog(cur.prog))
    }
    for(e in rev(prune.seq)[-1]){
      #sample rates at edge's descendant node...
      tmp.ndes<-get.ndes()
      if(tmp.ndes){
        tmp.v<-calc.v()
        tmp.r<-calc.r()
      }else{
        tmp.v<-elen[e]/2
        tmp.r<-dt.Rs[,e]*elen[e]+Rmu*nhgt[e]
      }
      #"quick" MCMC?
      #likelihood is simple enough--just the sum of squares of the raw seed + squared deviation between resulting integral and target divided by some small sd, I think...
      #maybe use a standard normal proposal distribution with sd of 0.1ish?
      #pseudo-MCMC actually seems remarkably quick to converge and results in much better approximation
      #able to use MUCH lower tolerances and higher number of iterations!!!
      tmp.nt<-get.nt()
      tmp.dt<-get.dt()
      log.dt<-log(tmp.dt)
      tmp.chol<-get.chol()
      tmp.mu<-get.mu()
      targ<-t.Rs[,e]
      tol2<-tol^2
      cur.seed<-matrix(seed[[e,1]],nsims,tmp.nt,byrow=TRUE)
      tmp.rs<-tmp.mu+sweep(cur.seed%*%tmp.chol,1,Rsig,'*',check.margin=FALSE)
      tmp.score<-(log(.rowSums(exp(tmp.rs),nsims,tmp.nt))+log.dt-targ)^2/tol2
      cur.lik<- -0.5*(.rowSums(cur.seed^2,nsims,tmp.nt)+tmp.score)
      inds<-which(tmp.score>1)
      n.inds<-length(inds)
      for(i in seq_len(max.iter)){
        prop.seed<-cur.seed[inds,,drop=FALSE]+rnorm(n.inds*tmp.nt,0,0.1)
        tmp.rs<-tmp.mu[inds,,drop=FALSE]+sweep(prop.seed%*%tmp.chol,1,Rsig[inds],'*',check.margin=FALSE)
        tmp.score<-(log(.rowSums(exp(tmp.rs),n.inds,tmp.nt))+log.dt-targ[inds])^2/tol2
        prop.lik<- -0.5*(.rowSums(prop.seed^2,n.inds,tmp.nt)+tmp.score)
        tmp<-tmp.score>1
        new.inds<-inds[tmp]
        n.inds<-sum(tmp)
        tmp[tmp]<-(prop.lik[tmp]-cur.lik[inds[tmp]])< -rexp(n.inds)
        tmp<-!tmp
        cur.seed[inds[tmp],]<-prop.seed[tmp,,drop=FALSE]
        cur.lik[inds[tmp]]<-prop.lik[tmp]
        inds<-new.inds
        if(!n.inds) break
      }
      tmp.rs<-tmp.mu+sweep(cur.seed%*%tmp.chol,1,Rsig,'*',check.margin=FALSE)
      seed[[e,1]][]<-as.vector(t(tmp.rs))
      scores[,e]<-abs(log(.rowSums(exp(tmp.rs),nsims,tmp.nt))+log.dt-targ)
      if(verbose){
        counter<-counter+1
        prop.prog<-floor(100*counter/tot)
        if(prop.prog>cur.prog){
          cur.prog<-prop.prog
          cat(.prog(cur.prog))
        }
      }
    }
  }
  
  #format output
  x<-seed
  attr(x,'ts')<-ts
  attr(x,'tree')<-tree
  attr(x,'maps')<-maps
  attr(x,'treeID')<-treeID
  traits<-colnames(Xsig2[[1]])
  attr(x,'traits')<-setNames(rep(1,length(traits)),traits)
  attr(x,'params')<-matrix(list(NULL,Xsig2,Ysig2,mu,lookup,list(fxn="as.contsimmap.evorates_fit")),
                           6,1,
                           dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  out<-.uncompress(x)
  #now for the million-dollar test...
  Ysig2<-if(any(dimnames(chains)[[2]]=='Y_sig2')) as.list(chains[,'Y_sig2',,drop=FALSE]) else NULL
  tmp<-fit[['call']][['trait.se']]
  not.nas<-!is.na(tmp)&!is.infinite(tmp)
  Ysig2<-c(list(Ysig2),setNames(as.list(tmp[not.nas]),rownames(tmp[not.nas,,drop=FALSE])))
  dat<-fit[['call']][['trait.data']]
  if(length(Rsig2)==1&Rsig2[1]==0&length(Rmu)==1&Rmu[1]==0){
    Xsig2<-as.matrix(lapply(R0,exp))
    form<-as.formula(paste0(colnames(dat),"~diffusion(Xsig2=Xsig2,Ysig2=Ysig2,trait.data=dat",if(verbose) ",verbose=TRUE",")"))
  }else{
    form<-as.formula(paste0(colnames(dat),"~diffusion(",colnames(dat),"_Xsig2='exp_R',Ysig2=Ysig2,trait.data=dat",if(verbose) ",verbose=TRUE",")"))
  }
  make.traits(out,exp_R~exp(R),form)
}
