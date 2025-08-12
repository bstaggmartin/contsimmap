.uncond.traversals<-function(prune.seq,anc,edge.mat,ntips,
                             maps,
                             X0,nobs,Xsig2,Ysig2,mu,lookup,
                             nts,seed,
                             Xsig2.mods=NULL,mu.mods=NULL,
                             verbose=FALSE){
  #same block of code as below, really-->undesirable redundancy
  edge.inds<-rownames(seed)
  tree.inds<-colnames(seed)
  ntrees<-length(tree.inds)
  traits<-rownames(X0)
  ntraits<-dim(seed[[1]])[2]
  diag.inds<-do.call(cbind,rep(list(seq_len(ntraits)),2))
  tmp<-c(0,seq_len(ntraits^2-1))
  diag.inds2<-tmp%%ntraits==tmp%/%ntraits
  sims.per.tree<-unlist(lapply(seed[1,],function(ii) dim(ii)[3]),use.names=FALSE)
  nsims<-sum(sims.per.tree)
  lookups.per.tree<-unlist(lapply(lookup,function(ii) nrow(ii[['table']])),use.names=FALSE)
  states<-colnames(Xsig2)
  if(!is.null(mu.mods)){
    mu.dims<-which(lengths(mu.mods)>0)
  }
  if(!is.null(Xsig2.mods)){
    Xsig2.dims<-which(lengths(Xsig2.mods)>0)
  }
  
  ####INITIALIZING ROOT VALUES####
  for(t in seq_len(ntrees)){
    seed[[1,t]][1,,]<-X0[,lookup[[t]][['nobs.X0']][,2],drop=FALSE]
  }
  
  ####INITIALIZING FUNCTIONS FOR PREORDER TRAVERSAL####
  get.nsim<-function(){
    sum(tmp.inds)
  }
  get.maps<-function(){
    maps[as.numeric(edge.inds[e]),as.numeric(tree.inds[t]),c('dts','state')]
  }
  get.nts<-function(){
    nts[e,t]
  }
  get.zs<-function(){
    aperm(sweep(seed[[e,t]],1,sqrt(tmp.maps[['dts']]),'*',check.margin=FALSE),c(1,3,2))
  }
  get.state.inds<-function(){
    out<-setNames(lapply(states,function(ii) tmp.maps[['state']]==ii),states)
    sums<-unlist(lapply(out,sum),use.names=FALSE)
    nz<-sums>0
    out<-out[nz]
    attr(out,'sums')<-setNames(sums[nz],names(out))
    out
  }
  get.inds<-function(){
    lookup[[t]][['inds']][l,]
  }
  #will have to revisit for scalars...
  trans.zs<-function(){
    matrix(tmp.z[tmp.state.inds[[s]],tmp.inds,,drop=FALSE],attr(tmp.state.inds,'sums')[s]*tmp.nsim,ntraits)%*%
      cholX[[lookup[[t]][['table']][l,2],s]]
  }
  if(is.null(mu.mods)){
    get.mu<-function(){
      sweep(aperm(array(unlist(mu[lookup[[t]][['table']][l,4],tmp.maps[['state']]],use.names=FALSE),c(ntraits,tmp.nts,tmp.nsim)),c(2,3,1)),
            1,tmp.maps[['dts']],'*',check.margin=FALSE)
    }
  }else{
    get.mu<-function(){
      tmp.mu<-aperm(array(unlist(mu[lookup[[t]][['table']][l,4],tmp.maps[['state']]],use.names=FALSE),c(ntraits,tmp.nts,tmp.nsim)),
                    c(2,3,1))
      for(i in mu.dims){
        tmp.mu[,,i]<-matrix(tmp.mu[,,i,drop=FALSE],tmp.nts,tmp.nsim)+mu.mods[[i]][[e,t]][,tmp.inds,drop=FALSE]
      }
      sweep(tmp.mu,1,tmp.maps[['dts']],'*',check.margin=FALSE)
    }
  }
  if(is.null(Xsig2.mods)){
    add.mu<-function(){
      tmp.z[,tmp.inds,,drop=FALSE]+tmp.mu
    }
  }else{
    add.mu<-function(){
      for(i in Xsig2.dims){
        tmp.z[,tmp.inds,i]<-matrix(tmp.z[,tmp.inds,i,drop=FALSE],tmp.nts,tmp.nsim)*Xsig2.mods[[i]][[e,t]][,tmp.inds,drop=FALSE]
      }
      tmp.z[,tmp.inds,,drop=FALSE]+tmp.mu
    }
  }
  accumulate.z<-function(){
    tmp.z<-aperm(tmp.z,c(1,3,2))
    tmp.z[1,,]<-tmp.z[1,,,drop=FALSE]+seed[[anc[[e]],t]][nts[anc[[e]],t],,,drop=FALSE]
    apply(tmp.z,c(2,3),cumsum)
  }
  
  ####PREORDER TRAVERSAL#####
  cholX<-Xsig2
  cholX[]<-lapply(Xsig2,.pseudo.chol,k=ntraits,diag.inds=diag.inds)
  if(verbose){
    counter<-cur.prog<-0
    tot<-length(prune.seq)-1
    cat("\n\nDoing preorder (root-to-tips) traversal...\n")
    cat(.prog(cur.prog))
  }
  for(e in rev(prune.seq)[-1]){
    for(t in seq_len(ntrees)){
      tmp.maps<-get.maps()
      tmp.nts<-get.nts()
      tmp.z<-get.zs()
      tmp.state.inds<-get.state.inds()
      for(l in seq_len(lookups.per.tree[t])){
        tmp.inds<-get.inds()
        tmp.nsim<-get.nsim()
        for(s in names(tmp.state.inds)){
          tmp.z[tmp.state.inds[[s]],tmp.inds,]<-trans.zs()
        }
        tmp.mu<-get.mu()
        tmp.z[,tmp.inds,]<-add.mu()
      }
      seed[[e,t]][]<-accumulate.z()
    }
    if(verbose){
      counter<-counter+1
      prop.prog<-floor(100*counter/tot)
      if(prop.prog>cur.prog){
        cur.prog<-prop.prog
        cat(.prog(cur.prog))
      }
    }
  }
  
  ####GENERATING TRAIT DATA####
  #having a lot of trouble coming up with indexing system for intraspecific trait data due to its multiple sources of raggedness...
  #I think the below works! Just have to find some way of splitting it all out efficiently...
  #generating observation codes...
  if(verbose) cat("\n\nSimulating trait data...\n")
  nobs.codes<-unlist(lapply(lookup,function(ii) ii[['nobs.X0']][,1]),use.names=FALSE)
  tot.nobs<-colSums(nobs)
  tot<-sum(tabulate(nobs.codes,nbins=length(tot.nobs))*tot.nobs)
  if(tot){
    #generating indicators for different tips...
    tmp.nms<-lapply(seq_len(ncol(nobs)),function(ii) rep(rownames(nobs),nobs[,ii,drop=FALSE]))
    nms<-unique(unlist(tmp.nms,use.names=FALSE))
    nms.edges<-match(nms,rownames(nobs)[c(ntips+1,edge.mat[as.numeric(edge.inds[-1]),2])])
    nms.inds<-setNames(lapply(nms,function(ii) lapply(tmp.nms,'==',ii)),nms)
    #generating Ysig2 codes....
    Ysig2.codes<-unlist(lapply(lookup,function(ii) ii[['table']][,3][ii[['matches']]]),use.names=FALSE)
    n.Ysig2.codes<-max(Ysig2.codes)
    Ysig2.inds<-Ysig2.inds.TF<-vector("list",n.Ysig2.codes)
    for(i in seq_len(n.Ysig2.codes)){
      tmp<-Ysig2.codes==i
      Ysig2.inds[[i]]<-which(tmp)
      Ysig2.inds.TF[[i]]<-rep(tmp,tot.nobs[nobs.codes])
    }
    #generating seed..
    obs.seed<-matrix(rnorm(tot*ntraits),
                     tot,ntraits)
    #helper functions
    get.nms.inds<-function(){
      unlist(nms.inds[[i]][nobs.codes],use.names=FALSE)
    }
    get.x<-function(){
      t(do.call(cbind,lapply(seq_len(ntrees),function(ii) matrix(seed[[e,ii]][nts[e,ii],,,drop=FALSE],ntraits,sims.per.tree[ii]))))
    }
    get.rep.vec<-function(){
      rep(Ysig2.inds[[j]],nobs[nms[i],nobs.codes[Ysig2.inds[[j]]]])
    }
    cholY<-Ysig2[,nms,drop=FALSE]
    cholY[]<-lapply(cholY,.pseudo.chol,k=ntraits,diag.inds=diag.inds)
    #main loop
    for(i in seq_along(nms)){
      e<-nms.edges[i]
      tmp.nms.inds<-get.nms.inds()
      tmp.x<-get.x()
      for(j in seq_len(n.Ysig2.codes)){
        tmp.inds<-tmp.nms.inds&Ysig2.inds.TF[[j]]
        obs.seed[tmp.inds,]<-tmp.x[get.rep.vec(),,drop=FALSE]+obs.seed[tmp.inds,,drop=FALSE]%*%cholY[[j,i]]
      }
    }
    #split up into list, rather than trying to do anything fancy with arrays
    trait.data<-vector('list',nsims)
    traits<-colnames(Xsig2[[1]])
    counter<-1
    for(i in seq_along(trait.data)){
      tmp.code<-nobs.codes[i]
      tmp.tot<-tot.nobs[tmp.code]
      if(tmp.tot){
        trait.data[[i]]<-matrix(obs.seed[counter:(counter+tmp.tot-1),,drop=FALSE],
                                nrow=tmp.tot,ncol=ntraits,
                                dimnames=list(tmp.nms[[tmp.code]],traits))
        counter<-counter+tmp.tot
      }else{
        trait.data[[i]]<-matrix(nrow=0,ncol=ntraits,
                                dimnames=list(NULL,traits))
      }
    }
  }else{
    trait.data<-rep(list(matrix(nrow=0,ncol=ntraits,dimnames=list(NULL,traits))),nsims)
  }
  
  if(verbose) cat("\n\nDone!")
  list(seed,trait.data)
}

#This actually might not be necessary--unconditional traversal is honestly pretty fast
# .uncond.traversals.uni<-function(prune.seq,anc,edge.mat,ntips,
#                                  maps,
#                                  X0,nobs,Xsig2,Ysig2,mu,lookup,
#                                  nts,seed,
#                                  Xsig2.mods=NULL,mu.mods=NULL,
#                                  verbose=FALSE){
#   
# }

.cond.traversals<-function(prune.seq,anc,des,ndes,
                           maps,
                           parsed.obs,parsed.mis,nobs,Xsig2,Ysig2,mu,lookup,
                           nts,NTS,t1s,seed,x,v,dx,dv,
                           Xsig2.mods=NULL,mu.mods=NULL,
                           verbose=FALSE,
                           inf.const=1e10){
  edge.inds<-rownames(x)
  tree.inds<-colnames(x)
  ntrees<-length(tree.inds)
  ntraits<-dim(x[[1]])[2]
  diag.inds<-do.call(cbind,rep(list(seq_len(ntraits)),2))
  tmp<-c(0,seq_len(ntraits^2-1))
  diag.inds2<-tmp%%ntraits==tmp%/%ntraits
  sims.per.tree<-unlist(lapply(x[1,],function(ii) dim(ii)[3]),use.names=FALSE)
  nsims<-sum(sims.per.tree)
  lookups.per.tree<-unlist(lapply(lookup,function(ii) nrow(ii[['table']])),use.names=FALSE)
  states<-colnames(Xsig2)
  if(!is.null(mu.mods)){
    mu.dims<-which(lengths(mu.mods)>0)
    base<-rep(FALSE,ntraits)
    foo<-function(i){
      base[i]<-TRUE
      base
    }
    mu.inds<-rep(list(NULL),ntraits)
    mu.inds[mu.dims]<-lapply(mu.dims,foo)
  }
  if(!is.null(Xsig2.mods)){
    Xsig2.dims<-which(lengths(Xsig2.mods)>0)
    base<-matrix(FALSE,ntraits,ntraits)
    foo1<-function(i){
      base[i,]<-TRUE
      base
    }
    foo2<-function(i){
      base[,i]<-TRUE
      base
    }
    Xsig2.inds2<-Xsig2.inds1<-rep(list(NULL),ntraits)
    Xsig2.inds1[Xsig2.dims]<-lapply(Xsig2.dims,foo1)
    Xsig2.inds2[Xsig2.dims]<-lapply(Xsig2.dims,foo2)
  }
  
  ####INITIALIZING FUNCTIONS FOR POSTORDER TRAVERSAL####
  #assumes potential scalars --> could be made more efficient/faster by assuming 1 covariance matrix per lookup
  get.pars<-function(){
    lookup[[t]][['table']][l,]
  }
  get.inds<-function(){
    lookup[[t]][['inds']][l,]
  }
  get.nsim<-function(){
    sum(tmp.inds)
  }
  get.nobs<-function(){
    nobs[e,tmp.pars[1]]
  }
  get.obs.x<-function(){
    parsed.obs[[e,tmp.pars[1]]]
  }
  get.obs.p<-function(){
    obs.p<-array(Ysig2[[tmp.pars[3],e]],c(ntraits,ntraits,tmp.nobs))
    obs.p[diag.inds2][parsed.mis[[e,tmp.pars[1]]]]<-Inf
    .solve(obs.p,tmp.nobs,ntraits,diag.inds)
  }
  get.x<-function(){
    out<-lapply(des[[e]],function(ee) matrix(x[[ee,t]][t1s[ee,t],,tmp.inds,drop=FALSE]-dx[[ee,t]][1,,tmp.inds,drop=FALSE],
                                             ntraits,tmp.nsim))
    if(!is.null(tmp.x)){
      out<-c(out,split(tmp.x,seq_len(tmp.nobs)))
    }
    out
  }
  get.p<-function(){
    out<-lapply(des[[e]],function(ee) .solve(array(v[[ee,t]][,,1,tmp.inds,drop=FALSE]+dv[[ee,t]][,,1,tmp.inds,drop=FALSE],
                                                   c(ntraits,ntraits,tmp.nsim)),
                                             tmp.nsim,ntraits,diag.inds))
    if(!is.null(tmp.p)){
      out<-c(out,split(tmp.p,rep(seq_len(tmp.nobs),each=ntraits^2)))
    }
    out
  }
  calc.v<-function(){
    .solve(Reduce('+',tmp.p),tmp.nsim,ntraits,diag.inds)
  }
  calc.x<-function(){
    tmp.p<-.resolve.infs.ls(tmp.p,tmp.nsim,tmp.nobs+tmp.ndes,ntraits,diag.inds2,inf.const=inf.const)
    .multAb(.solve(Reduce('+',tmp.p),tmp.nsim,ntraits,diag.inds,z2z=TRUE),
            Reduce('+',lapply(seq_len(tmp.nobs+tmp.ndes),function(ii) .multAb(tmp.p[[ii]],tmp.x[[ii]],tmp.nsim,ntraits))),
            tmp.nsim,ntraits)
  }
  calc.obs.v<-function(){
    .solve(.sum3d(tmp.p,tmp.nobs),1,ntraits,diag.inds)
  }
  calc.obs.x<-function(){
    tmp.p<-.resolve.infs(tmp.p,tmp.nobs,ntraits,diag.inds2,inf.const=inf.const)
    .colSums(.multbA(tmp.x,tmp.p,tmp.nobs,ntraits),tmp.nobs,ntraits)%*%matrix(.solve(.sum3d(tmp.p,tmp.nobs),1,ntraits,diag.inds,z2z=TRUE),ntraits,ntraits)
  }
  get.look<-function(){
    lookup[[t]][['matches']]
  }
  if(is.null(mu.mods)&is.null(Xsig2.mods)){
    get.maps<-function(){
      maps[as.numeric(edge.inds[e]),as.numeric(tree.inds[t]),c('coarse','inds')]
    }
    get.dts<-get.tpts<-function(){
      
    }
  }else{
    get.maps<-function(){
      maps[as.numeric(edge.inds[e]),as.numeric(tree.inds[t]),c('coarse','inds','dts')]
    }
    get.tpts<-function(){
      if(k>1){
        tmp<-tmp.maps[['inds']][k-1]+1
      }else{
        tmp<-1
      }
      if(tmp.maps[['inds']][k]>tmp){
        tmp:tmp.maps[['inds']][k]
      }else{
        tmp
      }
    }
    get.dts<-function(){
      tmp.maps[['dts']][tmp.tpts]/sum(tmp.maps[['dts']][tmp.tpts])
    }
  }
  get.NTS<-function(){
    NTS[e,t]
  }
  if(is.null(mu.mods)){
    get.dx<-function(){
      tmp.maps[['coarse']][k]*unlist(mu[lookup[[t]][['table']][,4],names(tmp.maps[['coarse']])[k]][tmp.look],use.names=FALSE)
    }
  }else{
    get.dx<-function(){
      tmp.mu<-unlist(mu[lookup[[t]][['table']][,4],names(tmp.maps[['coarse']])[k]][tmp.look],use.names=FALSE)
      for(i in mu.dims){
        tmp.mu[mu.inds[[i]]]<-tmp.mu[mu.inds[[i]]]+
          .colSums(mu.mods[[i]][[e,t]][tmp.tpts,,drop=FALSE]*tmp.dts,tmp.nts,sims.per.tree[t])
      }
      tmp.maps[['coarse']][k]*tmp.mu
    }
  }
  if(is.null(Xsig2.mods)){
    get.dv<-function(){
      tmp.maps[['coarse']][k]*unlist(Xsig2[lookup[[t]][['table']][,2],names(tmp.maps[['coarse']])[k]][tmp.look],use.names=FALSE)
    }
  }else{
    get.dv<-function(){
      tmp.Xsig2<-unlist(Xsig2[lookup[[t]][['table']][,2],names(tmp.maps[['coarse']])[k]][tmp.look],use.names=FALSE)
      for(i in Xsig2.dims){
        tmp.scalars<-rep(.colSums(Xsig2.mods[[i]][[e,t]][tmp.tpts,,drop=FALSE]*tmp.dts,tmp.nts,sims.per.tree[t]),each=ntraits)
        tmp.Xsig2[Xsig2.inds1[[i]]]<-tmp.Xsig2[Xsig2.inds1[[i]]]*tmp.scalars
        tmp.Xsig2[Xsig2.inds2[[i]]]<-tmp.Xsig2[Xsig2.inds2[[i]]]*tmp.scalars
      }
      tmp.maps[['coarse']][k]*tmp.Xsig2
    }
  }
  update.x<-function(){
    x[[e,t]][tmp.maps[['inds']][k],,,drop=FALSE]-dx[[e,t]][k,,,drop=FALSE]
  }
  update.v<-function(){
    v[[e,t]][,,k,,drop=FALSE]+dv[[e,t]][,,k,,drop=FALSE]
  }
  
  ####POSTORDER TRAVERSAL####
  if(verbose){
    counter<-cur.prog<-0
    tot<-length(prune.seq)
    cat("\n\nDoing postorder (tips-to-root) traversal...\n")
    cat(.prog(cur.prog))
  }
  for(e in prune.seq){
    not.root<-e>1
    tmp.ndes<-ndes[e]
    tmp.has.des<-tmp.ndes>0
    for(t in seq_len(ntrees)){
      for(l in seq_len(lookups.per.tree[t])){
        tmp.pars<-get.pars()
        tmp.inds<-get.inds()
        tmp.nsim<-get.nsim()
        tmp.nobs<-get.nobs()
        tmp.x<-NULL
        tmp.p<-NULL
        if(tmp.nobs>0){
          tmp.x<-get.obs.x()
          tmp.p<-get.obs.p()
        }
        if(tmp.has.des){
          tmp.x<-get.x()
          tmp.p<-get.p()
          v[[e,t]][,,NTS[e,t],tmp.inds]<-calc.v()
          x[[e,t]][nts[e,t],,tmp.inds]<-calc.x()
        }else if(tmp.nobs>0){
          v[[e,t]][,,NTS[e,t],tmp.inds]<-calc.obs.v()
          x[[e,t]][nts[e,t],,tmp.inds]<-calc.obs.x()
        }
      }
      if(not.root){
        tmp.look<-get.look()
        tmp.maps<-get.maps()
        for(k in get.NTS():1){
          tmp.tpts<-get.tpts()
          tmp.dts<-get.dts()
          tmp.nts<-length(tmp.tpts)
          dx[[e,t]][k,,]<-get.dx()
          dv[[e,t]][,,k,]<-get.dv()
          if(k>1){
            x[[e,t]][tmp.maps[['inds']][k-1],,]<-update.x()
            v[[e,t]][,,k-1,]<-update.v()
          }
        }
      }
    }
    if(verbose){
      counter<-counter+1
      prop.prog<-floor(100*counter/tot)
      if(prop.prog>cur.prog){
        cur.prog<-prop.prog
        cat(.prog(cur.prog))
      }
    }
  }
  
  ####SIMULATING ROOT VALUES####
  #what to do when elements of v[1,] are infinite?
  #default to sampling from REALLY wide normal distribution, I guess...
  #maybe allow users to specify how wide this distribution is down the line...
  #or just tell folks to set a prior on the root state accordingly, honestly
  tmp.v<-array(unlist(v[1,],use.names=FALSE),c(ntraits,ntraits,nsims))
  chol.v<-aperm(.chol(tmp.v,nsims,ntraits,diag.inds),
                c(2,1,3))
  zz<-matrix(unlist(seed[1,],use.names=FALSE),c(ntraits,nsims))
  xx<-matrix(unlist(x[1,],use.names=FALSE),c(ntraits,nsims))+.multAb(chol.v,zz,nsims,ntraits)
  infs<-is.infinite(tmp.v[diag.inds2])
  if(any(infs)){
    xx[infs]<-runif(sum(infs),-1e9,1e9)
  }
  counter<-1
  for(t in seq_len(ntrees)){
    x[[1,t]][1,,]<-xx[,counter:(counter+sims.per.tree[t]-1)]
    counter<-counter+sims.per.tree[t]
  }
  
  #NEED TO ADD MOD STUFF DOWN HERE, BUT THE ABOVE SHOULD WORK!
  ####INITIALIZING FUNCTIONS FOR PREORDER TRAVERSAL####
  get.nsim<-function(){
    sims.per.tree[t]
  }
  if(is.null(Xsig2.mods)){
    get.maps<-function(){
      maps[as.numeric(edge.inds[e]),as.numeric(tree.inds[t]),c('coarse','inds','bb.sds','bb.dts')]
    }
  }else{
    get.maps<-function(){
      maps[as.numeric(edge.inds[e]),as.numeric(tree.inds[t]),c('coarse','inds','dts')]
    }
  }
  get.x<-function(){
    if(k>1){
      tmp<-x[[e,t]][tmp.maps[['inds']][k-1],,,drop=FALSE]
    }else{
      tmp<-x[[anc[[e]],t]][nts[anc[[e]],t],,,drop=FALSE]
    }
    list(matrix(tmp+dx[[e,t]][k,,,drop=FALSE],ntraits,tmp.nsim),
         matrix(x[[e,t]][tmp.maps[['inds']][k],,,drop=FALSE],ntraits,tmp.nsim))
  }
  get.p<-function(){
    list(.solve(array(dv[[e,t]][,,k,,drop=FALSE],c(ntraits,ntraits,tmp.nsim)),tmp.nsim,ntraits,diag.inds),
         .solve(array(v[[e,t]][,,k,,drop=FALSE],c(ntraits,ntraits,tmp.nsim)),tmp.nsim,ntraits,diag.inds))
  }
  get.z<-function(){
    matrix(seed[[e,t]][tmp.maps[['inds']][k],,,drop=FALSE],ntraits,tmp.nsim)
  }
  #8/10/25: I finally get it--something goes wrong specifically when you have an exact/missing measurement for the same tip
  #the covariances get set to 0 and are thus not taken into account...
  #how to solve this??? some kind of limit?
  #I have to think about this...jeez
  #It's interesting...the covariance does truly seem to approach 0 as the appropriate entry --> Inf
  #I mean, I guess that makes sense...it's not like the z bit should include any variance...
  #So it's more about calculating the distribution's proper mean, no?
  #Of course! It ends up centered around 0! This is the issue!!!
  #Yeah, you just need a smarter way to resolve the infinities such that you don't simply eliminate the covariances...
  #But how???
  #A hacky approach is to just set any infinities to arbitrarily high numbers...but is there a better approach?
  
  # .multAb(.solve(tmp.p[[1]]+tmp.p[[2]],tmp.nsim,ntraits,diag.inds), #this part should theoretically be okay?
  #         .multAb(tmp.p[[1]],tmp.x[[1]],tmp.nsim,ntraits)+.multAb(tmp.p[[2]],tmp.x[[2]],tmp.nsim,ntraits),
  #         #^it's the above that's problematic, I think...
  #         #yeah, the issue is here--you would need a way to store Inf+a kinds of values...
  #         #THAT's the problem!
  #         #So, in the end, just replacing Inf's with arbitrarily large values may be the only way...
  #         #But like, ugh...is that really it?
  #         tmp.nsim,ntraits)
  # tmp1<-.solve(tmp.p[[1]]+tmp.p[[2]],tmp.nsim,ntraits,diag.inds)
  # tmp2<-.multAb(tmp.p[[1]],tmp.x[[1]],tmp.nsim,ntraits)+.multAb(tmp.p[[2]],tmp.x[[2]],tmp.nsim,ntraits)
  # #wait...
  # tmp1[tmp1==0]<- -tmp.x[[2]][1,1]
  # tmp2[is.infinite(tmp2)]<-1
  # .multAb(tmp1,tmp2,
  #         tmp.nsim,ntraits)
  # #nah, doesn't quite work...
  # #you seem to lose critical sign information...
  # tmp1[1,,1]*tmp2[,1] #you want it to be the exact tip value + 0...
  # tmp1[2,,1]*tmp2[,1] #you want it to be the exact tip value(*some kind of covariance?) + normal value...
  # #so just...ensure the 0*Infs cancel out to og value...somehow
  # tmp1[1,1,]<-1
  # tmp1[2,1,]<- -cov2cor(tmp.p[[1]][,,1])[2,1] #maybe? feels wrong...yeah, not quite right...
  # tmp2[1,]<-tmp.x[[2]][1,1]
  # tmp1[,,1]%*%tmp2[,1]
  # #this all does seem to suggest that there might be a way to modify the infinity resolution procedure to be more robust when there's covariance...
  # #I think the key is you have to divide out stuff...you can't just set things to 0...
  # #But I just can't figure out how to do this in a general way...
  # 
  # #wait a minute!
  # tmp.p<-old.p
  # # tmp.p[[1]][1,1,]<-0
  # tmp.p[[2]][1,1,]<-tmp.p[[1]][1,1,]
  # #oof, no, that turned out terribly...
  # #damn, no free lunch!
  # #probably will require just setting inf's to arbitrarily large value...sigh...
  # #all that work on resolving infinities for nothing! whatcha gonna do...
  # #it at least only affects the preorder traversal in cases where you have exact measurements, so far as I can discern...
  # #basically just have to make sure precision for Inf entry is really high compared to others!
  # #the only other thing I can think to do is set precision for non-Inf entries arbitrarily low...
  # #dividing whatever factor used through the off-diagonals as well!
  # #maybe that's the more consistent thing to do...
  # #Interestingly, that doesn't seem to work either!
  # #So no substitute...you gotta just...add infinity somehow!
  # #Ugh
  # tmp.p<-old.p
  # fac<-sqrt(1e-10/tmp.p[[1]][1,1,])
  # tmp.p[[1]][,1,]<-tmp.p[[1]][,1,]*rep(fac,each=2)
  # tmp.p[[1]][1,,]<-tmp.p[[1]][1,,]*fac
  # tmp.p[[2]][1,1,]<-1
  
  
  # (p1*x1+p2*x2)/(p1+p2) --> can this be rearranged in a better way?
  #(p1+p2)^-1*p1*x1 + (p1+p2)^-1*p2*x2
  #solve(p1+p2,p1*x1)+solve(p1+p2,p2*x2)
  #let's think... (p1+p2)y=p2*x2
  #p1*y+p2*y=p2*x2
  #p1*y+p2*y-p2*x2=0
  #p1*y+p2*(y-x2)=0
  #by extension, p1*(z-x1)+p2*z=0
  #p1*y=p2*(y-x2)
  #y=p1^-1*p2*(y-x2)
  
  #8/12:
  #I really tried to figure out another way
  #In the end, I can't really see a better way than to just replace Infs with arbitrarily high precision...
  #I just modified the resolve.infs functions to take care of this
  #Renders most of their arguments defunct, but keep things as is for convenience
  
  calc.x<-function(){
    chol.v<-aperm(.chol.solve(tmp.p[[1]]+tmp.p[[2]],tmp.nsim,ntraits,diag.inds),c(2,1,3))
    tmp.p<-.resolve.infs.ls(tmp.p,tmp.nsim,2,ntraits,diag.inds2,precedence=TRUE,inf.const=inf.const)
    .multAb(.solve(tmp.p[[1]]+tmp.p[[2]],tmp.nsim,ntraits,diag.inds),
            .multAb(tmp.p[[1]],tmp.x[[1]],tmp.nsim,ntraits)+.multAb(tmp.p[[2]],tmp.x[[2]],tmp.nsim,ntraits),
            tmp.nsim,ntraits)+
      .multAb(chol.v,tmp.z,tmp.nsim,ntraits)
  }
  get.tpts<-function(){
    if(k>1){
      tmp<-tmp.maps[['inds']][k-1]+1
    }else{
      tmp<-1
    }
    if(tmp.maps[['inds']][k]>tmp){
      tmp:(tmp.maps[['inds']][k]-1)
    }else{
      integer(0)
    }
  }
  if(is.null(Xsig2.mods)){
    get.zs<-function(){
      tmp.z<-matrix(aperm(sweep(seed[[e,t]][tmp.tpts,,,drop=FALSE],1,tmp.maps[['bb.sds']][tmp.tpts],'*',check.margin=FALSE),c(3,1,2)),
                    tmp.nts*tmp.nsim,ntraits)
      for(l in seq_len(lookups.per.tree[t])){
        tmp.inds<-lookup[[t]][['inds']][l,]
        tmp.z[tmp.inds,]<-tmp.z[tmp.inds,,drop=FALSE]%*%cholX[[lookup[[t]][['table']][l,2],names(tmp.maps[['coarse']])[k]]]
      }
      aperm(array(tmp.z,c(tmp.nsim,tmp.nts,ntraits)),c(2,3,1))
    }
  }else{
    #I wonder if this could be made more efficient...but seems to do job for now
    #Could make everything a lot smoother if rate matrix is scaled entirely, but this will only affect so many cases...
    get.zs<-function(){
      tmp.dim<-tmp.nts*tmp.nsim
      tmp.z<-matrix(aperm(seed[[e,t]][tmp.tpts,,,drop=FALSE],c(1,3,2)),
                    tmp.dim,ntraits)
      tmp.tmp.tpts<-c(tmp.tpts,tmp.tpts[tmp.nts]+1)
      # tmp.p<-array(1,c(ntraits,ntraits,tmp.nts+1,tmp.nsim))
      tmp.p<-aperm(array(tmp.maps[['dts']][tmp.tmp.tpts],c(tmp.nts+1,ntraits,ntraits,tmp.nsim)),c(2,3,1,4))
      for(i in Xsig2.dims){
        #maybe this extraction could be made more efficient somehow? Hard to say...
        tmp.scalars<-rep(as.vector(Xsig2.mods[[i]][[e,t]][tmp.tmp.tpts,,drop=FALSE]),each=ntraits)
        tmp.p[,i,,]<-tmp.p[,i,,,drop=FALSE]*tmp.scalars
        tmp.p[i,,,]<-tmp.p[i,,,,drop=FALSE]*tmp.scalars
      }
      tmp.p<-rep(list(tmp.p),2)
      for(i in (tmp.nts:1)[-tmp.nts]){
        tmp.p[[2]][,,i,]<-tmp.p[[2]][,,i+1,,drop=FALSE]+tmp.p[[2]][,,i,,drop=FALSE]
      }
      tmp.p[[1]]<-tmp.p[[1]][,,-(tmp.nts+1),,drop=FALSE]
      tmp.p[[2]]<-tmp.p[[2]][,,-1,,drop=FALSE]
      for(l in seq_len(lookups.per.tree[t])){
        tmp.inds<-lookup[[t]][['inds']][l,]
        tmp.invX<-as.vector(invX[[lookup[[t]][['table']][l,2],names(tmp.maps[['coarse']])[k]]])
        tmp.p[[1]][,,,tmp.inds]<-1/tmp.p[[1]][,,,tmp.inds,drop=FALSE]*tmp.invX
        tmp.p[[2]][,,,tmp.inds]<-1/tmp.p[[2]][,,,tmp.inds,drop=FALSE]*tmp.invX
      }
      tmp.p<-lapply(tmp.p,array,dim=c(ntraits,ntraits,tmp.dim))
      chol.v<-.chol.solve(tmp.p[[1]]+tmp.p[[2]],tmp.dim,ntraits,diag.inds)
      tmp.p<-.resolve.infs.ls(tmp.p,tmp.dim,2,ntraits,diag.inds2,precedence=TRUE,inf.const=inf.const)
      #tmp.z includes a lot more components if Xsig2 modifiers are present!
      #need precision matrix info for each time step!
      c(list(asplit(aperm(array(.multbA(tmp.z,chol.v,tmp.dim,ntraits),c(tmp.nts,tmp.nsim,ntraits)),c(3,1,2)),2)),
        lapply(c(tmp.p,list(.solve(tmp.p[[1]]+tmp.p[[2]],tmp.dim,ntraits,diag.inds))),
               function(ii) asplit(array(ii,c(ntraits,ntraits,tmp.nts,tmp.nsim)),3)))
    }
  }
  if(is.null(Xsig2.mods)){
    get.next.x<-function(){
      x[[e,t]][tmp.maps[['inds']][k],,,drop=FALSE]
    }
    get.cur<-function(){
      if(k==1){
        x[[anc[[e]],t]][nts[anc[[e]],t],,,drop=FALSE]
      }else{
        x[[e,t]][tmp.maps[['inds']][k-1],,,drop=FALSE]
      }
    }
    update.cur<-function(){
      cur+(tmp.x-cur)*tmp.maps[['bb.dts']][tmp.tpts][m]+tmp.z[m,,,drop=FALSE]
    }
  }else{
    get.next.x<-function(){
      matrix(x[[e,t]][tmp.maps[['inds']][k],,,drop=FALSE],ntraits,tmp.nsim)
    }
    get.cur<-function(){
      if(k==1){
        matrix(x[[anc[[e]],t]][nts[anc[[e]],t],,,drop=FALSE],ntraits,tmp.nsim)
      }else{
        matrix(x[[e,t]][tmp.maps[['inds']][k-1],,,drop=FALSE],ntraits,tmp.nsim)
      }
    }
    update.cur<-function(){
      .multAb(tmp.z[[4]][[m]],
              .multAb(tmp.z[[2]][[m]],cur,tmp.nsim,ntraits)+
                .multAb(tmp.z[[3]][[m]],tmp.x,tmp.nsim,ntraits),
              tmp.nsim,ntraits)+
        tmp.z[[1]][[m]]
    }
  }
  
  ####PREORDER TRAVERSAL####
  invX<-cholX<-Xsig2
  cholX[]<-lapply(Xsig2,.pseudo.chol,k=ntraits,diag.inds=diag.inds)
  invX[]<-lapply(Xsig2,.pseudo.solve,k=ntraits,diag.inds=diag.inds)
  if(verbose){
    counter<-cur.prog<-0
    tot<-length(prune.seq)-1
    cat("\n\nDoing preorder (root-to-tips) traversal...\n")
    cat(.prog(cur.prog))
  }
  for(e in rev(prune.seq)[-1]){
    for(t in seq_len(ntrees)){
      tmp.nsim<-get.nsim()
      tmp.maps<-get.maps()
      for(k in seq_len(get.NTS())){
        tmp.x<-get.x()
        tmp.p<-get.p()
        tmp.z<-get.z()
        x[[e,t]][tmp.maps[['inds']][k],,]<-calc.x()
        tmp.tpts<-get.tpts()
        tmp.nts<-length(tmp.tpts)
        if(tmp.nts){
          tmp.z<-get.zs()
          tmp.x<-get.next.x()
          cur<-get.cur()
          for(m in seq_len(tmp.nts)){
            x[[e,t]][tmp.tpts[m],,]<-cur<-update.cur()
          }
        }
      }
    }
    if(verbose){
      counter<-counter+1
      prop.prog<-floor(100*counter/tot)
      if(prop.prog>cur.prog){
        cur.prog<-prop.prog
        cat(.prog(cur.prog))
      }
    }
  }
  
  if(verbose) cat("\n\nDone!")
  x
}

# .cond.traversals.uni<-function(prune.seq,anc,des,ndes,
#                                maps,
#                                parsed.obs,parsed.mis,nobs,Xsig2,Ysig2,mu,lookup,
#                                nts,NTS,t1s,seed,x,v,dx,dv,
#                                Xsig2.mods=NULL,mu.mods=NULL,
#                                verbose=FALSE){
#   edge.inds<-rownames(x)
#   tree.inds<-colnames(x)
#   ntrees<-length(tree.inds)
#   sims.per.tree<-unlist(lapply(x[1,],function(ii) dim(ii)[3]),use.names=FALSE)
#   nsims<-sum(sims.per.tree)
#   states<-colnames(Xsig2)
#   sum.obs<-matrix(unlist(lapply(parsed.obs,sum),use.names=FALSE),
#                   nrow(parsed.obs),ncol(parsed.obs))
#   
#   #figure these guys out later...
#   # if(!is.null(mu.mods)){
#   #   mu.dims<-which(lengths(mu.mods)>0)
#   #   base<-rep(FALSE,ntraits)
#   #   foo<-function(i){
#   #     base[i]<-TRUE
#   #     base
#   #   }
#   #   mu.inds<-rep(list(NULL),ntraits)
#   #   mu.inds[mu.dims]<-lapply(mu.dims,foo)
#   # }
#   # if(!is.null(Xsig2.mods)){
#   #   Xsig2.dims<-which(lengths(Xsig2.mods)>0)
#   #   base<-matrix(FALSE,ntraits,ntraits)
#   #   foo1<-function(i){
#   #     base[i,]<-TRUE
#   #     base
#   #   }
#   #   foo2<-function(i){
#   #     base[,i]<-TRUE
#   #     base
#   #   }
#   #   Xsig2.inds2<-Xsig2.inds1<-rep(list(NULL),ntraits)
#   #   Xsig2.inds1[Xsig2.dims]<-lapply(Xsig2.dims,foo1)
#   #   Xsig2.inds2[Xsig2.dims]<-lapply(Xsig2.dims,foo2)
#   # }
#   
#   get.pars<-function(){
#     lookup[[t]][['table']]
#   }
#   get.nsim<-function(){
#     sims.per.tree[t]
#   }
#   get.matches<-function(){
#     lookup[[t]][['matches']]
#   }
#   get.obs.p<-function(){
#     tmp<-nobs[e,tmp.pars[,1]]/unlist(Ysig2[tmp.pars[,3],e],use.names=FALSE)
#     tmp[is.nan(tmp)]<-0
#     tmp
#   }
#   get.obs.x<-function(){
#     
#   }
#   
#   get.obs.v<-function(){
#     unlist(Ysig2[tmp.pars[,3],e],use.names=FALSE)
#   }
#   get.exact.inds<-function(){
#     tmp.obs.v==0
#   }
#   get.obs.p<-function(){
#     tmp<-nobs[e,tmp.pars[,1]]/unlist(Ysig2[tmp.pars[,3],e],use.names=FALSE)
#     tmp[is.nan(tmp)]<-0
#     tmp
#   }
#   get.obs.x<-function(){
#     infs<-is.infinite(tmp)
#     tmp[infs]<-nobs[e,]
#     tmp.obs.v[tmp.exact.inds]<-1
#     sum.obs[e,tmp.pars[,1]]/tmp.obs.v
#   }
#   
#   get.x<-function(){
#     out<-
#     out<-lapply(des[[e]],function(ee) matrix(x[[ee,t]][t1s[ee,t],,tmp.inds,drop=FALSE]-dx[[ee,t]][1,,tmp.inds,drop=FALSE],
#                                              ntraits,tmp.nsim))
#     if(!is.null(tmp.x)){
#       out<-c(out,split(tmp.x,seq_len(tmp.nobs)))
#     }
#     out
#   }
#   get.p<-function(){
#     out<-lapply(des[[e]],function(ee) .solve(array(v[[ee,t]][,,1,tmp.inds,drop=FALSE]+dv[[ee,t]][,,1,tmp.inds,drop=FALSE],
#                                                    c(ntraits,ntraits,tmp.nsim)),
#                                              tmp.nsim,ntraits,diag.inds))
#     if(!is.null(tmp.p)){
#       out<-c(out,split(tmp.p,rep(seq_len(tmp.nobs),each=ntraits^2)))
#     }
#     out
#   }
#   
#   for(e in prune.seq){
#     not.root<-e>1
#     tmp.ndes<-ndes[e]
#     tmp.has.des<-tmp.ndes>0
#     for(t in seq_len(ntrees)){
#       tmp.pars<-get.pars()
#       tmp.nsim<-get.nsim()
#       tmp.matches<-get.matches()
#       tmp.obs.v<-get.obs.v()
#       tmp.exact.inds<-get.exact.inds()
#       tmp.x<-get.obs.x()[tmp.matches]
#       tmp.p<-get.obs.p()[tmp.matches]
#       tmp.exact.inds<-tmp.exact.inds[tmp.matches]
#       if(tmp.has.des){
#         for(ee in des[[e]]){
#           
#         }
#       }
#       if(tmp.nobs>0){
#         tmp.x<-get.obs.x()
#         tmp.p<-get.obs.p()
#       }
#       if(tmp.ndes>0){
#         tmp.x<-get.x()
#         tmp.p<-get.p()
#         v[[e,t]][,,NTS[e,t],tmp.inds]<-calc.v()
#         x[[e,t]][nts[e,t],,tmp.inds]<-calc.x()
#       }else if(tmp.nobs>0){
#         v[[e,t]][,,NTS[e,t],tmp.inds]<-calc.obs.v()
#         x[[e,t]][nts[e,t],,tmp.inds]<-calc.obs.x()
#       }
#       if(not.root){
#         tmp.look<-get.look()
#         tmp.maps<-get.maps()
#         for(k in get.NTS():1){
#           tmp.tpts<-get.tpts()
#           tmp.dts<-get.dts()
#           tmp.nts<-length(tmp.tpts)
#           dx[[e,t]][k,,]<-get.dx()
#           dv[[e,t]][,,k,]<-get.dv()
#           if(k>1){
#             x[[e,t]][tmp.maps[['inds']][k-1],,]<-update.x()
#             v[[e,t]][,,k-1,]<-update.v()
#           }
#         }
#       }
#     }
#   }
# }