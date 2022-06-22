#' @export
prep.dat<-function(contsimmap,trait.data){
  
  ntrees<-max(contsimmap[['perm.inds']][,1])
  trees<-contsimmap[['tree']][contsimmap[['perm.inds']][,1]]
  nsims<-length(trees)
  states<-contsimmap:::.get.states(trees)
  S<-length(states)
  
  #basic info
  N<-length(trees[[1]]$tip.label)
  K<-ncol(trait.data)
  T<-nsims
  E<-sum(unlist(lapply(trees,function(ii) nrow(ii$edge)),use.names=FALSE))+T #include root edge...
  
  #trait info
  tip.labels<-trees[[1]]$tip.label
  obs<-lapply(tip.labels,function(ii) which(rownames(trait.data)==ii))
  n_obs<-lengths(obs)
  ind_obs<-unlist(obs,use.names=FALSE)
  Y<-t(trait.data[ind_obs,,drop=FALSE])
  
  #"flatten out" trees to get topological info
  #works for now, but need to modify to make sure tips remain correct even when ordering differs between trees
  #(just use tip labels for this)
  des<-vector('list',E)
  ind_root<-vector('integer',T)
  ind_tip<-matrix(0L,T,N)
  ind_prune<-vector('integer',E-N*T)
  counter<-0
  prune.counter<-0
  for(i in seq_along(trees)){
    tmp.des<-c(list(contsimmap:::root.edges(trees[[i]])),des.edges(trees[[i]]))
    n.tmp.des<-length(tmp.des)
    tmp.des<-lapply(tmp.des,function(ii) ii+counter+1)
    tmp.seq<-seq_len(n.tmp.des)+counter
    des[tmp.seq]<-tmp.des
    
    ind_root[i]<-counter+1
    
    ind_tip[i,]<-tip.edges(trees[[i]],include.names=TRUE)[tip.labels]+counter+1
    
    tmp.prune.seq<-c(reorder(trees[[i]],'postorder',TRUE)+counter+1,ind_root[i])
    tmp.prune.seq<-tmp.prune.seq[is.na(match(tmp.prune.seq,ind_tip[i,]))]
    n.tmp.prune<-length(tmp.prune.seq)
    ind_prune[seq_len(n.tmp.prune)+prune.counter]<-tmp.prune.seq
    
    counter<-counter+n.tmp.des
    prune.counter<-prune.counter+n.tmp.prune
  }
  n_des<-lengths(des)
  pos_des<-c(0,cumsum(n_des[-counter]))+1
  ind_des<-unlist(des,use.names=FALSE)
  
  #missing data stuff
  obs.inf.dim<-!is.na(Y)
  tip.inf.dim<-matrix(FALSE,K,N)
  nod.inf.dim<-matrix(FALSE,K,E)
  for(i in seq_len(N)){
    if(n_obs[i]){
      tmp.obs<-obs.inf.dim[,obs[[i]],drop=FALSE]
      tip.inf.dim[,i]<-Reduce('|',asplit(tmp.obs,2))
      nod.inf.dim[,ind_tip[,i]]<-tip.inf.dim[,i]
    }else{
      tip.inf.dim[,i]<-FALSE
      nod.inf.dim[,ind_tip[,i]]<-FALSE
    }
  }
  for(i in ind_prune){
    tmp.des<-nod.inf.dim[,des[[i]],drop=FALSE]
    nod.inf.dim[,i]<-Reduce('|',asplit(tmp.des,2))
  }
  obs.inf.dim[]<-as.numeric(obs.inf.dim)
  obs.codes<-do.call(paste,c(asplit(obs.inf.dim,1),list(sep='')))
  tip.inf.dim[]<-as.numeric(tip.inf.dim)
  tip.codes<-do.call(paste,c(asplit(tip.inf.dim,1),list(sep='')))
  nod.inf.dim[]<-as.numeric(nod.inf.dim)
  nod.codes<-do.call(paste,c(asplit(nod.inf.dim,1),list(sep='')))
  codes<-unique(c(obs.codes,tip.codes,nod.codes))
  C<-length(codes)
  code_obs<-match(obs.codes,codes)
  code_tip<-match(tip.codes,codes)
  code_nod<-match(nod.codes,codes)
  parsed.codes<-lapply(strsplit(codes,''),as.numeric)
  n_dim<-unlist(lapply(parsed.codes,sum),use.names=FALSE)
  pos_dim<-c(0,cumsum(n_dim)[-C])+1
  ind_dim<-unlist(lapply(parsed.codes,function(ii) which(ii==1)),use.names=FALSE)
  Y[is.na(Y)]<-0
  
  #increment stuff...
  #insert check for variables that don't match up!
  t<-s<-vector('list',ntrees)
  traits<-contsimmap[['traits']]
  z<-setNames(rep(list(vector('list',nsims)),length(traits)),traits)
  counter<-0
  for(i in seq_len(ntrees)){
    tmp.ord<-order(as.numeric(names(contsimmap[['x']])))
    tmp<-lapply(contsimmap[['maps']][tmp.ord],'[[',i)
    t[[i]]<-c(list(NULL),lapply(tmp,'[[','dts')) #NULL for root edge...
    s[[i]]<-lapply(tmp,'[[','state')
    for(j in seq_along(traits)){
      tmp<-matrix(unlist(lapply(contsimmap[['x']][tmp.ord],function(ii) asplit(ii[[i]][,j,,drop=FALSE],3)),
                         recursive=FALSE,
                         use.names=FALSE),
                  ncol=length(s[[i]]))
      if(j==1){
        tmp.n<-nrow(tmp)
        tmp.seq<-seq_len(tmp.n)
      }
      z[[j]][counter+tmp.seq]<-asplit(tmp,1)
    }
    counter<-tmp.n
  }
  t<-unlist(t[contsimmap[['perm.inds']][,1]],recursive=FALSE,use.names=FALSE)
  n_inc<-lengths(t)
  pos_inc<-c(0,cumsum(n_inc[-E]))+1
  t<-unlist(t,use.names=FALSE)
  s<-factor(unlist(s[contsimmap[['perm.inds']][,1]],use.names=FALSE),levels=states)
  z<-lapply(z,function(ii) unlist(ii[contsimmap[['perm.inds']][,3]],use.names=FALSE))
  z[['STATE']]<-s
  
  #final data list
  list('N'=N,'K'=K,'T'=T,'E'=E,'C'=C,
       'n_obs'=n_obs,'ind_obs'=ind_obs,'Y'=Y,
       'n_inc'=n_inc,'pos_inc'=pos_inc,'t'=t,
       'n_des'=n_des,'ind_des'=ind_des,'pos_des'=pos_des,
       'ind_tip'=ind_tip,'ind_root'=array(ind_root),'ind_prune'=ind_prune,
       'n_dim'=array(n_dim),'ind_dim'=array(ind_dim),'pos_dim'=array(pos_dim),
       'code_obs'=code_obs,'code_tip'=code_tip,'code_nod'=code_nod, 
       'z'=z)
  
}

#' @export
fit.model<-function(dat,beta=~1,alpha=~1,gamma=~1,
                    niter=100,ncores=1,naive.marg=TRUE,...){
  
  null.mat<-matrix(1,nrow=sum(dat[['n_inc']]),ncol=1)
  
  if(!is.null(alpha)){
    logit<-TRUE
    U<-model.matrix(alpha,data=dat[['z']])
    if(!length(U)){
      U<-null.mat
    }
    A<-ncol(U)
  }else{
    A<-0
    logit<-FALSE
  }
  
  W<-model.matrix(beta,data=dat[['z']])
  if(!length(W)){
    W<-null.mat
  }
  B<-ncol(W)
  
  if(!is.null(gamma)){
    Z<-model.matrix(gamma,data=dat[['z']])
    if(!length(Z)){
      Z<-null.mat
    }
  }else{
    Z<-null.mat[,-1,drop=FALSE]
  }
  G<-ncol(Z)
  
  if(logit){
    dat[['A']]<-A
    dat[['U']]<-U
  }
  dat[['B']]<-B
  dat[['W']]<-W
  dat[['G']]<-G
  dat[['Z']]<-Z
  dat[['z']]<-NULL
  
  out<-vector('list',niter)
  if(logit){
    if(naive.marg){
      mod<-stanmodels[['multirate_BM_logitlink_univar_nointravar_naive']]
    }else{
      mod<-stanmodels[['multirate_BM_logitlink_univar_nointravar']]
    }
  }else{
    if(naive.marg){
      mod<-stanmodels[['multirate_BM_loglink_univar_nointravar_naive']]
    }else{
      mod<-stanmodels[['multirate_BM_loglink_univar_nointravar']]
    }
  }
  try<-1
  cur<-list('return_code'=123,'value'=-Inf,'par'=rep(0,A+B+G))
  while(try<niter){
    prop<-optimizing(mod,dat=dat,...)
    if(prop[['return_code']]==0){
      if(prop[['value']]>cur[['value']]){
        cur<-prop
      }
    }
    try<-try+1
  }
  cur
  # if(ncores==1){
  #   # prog<-function(iter,tot,cur){
  #   #   prop<-floor(100*iter/tot)
  #   #   if(prop>cur){
  #   #     cur<-prop
  #   #     cat(paste0('|',paste(rep('=',cur),collapse=''),
  #   #                paste(rep(' ',100-cur),collapse=''),'| ',
  #   #                cur,'%',if(cur<100) '\r'))
  #   #   }
  #   #   cur
  #   # }
  #   # cur<-prog(0,niter,-1)
  #   for(i in seq_len(niter)){
  #     out[[i]]<-optimizing(mod,dat=dat,...)
  #     # cur<-prog(i,niter,cur)
  #   }
  #   cat('\n')
  # }else{
  #   tmp.args<-c(object=list(mod),dat=list(dat),list(...))
  #   cl<-parallel::makeCluster(ncores)
  #   #this seems faster than simple call to clusterExport
  #   tmp.file<-tempfile()
  #   saveRDS(tmp.args,tmp.file)
  #   parallel::clusterExport(cl,'tmp.file',envir=environment())
  #   parallel::clusterEvalQ(cl,{
  #     library(contsimmap)
  #     tmp.args<-readRDS(tmp.file)
  #     })
  #   out<-parallel::parLapply(cl,seq_len(niter),function(ii) do.call(optimizing,tmp.args))
  #   parallel::stopCluster(cl)
  # }
  # 
  # tmp<-lapply(out,'[[','value')
  # tmp[!lengths(tmp)]<-list(-Inf)
  # max.lik<-which.max(unlist(tmp,use.names=FALSE))
  # out[[max.lik]]
  
}