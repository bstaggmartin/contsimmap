rm(list=ls())
set.seed(123)
library(phytools)
library(contsimmap)
library(rstan)

##INSTEAD OF MAKING NULL IT'S OWN THING, JUST MAKE A NULL MODEL TO RUN!

prog<-function(iter,tot,cur){
  prop<-floor(100*iter/tot)
  if(prop>cur){
    cur<-prop
    cat(paste0('|',paste(rep('=',cur),collapse=''),
               paste(rep(' ',100-cur),collapse=''),'| ',
               cur,'%',if(cur<100) '\r'))
  }
  cur
}

####CONSTRUCTING SIMULATION SCENARIOS####

#sample size
ntax<-c(100)
#functional relationship type
func<-c('null','simple','sweetspot','threshold')
#related to focal trait or not?
foc<-c(FALSE,TRUE)
#noise or not?
noise<-c(FALSE,TRUE)
#number of replicates per scenario
nrep<-10
#make data.frame of scenarios
df<-expand.grid(ntax,func,foc,noise)
var.nms<-c('ntax','func','foc','noise')
colnames(df)<-var.nms
df<-df[!(!df$foc&func=='null'),] #eliminate null, non-focal simulations
nscen<-nrow(df)
tmp<-rep(seq_len(nscen),each=nrep)
df<-df[tmp,]
rownames(df)<-paste(tmp,seq_len(nrep),sep='_')
ntot<-nscen*nrep

####DOING SIMULATIONS####

#some things are constant across all simulations
#each simulation will be on a pure-birth tree with total height 1
#each simulation will consist of two "explanatory" traits, x and z, evolving via independent BM processes
#   these processes will always start at 0 and have rate 1
#   resolution of 1,000
#   for now, only 1 observation per species and no intraspecific variation
#   in all cases, x will be the observed explanatory trait
#each simulation will also have a single "response" trait, y, evolving with rate, r, dependent upon either x or z
#   when func is null: r = exp(e)*a
#   when func is simple: r = exp(e)*(a*exp(b*x)+c)
#   when func is sweetspot: r = exp(e)*(a*exp(-b/2*(x-d)^2)+c)
#   when func is threshold: r = exp(e)*(a/(1+exp(-b*(x-d)))+c)
#   if foc is TRUE, r will depend on x, otherwise it will depend on z
#   if noise is TRUE, e will be ~normal(0,0.5), otherwise it will be ~normal(0,0)
#   most helpful to tweak parameters individually; also note:
#     a must always be positive
#     b must always be positive in the "sweetspot" case
#     c must always be positive
#     d has no restrictions
#deciding on final parameter sets
xx<-seq(-2,2,length.out=100)
ee<-0
plot(exp(rnorm(100,0,ee))*(exp(1+1*xx)+exp(0))~xx,type='l',ylim=c(0,30))
lines(exp(rnorm(100,0,ee))*(0.5*(exp(1+1*xx)+exp(0))+0.5)~xx,lty=2)
lines(exp(rnorm(100,0,ee))*(exp(1+8*xx-8*xx^2)+exp(0))~xx,col='red')
lines(exp(rnorm(100,0,ee))*(0.5*(exp(1+8*xx-8*xx^2)+exp(0))+0.5)~xx,col='red',lty=2)
lines(exp(rnorm(100,0,ee))*(exp(3)*plogis(-2+6*xx)+exp(0))~xx,col='blue')
lines(exp(rnorm(100,0,ee))*(0.5*(exp(3)*plogis(-2+6*xx)+exp(0))+0.5)~xx,col='blue',lty=2)
abline(h=exp(2),col='gray',lty=2)
#these look good!
e_sig<-c(0,0.5)
zx<-c(quote(z),quote(x))
# base.formulae<-c('(exp(1+1*x)+1)','(exp(1+8*x-8*x^2)+1)','(14*plogis(-2+10*x)+1)')
# formulae<-setNames(c('r~exp(e)*5',paste0('r~exp(e)*',c('(0.5*',''),rep(base.formulae,each=2),c('+0.5)',''))),
#                    func)
# formulae<-lapply(formulae,str2lang)
formulae<-setNames(c(r~exp(e)*exp(2),
                     r~(exp(1+1*x)+exp(0)),
                     r~exp(e)*(exp(1+8*x-8*x^2)+exp(0)),
                     r~exp(e)*(exp(3)*plogis(-2+6*x)+exp(0))),
                   func)
df$maps<-df$tip_vals<-df$sims<-vector('list',ntot)
rep.seq<-seq_len(nrep)
cur<- -1
for(i in seq.int(1,ntot,nrep)){
  conds<-df[i,var.nms]
  env<-list('x'=zx[[conds$foc+1]])
  trees<-as.multiPhylo(pbtree(n=df[i,'ntax'],scale=1,nsim=nrep))
  inds<-i+rep.seq-1
  for(j in rep.seq){
    tmp<-sim.conthistory(trees[[j]],ntraits=2,traits=c('x','z'),
                         nsims=1,res=500)
    form<-formulae[[conds$func]]
    form<-do.call(substitute,list(expr=form,env=env))
    tmp<-make.traits(tmp,
                     e~rnorm(placeholder_n,0,e_sig[conds$noise+1]),
                     form,
                     log_r~log(r),
                     y~diffusion('r'))
    df$sims[[inds[j]]]<-tmp
    tmp<-get.tip.vals(tmp)[,c('x','y'),1]
    df$tip_vals[[inds[j]]]<-tmp
    tmp<-tmp[,'x',drop=FALSE]
    est.Xsig2<-mean(pic(tmp,trees[[j]])^2)
    df$maps[[inds[j]]]<-make.contsimmap(tmp,trees[[j]],Xsig2=est.Xsig2,Ysig2=0)
    cur<-prog(inds[j],ntot,cur)
  }
}

make.dat<-function(contsimmap,trait.data,beta=~1,alpha=NULL,gamma=~1){
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
  
  #increment stuff...
  foc.traits<-unique(unlist(lapply(c(beta,alpha,gamma),all.vars),use.names=FALSE))
  trait.inds<-match(foc.traits,contsimmap[['traits']])
  #insert check for variables that don't match up!
  t<-s<-vector('list',ntrees)
  z<-setNames(rep(list(vector('list',nsims)),length(foc.traits)),foc.traits)
  counter<-0
  for(i in seq_len(ntrees)){
    tmp.ord<-order(as.numeric(names(contsimmap[['x']])))
    tmp<-lapply(contsimmap[['maps']][tmp.ord],'[[',i)
    t[[i]]<-c(list(NULL),lapply(tmp,'[[','dts')) #NULL for root edge...
    s[[i]]<-lapply(tmp,'[[','state')
    for(j in seq_along(foc.traits)){
      tmp<-matrix(unlist(lapply(contsimmap[['x']][tmp.ord],function(ii) asplit(ii[[i]][,trait.inds[j],,drop=FALSE],3)),
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
  ones<-matrix(1,nrow=sum(n_inc),ncol=1)
  W<-model.matrix(beta,data=z)
  if(!length(W)){
    W<-ones
  }
  B<-ncol(W)
  if(!is.null(alpha)){
    logit<-TRUE
    U<-model.matrix(alpha,data=z)
    if(!length(U)){
      U<-ones
    }
    A<-ncol(U)
  }else{
    logit<-FALSE
  }
  Z<-model.matrix(gamma,data=z)
  if(!length(Z)){
    Z<-ones
  }
  G<-ncol(Z)
  
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
  
  #final data list
  out<-list('N'=N,'K'=K,'T'=T,'E'=E,'C'=C,'B'=B,'G'=G,
            'n_obs'=n_obs,'ind_obs'=ind_obs,'Y'=Y,
            'n_inc'=n_inc,'pos_inc'=pos_inc,'t'=t,'W'=W,'Z'=Z,
            'n_des'=n_des,'ind_des'=ind_des,'pos_des'=pos_des,
            'ind_tip'=ind_tip,'ind_root'=array(ind_root),'ind_prune'=ind_prune,
            'n_dim'=array(n_dim),'ind_dim'=array(ind_dim),'pos_dim'=array(pos_dim),
            'code_obs'=code_obs,'code_tip'=code_tip,'code_nod'=code_nod)
  if(logit){
    out[['A']]<-A
    out[['U']]<-U
  }
  out
}

dat<-make.dat(df$maps[[136]],df$tip_vals[[136]][,'y',drop=FALSE],beta=~x+I(x^2),alpha=~1)
log.stan.mod<-stan_model('~/../Documents/git/contsimmap/R_tests/multirate_BM_loglink_univar_nointravar.stan')
logit.stan.mod<-stan_model('~/../Documents/git/contsimmap/R_tests/multirate_BM_logitlink_univar_nointravar.stan')
optimizing(log.stan.mod,data=dat,verbose=TRUE)
