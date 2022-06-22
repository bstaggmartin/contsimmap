rm(list=ls())
set.seed(321)
library(phytools)
library(contsimmap)
library(parallel)

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

args<-as.numeric(commandArgs(trailingOnly=TRUE))

####CONSTRUCTING SIMULATION SCENARIOS####

#sample size
ntax<-args[1]
#functional relationship type
func<-c('null','simple','sweetspot','threshold')
#related to focal trait or not?
foc<-c(FALSE,TRUE)
#noise or not?
noise<-c(FALSE,TRUE)
#number of replicates per scenario
nrep<-args[2]
#make data.frame of scenarios
df<-expand.grid(ntax,func,foc,noise)
var.nms<-c('ntax','func','foc','noise')
colnames(df)<-var.nms
df<-df[!(!df$foc&df$func=='null'),] #eliminate null, non-focal simulations
nscen<-nrow(df)
tmp<-rep(seq_len(nscen),each=nrep)
df<-df[tmp,]
rownames(df)<-paste(tmp,seq_len(nrep),sep='_')
ntot<-nscen*nrep
#some things are constant across all simulations
#each simulation will be on a pure-birth tree with total height 1
#each simulation will consist of two "explanatory" traits, x and z, evolving via independent BM processes
#   these processes will always start at 0 and have rate 1
#   resolution of 500
#   for now, only 1 observation per species and no intraspecific variation
#   in all cases, x will be the observed explanatory trait
#each simulation will also have a single "response" trait, y, evolving with rate, r, dependent upon either x or z
#   when func is null: r = e*exp(b1)
#   when func is simple: r = e*(exp(b1+b2*x)+exp(g1))
#   when func is sweetspot: r = e*(exp(b1+b2*x+b3*x^2)+exp(g1))
#   when func is threshold: r = e*(exp(a1)*plogis(b1+b2*x)+exp(g1))
#   if foc is TRUE, r will depend on x, otherwise it will depend on z
#   if noise is TRUE, e will follow a "normalized" gamma process (i.e., a gamma process with mean 1 divided by elapsed time)
#     if I'm thinking about this correctly, this essentially means e will be the average rate due to gamma noise
#   most helpful to tweak parameters individually
#deciding on final parameter sets...
# hist(rgamma(1e6,10/10,10)/10/(1/10^2))
#can divide by time increment to "standardize" rate change implied by variance gamma, though the distribution is NOT linked to time
#by a simple scalar!
#let's turn this into a simpler function ("norm" denotes that it's normalized by time and set to have an average value 1):
# norm.gamma.proc2<-function(n,nu,t){
#   rgamma(n,t/nu,nu)*nu^2/t
# }
# #below is simpler and equivalent
# norm.gamma.proc<-function(n,nu,t){
#   rgamma(n,t/nu,t/nu)
# }
rngamproc<-function(n,nu,t){
  out<-vector('numeric',n)
  valid<-nu!=0&t!=0
  tmp<-rep(t,length.out=n)[valid]/rep(nu,length.out=n)[valid]
  out[valid]<-rgamma(sum(valid),tmp,tmp)
  out[!valid]<-1
  out
}
qngamproc<-function(x,nu,t){
  len<-max(unlist(lengths(list(x,nu,t)),use.names=FALSE))
  out<-vector('numeric',len)
  valid<-nu!=0&t!=0
  tmp<-rep(t,length.out=len)[valid]/rep(nu,length.out=len)[valid]
  out[valid]<-qgamma(rep(x,length.out=len)[valid],tmp,tmp)
  out[!valid]<-1
  out
}
xx<-seq(-2,2,length.out=100)
nu<-0.01
t<-0.01
tmp<-cbind(exp(2.3),
           exp(1.5+1*xx)+exp(-1),
           exp(3+0*xx-3*xx^2)+exp(-1),
           exp(3)*plogis(0+8*xx)+exp(-1))
matplot(x=xx,y=cbind(tmp,qngamproc(0.975,nu,t)*tmp,qngamproc(0.025,nu,t)*tmp),type='l',
        lty=rep(c(1,2),c(1,2)*ncol(tmp)),
        lwd=rep(c(4,2),c(1,2)*ncol(tmp)),
        col=rgb(c(0.5,0,1,0),
                c(0.5,0,0,0),
                c(0.5,0,0,1),
                rep(c(1,0.5),c(1,2)*ncol(tmp))),
        xlab='trait',ylab='rate')
#these look good!
e_nu<-c(0,0.01)
zx<-c(quote(z),quote(x))
# base.formulae<-c('(exp(1+1*x)+1)','(exp(1+8*x-8*x^2)+1)','(14*plogis(-2+10*x)+1)')
# formulae<-setNames(c('r~exp(e)*5',paste0('r~exp(e)*',c('(0.5*',''),rep(base.formulae,each=2),c('+0.5)',''))),
#                    func)
# formulae<-lapply(formulae,str2lang)
formulae<-setNames(c(r~e*exp(2.3),
                     r~e*(exp(1.5+1*x)+exp(-1)),
                     r~e*(exp(3+0*x-3*x^2)+exp(-1)),
                     r~e*(exp(3)*plogis(0+8*x)+exp(-1))),
                   func)

####DOING SIMULATIONS####

df$maps<-df$tip_vals<-df$sims<-vector('list',ntot)
rep.seq<-seq_len(nrep)
for(i in seq.int(1,ntot,nrep)){
  if(i==1){
    cur<-prog(0,ntot,-1)
  }
  conds<-df[i,var.nms]
  env<-list('x'=zx[[conds$foc+1]])
  trees<-as.multiPhylo(pbtree(n=df[i,'ntax'],scale=1,nsim=nrep))
  inds<-i+rep.seq-1
  for(j in rep.seq){
    tmp<-sim.conthistory(trees[[j]],ntraits=2,traits=c('x','z'),
                         Xsig2=4,
                         nsims=1,res=500)
    form<-formulae[[conds$func]]
    form<-do.call(substitute,list(expr=form,env=env))
    tmp<-make.traits(tmp,
                     e~rngamproc(n_PH,e_nu[conds$noise+1],dt_PH),
                     form,
                     log_r~log(r),
                     y~diffusion('r'))
    df$sims[[inds[j]]]<-tmp
    tmp<-get.tip.vals(tmp)[,c('x','y'),1]
    df$tip_vals[[inds[j]]]<-tmp
    tmp<-tmp[,'x',drop=FALSE]
    est.Xsig2<-mean(pic(tmp,trees[[j]])^2)
    tmp<-make.contsimmap(tmp,trees[[j]],Xsig2=est.Xsig2,Ysig2=0)
    tmp<-make.traits(tmp,d~diffusion(d=tmp[['params']][['Xsig2']],X0=tmp[['nodes']][['0']]))
    df$maps[[inds[j]]]<-tmp
    cur<-prog(inds[j],ntot,cur)
  }
}

####FITTING MODELS####

alpha.formulae<-list(NULL,NULL,NULL,~1,NULL,NULL,~1)
beta.formulae<-c(~1,~x,~x+I(x^2),~x,~d,~d+I(d^2),~d)
gamma.formulae<-list(NULL,~1,~1,~1,~1,~1,~1)
#use naive marginalization for non-dummy variables (makes NULL models more competitive)?
if(args[4]==0){
  naive.marg<-c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE)
}else{
  naive.marg<-rep(TRUE,7)
}
foo<-function(i){
  tmp<-as.list(df[i,c('maps','tip_vals')])
  tmp[['maps']]<-tmp[['maps']][[1]]
  tmp[['tip_vals']]<-tmp[['tip_vals']][[1]][,'y',drop=FALSE]
  tmp[['alpha']]<-alpha.formulae
  tmp[['beta']]<-beta.formulae
  tmp[['gamma']]<-gamma.formulae
  tmp[['naive.marg']]<-naive.marg
  tmp[['niter']]<-args[3]
  tmp
}
split.df<-lapply(seq_len(ntot),foo)
foo<-function(df.row){
  tmp.dat<-contsimmap::prep.dat(df.row[['maps']],df.row[['tip_vals']])
  lapply(seq_along(df.row[['alpha']]),function(j) contsimmap::fit.model(tmp.dat,
                                                                        beta=df.row[['beta']][[j]],
                                                                        alpha=df.row[['alpha']][[j]],
                                                                        gamma=df.row[['gamma']][[j]],
                                                                        naive.marg=df.row[['naive.marg']][[j]],
                                                                        niter=df.row[['niter']]))
}
cl<-makeCluster(14)
res<-parLapply(cl=cl,X=split.df,fun=foo)
stopCluster(cl)
res<-do.call(rbind,res)
colnames(res)<-c(func,paste0('d_',func[-1]))
df<-cbind(df,res)

####EXPORT RESULTS####

saveRDS(df,paste0('/mnt/gs18/scratch/users/mart2121/contsimmap_sim_study_',paste(ntax,collapse='_'),'_',Sys.Date()))
