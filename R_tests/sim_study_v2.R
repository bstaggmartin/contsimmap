rm(list=ls())
set.seed(123)
library(phytools)
library(contsimmap)
library(rstan)

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
#some things are constant across all simulations
#each simulation will be on a pure-birth tree with total height 1
#each simulation will consist of two "explanatory" traits, x and z, evolving via independent BM processes
#   these processes will always start at 0 and have rate 1
#   resolution of 1,000
#   for now, only 1 observation per species and no intraspecific variation
#   in all cases, x will be the observed explanatory trait
#each simulation will also have a single "response" trait, y, evolving with rate, r, dependent upon either x or z
#   when func is null: r = exp(e)*exp(b0)
#   when func is simple: r = exp(e)*(exp(b0+b1*x)+exp(g0))
#   when func is sweetspot: r = exp(e)*(exp(b0+b1*x+b2*x^2)+exp(g0))
#   when func is threshold: r = exp(e)*(exp(a0)*plogis(b0+b1*x)+exp(g0))
#   if foc is TRUE, r will depend on x, otherwise it will depend on z
#   if noise is TRUE, e will be ~normal(0,0.5), otherwise it will be ~normal(0,0)
#   most helpful to tweak parameters individually
#deciding on final parameter sets...
xx<-seq(-2,2,length.out=100)
ee<-0.5
tmp<-cbind(exp(2.3),
           exp(1+1*xx)+exp(0),
           exp(2+4*xx-4*xx^2)+exp(0),
           exp(3)*plogis(-1+6*xx)+exp(0))
matplot(x=xx,y=cbind(tmp,exp(2*ee)*tmp,exp(-2*ee)*tmp),type='l',
        lty=rep(c(1,2),c(1,2)*ncol(tmp)),
        lwd=rep(c(4,2),c(1,2)*ncol(tmp)),
        col=rgb(c(0.5,0,1,0),
                c(0.5,0,0,0),
                c(0.5,0,0,1),
                rep(c(1,0.5),c(1,2)*ncol(tmp))),
        xlab='trait',ylab='rate')
# tmp<-0.25*tmp+0.75 #"weak" versions?
# matplot(x=xx,y=cbind(tmp,exp(2*ee)*tmp,exp(-2*ee)*tmp),type='l',
#         lty=rep(c(1,2),c(1,2)*ncol(tmp)),
#         lwd=rep(c(4,2),c(1,2)*ncol(tmp)),
#         col=rgb(c(0.5,0,1,0),
#                 c(0.5,0,0,0),
#                 c(0.5,0,0,1),
#                 rep(c(1,0.5),c(1,2)*ncol(tmp))),
#         xlab='trait',ylab='rate')
#these look good!
e_sig<-c(0,0.5)
zx<-c(quote(z),quote(x))
# base.formulae<-c('(exp(1+1*x)+1)','(exp(1+8*x-8*x^2)+1)','(14*plogis(-2+10*x)+1)')
# formulae<-setNames(c('r~exp(e)*5',paste0('r~exp(e)*',c('(0.5*',''),rep(base.formulae,each=2),c('+0.5)',''))),
#                    func)
# formulae<-lapply(formulae,str2lang)
formulae<-setNames(c(r~exp(e)*exp(2.3),
                     r~(exp(1+1*x)+exp(0)),
                     r~exp(e)*(exp(1+4*x-4*x^2)+exp(0)),
                     r~exp(e)*(exp(3)*plogis(-1+6*x)+exp(0))),
                   func)

####DOING SIMULATIONS####

df$maps<-df$tip_vals<-df$sims<-vector('list',ntot)
rep.seq<-seq_len(nrep)
cur<-prog(0,ntot,-1)
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

df[,func]<-rep(list(vector('list',ntot)),length(func))
alpha.formulae<-setNames(list(NULL,NULL,NULL,~1),func)
beta.formulae<-setNames(c(~1,~x,~x+I(x^2),~x),func)
gamma.formulae<-setNames(list(NULL,~1,~1,~1),func)
for(i in seq_len(ntot)){
  tmp.dat<-prep.dat(df$maps[[i]],df$tip_vals[[i]][,'y',drop=FALSE])
  cat('Begining model fits for simulation ',i,' out of ',ntot,'...\n',sep='')
  for(j in func){
    df[[j]][[i]]<-fit.model(tmp.dat,
                            beta=beta.formulae[[j]],
                            alpha=alpha.formulae[[j]],
                            gamma=gamma.formulae[[j]],
                            niter=10,ncores=4)
    cat(j,' done\r',sep='')
  }
  cat('\n')
}
