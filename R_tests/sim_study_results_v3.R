rm(list=ls())
library(contsimmap)

df<-rbind(readRDS('D:/contsimmap_sim_study/contsimmap_sim_study_50_2022-06-13'),
          readRDS('D:/contsimmap_sim_study/contsimmap_sim_study_100_2022-06-13'),
          readRDS('D:/contsimmap_sim_study/contsimmap_sim_study_200_2022-06-13'))
func<-colnames(df)[-seq_len(which(colnames(df)=='maps'))]

codes<-do.call(cbind,lapply(df[,func],function(ii) unlist(lapply(ii,'[[','return_code'),use.names=FALSE))) #no failed fits now!
lnL<-do.call(cbind,lapply(df[,func],function(ii) unlist(lapply(ii,'[[','value'),use.names=FALSE)))
hist(lnL) #shit, but one fit still got really wonky
length(lnL)
sum(lnL<0)
lnL[lnL>0] #just three went "off the rails"
inds<-which(lnL>100,arr.ind=TRUE)
df[inds[,1],func[inds[,2]]] #all threshold models--probably some sort of overflow issue
df[inds[,1],seq_len(which(colnames(df)=='sims')-1)] #all 50 tip trees!
#replace these with -Inf lnL
lnL[lnL>100]<- -Inf
hist(lnL) #much better
npars<-unlist(lapply(df[1,func],function(ii) length(ii[[1]][['par']])),use.names=FALSE)
AIC<- -2*sweep(lnL,2,npars,'-')
cor<-sweep(matrix(df$ntax-1,nrow=nrow(df),ncol=length(npars)),2,npars,'-')
cor<-sweep(1/cor,2,2*npars*(npars+1),'*')
AICc<-AIC+cor
dAIC<-t(apply(AICc,1,function(ii) ii-min(ii,na.rm=TRUE)))
AIC.wgts<-exp(-0.5*dAIC)
AIC.wgts[is.na(AIC.wgts)]<-0
AIC.wgts<-AIC.wgts/rowSums(AIC.wgts)
foo<-function(i){
  if(is.logical(df[[i]])){
    out<-as.factor(df[[i]])
    levels(out)<-paste0(c('!',''),i)
    out
  }else{
    as.factor(df[[i]])
  }
}
tmp<-colnames(df)[seq_len(which(colnames(df)=='sims')-1)]
facs<-lapply(setNames(tmp,tmp),foo)

summarize.wgts<-function(wgts,fac,models=NULL,subset=NULL){
  if(!is.null(models)){
    wgts<-do.call(cbind,lapply(models,function(ii) rowSums(wgts[,ii,drop=FALSE])))
  }
  tmp.facs<-facs
  if(!is.null(subset)){
    wgts<-wgts[subset,,drop=FALSE]
    for(i in seq_along(tmp.facs)){
      tmp.facs[[i]]<-tmp.facs[[i]][subset]
    }
  }
  new.fac<-fac
  nms<-unlist(lapply(fac,function(ii) as.character(ii[2])),use.names=FALSE)
  tmp<-lapply(fac,'[[',3)
  for(i in seq_along(fac)){
    new.fac[[i]]<-eval(tmp[[i]],envir=tmp.facs)
  }
  names(new.fac)<-nms
  out<-lapply(asplit(wgts,2),function(ii) tapply(ii,new.fac,mean))
  nms<-c(dimnames(out[[1]]),'model'=list(names(out)))
  ndims<-length(nms)
  out<-aperm(array(unlist(out,use.names=FALSE),lengths(nms),nms),c(ndims,seq_len(ndims-1)))
  sweep(out,seq_len(ndims)[-1],apply(out,-1,sum),'/')
}
nice.image<-function(x,xlab=names(dimnames(x))[1],ylab=names(dimnames(x))[2],...){
  dims<-dim(x)
  image(x[,rev(seq_len(dims[2])),drop=FALSE],...,zlim=c(0,1),
        xaxt='n',xlab=xlab,
        yaxt='n',ylab=ylab)
  axis(1,seq(0,1,length.out=dims[1]),dimnames(x)[[1]],tick=FALSE)
  axis(2,seq(0,1,length.out=dims[2]),rev(dimnames(x)[[2]]),tick=FALSE)
}

#foc only, noise marginalized, dummy models expanded
tmp<-summarize.wgts(AIC.wgts,
                    list(func~func,
                         size~ntax),
                    NULL,
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  nice.image(tmp[,,i],col=hcl.colors(palette='Oslo',100),main=i,
             xlab='fitted model',ylab='simulated model')
}
#foc only, dummy models expanded
tmp<-summarize.wgts(AIC.wgts,
                    list(func~func,
                         size~ntax,
                         noise~noise),
                    NULL,
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  for(j in dimnames(tmp)[['noise']]){
    nice.image(tmp[,,i,j],col=hcl.colors(palette='Oslo',100),main=paste(i,j,sep='; '),
               xlab='fitted model',ylab='simulated model')
  }
}
#foc only, noise marginalized, dummy models collapsed
tmp<-summarize.wgts(AIC.wgts,
                    list(func~func,
                         size~ntax),
                    list('null'=c('null','d_simple','d_sweetspot','d_threshold'),
                         'simple'='simple',
                         'sweetspot'='sweetspot',
                         'threshold'='threshold'),
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  nice.image(tmp[,,i],col=hcl.colors(palette='Oslo',100),main=i,
             xlab='fitted model',ylab='simulated model')
}
#foc only, dummy models collapsed
tmp<-summarize.wgts(AIC.wgts,
                    list(func~func,
                         size~ntax,
                         noise~noise),
                    list('null'=c('null','d_simple','d_sweetspot','d_threshold'),
                         'simple'='simple',
                         'sweetspot'='sweetspot',
                         'threshold'='threshold'),
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  for(j in dimnames(tmp)[['noise']]){
    nice.image(tmp[,,i,j],col=hcl.colors(palette='Oslo',100),main=paste(i,j,sep='; '),
               xlab='fitted model',ylab='simulated model')
  }
}
#foc only, noise marginalized, dummy models excluded
tmp<-summarize.wgts(AIC.wgts,
                    list(func~func,
                         size~ntax),
                    list('null'='null',
                         'simple'='simple',
                         'sweetspot'='sweetspot',
                         'threshold'='threshold'),
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  nice.image(tmp[,,i],col=hcl.colors(palette='Oslo',100),main=i,
             xlab='fitted model',ylab='simulated model')
}
#foc only, dummy models excluded
tmp<-summarize.wgts(AIC.wgts,
                    list(func~func,
                         size~ntax,
                         noise~noise),
                    list('null'='null',
                         'simple'='simple',
                         'sweetspot'='sweetspot',
                         'threshold'='threshold'),
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  for(j in dimnames(tmp)[['noise']]){
    nice.image(tmp[,,i,j],col=hcl.colors(palette='Oslo',100),main=paste(i,j,sep='; '),
               xlab='fitted model',ylab='simulated model')
  }
}
#relationship or not?, foc only, noise marginalized, dummy models excluded
tmp<-summarize.wgts(AIC.wgts,
                    list(func~c('null','not null')[as.numeric(func=='null')+1],
                         size~ntax),
                    list('null'='null',
                         'not null'=c('simple','sweetspot','threshold')),
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  nice.image(tmp[,,i],col=hcl.colors(palette='Oslo',100),main=i,
             xlab='fitted model',ylab='simulated model')
}
#relationship or not?, foc only, dummy models excluded
tmp<-summarize.wgts(AIC.wgts,
                    list(func~c('null','not null')[as.numeric(func=='null')+1],
                         size~ntax,
                         noise~noise),
                    list('null'='null',
                         'not null'=c('simple','sweetspot','threshold')),
                    df$foc)
for(i in dimnames(tmp)[['size']]){
  for(j in dimnames(tmp)[['noise']]){
    nice.image(tmp[,,i,j],col=hcl.colors(palette='Oslo',100),main=paste(i,j,sep='; '),
               xlab='fitted model',ylab='simulated model')
  }
}
#relationship or not?, foc/noise marginalized
tmp<-summarize.wgts(AIC.wgts,
                    list(func~c('null','not null')[as.numeric(df$func=='null'|!df$foc)+1],
                         size~ntax),
                    list('null'=c('null','d_simple','d_sweetspot','d_threshold'),
                         'not null'=c('simple','sweetspot','threshold')),
                    NULL)
for(i in dimnames(tmp)[['size']]){
  nice.image(tmp[,,i],col=hcl.colors(palette='Oslo',100),main=i,
             xlab='fitted model',ylab='simulated model')
}
#relationship or not?, foc marginalized
tmp<-summarize.wgts(AIC.wgts,
                    list(func~c('null','not null')[as.numeric(df$func=='null'|!df$foc)+1],
                         size~ntax,
                         noise~noise),
                    list('null'=c('null','d_simple','d_sweetspot','d_threshold'),
                         'not null'=c('simple','sweetspot','threshold')),
                    NULL)
for(i in dimnames(tmp)[['size']]){
  for(j in dimnames(tmp)[['noise']]){
    nice.image(tmp[,,i,j],col=hcl.colors(palette='Oslo',100),main=paste(i,j,sep='; '),
               xlab='fitted model',ylab='simulated model')
  }
}

#AIC (even AICc) is biased towards more complex models, even when there's no noise
#the dummy models, however, mitigate this bias almost completely (at the cost of lesser power)
foo<-function(x){
  tmp<-quantile(contsimmap:::.extract.traits(x,foc.traits='x',lhs=TRUE)[['x']],prob=c(0.1,0.9))
  seq(tmp[1],tmp[2],length.out=100)
}
xxs<-lapply(df[['maps']],foo)
pars<-lapply(c('null','simple','sweetspot','threshold'),function(ii) lapply(df[[ii]],'[[','par'))
formulae<-c(~rep(exp(par[1]),length(x)),
            ~exp(par[1]+par[2]*x)+exp(par[3]),
            ~exp(par[1]+par[2]*x+par[3]*x^2)+exp(par[4]),
            ~exp(par[1])*plogis(par[2]+par[3]*x)+exp(par[4]))
formulae<-lapply(formulae,'[[',2)
preds<-pars
avg.preds<-matrix(NA,length(xxs[[1]]),length(xxs))
tmp.wgts<-AIC.wgts
tmp.wgts[,'null']<-tmp.wgts[,'null']+rowSums(tmp.wgts[,grepl('^d_',colnames(tmp.wgts))])
tmp.wgts<-tmp.wgts[,!grepl('^d_',colnames(tmp.wgts))]
for(i in seq_along(xxs)){
  for(j in seq_along(pars)){
    tmp<-eval(formulae[[j]],list(par=pars[[j]][[i]],x=xxs[[i]]))
    tmp[is.na(tmp)]<-0 #fine since NA pars associated with 0 AIC weight anyways
    preds[[j]][[i]]<-tmp
  }
  avg.preds[,i]<-Reduce('+',lapply(seq_along(pars),function(jj) tmp.wgts[i,jj]*preds[[jj]][[i]]))
}
xx<-do.call(cbind,xxs)

inds<-df$func=='null'|!df$foc
matplot(x=xx[,inds],y=log(avg.preds[,inds]),lwd=4,lty=1,type='l',
        col=hcl.colors(3,palette='sunset',rev=TRUE,alpha=0.2)[as.factor(df$ntax[inds])],
        xlab='trait value',ylab='ln(rate)')
xlim<-par()$usr[1:2]
xlim<-xlim+c(1,-1)*0.02*diff(xlim)
lines(x=xlim,y=rep(2.3,2),lwd=6,lty=2,lend=2)

inds<-df$func=='simple'&df$foc
matplot(x=xx[,inds],y=log(avg.preds[,inds]),lwd=4,lty=1,type='l',
        col=hcl.colors(3,palette='sunset',rev=TRUE,alpha=0.2)[as.factor(df$ntax[inds])],
        xlab='trait value',ylab='ln(rate)')
xlim<-par()$usr[1:2]
xlim<-xlim+c(1,-1)*0.02*diff(xlim)
tmp.x<-seq(xlim[1],xlim[2],length.out=100)
lines(x=tmp.x,y=log(exp(1.5+1*tmp.x)+exp(-1)),lwd=6,lty=2,lend=2)

inds<-df$func=='sweetspot'&df$foc
matplot(x=xx[,inds],y=log(avg.preds[,inds]),lwd=4,lty=1,type='l',
        col=hcl.colors(3,palette='sunset',rev=TRUE,alpha=0.2)[as.factor(df$ntax[inds])],
        xlab='trait value',ylab='ln(rate)')
xlim<-par()$usr[1:2]
xlim<-xlim+c(1,-1)*0.02*diff(xlim)
tmp.x<-seq(xlim[1],xlim[2],length.out=100)
lines(x=tmp.x,y=log(exp(3+0*tmp.x-3*tmp.x^2)+exp(-1)),lwd=6,lty=2,lend=2)

inds<-df$func=='threshold'&df$foc
matplot(x=xx[,inds],y=log(avg.preds[,inds]),lwd=4,lty=1,type='l',
        col=hcl.colors(3,palette='sunset',rev=TRUE,alpha=0.2)[as.factor(df$ntax[inds])],
        xlab='trait value',ylab='ln(rate)')
xlim<-par()$usr[1:2]
xlim<-xlim+c(1,-1)*0.02*diff(xlim)
tmp.x<-seq(xlim[1],xlim[2],length.out=100)
lines(x=tmp.x,y=log(exp(3)*plogis(0+8*tmp.x)+exp(-1)),lwd=6,lty=2,lend=2)

tmp.x<-xx[,df$foc]
true.y<-tmp.y<-avg.preds[,df$foc]
for(i in seq_len(ncol(tmp.x))){
  true.y[,i]<-switch(as.character(df$func)[df$foc][i],
                     'null'=exp(2.3),
                     'simple'=exp(1.5+1*tmp.x[,i])+exp(-1),
                     'sweetspot'=exp(3+0*tmp.x[,i]-3*tmp.x[,i]^2)+exp(-1),
                     'threshold'=exp(3)*plogis(0+8*tmp.x[,i])+exp(-1))
}
matplot(x=log(true.y),y=log(tmp.y),lwd=4,lty=1,type='l',
        col=hcl.colors(3,palette='sunset',rev=TRUE,alpha=0.1)[as.factor(df$ntax[df$foc])],
        xlab='ln(true rate)',ylab='ln(fitted rate)')
abline(0,1,lty=2,lwd=3)
cor(as.vector(log(tmp.y)),as.vector(log(true.y))) #nice
