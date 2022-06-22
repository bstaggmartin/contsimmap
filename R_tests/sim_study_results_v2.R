rm(list=ls())
library(contsimmap)

df<-rbind(readRDS('D:/contsimmap_sim_study/contsimmap_sim_study_50_2022-06-11'),
          readRDS('D:/contsimmap_sim_study/contsimmap_sim_study_100_2022-06-11'),
          readRDS('D:/contsimmap_sim_study/contsimmap_sim_study_200_2022-06-11'))
func<-colnames(df)[-seq_len(which(colnames(df)=='maps'))]

codes<-do.call(cbind,lapply(df[,func],function(ii) unlist(lapply(ii,'[[','return_code'),use.names=FALSE)))
lnL<-do.call(cbind,lapply(df[,func],function(ii) unlist(lapply(ii,'[[','value'),use.names=FALSE)))
hist(lnL) #one fit got really wonky
codes[which.max(lnL)%%nrow(lnL),] #shoot, the problem is you're using maxes, you should make sure return_codes are 0 in next iteration
which(codes!=0,arr.ind = TRUE)
#interestingly, it's the threshold models that seem to mess up most often...
#might be worth switching back to stan's inv_logit function?
df[,-seq_len(which(colnames(df)=='maps'))][which(codes!=0,arr.ind = TRUE)]<-
  lapply(df[,-seq_len(which(colnames(df)=='maps'))][which(codes!=0,arr.ind = TRUE)],
         function(ii) list('return_code'=ii[['return_code']],'value'=NA,'par'=rep(NA,length(ii[['par']]))))
lnL<-do.call(cbind,lapply(df[,func],function(ii) unlist(lapply(ii,'[[','value'),use.names=FALSE)))
hist(lnL) #seems more promising
AIC<- -2*sweep(lnL,2,
               unlist(lapply(df[1,func],function(ii) length(ii[[1]][['par']])),use.names=FALSE))
dAIC<-t(apply(AIC,1,function(ii) ii-min(ii,na.rm=TRUE)))
AIC.wgts<-exp(-0.5*dAIC)
AIC.wgts[is.na(AIC.wgts)]<-0
AIC.wgts<-AIC.wgts/rowSums(AIC.wgts)

full.fac<-Reduce(':',lapply(df[,seq_len(which(colnames(df)=='sims')-2)],as.factor))
sum.AIC.wgts<-apply(AIC.wgts,2,function(ii) tapply(ii,full.fac,mean,na.rm=TRUE))
sum.AIC.wgts<-sum.AIC.wgts[,order(grepl('null|^d_',func),decreasing=TRUE)]
barplot(t(sum.AIC.wgts),col=hcl.colors(4)[c(rep(1,4),2:4)],
        las=2,cex.names=0.6,
        legend.text=colnames(sum.AIC.wgts),args.legend=list(bty='n'))

AIC.wgts2<-AIC.wgts[,1:4]
AIC.wgts2<-AIC.wgts2/rowSums(AIC.wgts2)
sum.AIC.wgts2<-apply(AIC.wgts2,2,function(ii) tapply(ii,full.fac,mean,na.rm=TRUE))
barplot(t(sum.AIC.wgts2),col=hcl.colors(4),
        las=2,cex.names=0.6,
        legend.text=colnames(sum.AIC.wgts2),args.legend=list(bty='n'))

tmp.AIC.wgts<-AIC.wgts[,order(grepl('null|^d_',func),decreasing=TRUE)]
barplot(t(tmp.AIC.wgts),col=hcl.colors(4)[c(rep(1,4),2:4)],border=NA,space=0)

#predicted relationships
xx<-seq(-3,3,length.out=1000)
tmp.AIC.wgts<-AIC.wgts[,!grepl('^d_',func)]
tmp.AIC.wgts[,'null']<-tmp.AIC.wgts[,'null']+rowSums(AIC.wgts[,grepl('^d_',func)])
null.pars<-lapply(df$null,'[[','par')
simple.pars<-lapply(df$simple,'[[','par')
sweetspot.pars<-lapply(df$sweetspot,'[[','par')
threshold.pars<-lapply(df$threshold,'[[','par')
null.preds<-do.call(cbind,lapply(null.pars,function(ii) exp(rep(ii,length(xx)))))
null.preds[is.na(null.preds)]<-0
simple.preds<-do.call(cbind,lapply(simple.pars,function(ii) exp(ii[1]+ii[2]*xx)+exp(ii[3])))
simple.preds[is.na(simple.preds)]<-0
sweetspot.preds<-do.call(cbind,lapply(sweetspot.pars,function(ii) exp(ii[1]+ii[2]*xx+ii[3]*xx^2)+exp(ii[4])))
sweetspot.preds[is.na(sweetspot.preds)]<-0
threshold.preds<-do.call(cbind,lapply(threshold.pars,function(ii) exp(ii[1])*plogis(ii[2]+ii[3]*xx)+exp(ii[4])))
threshold.preds[is.na(threshold.preds)]<-0
avg.preds<-sweep(null.preds,2,tmp.AIC.wgts[,'null'],'*')+
  sweep(simple.preds,2,tmp.AIC.wgts[,'simple'],'*')+
  sweep(sweetspot.preds,2,tmp.AIC.wgts[,'sweetspot'],'*')+
  sweep(threshold.preds,2,tmp.AIC.wgts[,'threshold'],'*')

matplot(x=xx,log(avg.preds[,df$func=='null'&df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(avg.preds[,df$ntax==50&(df$func=='null'&df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==100&(df$func=='null'&df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==200&(df$func=='null'&df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(rep(2.3,length(xx))~xx,lwd=8)

matplot(x=xx,log(avg.preds[,df$func=='null'|!df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(avg.preds[,df$ntax==50&(df$func=='null'|!df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==100&(df$func=='null'|!df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==200&(df$func=='null'|!df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(rep(2.3,length(xx))~xx,lwd=8)

matplot(x=xx,log(avg.preds[,df$func=='simple'&df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(avg.preds[,df$ntax==50&df$func=='simple'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==100&df$func=='simple'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==200&df$func=='simple'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(log(exp(1.5+1*xx)+exp(-1))~xx,lwd=8)

matplot(x=xx,log(avg.preds[,df$func=='sweetspot'&df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(avg.preds[,df$ntax==50&df$func=='sweetspot'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==100&df$func=='sweetspot'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(avg.preds[,df$ntax==200&df$func=='sweetspot'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(log(exp(3+0*xx-3*xx^2)+exp(-1))~xx,lwd=8)

matplot(x=xx,log(avg.preds[,df$func=='threshold'&df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(avg.preds[,df$ntax==50&df$func=='threshold'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0.5,0.5,0.5,0.5))
tmp<-apply(log(avg.preds[,df$ntax==100&df$func=='threshold'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0.25,0.25,0.25,0.5))
tmp<-apply(log(avg.preds[,df$ntax==200&df$func=='threshold'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,0,0.5))
lines(log(exp(3)*plogis(0+8*xx)+exp(-1))~xx,lwd=8)

#thinking in terms of thresholds rather than weights to better parse things out...
tmp.dAIC<-dAIC
tmp.dAIC[,'null']<-do.call(pmin,c(asplit(dAIC[,grepl('null|^d_',func)],2),na.rm=TRUE))
tmp.dAIC<-tmp.dAIC[,!grepl('^d_',func)]
sig.dAIC<-tmp.dAIC<2
sig.dAIC[is.na(sig.dAIC)]<-FALSE
props<-do.call(rbind,tapply(asplit(sig.dAIC,1),full.fac,function(ii) Reduce('+',ii)/length(ii)))
npars<-c(1,3,4,4)
foo<-function(i){
  tmp<-npars[sig.dAIC[i,]]
  tmp.min<-min(tmp)
  out<-which(npars==tmp.min)
  if(length(out)>1){
    out<-out[which.min(tmp.dAIC[i,out])]
  }
  out
}
best.mod<-unlist(lapply(seq_len(nrow(tmp.dAIC)),foo),use.names=FALSE)
plot(as.numeric(as.factor(best.mod)),type='l')
preds<-array(c(null.preds,simple.preds,sweetspot.preds,threshold.preds),c(dim(null.preds),4))
best.preds<-matrix(preds[cbind(seq_len(dim(preds)[1]),rep(seq_len(dim(preds)[2]),each=dim(preds)[1]),rep(best.mod,each=dim(preds)[1]))],
                   dim(preds)[1])

matplot(x=xx,log(best.preds[,df$func=='null'|!df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(best.preds[,df$ntax==50&(df$func=='null'|!df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==100&(df$func=='null'|!df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==200&(df$func=='null'|!df$foc)]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(rep(2.3,length(xx))~xx,lwd=8)

matplot(x=xx,log(best.preds[,df$func=='simple'&df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(best.preds[,df$ntax==50&df$func=='simple'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==100&df$func=='simple'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==200&df$func=='simple'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(log(exp(1.5+1*xx)+exp(-1))~xx,lwd=8)

matplot(x=xx,log(best.preds[,df$func=='sweetspot'&df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(best.preds[,df$ntax==50&df$func=='sweetspot'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==100&df$func=='sweetspot'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==200&df$func=='sweetspot'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(log(exp(3+0*xx-3*xx^2)+exp(-1))~xx,lwd=8)

matplot(x=xx,log(best.preds[,df$func=='threshold'&df$foc]),type='l',col=rgb(0,0,0,0),lty=1,lwd=3)
tmp<-apply(log(best.preds[,df$ntax==50&df$func=='threshold'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(1,0,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==100&df$func=='threshold'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,1,0,0.1))
tmp<-apply(log(best.preds[,df$ntax==200&df$func=='threshold'&df$foc]),1,quantile,prob=c(0.025,0.975))
polygon(c(tmp[1,],rev(tmp[2,]))~c(xx,rev(xx)),border=NA,col=rgb(0,0,1,0.1))
lines(log(exp(3)*plogis(0+8*xx)+exp(-1))~xx,lwd=8)
#pretty much equivalent to weighted results...

tmp<-props[grepl('TRUE',rownames(props)),]
image(t(tmp[4:1,]),col=hcl.colors(palette='Oslo',20))
image(t(tmp[4:1+4,]),col=hcl.colors(palette='Oslo',20))
image(t(tmp[4:1+8,]),col=hcl.colors(palette='Oslo',20))

props2<-do.call(rbind,tapply(best.mod,full.fac,tabulate,nbin=4))/20
tmp<-props2[grepl('TRUE',rownames(props2)),]
image(t(tmp[4:1,]),col=hcl.colors(palette='Oslo',20))
image(t(tmp[4:1+4,]),col=hcl.colors(palette='Oslo',20))
image(t(tmp[4:1+8,]),col=hcl.colors(palette='Oslo',20))

tmp.AIC.wgts<-AIC.wgts[,!grepl('^d_',func)]
tmp.AIC.wgts[,'null']<-tmp.AIC.wgts[,'null']+rowSums(AIC.wgts[,grepl('^d_',func)])
sum.AIC.wgts<-apply(tmp.AIC.wgts,2,function(ii) tapply(ii,full.fac,mean,na.rm=TRUE))
tmp<-sum.AIC.wgts[grepl('TRUE',rownames(sum.AIC.wgts)),]
image(t(tmp[4:1,]),col=hcl.colors(palette='Oslo',20))
image(t(tmp[4:1+4,]),col=hcl.colors(palette='Oslo',20))
image(t(tmp[4:1+8,]),col=hcl.colors(palette='Oslo',20))
