library(contsimmap)
df3060<-readRDS('D:/contsimmap_sim_study/sim_study_30_60_2022-06-10')
df120<-readRDS('D:/contsimmap_sim_study/contsimmap_sim_study_2022-06-10')
df<-rbind(df3060,df120)
tmp<-which(colnames(df)=='maps')
func<-colnames(df)[-seq_len(tmp)]
tmp<-do.call(cbind,lapply(df[,func],function(ii) unlist(lapply(ii,'[[','value'),use.names=FALSE)))
colnames(tmp)<-paste0(func,'_lnL')
df<-cbind(df,tmp)
tmp<-sweep(tmp,2,
           unlist(lapply(df[1,func],function(ii) length(ii[[1]][['par']])),use.names=FALSE))
tmp<- -2*tmp
colnames(tmp)<-paste0(func,'_AIC')
df<-cbind(df,tmp)
tmp<-t(apply(tmp,1,function(ii) ii-min(ii)))
colnames(tmp)<-paste0(func,'_dAIC')
df<-cbind(df,tmp)
fac<-df$func:as.factor(df$foc)
for(i in func){
  resp<-df[[paste0(i,'_dAIC')]]
  bw<-bw.nrd0(resp)
  dens<-tapply(resp,fac,density,bw=bw,from=0)
  plot(NA,xlim=c(0,length(dens)+1),ylim=c(0,min(100,max(unlist(lapply(dens,'[[','x'),use.names=FALSE)))),
       xaxt='n',xlab='',ylab=bquote(Delta*"AIC for "~.(i)~~"model"))
  axis(1,seq_along(dens),names(dens),las=2)
  for(j in seq_along(dens)){
    tmp.dens<-dens[[j]]
    if(!is.null(tmp.dens)){
      tmp.max<-max(tmp.dens$y)
      polygon(c(tmp.dens$x,rev(tmp.dens$x))~j+rep(0.25/tmp.max*c(-1,1),each=length(tmp.dens$y))*c(tmp.dens$y,rev(tmp.dens$y)),
              col='gray',border=NA)
      tmp.resp<-resp[fac==names(dens)[j]]
      points(x=rep(j,length(tmp.resp)),y=tmp.resp)
    }
  }
  abline(h=2,col='red')
}

#AIC weights would probably be the best way to do this...
tmp<-exp(-0.5*tmp)
tmp<-sweep(tmp,1,rowSums(tmp),'/')
ranked<-t(tmp[df$foc,])
inds<-cbind(as.numeric(df$func[df$foc]),seq_len(ncol(ranked)))
ranked<-ranked[,order(df$func[df$foc],ranked[inds],decreasing=TRUE)]
inds<-order(sort(df$func[df$foc],decreasing=TRUE))
ranked<-ranked[,inds]
inds<-seq_len(ncol(ranked))
reps<-rep(1,length(inds))
reps[cumsum(table(df$func[df$foc]))]<-4
inds<-rep(inds,reps)
inds[cumsum(table(df$func[df$foc])+3)-rep(seq_len(3)-1,each=length(levels(df$func[df$foc])))]<-NA
pos<-barplot(ranked[,inds],border=NA,col=hcl.colors(7),space=0,legend.text=func,
             args.legend=list(bty='n'))
pos<-c(0,cumsum(table(df$func[df$foc])[-length(levels(df$func[df$foc]))]+3))+
  cumsum(table(df$func[df$foc])+3)-3
pos<-pos/2
axis(1,pos,c('null','simple','sweetspot','threshold'),tick=FALSE)

test<-t(tmp)
barplot(test,border=NA,col=hcl.colors(7),space=0,legend.text=func,args.legend=list(bty='n'),
        names.arg=paste(df$func,df$foc,sep=';'),las=2)
barplot(test,border=NA,col=hcl.colors(5)[c(1,3,4,5,rep(2,3))],space=0,legend.text=func,args.legend=list(bty='n'),
        names.arg=paste(df$func,df$foc,sep=';'),las=2)
test2<-test[seq_len(4),]
test2[1,]<-test[1,]+colSums(test[-seq_len(4),])
barplot(test2,border=NA,col=hcl.colors(4),space=0,legend.text=func[seq_len(4)],args.legend=list(bty='n'),
        names.arg=paste(df$func,df$foc,sep=';'),las=2)
test3<-apply(test2,1,function(ii) tapply(ii,as.factor(df$ntax):as.factor(df$foc):as.factor(df$func):as.factor(df$noise),mean))
barplot(t(test3),col=hcl.colors(4),space=0,legend.text=func[seq_len(4)],args.legend=list(bty='n'),las=2,cex.names=0.6)

test4<-apply(test,1,function(ii) tapply(ii,as.factor(df$ntax):as.factor(df$foc):as.factor(df$func):as.factor(df$noise),mean))
barplot(t(test4),col=hcl.colors(7),space=0,legend.text=func,args.legend=list(bty='n'))

#adequacy?
pics<-lapply(seq_len(nrow(df)),function(ii) pic(df$tip_vals[[ii]][,'y'],df$sims[[ii]][['tree']][[1]]))
sim_pics<-sim_tip_vals<-tmp<-df$maps
#may want to make a "fast_sim" function or something that drops down to res 1 and only returns tips...
for(i in seq_along(tmp)){
  pars<-lapply(df[i,func],function(ii) ii[[1]]$par)
  tmp[[i]]<-make.traits(tmp[[i]],
                        simple_r~exp(pars[[2]][1]+pars[[2]][2]*x)+exp(pars[[2]][3]),
                        sweetspot_r~exp(pars[[3]][1]+pars[[3]][2]*x+pars[[3]][3]*x^2)+exp(pars[[3]][4]),
                        threshold_r~exp(pars[[4]][1])*plogis(pars[[4]][2]+pars[[4]][3]*x)+exp(pars[[4]][4]),
                        c(null_y,
                          simple_y,
                          sweetspot_y,
                          threshold_y)~
                          diffusion(null_y=exp(pars[[1]]),
                                    simple_y='simple_r',
                                    sweetspot_y='sweetspot_r',
                                    threshold_y='threshold_r'))
  sim_tip_vals[[i]]<-get.tip.vals(tmp[[i]])[,c('null_y','simple_y','sweetspot_y','threshold_y'),,drop=FALSE]
  sim_pics[[i]]<-apply(sim_tip_vals[[i]],c(2,3),function(ii) pic(ii,tmp[[i]][['tree']][[1]]))
  cat(i,'\n')
}
smoothScatter(log(as.vector(sim_pics[[8]][,3,])^2)~rep(log(pics[[8]]^2),100),nbin=c(50,100))
cor(log(as.vector(sim_pics[[25]][,3,])^2),rep(log(pics[[25]]^2),100))
test<-do.call(rbind,lapply(asplit(sim_pics[[25]],3),function(ii) apply(ii,2,function(jj) cor(log(jj^2),log(pics[[25]]^2)))))
hist(test[,2])

test2<-do.call(rbind,lapply(asplit(sim_pics[[42]],3),function(ii) apply(ii,2,function(jj) cor(log(jj^2),log(pics[[42]]^2)))))
hist(test2[,4])
cor(log(as.vector(sim_pics[[1]][,3,])^2),rep(log(pics[[1]]^2),100))
#there's definitely something here, I think-->true variables tend to yield correlations on the order of 0.05 to 0.2
#false one tend to remain <0.05 so far...
#Wow! discrimination based on distribution of correlation coefficients seems to work very well!
#let's scale it up
test<-matrix(FALSE,length(sim_pics),4)


#let's look at parameter estimates...
xx<-seq(-2,2,length.out=100)

null<-lapply(df$null[df$foc&df$func=='null'],'[[','par')
null.tmp<-do.call(cbind,lapply(null,function(ii) exp(rep(ii,length(xx)))))
matplot(x=xx,null.tmp,lty=1,type='l',col=rgb(0,0,0,0.5))
lines(rep(exp(2.3),length(xx))~xx,lwd=4)

simple<-lapply(df$simple[df$foc&df$func=='simple'],'[[','par')
simple.tmp<-do.call(cbind,lapply(simple,function(ii) exp(ii[1]+ii[2]*xx)+exp(ii[3])))
matplot(x=xx,simple.tmp,lty=1,type='l',col=rgb(0,0,0,0.5),ylim=c(0,100))
lines(exp(1+1*xx)+exp(-2)~xx,lwd=4)

sweetspot<-lapply(df$sweetspot[df$foc&df$func=='sweetspot'&df$ntax==120],'[[','par')
sweetspot.tmp<-do.call(cbind,lapply(sweetspot,function(ii) exp(ii[1]+ii[2]*xx+ii[3]*xx^2)+exp(ii[4])))
matplot(x=xx,sweetspot.tmp,lty=1,type='l',col=rgb(0,0,0,0.5),ylim=c(0,30))
lines(exp(2.5+2*xx-2*xx^2)+exp(-2)~xx,lwd=4)

threshold<-lapply(df$threshold[df$foc&df$func=='threshold'&df$ntax==120],'[[','par')
threshold.tmp<-do.call(cbind,lapply(threshold,function(ii) exp(ii[1])*plogis(ii[2]+ii[3]*xx)+exp(ii[4])))
matplot(x=xx,threshold.tmp,lty=1,type='l',col=rgb(0,0,0,0.5))
lines(exp(3)*plogis(-1+4*xx)+exp(-2)~xx,lwd=4)

#looks like it's doing the job!

tmp<-tmp<2
colnames(tmp)<-paste0(func,'_AIC_test')
df<-cbind(df,tmp)
table(df$null_AIC_test,func=df$func,foc=df$foc,noise=df$noise)
table(df$simple_AIC_test,func=df$func,foc=df$foc,noise=df$noise)
table(df$sweetspot_AIC_test,func=df$func,foc=df$foc,noise=df$noise)
table(df$threshold_AIC_test,func=df$func,foc=df$foc,noise=df$noise)
#proposed updates: variance gamma process for noise
#                  also re-jigger the functions to make them more distinguishable?
tmp<-do.call(cbind,lapply(df[,func],function(ii) do.call(rbind,lapply(ii,'[[','par'))))
plot(tmp[df$func=='threshold'&df$foc,'alpha[1]'])
abline(h=3)
plot(tmp[df$func=='sweetspot'&df$foc,'beta[3]'],ylim=c(-10,0))
abline(h=-4)
plot(tmp[df$func=='simple'&df$foc,'beta[2]'])
abline(h=1)
plot(tmp[df$func=='simple'&df$foc,'beta[1]'])
abline(h=2)


AICs<-df[,20:23]
AICs<-exp(-0.5*AICs)
AICs<-sweep(AICs,1,rowSums(AICs),'/')
AICs<-t(AICs[!df$foc,])
blah<-2
AICs<-AICs[,order(AICs[blah,])]
barplot(AICs[c(blah,seq_len(4)[-blah]),],border=NA,col=hcl.colors(4)[c(blah,seq_len(4)[-blah])],space=0)
