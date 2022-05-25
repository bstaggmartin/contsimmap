##12/6 after a BUNCH of experiments, I have determined that the best strategy indeed seems to be simply be using the native R sorting function...
##Some helpful nuggets in here, though, including a multi-append function

x<-seq(0,100,1)
values<-sort(c(runif(5,0,100),100,50.5,50.6))
after<-findInterval(values,x,all.inside=TRUE)
values<-c(split(values,after),list(NULL))
after<-unique(after)
#need to make sure after is greater than 1 and less than length(x)!
mappend<-function(x,values,after,presorted=FALSE){
  values<-rep(values,length.out=length(after))
  pre<-after==0
  npre<-sum(pre)
  if(npre){
    x<-c(values[pre],x)
    values<-values[!pre]
    after<-after[!pre]
    after<-after+npre
  }
  reps<-tabulate(after,length(x))+1
  if(!presorted){
    after<-sort.int(after,index.return=TRUE)
    values<-values[after$ix]
    after<-after$x
  }
  after<-after+seq_along(after)
  x<-rep(x,reps)
  x[after]<-values
  x
}
microbenchmark::microbenchmark(mappend(vec,'here!',c(0,3,5,10,40),presorted=TRUE)) #shaved off 10 microseconds with 'rep and replace' method above
#crazy fast if you assume after is presorted

mappend(x,values,after)
microbenchmark::microbenchmark(sort(c(values,x),method='quick'))
microbenchmark::microbenchmark(findInterval(values,x)) #much faster
microbenchmark::microbenchmark(mappend(x,values,findInterval(values,x,all.inside=TRUE))) #only marginally faster then quick sort--but I'll take it!
#preserves ties and all that; not sure about names though
names(values)<-letters[1:5]
mappend(x,values,findInterval(values,x)) #preserves names!
#basically the same speed as quick sort but stable, name-preserving sorting
#quick partial sort
#assumes values and x are already sorted and that no values is below or above all entries in x
qpartsort<-function(values,x){
  nx<-length(x)
  after<-findInterval(values,x)
  values<-c(split(values,after),list(NULL))
  after<-unique(after)
  inds<-matrix(c(1,after+1,after,nx),ncol=2)
  #makes sure ties for last entry are sorted correctly
  ninds<-nrow(inds)
  if(after[ninds-1]==nx){
    inds<-inds[-ninds,,drop=FALSE]
    values<-values[-ninds]
  }
  do.call(c,
          lapply(seq_along(values),
                 function(i) 
                   c(x[seq.int(inds[i,1],inds[i,2])],values[[i]])
                 )
          )
}
#might as well assume that single last tie is there for your case...
qpartsort<-function(values,x){
  nvals<-length(values)
  x<-c(x,values[nvals])
  values<-values[-nvals]
  after<-findInterval(values,x)
  tmp<-unique(after)
  inds<-matrix(c(1,tmp+1,tmp,length(x)),ncol=2)
  do.call(c,
          lapply(seq_len(length(tmp)+1),
                 function(i) 
                   c(x[seq.int(inds[i,1],inds[i,2])],values[after==inds[i,2]])
          )
  )
}
microbenchmark::microbenchmark(sort(c(sort(c(runif(5),1)),seq(0,1,0.001))))
microbenchmark::microbenchmark(qpartsort(sort(c(runif(5),1)),seq(0,1,0.001))) #a bit slower here--must be the split and unique operations...
#still I think the modest speed loss is worth tie preserving--it's nice to know where the ties will land!
qpartsort(setNames(c(99.1,99.2,100),c('here!','here!','here!')),seq(0,100))
#getting rid of split made it a better rival, but it still is outpaced by a normal sort
#one last strategy
qpartsort<-function(values,x){
  out<-c(x,values)
  tmp<-findInterval(out,x)
  out[sort.int(tmp,index.return=TRUE)$ix]
}
microbenchmark::microbenchmark(sort(c(sort(c(runif(7),1)),seq(0,1,0.01))))
microbenchmark::microbenchmark(qpartsort(sort(c(runif(7),1)),seq(0,1,0.01))) #speedup with only a small handful of additional elements...
qpartsort(setNames(c(99.1,99.2,100),c('here!','here!','here!')),seq(0,100))
qpartsort<-function(values,x){
  sort(c(x,values))
}
qpartsort(setNames(c(99.1,99.2,100),c('here!','here!','here!')),seq(0,100))

#one last thing--partial sorting while preserving names?
#only speeds up when length of x is modest...
qpartsort<-function(values,x){
  out<-sort(c(values,x),partial=seq_along(values))
  out
}
#nonetheless, let's try...
qpartsort<-function(values,x){
  after<-findInterval(values,x)
  out<-sort(c(values,x),partial=seq_along(values)) #this doesn't even work! consistently places first one 1 index right of where it should be...
  names(out)[after+seq_along(after)]<-names(values)
  out
}
#yep, sort's the fastest--partial might be faster if it worked
#even then, it's only faster under pretty restrictive conditions...

x<-seq(0,1,length.out=100)
values<-sort(setNames(runif(5),rbinom(5,1,0.5)))
microbenchmark::microbenchmark({out<-sort(c(x,values),index.return=TRUE);out$x>length(x)})
microbenchmark::microbenchmark({out<-sort(c(x,values));nzchar(names(out))}) #faster
out<-sort(c(x,values));nzchar(names(out))
