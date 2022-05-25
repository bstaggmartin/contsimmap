#reordering test: how to order anc to ensure ancestral come before descending edges?
tree<-pbtree(n=1000)
tree<-reorder(tree,'pruning')
anc<-evorates:::anc.edges(tree)
des<-evorates:::des.edges(tree)
anc[lengths(anc)==0]<-0
anc<-unlist(anc)
out<-rep(FALSE,length(anc))
test.loop<-function(seq){
  for(i in seq){
    anc.ind<-anc[i]
    if(!anc.ind){
      out[i]<-TRUE
    }else{
      if(out[anc.ind]){
        out[i]<-TRUE
      }
    }
  }
  all(out)
}
test.loop(1:length(anc)) #doesn't work of course
test.loop(order(anc)) #also doesn't work
debug(test.loop)
test.loop(order(anc))
rbind(1:length(anc),anc)
#you have to make sure the above row's numbers come AFTER the below row's numbers, but you have to keep the couplings...
tmp<-cbind(anc,1:length(anc)) #now left col's numbers need to occur BEFORE right col's
tmp[order(tmp[,1],tmp[,2]),] #no...maybe a for loop?
#overly-complicated--I think this might work? Find which edge corresponds to the descendant of current edge, and make sure thing are ordered
#such that these edges come after the other? I can't make sense of this, but this intutively seems to make sense...
tmp[order(match(tmp[,2],tmp[,1])),]
#nope, definitely not...
#so left is index of ancestral edge, right is index of descendant edge
#you want to ensure ancestors come before descendants...
#maybe start out smaller... set n from 50 to 5
#let's try for loop again
ord<-order(anc)
ord<-list(NULL,ord)
while(length(ord[[2]])){
  foc<-ord[[2]][[1]]
  inds1<-which(anc==foc)
  ord[[1]]<-c(ord[[1]],foc)
  ord[[2]]<-ord[[2]][-1]
  if(length(inds1)){
    inds2<-ord[[1]]%in%inds1
    inds3<-ord[[2]]%in%inds1
    ord[[1]]<-c(ord[[1]][!inds2],inds1)
    ord[[2]]<-ord[[2]][!inds3]
  }
}
test.loop(ord) #works! not particularly efficient, but it works!
#unfortunately the while loop, while clever, is slower than your naive implementation...
ord<-order(anc)
#this was the fastest I could come up with
for(i in seq_along(ord)){
  inds1<-which(anc==ord[i])
  tmp.len<-length(inds1)
  if(tmp.len){
    inds2<-which(ord%in%inds1)
    if(inds2[1]>i) inds2<-inds2+tmp.len
    ord<-append(ord,inds1,i)[-inds2]
  }
}
#oh duh
ord<-reorder(tree,index.only=TRUE)
test.loop(ord)
microbenchmark::microbenchmark(reorder(tree,index.only=TRUE))
#awks, by far the best
#oh damn, but this requires having the tree's edges subset...I think the slower function might be the better choice after all

ord<-order(anc)
#this was the fastest I could come up with
for(i in seq_along(ord)){
  inds<-which(anc==ord[i])
  if(length(inds)){
    pos<-match(inds,ord)
    after<-i-sum(pos<i)
    ord<-append(ord[-pos],inds,after)
  }
}
ord
test.loop(ord)

#using descendants alg should work I think...
des<-unlist(lapply(seq_along(des),function(ii) if(length(des[[ii]])) c(ii,des[[ii]])))
test.loop(des)

#nope--need to revisit how I get des...
len<-length(anc)
tmp.out<-rep(list(integer(0)),len)
inds<-!!anc
tmp.anc<-anc[inds]
des<-split((1:len)[inds],tmp.anc)
# nms<-names(des)
# des<-unlist(lapply(nms,function(ii) c(as.numeric(ii),des[[ii]])))
test.loop(des) #damn, what's wrong?
#doesn't order the ancestral edges correctly, only the descendants immediately following...
len<-length(des)
counter<-len
for(i in rev(seq_along(des))){
  tmp.des<-des[counter]
  des<-des[-counter]
  tmp.anc<-as.numeric(names(tmp.des))
  pos<-match(tmp.anc,unlist(des,use.names=FALSE)) #wow! setting use.names to FALSE speeds things up a lot...
  if(is.na(pos)){
    append(des,tmp.des,0)
  }else{
    pos<-as.numeric(rep(names(des),lengths(des))[pos])
    append(des,tmp.des,pos)
  }
}
#this is actually just as complicated if not more...

#takes all subtrees grabs 1st one, tacks on any descending subtrees to right if available, tacks on 1st subtree in grab bag is no descendants available
#rinse and repeat...
names(des)<-seq_along(des)
des2<-des[lengths(des)>0]
grabbag<-des2
out<-des2[1]
grabbag<-grabbag[-1]
foc<-1
while(length(grabbag)){
  tmp.des<-out[[foc]]
  inds<-match(tmp.des,as.numeric(names(grabbag)))
  inds<-inds[!is.na(inds)]
  if(length(inds)){
    out<-append(out,grabbag[inds],foc)
    grabbag<-grabbag[-inds]
    foc<-foc+1
  }else{
    out<-append(out,grabbag[1],0)
    grabbag<-grabbag[-1]
    foc<-1
  }
}
out<-unlist(out,use.names=FALSE)
tmp<-as.numeric(names(des2))
tmp<-tmp[is.na(match(tmp,out))]
test.loop(c(tmp,out))
#yeah, I think there's some way to speed up the algorithm here...but it's trickier than I anticipated...
#the above is a neat one, but still slower than current algorithm