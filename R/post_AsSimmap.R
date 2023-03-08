#convert contSimmap to simmap via thresholding
#for simplicity, let's make it only do 1 trait at a time for now...
#... is extra arguments to pass to threshold()
#now should work, even when multiple thresholds are crossed in a single time interval
#' @export
as.simmap<-function(contsimmap,trait=1,...){
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(.stored.nodes(contsimmap))==1){
      stop('The as.simmap() function will support single subtrees in the future, but not yet')
    }else{
      stop('The as.simmap() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  
  #make appropriate formula
  #no NA traits for now!
  if(is.character(trait)){
    trait<-pmatch(trait,dimnames(contsimmap)[[2]])
  }
  trait<-dimnames(contsimmap)[[2]][trait]
  trait<-trait[!is.na(trait)]
  if(!length(trait)) trait<-dimnames(contsimmap)[[2]] #need to account for possibility of NA trait dimension eventually...
  trait<-trait[1]
  contsimmap<-contsimmap[seq_len(Nedge(contsimmap)),trait,]
  ex.args<-paste(paste(names(list(...)),list(...),sep='='),collapse=',')
  form<-as.formula(paste0('new_',trait,'~threshold(',trait,',',ex.args,')'))
  contsimmap<-make.traits(contsimmap,form)
  tmp<-attr(contsimmap,'params')[['call_info',attr(contsimmap,'traits')[paste0('new_',trait)]]]
  brk<-tmp[['breaks']]
  nms<-tmp[['state_names']]
  
  
  
  ##THIS PART COULD BE BETTER CLEANED UP, BUT IT SEEMS TO WORK WELL AND FAST FOR NOW##
  #main part--getting maps element
  #could make a custom method for more efficiently converting edges to the granular format in the case of plot and as.simmap, for sure
  #don't need to keep replicating the time and state vectors necessarily, for example!
  #initialization stuff
  lens<-.get.ns(contsimmap,'nts',uncompress=TRUE)[-1,,drop=FALSE]+1
  tmp<-contsimmap[,1,]
  tmp<-tmp[seq_along(tmp)]
  values<-unlist(lapply(tmp,'[[','values'),use.names=FALSE)
  ts<-unlist(lapply(tmp,'[[','ts'),use.names=FALSE)
  tmp<-contsimmap[,2,]
  tmp<-tmp[seq_along(tmp)]
  states<-unlist(lapply(tmp,'[[','values'),use.names=FALSE)
  #locating endpoints of shifts across all edges/sims
  #only incl indices not corresponding to beginning of edge
  incl<-rep(TRUE,sum(lens))
  incl[cumsum(lens[-length(lens)])+1]<-FALSE
  ends<-which(c(FALSE,states[-1]!=states[-length(states)])&incl)
  starts<-ends-1
  x0<-values[starts]
  x1<-values[ends]
  t0<-ts[starts]
  t1<-ts[ends]
  s0<-states[starts]
  s1<-states[ends]
  #calculating sign/number of threshold crosses in each shift
  cross<-s1-s0
  polarity<-sign(cross)
  ncross<-abs(cross)
  #taking care of multiple crosses
  #replicate s0 accordingly
  #renders s1 useless since each cross is now only 1 unit
  s0<-rep(s0,ncross)
  nncross<-ncross>1
  tmp<-which(nncross)
  tmp.seqs<-unlist(lapply(tmp,function(ii) polarity[ii]*seq_len(ncross[ii]-1)),use.names=FALSE)
  tmp<-rep(tmp,ncross[nncross]-1)
  tmp<-tmp+seq_along(tmp)
  s0[tmp]<-s0[tmp]+tmp.seqs
  #replicate x0, x1, t0, t1, polarity accordingly
  x0<-rep(x0,ncross)
  x1<-rep(x1,ncross)
  t0<-rep(t0,ncross)
  t1<-rep(t1,ncross)
  polarity<-rep(polarity,ncross)
  #get actual thresholds and shift pts...
  thresholds<-brk[s0-(polarity<0)]
  shift.ts<-(thresholds-x0)/(x1-x0)*(t1-t0)+t0
  shift.states<-nms[s0]
  #now the trick is to resplit...
  #recalculate lens matrix...
  tmp<-findInterval(rep(ends,ncross),c(0,cumsum(lens)),left.open=TRUE)
  lens<-lens+tabulate(tmp,nbins=length(lens))
  #replicate states/times accordingly...
  #add a flag to states for the rle step...
  nreps<-rep(1,length(states))
  nreps[starts]<-ncross+1
  #need to insert shift times and states in the right places--tricky!
  #maybe this is the most efficient way? Hard to say...
  tmp<-list(lengths=as.vector(rbind(1,nreps-1)),
            values=rep(c(FALSE,TRUE),length(states)))
  tmp<-inverse.rle(tmp)
  states<-rep(nms[states],nreps)
  states[tmp]<-shift.states
  ts<-rep(ts,nreps)
  ts[tmp]<-shift.ts
  #add a flag for edge index so rle looks for runs of BOTH states and edge indices
  states<-paste0(states,'@@@',rep(seq_along(lens),lens))
  states<-rle(states)
  #form the unlisted maps along with their edge indices
  splits<-strsplit(states[['values']],'@@@')
  maps<-setNames(ts[cumsum(states[['lengths']])],
                 unlist(lapply(splits,'[[',1),use.names=FALSE))
  splits<-as.numeric(unlist(lapply(splits,'[[',2),use.names=FALSE))
  #convert from time points to time intervals
  #first subtract time of node ancestral to edges
  maps<-maps-edge.ranges(contsimmap)[(splits-1)%%Nedge(contsimmap)+1,1]
  #then subtract the preceding times--but exclude the first entry for each edge!
  tmp<-rep(0,length(splits))
  inds<-which(c(FALSE,splits[-1]==splits[-length(splits)]))
  indsm1<-inds-1
  tmp[inds]<-maps[indsm1]
  maps<-maps-tmp
  #form the extras
  #form mapped.edge
  edges<-(splits-1)%%Nedge(contsimmap)+1
  sims<-(splits-1)%/%Nedge(contsimmap)+1
  mapped.edge<-tapply(maps,list(edges,names(maps),sims),sum)
  mapped.edge[is.na(mapped.edge)]<-0
  edge.mat<-attr(contsimmap,'tree')[[1]][['edge']]
  dimnms<-list(paste0(edge.mat[,1],edge.mat[,2],sep=","),nms)
  #form node.states
  inds<-c(splits[-1]!=splits[-length(splits)],TRUE)
  node.states<-matrix(names(maps)[inds],Nedge(contsimmap),dim(contsimmap)[3])
  #need to get first state in each root edge
  #should be fine since contsimmap construction ENSURES ancestral node for each root edge will be identical
  root.edge<-root.edges(attr(contsimmap,'tree')[[1]])[1]
  inds<-edges==root.edge&c(TRUE,splits[-1]!=splits[-length(splits)])
  root.states<-names(maps)[inds]
  tmp<-match(edge.mat[,1],edge.mat[,2])
  node.states<-aperm(array(c(node.states[tmp,,drop=FALSE],node.states),c(dim(node.states),2)),c(1,3,2))
  inds<-is.na(tmp)
  n.inds<-sum(inds)
  node.states[inds,1,]<-rep(root.states,each=n.inds)
  #form states
  states<-matrix(node.states[tip.edges(contsimmap),2,],Ntip(contsimmap),dim(contsimmap)[3])
  rownames(states)<-attr(contsimmap,'tree')[[1]][['tip.label']]
  #finally split into the final maps
  maps<-matrix(split(maps,splits),Nedge(contsimmap),dim(contsimmap)[3])
  
  
  
  #form output
  out<-attr(contsimmap,'tree')[[1]]
  class(out)<-c('simmap','phylo')
  out[['maps']]<-maps[,1]
  nedges<-Nedge(contsimmap)
  nstates<-dim(mapped.edge)[2]
  out[['mapped.edge']]<-out[['node.states']]<-out[['states']]<-out[['Q']]<-out[['logL']]<-NULL
  out<-rep(list(out),dim(contsimmap)[3])
  for(i in seq_along(out)){
    out[[i]][['maps']]<-maps[,i]
    out[[i]][['mapped.edge']]<-matrix(mapped.edge[,,i,drop=FALSE],c(nedges,nstates),dimnames=dimnms)
    out[[i]][['node.states']]<-matrix(node.states[,,i,drop=FALSE],c(nedges,2))
    out[[i]][['states']]<-states[,i]
  }
  class(out)<-c('multiSimmap','multiPhylo')
  out
}