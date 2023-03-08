#function to put in make.traits formulae to threshold traits!
#' @export
threshold<-function(x,breaks=NULL,nbreaks=10,
                    state.names=NULL,numeric.names=FALSE,
                    ...){
  #nbreaks processing
  if(is.null(breaks)){
    if(is.null(nbreaks)){
      stop('PLACEHOLDER ABOUT NO BREAKS SPECIFIED')
    }
    ran<-range(x,na.rm=TRUE)
    breaks<-seq(ran[1],ran[2],length.out=nbreaks+2)[-c(1,nbreaks+2)]
  }
  #naming
  term.nm<-NA
  if(!is.null(state.names)){
    nnms<-length(state.names)
    nbreaks<-length(breaks)
    if(nnms==nbreaks){
      names(breaks)<-state.names
    }else{
      if(nnms==nbreaks+1){
        names(breaks)<-state.names[-nnms]
        term.nm<-state.names[nnms]
      }else{
        stop('PLACEHOLDER ABOUT STATE.NAMES BEING OF WRONG LENGTH')
      }
    }
  }
  nms<-names(breaks)
  if(is.null(nms)){
    nms<-rep(NA,length(breaks))
  }
  ord<-order(breaks)
  breaks<-breaks[ord]
  nms<-nms[ord]
  nms<-c(nms,term.nm) #for last, always unnamed, group
  prob.nms<-!nzchar(nms)|is.na(nms)
  len<-length(nms)
  if(any(prob.nms)){
    if(numeric.names){
      nn<-ceiling(log(len,10))
      tmp<-do.call(expand.grid,rep(list(as.character(c(0,seq_len(9)))),nn))
    }else{
      nn<-ceiling(log(len,26))
      tmp<-do.call(expand.grid,rep(list(LETTERS),nn))
    }
    tmp<-as.matrix(tmp)[seq_len(len),,drop=FALSE]
    nms[prob.nms]<-do.call(paste,c(rev(asplit(tmp,2)),sep=''))[prob.nms]
  }
  #output
  out<-findInterval(x,breaks,...)+1
  attr(out,"state_names")<-nms
  attr(out,"breaks")<-breaks
  out
}