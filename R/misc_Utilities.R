#gets treeIDs
.get.treeID<-function(contsimmap){
  nms<-dimnames(contsimmap)[[3]]
  as.numeric(substr(nms,5,regexpr('_',nms)-1))
}

#gets specific elements of maps attribute
#can "uncompress" maps to match indexing of contsimmap object
.get.maps<-function(contsimmap,element=c("coarse","ts","dts","bb.dts","bb.sds","state","incl","inds"),uncompress=FALSE){
  element<-pmatch(element[1],c("coarse","ts","dts","bb.dts","bb.sds","state","incl","inds"))
  out<-attr(contsimmap,'maps')[,,element,drop=FALSE]
  out<-matrix(out,dim(out)[1],dim(out)[2])
  if(uncompress){
    tmp<-dimnames(contsimmap)[[1]]
    nodes.ind<-grepl('N',tmp)
    roots.ind<-tmp=='N0'
    roots.value<-NULL
    if(any(roots.ind)){
      if(element==6){
        #implicitly picks first edge and sticks with that--might cause some issues in weird edge cases
        #for example, if descending edges from node have different initial states
        #that's impossible under a normal markov model, but in the case of regime mapping things could get weird...
        if(inherits(contsimmap,"summarized_contsimmap")){
          roots.value<-rep(lapply(out[root.edges(contsimmap)[1],],function(ii) ii[1,,drop=FALSE]),each=sum(roots.ind))
        }else{
          roots.value<-rep(lapply(out[root.edges(contsimmap)[1],],'[[',1),each=sum(roots.ind))
        }
      }else if(element==7){
        roots.value<-list(TRUE)
      }else if(element==8){
        roots.value<-list(1)
      }else{
        roots.value<-list(0)
      }
    }
    tmp[nodes.ind]<-gsub('N','',tmp[nodes.ind])
    out<-out[match(tmp,seq_len(Nedge(contsimmap))),,drop=FALSE]
    out[nodes.ind,]<-lapply(out[nodes.ind,],function(ii) ii[length(ii)])
    if(!is.null(roots.value)){
      out[roots.ind,]<-roots.value
    }
    out<-out[,.get.treeID(contsimmap),drop=FALSE]
  }
  out
}

#gets either nts (number of timepoints per edge) or NTs (number of critical timepoints per edge)
#can "uncompress" to match indexing of contsimmap object
.get.ns<-function(contsimmap,type=c("nts","NTs"),uncompress=FALSE){
  type<-pmatch(type[1],c("nts","NTs"))
  out<-lengths(.get.maps(contsimmap,c("ts","coarse")[type]))
  if(uncompress){
    tmp<-dimnames(contsimmap)[[1]]
    matches<-match(tmp,seq_len(Nedge(contsimmap)))
    treeID<-.get.treeID(contsimmap)
    out<-out[matches,treeID,drop=FALSE]
    out[is.na(matches),]<-1
    out[is.na(tmp),]<-0
    out[,is.na(treeID)]<-0
  }
  out
}

.stored.nodes<-function(contsimmap){
  as.numeric(gsub('N','',unique(grep('N',dimnames(contsimmap)[[1]],value=TRUE))))
}

.report.names<-function(prefix=NULL,nms,suffix=NULL,printlen=Inf,combine='and'){
  nn<-length(nms)
  if(nn==1){
    paste(prefix[1],
          nms,
          suffix[1])
  }else if(nn>printlen){
    paste(prefix[2],
          paste(c(nms[seq_len(printlen)],'...'),collapse=', '),
          suffix[2])
  }else{
    if(nn==2){
      paste(prefix[2],
            paste(nms,collapse=paste0(' ',combine,' ')),
            suffix[2])
    }else{
      paste0(prefix[2],' ',
             paste(nms[-nn],collapse=', '),', ',combine,' ',nms[nn],
             ' ',suffix[2])
    }
  }
}

.prog<-function(cur){
  paste0('|',paste(rep('=',cur),collapse=''),paste(rep(' ',100-cur),collapse=''),'| ',cur,'%',if(cur<100) '\r')
}
