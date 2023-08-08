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

#could probably make more use of this function...
.stored.nodes<-function(contsimmap){
  as.numeric(gsub('N','',unique(grep('N',dimnames(contsimmap)[[1]],value=TRUE))))
}

.report.names<-function(prefix=NULL,nms,suffix=NULL,printlen=Inf,combine='and'){
  nn<-length(nms)
  if(nn==1){
    out<-paste(prefix[1],
               nms,
               suffix[1])
  }else if(nn>printlen){
    out<-paste(prefix[2],
               paste(c(nms[seq_len(printlen)],'...'),collapse=', '),
               suffix[2])
  }else{
    if(nn==2){
      out<-paste(prefix[2],
                 paste(nms,collapse=paste0(' ',combine,' ')),
                 suffix[2])
    }else{
      out<-paste0(prefix[2],' ',
                  paste(nms[-nn],collapse=', '),', ',combine,' ',nms[nn],
                  ' ',suffix[2])
    }
  }
  trimws(out)
}


.nice.array.print<-function(x,def.nm,nrows,ncols,nslices,...){
  if(!is.list(x)){
    if(length(dim(x))<3){
      x<-list(x)
    }else{
      x<-asplit(x,3)
    }
  }else{
    if(length(dim(x))==2){
      x<-asplit(x,1)
    }
  }
  tmp<-length(x)-nslices
  if(tmp>0){
    x<-x[1:nslices]
    slice.message<-paste0("\n\n\t\t [ omitted ",
                          tmp,
                          if(tmp>1) " slices" else " slice",
                          " ] ")
  }else{
    slice.message<-NULL
  }
  nms<-names(x)
  if(is.null(nms)){
    nms<-rep(NA,length(x))
  }
  nas<-which(is.na(nms))
  nms<-paste0("$",nms)
  nms[nas]<-paste0(def.nm,"[[",nas,']]')
  for(i in seq_along(x)){
    tmp<-x[[i]]
    twod<-length(dim(tmp))==2
    tmp.messages<-NULL
    tmp.rows<-NROW(tmp)-nrows
    if(tmp.rows>0){
      tmp<-if(twod) tmp[1:nrows,,drop=FALSE] else tmp[1:nrows,drop=FALSE]
      tmp.messages<-c(tmp.messages,paste0(tmp.rows,if(tmp.rows>1) " rows" else " row"))
    }
    if(twod){
      tmp.cols<-ncol(tmp)-ncols
      if(tmp.cols>0){
        tmp<-tmp[,1:ncols,drop=FALSE]
        tmp.messages<-c(tmp.messages,paste0(tmp.cols,if(tmp.cols>1) " columns" else " column"))
      }
    }
    if(length(tmp.messages)){
      rowcol.message<-paste0("\t\t [ omitted ",.report.names(nms=tmp.messages)," ] ")
    }else{
      rowcol.message<-NULL
    }
    if(length(x)>1){
      cat("\n\n\t\t",nms[i],"\n",sep="")
    }else{
      cat("\n\n")
    }
    writeLines(paste0("\t\t",capture.output(print(tmp,max=1e8,...))))
    cat(rowcol.message)
  }
  cat(slice.message)
}

.prog<-function(cur){
  paste0('|',paste(rep('=',cur),collapse=''),paste(rep(' ',100-cur),collapse=''),'| ',cur,'%',if(cur<100) '\r')
}
