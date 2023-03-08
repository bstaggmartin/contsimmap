
####INDEXING####

.get.nms<-function(i,nms,multi.index=TRUE){
  if(length(dim(i))==2&NCOL(i)==3&!is.logical(i)){
    if(is.character(i)){
      tmp<-lapply(seq_len(3),function(ii) match(i[,ii],nms[[ii]]))
    }else{
      mode(i)<-'numeric'
      tmp<-lapply(seq_len(3),function(ii) i[,ii])
    }
  }else{
    dims<-lengths(nms)
    if(is.logical(i)){
      if(multi.index) i<-rep(i,length.out=prod(dims))
      i[!is.na(i)]<-which(i)
    }else{
      i<-as.numeric(i)
      if(min(i,na.rm=TRUE)<0){
        i<-seq_len(prod(dims))[i]
      }
    }
    i<-i[i!=0]
    tmp<-list((i-1)%%dims[1]+1,
              ((i-1)%/%dims[1])%%dims[2]+1,
              (i-1)%/%prod(dims[c(1,2)])+1)
  }
  lapply(seq_len(3),function(ii) nms[[ii]][tmp[[ii]]])
}

.sub.anc.edges<-function(phy,sub=NULL){
  if(is.null(sub)){
    match(phy[['edge']][,1],phy[['edge']][,2],nomatch=0L)
  }else{
    match(phy[['edge']][sub,1],phy[['edge']][,2],nomatch=0L)
  }
}

.check.seq<-function(dp.x,dp.i,only.rows=FALSE){
  if(only.rows){
    patterns<-paste0(c('dim','(nrow|NROW)'),
                     '\\(',dp.x,'\\)')
    replacements<-c('dims','dims[1]')
  }else{
    patterns<-paste0(c('seq_along','length','dim','(nrow|NROW)'),
                     '\\(',dp.x,'\\)')
    replacements<-c('seq_len(prod(dims))','prod(dims)','dims','dims[1]')
  }
  checks<-unlist(lapply(patterns,grepl,x=dp.i),use.names=FALSE)
  if(any(checks)){
    for(i in which(checks)){
      dp.i<-gsub(patterns[i],replacements[i],dp.i)
    }
    dp.i
  }else{
    NULL
  }
}

#various attribute information could be dropped for efficiency/space reasons, but I haven't looked into this too deeply yet
#now deals with NAs like a champ in the cases tested so far...
#Important difference from base R behavior --> NULLs treated like missing arguments!
#' @export
`[.contsimmap`<-function(x,i,j,k){
  if(missing(i)) i<-NULL
  if(missing(j)) j<-NULL
  if(missing(k)) k<-NULL
  ijk<-list(i,j,k)
  nodes<-grepl('N',dimnames(x)[[1]])
  if(nargs()==2){
    if(is.null(ijk[[1]])) x else{
      dims<-dim(x)
      dims[1]<-dims[1]-sum(nodes)
      tmp<-.check.seq(deparse(substitute(x)),deparse(substitute(i)))
      if(!is.null(tmp)) ijk[[1]]<-eval(str2lang(tmp))
      tmp.out<-unclass(x)[!nodes,,,drop=FALSE][ijk[[1]]]
      nms<-dimnames(x)
      nms[[1]]<-nms[[1]][!nodes]
      ijk<-.get.nms(ijk[[1]],nms)
      anc<-as.character(.sub.anc.edges(attr(x,'tree')[[1]],as.numeric(ijk[[1]])))
      anc[is.na(ijk[[1]])]<-NA
      anc.inds<-tmp.anc<-unique(anc)
      ts<-attr(x,'ts')[tmp.anc]
      tmp<-tmp.anc=='0'
      ts[tmp]<-0
      names(ts)[tmp]<-'0'
      lens<-lengths(ts)
      tmp<-which(lens>1)
      ts[tmp]<-lapply(tmp,function(ii) ts[[ii]][lens[ii]])
      ts<-if(length(ts)) unlist(ts) else numeric(0)
      probs<-paste0('N',tmp.anc)%in%dimnames(x)[[1]][nodes]
      anc.inds[probs]<-paste0('N',tmp.anc[probs])
      j<-unique(ijk[[2]])
      j<-j[!is.na(j)]
      k<-unique(ijk[[3]])
      k<-k[!is.na(k)]
      nodes<-unclass(x)[anc.inds[!is.na(anc.inds)],j,k,drop=FALSE]
      dimnames(nodes)[[1]]<-tmp.anc[!is.na(anc.inds)]
      lens<-lengths(nodes)
      tmp<-which(lens>1)
      nodes[tmp]<-lapply(tmp,function(ii) nodes[[ii]][lens[ii]])
      nodes<-array(if(length(nodes)) unlist(nodes,use.names=FALSE) else numeric(0),
                   dim(nodes),dimnames(nodes))
      treeID<-as.numeric(substr(ijk[[3]],5,regexpr('_',ijk[[3]])-1))
      out<-list(values=NULL,ts=NULL,states=NULL)
      attr(out,'info')<-setNames(vector('character',3),c('edge','trait','sim'))
      out<-rep(list(out),length(tmp.out))
      summarized.flag<-inherits(x,"summarized_contsimmap")
      for(l in seq_along(tmp.out)){
        tmp.ijk<-unlist(lapply(ijk,'[',l),use.names=FALSE)
        if(!is.na(tmp.ijk[1])){
          if(is.na(tmp.ijk[3])){
            out[[l]][['ts']]<-attr(x,'ts')[[tmp.ijk[1]]]
          }else{
            out[[l]][c('ts','states')]<-attr(x,'maps')[tmp.ijk[1],treeID[l],c('ts','state')]
            out[[l]][['ts']]<-unname(c(ts[anc[l]],out[[l]][['ts']]))
            if(summarized.flag){
              out[[l]][['states']]<-rbind(out[[l]][['states']][1,,drop=FALSE],out[[l]][['states']])
            }else{
              out[[l]][['states']]<-c(out[[l]][['states']][1],out[[l]][['states']])
            }
            
          }
          if(any(is.na(tmp.ijk))){
            out[[l]][['values']]<-rep(NA,length(out[[l]][['ts']]))
          }else{
            out[[l]][['values']]<-c(nodes[anc[l],tmp.ijk[2],tmp.ijk[3]],tmp.out[[l]])
          }
        }
        attr(out[[l]],'info')[]<-tmp.ijk
      }
      out
    }
  }else{
    unmod<-unlist(lapply(ijk,is.null),use.names=FALSE)
    if(all(unmod)) x else{
      summarized.flag<-inherits(x,"summarized_contsimmap")
      dims<-dim(x)
      ijk[unmod]<-lapply(which(unmod),function(ii) seq_len(dims[ii]))
      dims[1]<-dims[1]-sum(nodes)
      if(unmod[1]) out<-unclass(x)[ijk[[1]],ijk[[2]],ijk[[3]],drop=FALSE] else{
        tmp<-.check.seq(deparse(substitute(x)),deparse(substitute(i)),only.rows=TRUE)
        if(!is.null(tmp)) ijk[[1]]<-eval(str2lang(tmp))
        out<-unclass(x)[!nodes,,,drop=FALSE][ijk[[1]],ijk[[2]],ijk[[3]],drop=FALSE]
        new.edges<-as.numeric(dimnames(out)[[1]])
        anc<-unique(.sub.anc.edges(attr(x,'tree')[[1]],new.edges[!is.na(new.edges)]))
        nodes<-dimnames(x)[[1]][nodes]
        nodes.inds<-as.numeric(gsub('N','',nodes))
        tmp<-nodes.inds%in%anc
        nodes<-nodes[tmp]
        nodes.inds<-nodes.inds[tmp]
        exist<-c(nodes.inds,new.edges)
        probs<-as.character(anc[!(anc%in%exist)])
        out.nodes<-aperm(unclass(x)[c(nodes,probs),ijk[[2]],ijk[[3]],drop=FALSE],c(2,3,1))
        lens<-lengths(out.nodes)
        tmp<-which(lens>1)
        out.nodes[tmp]<-lapply(tmp,function(ii) out.nodes[[ii]][lens[ii]])
        if(length(probs)){
          dimnames(out.nodes)[[3]]<-c(nodes,paste0('N',probs))
        }
        out<-aperm(
          array(c(out.nodes,aperm(out,c(2,3,1))),
                c(dim(out)[2:3],dim(out.nodes)[3]+dim(out)[1]),
                c(dimnames(out)[2:3],edge=list(c(dimnames(out.nodes)[[3]],dimnames(out)[[1]])))),
          c(3,1,2))
      }
      for(l in c('ts','tree','maps','traits','params')){
        attr(out,l)<-attr(x,l)
      }
      #just reconfigure params and traits accordingly
      #you technically don't need to name the traits object I realized...but I'll probably just keep things as is for now
      if(!unmod[2]){
        old.traits<-attr(out,'traits')
        new.traits<-unique(dimnames(out)[[2]])
        new.traits<-setNames(old.traits[new.traits],new.traits)
        traits.levs<-sort(unique(new.traits))
        new.params<-attr(out,'params')[,traits.levs,drop=FALSE]
        new.traits<-setNames(match(new.traits,traits.levs),names(new.traits))
        attr(out,'traits')<-new.traits
        attr(out,'params')<-new.params
      }
      if(!unmod[3]&!summarized.flag){
        nms<-dimnames(out)[[3]]
        matches<-match(nms,dimnames(x)[[3]])
        for(l in seq_along(attr(out,'params')['lookup',])){
          if(!is.null(attr(out,'params')[['lookup',l]])){
            attr(out,'params')[['lookup',l]]<-attr(out,'params')[['lookup',l]][matches,,drop=FALSE]
          }
        }
        treeID<-as.numeric(substr(nms,5,regexpr('_',nms)-1))
        for(l in unique(treeID)){
          inds<-which(treeID==l)
          nms[inds]<-paste0('tree',treeID[inds],'_sim',seq_len(length(inds)))
        }
        dimnames(out)[[3]]<-nms
      }
      if(summarized.flag){
        class(out)<-c('contsimmap','summarized_contsimmap')
      }else{
        class(out)<-'contsimmap'
      }
      out
    }
  }
}

#' @export
`[[.contsimmap`<-function(x,i,j,k){
  nodes<-grepl('N',dimnames(x)[[1]])
  out<-list('values'=NULL,
            'ts'=NULL,
            'states'=NULL)
  dims<-dim(x)
  dims[1]<-dims[1]-sum(nodes)
  if(nargs()==2){
    tmp<-.check.seq(deparse(substitute(x)),deparse(substitute(i)))
    if(!is.null(tmp)) i<-eval(str2lang(tmp))
    tmp.values<-unclass(x)[!nodes,,,drop=FALSE][[i]]
    nms<-dimnames(x)
    nms[[1]]<-nms[[1]][!nodes]
    ijk<-setNames(unlist(.get.nms(i,nms,multi.index=FALSE),use.names=FALSE),
                  c('edge','trait','sim'))
  }else{
    tmp<-.check.seq(deparse(substitute(x)),deparse(substitute(i)),only.rows=TRUE)
    if(!is.null(tmp)) i<-eval(str2lang(tmp))
    tmp.values<-unclass(x)[!nodes,,,drop=FALSE][[i,j,k]]
    ijk<-unlist(dimnames(unclass(x)[!nodes,,,drop=FALSE][i,j,k,drop=FALSE]))
  }
  if(!is.na(ijk[1])){
    if(is.na(ijk[3])){
      out[['ts']]<-attr(x,'ts')[[ijk[1]]]
    }else{
      out[c('ts','states')]<-attr(x,'maps')[ijk[1],as.numeric(substr(ijk[3],5,regexpr('_',ijk[3])-1)),c('ts','state')]
      anc<-.sub.anc.edges(attr(x,'tree')[[1]],as.numeric(ijk[1]))
      out[['ts']]<-unname(c(if(anc==0) 0 else attr(x,'ts')[[anc]][length(attr(x,'ts')[[anc]])],out[['ts']]))
      if(inherits(x,"summarized_contsimmap")){
        out[['states']]<-rbind(out[['states']][1,,drop=FALSE],out[['states']])
      }else{
        out[['states']]<-c(out[['states']][1],out[['states']])
      }
    }
    if(any(is.na(ijk))){
      out[['values']]<-rep(NA,length(out[['ts']]))
    }else{
      tmp<-unclass(x)[[grep(paste0('N',anc,'$|^',anc,'$'),dimnames(x)[[1]]),ijk[2],ijk[3]]]
      out[['values']]<-c(tmp[length(tmp)],tmp.values)
    }
  }
  attr(out,'info')<-ijk
  out
}

#' @export
`[<-.contsimmap`<-function(x,i,j,k,value){
  stop("There is no index assignment method for contsimmaps yet! Doing this would really break the code at this point; use make.traits() instead.")
}

#' @export
`[[<-.contsimmap`<-function(x,i,j,k,value){
  stop("There is no index assignment method for contsimmaps yet! Doing this would really break the code at this point; use make.traits() instead.")
}

####EXTRACTION####

.phy.method.contsimmap<-function(FUN,contsimmap,...){
  do.call(FUN,c(attr(contsimmap,'tree')[1],...))
}

#' @export
Ntip.contsimmap<-function(contsimmap){
  .phy.method.contsimmap(Ntip.phylo,contsimmap)
}

#' @export
Nedge.contsimmap<-function(contsimmap){
  .phy.method.contsimmap(Nedge.phylo,contsimmap)
}

#' @export
Nnode.contsimmap<-function(contsimmap,internal.only=TRUE){
  .phy.method.contsimmap(Nnode.phylo,contsimmap,internal.only)
}

#' @export
#' @method edge.ranges contsimmap
edge.ranges.contsimmap<-function(contsimmap){
  .phy.method.contsimmap(.edge.ranges,contsimmap)
}

#' @export
#' @method anc.edges contsimmap
anc.edges.contsimmap<-function(contsimmap){
  .phy.method.contsimmap(.anc.edges,contsimmap)
}

#' @export
#' @method des.edges contsimmap
des.edges.contsimmap<-function(contsimmap){
  .phy.method.contsimmap(.des.edges,contsimmap)
}

#' @export
#' @method sis.edges contsimmap
sis.edges.contsimmap<-function(contsimmap){
  .phy.method.contsimmap(.sis.edges,contsimmap)
}

#' @export
#' @method root.edges contsimmap
root.edges.contsimmap<-function(contsimmap){
  .phy.method.contsimmap(.root.edges,contsimmap)
}

#' @export
#' @method tip.edges contsimmap
tip.edges.contsimmap<-function(contsimmap,include.names=TRUE){
  .phy.method.contsimmap(.tip.edges,contsimmap,include.names)
}

#' @export
edge.inds<-function(contsimmap,unique.only=FALSE){
  out<-dimnames(contsimmap)[[1]]
  out<-as.numeric(out[!(grepl('N',out)|is.na(out))])
  if(unique.only) unique(out) else out
}

#' @export
nedge<-function(contsimmap,unique.only=TRUE){
  length(edge.inds(contsimmap,unique.only))
}

#' @export
tip.inds<-function(contsimmap,unique.only=FALSE){
  out<-attr(contsimmap,'tree')[[1]][['edge']][edge.inds(contsimmap,unique.only),2]
  out[out>Ntip(contsimmap)]<-NA
  out
}

#' @export
ntip<-function(contsimmap,unique.only=TRUE){
  sum(!is.na(tip.inds(contsimmap,unique.only)))
}

#' @export
node.inds<-function(contsimmap,unique.only=FALSE){
  attr(contsimmap,'tree')[[1]][['edge']][edge.inds(contsimmap,unique.only),,drop=FALSE]
}

#' @export
nnode<-function(contsimmap,unique.only=TRUE,internal.only=TRUE){
  out<-unique(as.vector(node.inds(contsimmap,unique.only)))
  if(internal.only) sum(out>Ntip(contsimmap)) else length(out)
}

####PRINT/SUMMARY####

#' @export
print.contsimmap<-function(contsimmap,printlen=6,...){
  traits<-dimnames(contsimmap)[[2]]
  nsims<-dim(contsimmap)[3]
  tree<-if(nsims==1) 'tree' else 'trees'
  cat(nsims,'phylogenetic',tree,'with',length(traits),'mapped continuous',
      if(length(traits)) .report.names(c('character:','characters:'),traits,printlen=printlen) else 'characters')
  edges<-edge.inds(contsimmap)
  flag<-FALSE
  if(length(edges)!=Nedge(contsimmap)){
    flag<-TRUE
  }
  #I don't think the below is necessary
  # else if(any(edges!=seq_len(Nedge(contsimmap)))){
  #   flag<-TRUE
  # }
  if(flag) cat('\t[ subsetted to ',if(length(edges)) .report.names(c('edge','edges'),nms=edges,printlen=printlen) else 'no edges ',']',sep='')
  invisible(flag)
}

#' @export
summary.contsimmap<-function(contsimmap,printlen=6,nrows=printlen,ncols=printlen,nslices=printlen,...){
  flag<-print(contsimmap,printlen=printlen)
  report.tree<-if(dim(contsimmap)[3]==1) 'Tree has' else 'Trees have'
  tree<-attr(contsimmap,'tree')[[1]]
  report.nodes<-if(tree[['Nnode']]==1) 'node' else 'nodes'
  cat('\n\n',report.tree,' ',Ntip(contsimmap),' tips and ',tree[['Nnode']],' internal ',report.nodes,sep='')
  if(flag) cat('\t[ ',ntip(contsimmap),' and ',nnode(contsimmap),', respectively, after accounting for subsetting ]',sep='')
  cat("\n\nTip labels:\n ",.report.names(nms=tree[['tip.label']],printlen=printlen))
  if(ntip(contsimmap)<Ntip(contsimmap)){
    tips<-tip.inds(contsimmap,unique.only=TRUE)
    tips<-tree[['tip.label']][sort(tips[!is.na(tips)])]
    cat('\t[ subsetted to',if(length(tips)) .report.names(nms=tips,printlen=printlen) else ' no tips ',']',sep='')
  }
  cat('\n\n',.report.names(c("Node label:\n ","Node labels:\n "),nms=tree[['node.label']],printlen=printlen),sep='')
  if(nnode(contsimmap)<Nnode(contsimmap)){
    nodes<-sort(unique(as.vector(node.inds(contsimmap,unique.only=TRUE))))
    nodes<-tree[['node.label']][nodes[nodes>Ntip(contsimmap)]-Ntip(contsimmap)]
    cat('\t[ subsetted to',if(length(nodes)) .report.names(nms=nodes,printlen=printlen) else ' no nodes ',']',sep='')
  }
  
  # 
  # 
  # 
  # if(!is.null(contsimmap[['trait.data']])){
  #   cat('\nIncludes trait data:\n\n')
  #   dat<-contsimmap[['trait.data']]
  #   dims<-dim(dat)
  #   report.dat<-head(contsimmap[['trait.data']],c(nrows,ncols,nslices))
  #   report.dims<-dim(report.dat)
  #   nmiss<-dims-report.dims
  #   message<-vector('character',3)
  #   for(i in seq_len(3)){
  #     if(nmiss[i]){
  #       if(i==1){
  #         message[i]<-paste(nmiss[i],
  #                           if(nmiss[i]==1) 'row' else 'rows')
  #       }else if(i==2){
  #         message[i]<-paste(nmiss[i],
  #                           if(nmiss[i]==1) 'column' else 'columns')
  #       }else{
  #         message[i]<-paste(nmiss[i],
  #                           if(nmiss[i]==1) 'matrix slice' else 'matrix slices')
  #       }
  #     }
  #   }
  #   print(report.dat,...)
  #   check<-nzchar(message)
  #   if(any(check)){
  #     message<-.report.names(nms=message[check])
  #     cat('\n[ omitted',message,']\n',sep='')
  #   }
  # }
  # Xsig2<-attr(contsimmap,'params')[['Xsig2',1]]
  # nstates<-ncol(Xsig2)
  # if(nstates>1){
  #   cat('\nIncludes discrete characters:\n  ',
  #       .report.names(nms=names(Xsig2),printlen=printlen),
  #       '\n\nWith rate matrices:\n\n',sep='')
  #   if(nstates>printlen){
  #     print(Xsig2[seq_len(printlen)],...)
  #     cat('\n...\n')
  #   }else{
  #     print(Xsig2,...)
  #   }
  # }else{
  #   cat('\nRate matrix:\n')
  #   print(unname(Xsig2),...)
  # }
}
