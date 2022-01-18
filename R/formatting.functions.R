.report.names<-function(prefix=NULL,nms,suffix=NULL,printlen=Inf){
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
            paste(nms,collapse=' and '),
            suffix[2])
    }else{
      paste0(prefix[2],' ',
             paste(nms[-nn],collapse=', '),', and ',nms[nn],
             ' ',suffix[2])
    }
  }
}

#Unfortunately, Ntip and Nedge don't accept additional arguments, so I instead split all these functions in Nxxx,
#which returns the FULL number of some quantity, and nxxx, which return the subsetted number of some quantity

#' @export
Ntip.contsimmap<-function(contsimmap){
  length(contsimmap[['tree']][[1]][['tip.label']])
}

#' @export
ntip<-function(contsimmap){
  ntips<-length(contsimmap[['tree']][[1]][['tip.label']])
  edge.mat<-contsimmap[['tree']][[1]][['edge']][as.numeric(names(contsimmap[['x']])),,drop=FALSE]
  sum(edge.mat[,2]<(ntips+1))
}

#' @export
Nedge.contsimmap<-function(contsimmap){
  nrow(contsimmap[['tree']][[1]][['edge']])
}

#' @export
nedge<-function(contsimmap){
  length(contsimmap[['x']])
}

#' @export
Nnode.contsimmap<-function(contsimmap,full=TRUE,internal.only=TRUE){
  if(internal.only){
    ntips<-0
  }else{
    ntips<-length(contsimmap[['tree']][[1]][['tip.label']])
  }
  if(full){
    contsimmap[['tree']][[1]][['Nnode']]+ntips
  }else{
    edge.mat<-contsimmap[['tree']][[1]][['edge']][as.numeric(names(contsimmap[['x']])),,drop=FALSE]
    length(unique(edge.mat[edge.mat>ntips]))
  }
}

#' @export
nnode<-function(contsimmap,internal.only=TRUE){
  Nnode(contsimmap,full=FALSE)
}

#' @export
ntrait<-function(contsimmap){
  length(contsimmap[['traits']])
}

#' @export
nsim<-function(contsimmap){
  nrow(contsimmap[['perm.inds']])
}

#about 100 milliseconds for 100 tips, 100 sims, 500 res, and 10 trees...might be better to pre-store this array...
#make it so that this can select specific simulations/traits/edges
#' @export
format.contsimmap<-function(contsimmap,edges=NULL,sims=NULL,traits=NULL){
  contsimmap<-subset(contsimmap,edges,sims,traits)
  ts<-contsimmap[['ts']]
  nts<-lengths(ts)
  t.nms<-format(unlist(ts,use.names=FALSE),digits=2)
  t.nms<-split(t.nms,rep(seq_along(nts),nts))
  ntraits<-ntrait(contsimmap)
  nsims<-nsim(contsimmap)
  nms<-list('time'=NULL,'trait'=contsimmap[['traits']],'sim'=NULL)
  foo<-function(ii){
    nms[[1]]<-t.nms[[ii]]
    tmp<-array(dim=c(nts[ii],ntraits,nsims),dimnames=nms)
  }
  out<-setNames(lapply(seq_along(nts),foo),names(ts))
  ord<-.get.cladewise.ord(contsimmap[['anc']])
  tree.seq<-seq_along(contsimmap[['tree']])
  #instead use  pre-stored anc vector and nodes list
  #will need to re-sort seq inds
  for(i in ord){
    x<-contsimmap[['x']][[i]]
    maps<-contsimmap[['maps']][[i]]
    incl<-lapply(maps,'[[','incl')
    x<-lapply(tree.seq,function(ii) x[[ii]][incl[[ii]],,,drop=FALSE])
    out[[i]][-1,,]<-unlist(x)
    out[[i]]<-out[[i]][,,contsimmap[['perm.inds']][,3],drop=FALSE]
    #add initial x
    #(is already permuted correctly)
    anc.ind<-contsimmap[['anc']][i]
    anc<-out[[anc.ind]]
    if(is.null(anc)){
      out[[i]][1,,]<-contsimmap[['nodes']][[anc.ind]]
    }else{
      out[[i]][1,,]<-anc[nts[anc.ind],,]
    }
  }
  out
}

#"bare-bones" print and summary functions
#might be good to have some functions for better extracting specific info from objects
#(specific traits, simulations, time points, etc.)
#could call these in plotting method down the road
#' @export
print.contsimmap<-function(contsimmap,printlen=6,...){
  ntraits<-ntrait(contsimmap)
  nsims<-nsim(contsimmap)
  tree<-if(nsims==1) 'tree' else 'trees'
  cat(nsims,'phylogenetic',tree,'with',ntraits,'mapped continuous',
      .report.names(c('character:','characters:'),contsimmap[['traits']],printlen=printlen))
  if(nedge(contsimmap)<Nedge(contsimmap)){
    edges<-names(contsimmap[['ts']])
    cat('(subsetted to ',.report.names(c('edge','edges'),nms=edges,printlen=printlen),')',sep='')
  }
}

#' @export
summary.contsimmap<-function(contsimmap,printlen=6,nrows=printlen,ncols=printlen,nslices=printlen,...){
  print(contsimmap)
  tree<-contsimmap[['tree']]
  report.tree<-if(length(tree)==1) 'Tree has' else 'Trees have'
  tree<-tree[[1]]
  ntips<-length(tree$tip.label)
  nnodes<-tree$Nnode
  report.nodes<-if(nnodes==1) 'node' else 'nodes'
  cat('\n\n',report.tree,' ',ntips,' tips and ',nnodes,' internal ',report.nodes,'\n',sep='')
  cat("\nTip labels:\n  ",.report.names(nms=tree$tip.label,printlen=printlen),'\n')
  if(!is.null(tree$nodel.label)){
    cat('\n',.report.names(c("Node label:\n  ","Node labels:\n  "),nms=tree$node.label,printlen=printlen),'\n')
  }
  if(!is.null(contsimmap[['trait.data']])){
    cat('\nIncludes trait data:\n\n')
    dat<-contsimmap[['trait.data']]
    dims<-dim(dat)
    report.dat<-head(contsimmap[['trait.data']],c(nrows,ncols,nslices))
    report.dims<-dim(report.dat)
    nmiss<-dims-report.dims
    message<-vector('character',3)
    for(i in seq_len(3)){
      if(nmiss[i]){
        if(i==1){
          message[i]<-paste(nmiss[i],
                            if(nmiss[i]==1) 'row' else 'rows')
        }else if(i==2){
          message[i]<-paste(nmiss[i],
                            if(nmiss[i]==1) 'column' else 'columns')
        }else{
          message[i]<-paste(nmiss[i],
                            if(nmiss[i]==1) 'matrix slice' else 'matrix slices')
        }
      }
    }
    print(report.dat,...)
    check<-nzchar(message)
    if(any(check)){
      message<-.report.names(nms=message[check])
      cat('\n[ omitted',message,']\n',sep='')
    }
  }
  Xsig2<-contsimmap[['params']][['Xsig2']]
  nstates<-length(Xsig2)
  if(nstates>1){
    cat('\nIncludes discrete characters:\n  ',
        .report.names(nms=names(Xsig2),printlen=printlen),
        '\n\nWith rate matrices:\n\n',sep='')
    if(nstates>printlen){
      print(Xsig2[seq_len(printlen)],...)
      cat('\n...\n')
    }else{
      print(Xsig2,...)
    }
  }else{
    cat('\nRate matrix:\n')
    print(unname(Xsig2),...)
  }
}
