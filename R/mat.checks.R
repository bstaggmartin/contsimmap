#' recycling explanation: Generally speaking, named entries of input will be matched up to specified trait/state names, then any entries with either no
#' or unmatched names are recycled in their given order to form an output of appropriate dimensions. When there are unmatched or unnamed entries, missing
#' entries are filled in with default values: 1 for variances and 0 for covariances (i.e., the identity matrix).
#' 
#' In the case of matrices specifically, the function first attempts to "symmetrize" the input. For vector inputs, this creates a diagonal matrix. For
#' matrix inputs, rows/columns with missing names are labelled either with their corresponding column/row names (if they exist) or a numerical index.
#' Extra rows/columns of NAs are tacked on to the matrix so each row has a corresponding column and vice versa, the matrix is rearranged such that its row
#' and column names are identical, and all NA entries i,j are replaced with non-NA entries j,i. From here, the output is formed from the symmetrized input
#' as described above. Unmatched rows/columns of the input matrix are recycled into a block-diagonal matrix, and any unspecified covariances remaining are
#' set to 0.

.fix.trait.names<-function(k,nms){
  nms<-c(nms,rep(NA,k-length(nms)))
  nas<-which(is.na(nms)|!nzchar(nms))
  nms[nas]<-paste0('trait_',nas)
  nms
}

.symmetrize.mat<-function(mat,nice.name,Ysig2.flag){
  ndims<-length(dim(mat))
  if(ndims>2){
    stop(nice.name,' has more than two dimensions; multiple covariance matrices should be stored in lists, not arrays')
  }else if(ndims<2){
    if(length(mat)>1){
      nms<-names(mat)
      mat<-diag(mat)
      rownames(mat)<-colnames(mat)<-nms
    }else{
      mat<-matrix(mat,1,1,dimnames=rep(list(names(mat))))
    }
  }
  dims<-dim(mat)
  #get names
  nms<-dimnames(mat)
  if(!length(nms)) nms<-rep(list(NULL),2)
  for(i in c(1,2)){
    if(is.null(nms[[i]])) nms[[i]]<-rep('',dims[i])
  }
  minseq<-seq_len(min(dims))
  offset<-0
  for(i in c(1,2)){
    inds<-!nzchar(nms[[i]])
    other<-nms[[c(2,1)[i]]]
    other.inds<-nzchar(other)
    #fill in with names from other dimension, if they exist
    nms[[i]][minseq][inds[minseq]&other.inds[minseq]]<-other[other.inds[minseq]]
    inds<-!nzchar(nms[[i]])
    nms[[i]][inds]<-paste0('TEMPORARY_',seq_len(sum(inds))+offset)
    offset<-sum(inds)
  }
  dimnames(mat)<-nms
  #append rows corresponding to cols
  missing.rows<-!(nms[[2]]%in%nms[[1]])
  mat<-do.call(rbind,c(list(mat),rep(NA,sum(missing.rows))))
  rownames(mat)<-c(nms[[1]],nms[[2]][missing.rows])
  nms<-dimnames(mat)
  #append cols corresponding to rows
  missing.cols<-!(nms[[1]]%in%nms[[2]])
  mat<-do.call(cbind,c(list(mat),rep(NA,sum(missing.cols))))
  colnames(mat)<-c(nms[[2]],nms[[1]][missing.cols])
  #sort
  nms<-sort(nms[[1]])
  mat<-mat[nms,nms,drop=FALSE]
  #reflect
  lo<-lower.tri(mat)
  up<-upper.tri(mat)
  tmat<-t(mat)
  inds<-is.na(mat)
  inds<-inds&t(!inds)
  mat[lo&inds]<-tmat[lo&inds]
  mat[up&inds]<-tmat[up&inds]
  #fill in the rest with defaults
  diagonal<-diag(mat)
  probs<-is.na(diagonal)
  if(any(probs)){
    warning('PLACEHOLDER: diagonal fill-ins')
  }
  diagonal[probs]<-if(Ysig2.flag) 0 else 1
  diag(mat)<-diagonal
  probs<-is.na(mat)
  if(any(probs)){
    warning('PLACEHOLDER: off-diagonal fill-ins')
  }
  mat[probs]<-0
  mat
}

.fix.mat<-function(mat,nice.name,k,traits,Ysig2.flag){
  mat<-.symmetrize.mat(mat,nice.name,Ysig2.flag)
  #form output
  out<-matrix(NA,k,k,dimnames=rep(list(traits),2))
  nms<-dimnames(mat)[[1]]
  prob.nms<-grepl('^TEMPORARY_\\d+$',nms)
  if(all(prob.nms)){
    has.prob.nms<-TRUE
    quiet.flag<-TRUE
  }else{
    prob.nms<-prob.nms|!(nms%in%traits)
    has.prob.nms<-any(prob.nms)
    quiet.flag<-FALSE
  }
  if(any(!prob.nms)){
    out[nms[!prob.nms],nms[!prob.nms]]<-mat[!prob.nms,!prob.nms]
  }
  diagonal<-diag(out)
  prob.traits<-is.na(diagonal)
  if(any(prob.traits)){
    if(has.prob.nms){
      #take elements with no matching names and fill in
      nrem<-sum(prob.traits)
      tmp<-kronecker(diag(ceiling(nrem/sum(prob.nms))),mat[prob.nms,prob.nms])
      n<-nrow(tmp)
      incl.inds<-seq_len(n)<=nrem
      tmp<-tmp[incl.inds,incl.inds]
      out[prob.traits,prob.traits]<-tmp
      message<-paste0("; missing entries of ",nice.name," recycled from other entries of ",nice.name)
    }else{
      #insert diagonal of 1s
      diagonal[prob.traits]<-if(Ysig2.flag) 0 else 1
      diag(out)<-diagonal
      message<-paste0("; missing entries of ",nice.name," filled in with",
                      if(Ysig2.flag) "0s" else "entries of identity matrix")
    }
    out[is.na(out)]<-0
    if(!quiet.flag){
      warning(.report.names(c('Trait','Traits'),traits[prob.traits],c('has','have')),
              ' no associated entries in ',
              nice.name,message,
              immediate.=TRUE)
    }
  }
  if(!isSymmetric(mat)|any(eigen(mat)$values<0)){
    stop('Failed to create proper covariance matrix from ',nice.name)
  }
  out
}

#may want a warning for empty state/trait names and trait names of pattern 'TEMPORARY_\\d+'
.fix.mat.list<-function(l,ntraits,traits,nstates,states){
  nice.name<-deparse(substitute(l))
  Ysig2<-nice.name=='Ysig2'
  traits<-.fix.trait.names(ntraits,traits)
  if(!is.list(l)){
    l<-list(l)
  }
  nms<-names(l)
  n<-length(l)
  if(!length(nms)){
    names(l)<-rep('',n)
    nms<-names(l)
    prob.nms<-rep(TRUE,n)
    has.prob.nms<-TRUE
    quiet.flag<-TRUE
  }else{
    prob.nms<-!(nms%in%states)
    has.prob.nms<-any(prob.nms)
    quiet.flag<-FALSE
  }
  for(i in seq_len(n)){
    nm<-nms[i]
    if(nzchar(nm)){
      report.i<-paste0("'",nm,"'")
    }else{
      report.i<-i
    }
    l[[i]]<-.fix.mat(l[[i]],paste0(nice.name,'[[',report.i,']]'),
                     ntraits,traits,
                     Ysig2)
  }
  out<-setNames(vector('list',nstates),states)
  if(any(!prob.nms)){
    out[nms[!prob.nms]]<-l[!prob.nms]
  }
  prob.states<-!lengths(out)
  if(any(prob.states)){
    if(has.prob.nms){
      #take Xsig2/Ysig2 elements with no matching names and fill in
      tmp<-l[prob.nms]
      message<-paste0("; missing ",nice.name,"'s recycled from elements of ",nice.name," list")
    }else{
      #take default Xsig2/Ysig2 element and fill in
      dimnms<-list(traits,traits)
      tmp<-list(diag(ntraits),dimnames=dimnms)
      if(Ysig2) tmp<-0*tmp
      message<-paste0("; missing ",nice.name,"'s filled in with",
                      if(Ysig2) "0s" else "identity matrices")
    }
    out[prob.states]<-rep(tmp,length.out=sum(prob.states))
    if(!quiet.flag){
      prefix<-if(Ysig2) c('Tip','Tips') else c('State','States')
      warning(.report.names(prefix,states[prob.states],c('has','have')),
              ' no associated elements in ',
              nice.name,message,
              immediate.=TRUE)
    }
  }
  out
}

.fix.nobs<-function(nobs,ntips,tips){
  nms<-names(nobs)
  n<-length(nobs)
  if(!length(nms)){
    names(nobs)<-rep('',n)
    nms<-names(nobs)
    prob.nms<-rep(TRUE,n)
    has.prob.nms<-TRUE
    quiet.flag<-TRUE
  }else{
    prob.nms<-!(nms%in%tips)
    has.prob.nms<-any(prob.nms)
    quiet.flag<-FALSE
  }
  out<-setNames(rep(NA,ntips),tips)
  if(any(!prob.nms)){
    out[nms[!prob.nms]]<-nobs[!prob.nms]
  }
  prob.tips<-is.na(out)
  if(any(prob.tips)){
    if(has.prob.nms){
      tmp<-nobs[prob.nms]
      message<-paste0("; missing entries of nobs recycled from other entries of nobs")
    }else{
      tmp<-1
      message<-paste0("; missing entries of nobs filled with 1s")
    }
    out[prob.tips]<-rep(tmp,length.out=sum(prob.tips))
    if(!quiet.flag){
      warning(.report.names(c('Tip','Tips'),tips[prob.tips],c('has','have')),
              ' no associated entries in nobs',message,
              immediate.=TRUE)
    }
  }
  out
}