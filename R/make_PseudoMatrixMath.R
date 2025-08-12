#shaved off ~4 microseconds off these functions...whoopee?

#Invert matrices including infinite and zero entries on diagonal
#Basically, 0s become Infs and Infs become 0s along diagonals, both have no covariance with other dimensions
.pseudo.solve<-function(mat,k,diag.inds,z2z=FALSE){
  diagonal<-mat[diag.inds]
  infs<-is.infinite(diagonal)
  zeros<-!diagonal
  inds<-!(infs|zeros)
  sum.inds<-sum(inds)
  if(sum.inds==k){
    solve(mat)
  }else{
    out<-mat
    out[]<-0
    if(!z2z){
      zeros<-diag.inds[zeros,,drop=FALSE]
      out[zeros]<-Inf
    }
    if(sum.inds){
      out[inds,inds]<-solve(mat[inds,inds])
    }
    out
  }
}

#Take cholesky composition of matrix including zero entries on diagonal (and infinites, I guess?)
#Returns expected result: t(out)%*%out returns original matrix (EXCEPT with infinites, but this shouldn't be
##issue in practice since they should never occur during Simmap construction)
#Keeping infinites cause proper propagation of Infs along diagonal, but also propagates NaNs during matrix
##multiplication, so should be avoided here if it can be helped
#The nice thing with the current construction is that both 0s and Infs cause diagonal entries to be 0 with 0
##correlation with other variables--this means that, in the context of multiplying uncorrelated variables to be
##correlated, dimensions with 0 and infinite variance will stay the same and not affect other dimensions; this
##also works in the case of infinite variance since we want these means to stay at 0
#Note that there are multiple decompostions for a matrices with 0 on diagonal in general (which is why normal
##Cholesky decompostion doesn't work); I avoid this problem here by making sure dimensions with 0 on the diagonal
##are filled up with 0s, ensuring a unique solution (well, unique according to what I read on Wikipedia)
.pseudo.chol<-function(mat,k,diag.inds){
  diagonal<-mat[diag.inds]
  inds<-is.infinite(diagonal)|!diagonal
  inds<-!inds
  sum.inds<-sum(inds)
  if(sum.inds==k){
    chol(mat)
  }else{
    out<-mat
    out[]<-0
    if(sum.inds){
      out[inds,inds]<-chol(mat[inds,inds])
    }
    out
  }
}

.pseudo.chol.solve<-function(mat,k,diag.inds){
  diagonal<-mat[diag.inds]
  inds<-is.infinite(diagonal)|!diagonal
  inds<-!inds
  sum.inds<-sum(inds)
  if(sum.inds==k){
    chol(solve(mat))
  }else{
    out<-mat
    out[]<-0
    if(sum.inds){
      out[inds,inds]<-chol(solve(mat[inds,inds]))
    }
    out
  }
}

.solve<-function(arr,dim,k,diag.inds,z2z=FALSE){
  for(i in seq_len(dim)){
    arr[,,i]<-.pseudo.solve(matrix(arr[,,i,drop=FALSE],k,k),k,diag.inds,z2z)
  }
  arr
}

.chol<-function(arr,dim,k,diag.inds){
  for(i in seq_len(dim)){
    arr[,,i]<-.pseudo.chol(matrix(arr[,,i,drop=FALSE],k,k),k,diag.inds)
  }
  arr
}

.chol.solve<-function(arr,dim,k,diag.inds){
  for(i in seq_len(dim)){
    arr[,,i]<-.pseudo.chol.solve(matrix(arr[,,i,drop=FALSE],k,k),k,diag.inds)
  }
  arr
}

.multAb<-function(arr1,arr2,dim,k){
  if(is.null(dim(arr2))){
    do.call(cbind,lapply(seq_len(dim),function(ii) matrix(arr1[,,ii,drop=FALSE],k,k)%*%arr2))
  }else{
    for(i in seq_len(dim)){
      arr2[,i]<-matrix(arr1[,,i,drop=FALSE],k,k)%*%arr2[,i,drop=FALSE]
    }
    arr2
  }
}

.multbA<-function(arr1,arr2,dim,k){
  for(i in seq_len(dim)){
    arr1[i,]<-matrix(arr1[i,,drop=FALSE],1,k)%*%matrix(arr2[,,i,drop=FALSE],k,k)
  }
  arr1
}

.sum3d<-function(arr,dim){
  if(dim>1){
    for(i in seq_len(dim-1)){
      arr[,,1]<-arr[,,1,drop=FALSE]+arr[,,i+1,drop=FALSE]
    }
    arr<-arr[,,1,drop=FALSE]
  }
  arr
}

.resolve.infs.ls<-function(arr.ls,dim1,dim2,k,diag.inds,precedence=FALSE,inf.const=1e10){
  # if(precedence){
  #   ones<-is.infinite(arr.ls[[1]][diag.inds])
  #   if(any(ones)){
  #     zeros<-array(rep(ones,each=k),c(k,k,dim1))
  #     zeros<-zeros|aperm(zeros,c(2,1,3))
  #     for(i in seq_len(dim2)){
  #       arr.ls[[i]][zeros]<-0
  #     }
  #     arr.ls[[1]][diag.inds][ones]<-1
  #   }
  # }
  # ones<-lapply(arr.ls,function(ii) is.infinite(ii[diag.inds]))
  # zeros<-array(rep(Reduce('|',ones),each=k),c(k,k,dim1))
  # if(any(zeros)){
  #   zeros<-zeros|aperm(zeros,c(2,1,3))
  #   for(i in seq_len(dim2)){
  #     arr.ls[[i]][zeros]<-0
  #     arr.ls[[i]][diag.inds][ones[[i]]]<-1
  #   }
  # }
  for(i in seq_len(dim2)){
    arr.ls[[i]][is.infinite(arr.ls[[i]])]<-inf.const
  }
  nulls<-unlist(lapply(arr.ls,function(ii) is.null(dim(ii))),use.names=FALSE)
  arr.ls[nulls]<-lapply(arr.ls[nulls],function(ii) array(ii,c(k,k,dim1)))
  arr.ls
}

.resolve.infs<-function(arr,dim,k,diag.inds,inf.const=1e10){
  # ones<-matrix(is.infinite(arr[diag.inds]),k,dim)
  # zeros<-apply(ones,1,any)
  # zeros<-zeros|rep(zeros,each=k)
  # arr[zeros]<-0
  # arr[diag.inds][ones]<-1
  arr[is.infinite(arr)]<-inf.const
  arr
}
