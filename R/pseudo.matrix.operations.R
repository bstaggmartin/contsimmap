#shaved off ~4 microseconds off these functions...whoopee?

#Invert matrices including infinite and zero entries on diagonal
#Basically, 0s become Infs and Infs become 0s along diagonals, both have no covariance with other dimensions
.pseudo.solve<-function(mat){
  mat<-as.matrix(mat)
  n<-nrow(mat)
  nn<-seq_len(n)
  tmp<-cbind(nn,nn)
  diagonal<-mat[tmp]
  infs<-is.infinite(diagonal)
  zeros<-!diagonal
  inds<-!(infs|zeros)
  sum.inds<-sum(inds)
  if(sum.inds==n){
    solve(mat)
  }else{
    out<-mat
    out[]<-0
    zeros<-tmp[zeros,,drop=FALSE]
    out[zeros]<-Inf
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
.pseudo.chol<-function(mat){
  mat<-as.matrix(mat)
  diagonal<-diag(mat)
  inds<-is.infinite(diagonal)|!diagonal
  inds<-!inds
  sum.inds<-sum(inds)
  if(sum.inds==nrow(mat)){
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