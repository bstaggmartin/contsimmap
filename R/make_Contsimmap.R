#2/13/23 Parameter formatting:
# - Generally expects list arrays with each column corresponding to a state (mu, Xsig2), tip (Ysig2, nobs), or trait (X0)
#    - Also accepts matrices in the case of X0/nobs
#    - Also accepts lists of lists, where each sublist represents a column
# - Unnamed elements readily recycled to fill missing portions of inputs, but if no unnamed elements are found, will be
#   filled in with defaults (identity for Xsig2, 0s for mu/Ysig2/nobs/X0)
#    - Entries of nobs corresponding to internal nodes are never derived from recycling and always default to 0 if unspecified

#' Simulate a continuous stochastic character map
#' 
#' This function creates a continuous stochastic character map (class "\code{contsimmap}") by simulating and
#' mapping evolving continuous characters on a given phylogeny.
#' 
#' @param tree A phylogeny or list of phylogenies with or without mapped discrete characters (classes 
#' "\code{phylo}", "\code{multiPhylo}", "\code{simmap}", or "\code{multiSimmap}"). Lists of phylogenies
#' with differing discrete character histories are allowed, but lists with \emph{differing topologies
#' are currently not supported}! I plan to implement "\code{multiContsimmap}"-type objects for this
#' purpose in the future.
#' @param ntraits The number traits to simulate.
#' @param traits A character vector specifying the names of each trait. The \code{<i>}th unnamed trait 
#' defaults to "\code{trait_<i>}".
#' @param nobs A numeric vector specifying the number of observations to simulate per tip/node. Entries
#' are assigned to tips/nodes if named according to \code{tree$tip.label}/\code{tree$node.label} (node labels default to 
#' their numeric index if not provided in \code{tree$node.label}). Unnamed/unassigned entries are recycled
#' to "fill in" values for any tips lacking input (\emph{but not internal nodes}) in order of increasing node index.
#' If no unnamed/unassigned entries are available, tips and internal nodes default to 1 and 0 observations, respectively.
#' To specify differing numbers of observations for each phylogeny/simulation, format \code{nobs} instead as a list of vectors
#' or matrix with each column/column name corresponding to different nodes (see section for more details recycling behavior
#' across phylogenies/simulations).
#' @param nsims The number of simulations to perform. Each simulation will correspond to
#' a separate phylogeny in \code{tree}, which is recycled as needed. For example, providing a
#' "\code{multiSimmap}" object of length 3 for \code{tree} and specifying \code{nsims=7} will result
#' in 7 simulations on phylogenies 1, 2, 3, 1, 2, 3, and 1, in that order.
#' @param res Controls the approximate number of timepoints at which to sample trait values across the
#' entire height of the phylogeny (i.e., from root to last-surviving tip). Higher values result
#' in more densely-sampled character histories but take longer to simulate and use more computer memory.
#' @param X0 A numeric vector specifying the starting trait values at the root of the phylogeny.
#' Entries are assigned to traits if named according to \code{traits}. Unnamed/unassigned entries are
#' recycled to "fill in" values for any traits lacking input in the same order as given in \code{traits}.
#' If no unnamed/unassigned entries are available, defaults to 0. To specify differing starting
#' trait values for each phylogeny/simulation, format \code{X0} instead as a list of vectors or matrix
#' with each column/column name corresponding to different traits (see section for more details recycling behavior
#' across phylogenies/simulations).
#' @param Xsig2 A list of numeric matrices/vectors specifying the evolutionary rate matrices for each discrete
#' character state mapped onto the phylogeny. Evolutionary rate matrices are symmetric matrices with
#' diagonal and off-diagonal entries corresponding to evolutionary rates and covariances, respectively. Vectors are
#' assumed to be diagonal matrices (i.e., no evolutionary covariance among traits). Matrix/vector entries are assigned
#' to pairs of traits based on associated row/column names (plain names in the case of vectors), which are matched to
#' \code{trait}. When possible, unspecified matrix entries (\code{NA} entries) and row/column names are automatically set to make
#' inputted matrices symmetric. Any remaining unnamed/unassigned rows/columns are recycled as a block-diagonal matrix to
#' "fill in" values for traits completely lacking input in the same order as given in \code{traits}. All remaining
#' unspecified rates and covariances default to 1 and 0, respectively. Specifying an invalid variance-covariance matrices
#' (either due to asymmetry or not being positive semidefinite) result in an error. List entries are assigned to discrete
#' character states if named according to the state names given in \code{tree}. Any unnamed/unassigned list entries are recycled to
#' "fill in" values for states lacking input in alphabetical order. If no unnamed/unassigned list entries are available,
#' defaults to identity matrix (i.e., all rates of 1 with no covariance). To specify differing evolutionary rate matrices
#' for each phylogeny/simulation, format \code{Xsig2} instead as a list of lists or matrix-shaped list with each column/column
#' name corresponding to different discrete character states (see section for more details recycling behavior
#' across phylogenies/simulations).
#' @param Ysig2 A list of numeric matrices/vectors specifying the intraspecifc and/or measurement error for each tip/node in
#' the phylogeny. These are symmetric matrices with diagonal and off-diagonal entries corresponding to the variances and covariances,
#' respectively, of trait measurements for a particular tip/node. Vectors are assumed to be diagonal matrices (i.e., no covariance among
#' trait measurements within a given node). Matrix/vector entries are assigned to pairs of traits based on associated row/column names
#' (plain names in the case of vectors), which are matched to \code{trait}. When possible, unspecified matrix entries (\code{NA} entries)
#' and row/column names are automatically set to make inputted matrices symmetric. Any remaining unnamed/unassigned rows/columns are recycled
#' as a block-diagonal matrix to "fill in" values for traits completely lacking input in the same order as given in \code{traits}. All remaining
#' unspecified matrix entries default to 0. Specifying an invalid variance-covariance matrices
#' (either due to asymmetry or not being positive semidefinite) result in an error. List entries are assigned to tips/nodes if named
#' according to \code{tree$tip.label}/\code{tree$node.label} (node labels default to 
#' their numeric index if not provided in \code{tree$node.label}). Any unnamed/unassigned list entries are recycled to
#' "fill in" values for nodes/tips lacking input in order of increasing node index. If no unnamed/unassigned list entries are available,
#' defaults to matrix of 0s. To specify differing intraspecific/measurement errors
#' for each phylogeny/simulation, format \code{Xsig2} instead as a list of lists or matrix-shaped list with each column/column
#' name corresponding to different tips/nodes (see section for more details recycling behavior
#' across phylogenies/simulations).
#' @param mu
#' @param verbose
#' 
#' @section Recycling Behavior:
#' Each parameter input system (\code{nobs}, \code{X0}, \code{Xsig2}, \code{Ysig2}, and \code{mu})
#' has certain idiosyncrasies, but they all follow a similar philosophy. There are three steps:
#' \enumerate{
#'  \item{Match labeled inputs to appropriate nodes/traits/states.}
#'  \item{"Fill in" for nodes/traits/states lacking inputs with unlabeled inputs if they exist 
#'    (recycling the unlabeled inputs as needed) and default values otherwise.}
#'  \item{If multiple parameter values are provided, different parameter values will be applied to
#'    each phylogeny in \code{tree}, recycling as necessary. This creates some interesting quirks 
#'    when the  length of different parameter values exceeds the number of phylogenies. For example,
#'    if \code{tree} contains 3 phylogenies and 4 lists of matrices are specified for \code{Xsig2}, then
#'    simulations on the 2nd and 3rd phylogenies will use the 2nd and 3rd \code{Xsig2} lists, respectively,
#'    while simulations on the 1st phylogeny will alternate between the 1st and 4th \code{Xsig2} lists.}
#' }
#' I tried to make warning messages informative enough to make potentially unexpected recycling behaviors
#' apparent. Users specifying complicated simulations involving varying parameter values can always
#' double-check their inputs were interpreted correctly using \code{get.param.info()}.
#' 
#' @export
sim.conthistory<-function(tree,ntraits=1,traits=paste0('trait_',seq_len(ntraits)),nobs=NULL,
                          nsims=100,res=100,
                          X0=NULL,Xsig2=NULL,Ysig2=NULL,mu=NULL,
                          verbose=FALSE){
  #takes care of all the necessary formatting work...
  list2env(.prep.contsimmap(tree,nsims,res,X0,Xsig2,Ysig2,mu,conditional=FALSE,ntraits,traits,nobs),envir=environment())
  
  #at some point want to generalize to have multiple nobs per call to simulate differently-structured trait.data...
  #a little confusing at this point since the trait.data lookup element refers to X0 here, but trait.data for make.contsimmap
  #but that can be resolved manually for now
  #^2/14/23: need to look up if the above still holds
  #Can now handle multiple nobs, but things are still a bit confusing to me with how lookups are handled
  #Should eventually become more concrete/organized as code matures
  tmp<-.uncond.traversals(prune.seq,anc,tree[[1]][['edge']],Ntip(tree[[1]]),
                          maps,
                          X0,nobs,Xsig2,Ysig2[,Yperm,drop=FALSE],mu,lookup,
                          nts,seed,
                          verbose=verbose)
  
  #format output
  x<-tmp[[1]]
  trait.data<-tmp[[2]]
  attr(x,'ts')<-ts
  attr(x,'tree')<-tree
  attr(x,'maps')<-maps
  attr(x,'treeID')<-treeID
  traits<-colnames(Xsig2[[1]])
  attr(x,'traits')<-setNames(rep(1,length(traits)),traits)
  #need to first reconfigure lookup to refer to include newly simulated trait.data...
  counter<-1
  for(i in seq_along(lookup)){
    tmp.n<-length(lookup[[i]][['matches']])
    lookup[[i]][['table']]<-lookup[[i]][['table']][lookup[[i]][['matches']],,drop=FALSE]
    lookup[[i]][['matches']]<-seq_len(tmp.n)
    lookup[[i]][['table']][,1]<-counter:(counter+tmp.n-1)
    counter<-counter+tmp.n
  }
  attr(x,'params')<-matrix(list(trait.data,Xsig2,Ysig2,mu,lookup,list(fxn="sim.conthistory")),
                           6,1,
                           dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  .uncompress(x)
}

#may want to make some function for auto-estimation of Xsig2/Ysig2
#' @export
make.contsimmap<-function(tree,trait.data,
                          nsims=100,res=100,
                          Xsig2=NULL,Ysig2=NULL,mu=NULL,
                          verbose=FALSE){
  #takes care of all the necessary formatting work...
  list2env(.prep.contsimmap(tree,nsims,res,trait.data,Xsig2,Ysig2,mu,TRUE),envir=environment())
  
  #does the post/preorder traversals, returns final simulated trait data
  x<-.cond.traversals(prune.seq,anc,des,ndes,
                      maps,
                      parsed.obs,parsed.mis,nobs,Xsig2,Ysig2[,Yperm,drop=FALSE],mu,lookup,
                      nts,NTS,t1s,seed,x,v,dx,dv,
                      verbose=verbose)
  
  #format output
  attr(x,'ts')<-ts
  attr(x,'tree')<-tree
  attr(x,'maps')<-maps
  attr(x,'treeID')<-treeID
  traits<-colnames(Xsig2[[1]])
  attr(x,'traits')<-setNames(rep(1,length(traits)),traits)
  attr(x,'params')<-matrix(list(trait.data,Xsig2,Ysig2,mu,lookup,list(fxn="make.contsimmap")),
                           6,1,
                           dimnames=list(c('trait.data','Xsig2','Ysig2','mu','lookup','call_info'),NULL))
  .uncompress(x)
}