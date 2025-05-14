#4/4 updates:
#function working for stuff with mapped traits--yay!
#univariate version works and is quite fast!
##  still need to test univariate version with mapped traits-->just did, broken :( But probs very close!
#as far as I can tell, everything works now!
#and pretty speedy in the univariate case so far as I can discern...
#feels a bit sluggish, as you were worried about
##  may not be a big issue for now
##  some potential shortcuts to explore:
##    include non-mapped versions? Maybe... (low-med priority)
##    include single ID version? I don't think this would save much time, but perhaps worth considering... (lower priority)
##    port likelihood function to C (ugh, I hope not --> lowest priority)
#4/12 update:
#function now works with positive semi-definite matrices in the multivariate case! Yay!
#so you can set Yvar elements to 0!
#still need to work on univariate case
#ALSO, note that it does the contsimmap approximation, where things are just averaged with 0 variance if you have multiple conflicting observations
#also makes it work with 0-length branches!
#One problem I realized here is what I'm terming "dimension-dropping"--in some cases, the likelihood might be arbitrarily increased by setting variance elements to 0 such that measurements are lumped together
#As the number of measurements decreases, the multivariate normal distribution will lose dimensions and the likelihood density will increase
#For safety I think I will just make it throw an error whenever multiple observations are averaged--this can be done by checking if the corresponding diagonal element of P!=1
#Might add an option to turn this on or off, but I think it's best to keep this on for now
#4/17 update:
#function now works with 0 variance in the univariate case too! Woo!
#4/19 update:
#added some more checks for variance components
# less than 0, NA, or inf values now result in -Inf likelihoods
# Also now round values less than a new argument, vartol, to 0--I think this is mainly what was causing "runaway" likelihoods...
# You could get values as low as ~5e-324 (e^-746ish) otherwise, at least on your machine, which was probably really fucking with things
#9/7 update:
#wow! Been a while, lots of changes
#Yeah, implementing vartol seemed to fix the runaway likelihoods--I don't see them anymore
#There was a math error with how likelihoods were weighted under the nuisance prior (downweighted because I added nuisance weighting ON TOP OF averaging)
# now fixed!
#Added option to do nuisance priors over root as well
#Working on implementing a way to subset which trees to calculate likelihoods for--annoying, but it will allow for (relatively) fast gradient calculations
# (because most of the likelihood at a point is conditional on just a few contsimmaps)
# I'm currently at the point where the function only calculates necessary mapped xvars/yvars/mus and returns them in smaller arrays
# still have to work out mechanisms to update other indexing to reflect the rearrangement of tree indices implied by this procedure...



#' Make a likelihood function with parameters dependent on continuous stochastic
#' character mapped variables
#' 
#' This function uses continuous stochastic character maps (class "\code{contsimmap}")
#' to generate a function that calculates the likelihood of a Brownian Motion model with
#' parameters (e.g., rates, trends) that depend on mapped variables.
#' 
#' 
#' 
#' @param tree A phylogeny or list of phylogenies with or without mapped 
#' continuous/discrete characters (classes "\code{contsimmap}", "\code{phylo}",
#' "\code{multiPhylo}", "\code{simmap}", or "\code{multiSimmap}"). Lists of
#' phylogenies with differing discrete character histories are allowed, but
#' lists with \emph{differing topologies are currently not supported}! I plan to 
#' implement "\code{multiContsimmap}"-type objects for this purpose in the 
#' future.
#' 
#' @param trait.data A numeric matrix/vector or data frame specifying
#' observed trait measurements to be modeled. Generally,
#' rows should correspond to different observations while columns correspond 
#' to different traits (vectors can only supply data for a single trait). Any 
#' tip/node in the phylogeny may be associated with 0 or more observations; use 
#' \code{NA} to specify missing trait measurements in the case of partial 
#' observations (i.e., only a subset of traits were measured in a given 
#' observation). 
#' 
#' To assign observations to tips/nodes, the data \emph{must} be
#' labelled in some way according to 
#' \code{tree$tip.label}/\code{tree$node.label} (node labels default to their 
#' numeric index if not provided in \code{tree$node.label}). These labels must
#' be provided as names for vectors, rownames for matrices, and 
#' a column of factor/character data for data frames. In the case of data 
#' frames, if multiple columns with factor/character data are found, the column 
#' with the most matches to \code{tree$tip.label}/\code{tree$node.label} is 
#' automatically chosen (any rownames are converted to a column before this 
#' step). Trait names are simply taken from the column names of matrices/data 
#' frames, and default to "\code{trait_<i>}" for the \code{<i>}th column of 
#' trait data if not provided. 
#' 
#' To specify different trait data for different phylogenies included in \code{tree}, format \code{trait.data} instead as a list of 
#' vectors/matrices/data frames (see Recycling Behavior section for further 
#' explanation of parameter/formula recycling across phylogenies/individual continuous
#' stochastic character maps). 
#' 
#' @param Xvar
#' @param Xcor
#' @param mu
#' @param Yvar
#' @param Ycor
#' @param nsim
#' @param root.nuisance.prior
#' @param wgts,wgt.by.nobs
#' @param tree.prior.exp
#' @param tree.balance.prior
#' 
#' @param vartol Should probably leave this alone unless you really know what 
#' you're doing! If an estimated variance parameter (i.e., evolutionary rates and
#' tip errors) falls below \code{vartol}, it is rounded down to 0. I picked the
#' square root of machine epsilon as a default through personal experimentation.
#' 
#' In my experience with the outputs of this function,
#' dividing by or inverting matrices with extremely small variances can result in 
#' accumulating numerical errors during likelihood calculations, such
#' that you generally get more accurate results just rounding small variances to
#' 0 (the functions use "pseudoinversion" to be robust to this; see Hassler et al.
#' 2020). Unlike the algorithm used to generate continuous stochastic character maps,
#' outputted likelihood function explicitly check for contradictory trait data (i.e.,
#' two different observations coming from a distribution inferred to have 0 variance).
#' In this case, the outputted likelihood function returns a likelihood of 0 (or
#' log-likelihood of \code{-Inf}).
#' 
#' @param ... Currenty unused; kept for backwards compatibility with earlier
#' versions of function.
#'
#' @export
make.lik.fun<-function(tree,trait.data,
                       Xvar=NULL,Xcor=NULL,
                       mu=NULL,
                       Yvar=NULL,Ycor=NULL,
                       nsim=NULL,
                       root.nuisance.prior=FALSE,
                       wgts=1,wgt.by.nobs=TRUE,
                       outlier.cutoff=0,winsorize=TRUE,
                       tree.prior.exp=1,tree.balance.prior=FALSE,
                       vartol=sqrt(.Machine$double.eps),
                       ...){
  
  ####INITIAL INPUT PROCESSING####
  
  if(hasArg(nuisance.prior)){
    if(list(...)[["nuisance.prior"]]){
      tree.prior.exp<-2
    }else{
      tree.prior.exp<-1
    }
  }
  if(hasArg(tree.nuisance.prior)){
    if(list(...)[["tree.nuisance.prior"]]){
      tree.prior.exp<-2
    }else{
      tree.prior.exp<-1
    }
  }
  
  if(outlier.cutoff>0){
    if(outlier.cutoff>1){
      stop("outlier cutoff cannot exceed 100% (and probably shouldn't even come close to that)")
    }
    quant.cutoffs<-c(outlier.cutoff/2,1-outlier.cutoff/2)
    outlier.cutoff<-TRUE
  }else{
    outlier.cutoff<-FALSE
  }
  #little helper flag for later...
  LL.ind.flag<-FALSE
  
  if(inherits(tree,"contsimmap")){
    if(nedge(tree)<Nedge(tree)){
      if(length(.stored.nodes(tree))==1){
        stop('The make.lik.fun() function will support single subtrees in the future, but not yet')
      }else{
        stop('The make.lik.fun() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
      }
    }
    nsim<-dim(tree)[3]
    tree.info<-.contsimmap2treeinfo(tree)
    contsimmap.traits<-dimnames(tree)[[2]]
  }else{
    if(is.null(nsim)){
      if(inherits(tree,"phylo")){
        nsim<-1
      }else if(inherits(tree,"multiPhylo")){
        nsim<-length(tree)
      }
    }
    tree.info<-.get.tree.info(tree,nsim)
    contsimmap.traits<-NULL
  }
  topo.info<-.get.topo.info(tree.info[['tree']][[1]])
  #recycling trait.data appropriately--cannibalized from .get.lookup in make_Utilites.R
  #Works if treeIDs have been permuted I believe (well, maybe it would break if any treeID is represented by 0 simulations)
  #Except that parameter recycling responds to the "tree" object is arranged, NOT how treeID is permuted
  #For example, if Xsig2=(column) list(A,B), A will go for tree1 simulations and B for tree2 simulations, EVEN IF tree2 sims come before tree1 sims in contsimmap!
  #Share this quirk with diffusion, and may want to revisit later
  trait.info<-.proc.trait.data(trait.data,
                               tree.info[['tips']],tree.info[['nodes']],
                               tree.info[['edges']],tree.info[['ntrees']],tree.info[['treeID']],
                               nsim)
  ntraits<-length(trait.info[['traits']])
  mult.traits<-ntraits>1
  wgts<-log(rep(wgts,length.out=nsim))
  if(wgt.by.nobs){
    tmp.nobs<-apply(trait.info[["parsed.mis"]],2,function(ii) sum(unlist(ii,use.names=FALSE)))
    wgts<-wgts-tmp.nobs*log(2*pi)/2
  }
  wgts<-wgts-max(wgts)
  wgts<-wgts-log(sum(exp(wgts)))
  
  ####GETTING PARAMETER FORMULAE####
  
  #need to figure out how to deal with stored variables in parent environment, I realized...
  #I guess see if any variables exist in parent R environment and don't treat them as variables if so...
  #Breaks if formula is provided as non-formula object is expression, which is undesirable and should be simple to fix!
  Xvar<-.fix.par.formulae(Xvar,trait.info[['traits']],tree.info[['states']],contsimmap.traits)
  Mu<-.fix.par.formulae(mu,trait.info[['traits']],tree.info[['states']],constimmap.traits)
  Yvar<-.fix.par.formulae(Yvar,trait.info[['traits']],c(tree.info[['tips']],tree.info[['nodes']]),contsimmap.traits)
  #fix ordering-->from nodewise to edgewise
  #This does NOT seem to be working as expected somehow...
  Yperm<-c(length(tree.info[["tips"]])+1,tree.info[["edges"]][,2])
  attr(Yvar,"inds")<-attr(Yvar,"inds")[Yperm]
  if(ntraits>1){
    Xcor<-.fix.par.formulae(Xcor,trait.info[['traits']],tree.info[['states']],contsimmap.traits)
    Ycor<-.fix.par.formulae(Ycor,trait.info[['traits']],c(tree.info[['tips']],tree.info[['nodes']]),contsimmap.traits)
    attr(Ycor,"inds")<-attr(Ycor,"inds")[Yperm]
  }else{
    Ycor<-Xcor<-NULL
  }
  
  #oooh I think I got it
  #just store that the unlisted trait as normal
  #BUT also store the indices that need to be grabbed for any particular formula...
  #Then modify formulae by adding [[indices]] to the end!
  
  #Yvar seems pointlessly complicated--why store a matrix for each node? Just store a vector!
  #Need to look into this more--there might be a reason you went for this approach, but I think it was just convenience...
  #Will simplify things a bit for bringing in tree indexing in likelihood function
  #NEVER MIND! Keep in mind j refers to formulae, and formulae can refer to multiple states/nodes...
  
  tmp<-list(Xvar,Xcor,Mu,Yvar,Ycor)
  inds<-setNames(lapply(tmp,attr,"inds"),
                 c("Xvar","Xcor","mu","Yvar","Ycor"))
  vars.per.formula<-lapply(tmp,function(ii) lapply(ii,function(jj) 
    if(is.call(jj)) all.vars(jj) else lapply(jj,all.vars))) #had to modify to prevent treating make.cor as variable
  par.nms<-unique(unlist(vars.per.formula,use.names=FALSE))
  conmaps<-list()
  conmaps[["RES_mapped_mu"]]<-rep(list(rep(FALSE,ntraits)),length(Mu))
  conmaps[["RES_mapped_xvar"]]<-rep(list(rep(FALSE,ntraits)),length(Xvar))
  conmaps[["RES_mapped_yvar"]]<-rep(list(rep(FALSE,ntraits)),length(Yvar))
  conrefs<-par.nms%in%contsimmap.traits
  if(any(conrefs)){
    conmaps.flag<-TRUE
    
    #could shut these off depending on whether conmaps apply to xvar/mu and/or yvar
    #have to get obs stuff for yvar indexing...
    has.obs<-trait.info[['parsed.mis']]
    nedges<-nrow(has.obs)
    ntrait.data<-ncol(has.obs)
    has.any.obs<-lengths(has.obs)>0
    has.obs[!has.any.obs]<-list(rep(FALSE,ntraits))
    has.obs[has.any.obs]<-lapply(has.obs[has.any.obs],apply,1,any)
    has.obs<-array(unlist(has.obs,use.names=FALSE),c(ntraits,nedges,ntrait.data))[,,trait.info[["traitID"]],drop=FALSE]
    
    #have to get state des stuff for xvar/mu indexing...
    states<-.get.maps(tree,element="state")
    #add in root state! There should always be only 1
    states<-rbind(lapply(states[root.edges(tree)[1],],'[[',1),states)
    states<-unlist(states[,tree.info[["treeID"]],drop=FALSE],use.names=FALSE)
    states<-lapply(tree.info[["states"]],"==",states)
    has.des<-has.obs
    des<-topo.info[["des"]]
    for(e in topo.info[["prune.seq"]]){
      has.des[,e,]<-matrix(has.obs[,e,,drop=FALSE],ntraits,nsim)|apply(has.des[,des[[e]],,drop=FALSE],c(1,3),any)
    }
    nts<-rbind(1,.get.ns(tree))[,tree.info[["treeID"]],drop=FALSE]
    has.des<-lapply(seq_len(ntraits),function(ii) rep(has.des[ii,,,drop=FALSE],nts))
    trees<-rep(seq_along(tree.info[["treeID"]]),colSums(nts))
    conmaps[["RES_trees"]]<-trees
    
    #setting up conmaps list...
    tmp.endpts<-cumsum(nts)
    conmaps[["RES_nodes_trees"]]<-trees[tmp.endpts]
    conmaps[["RES_nodes_inds"]]<-tmp.endpts
    tmp.seq<-seq_along(nts)
    conmaps[["RES_splits"]]<-factor(rep(tmp.seq,nts),levels=tmp.seq)
    len.conmaps<-length(conmaps[["RES_splits"]])
    conmaps[["RES_mu_holder"]]<-numeric(len.conmaps)
    conmaps[["RES_xvar_holder"]]<-rep(1,len.conmaps)
    conmaps[["RES_yvar_holder"]]<-matrix(1,nedges,nsim)
    vars.per.formula<-vars.per.formula[c(1,3,4)]
    conrefs.nms<-par.nms[conrefs]
    for(i in seq_len(3)){
      for(j in seq_along(vars.per.formula[[i]])){
        for(k in seq_along(vars.per.formula[[i]][[j]])){
          tmp.conrefs<-unique(conrefs.nms[conrefs.nms%in%vars.per.formula[[i]][[j]][[k]]])
          if(length(tmp.conrefs)){
            if(i==3){
              conmaps[["RES_mapped_yvar"]][[j]][k]<-TRUE
              tmp.nm<-paste0("RES_Yvar_",j,"_",k,"_inds")
              conmaps[[tmp.nm]]<-
                as.vector(has.obs[k,inds[["Yvar"]]==j,])
              Yvar[[j]][[k]]<-do.call(substitute,
                                      list(expr=Yvar[[j]][[k]],
                                           env=setNames(lapply(paste0(tmp.conrefs,"[RES_nodes_inds[",tmp.nm,"]]"),str2lang),tmp.conrefs)))
            }else{
              tmp.nm<-paste0("RES_",c("Xvar","mu")[i],"_",j,"_",k,"_inds")
              conmaps[[tmp.nm]]<-
                Reduce("|",states[inds[[c("Xvar","mu")[i]]]==j])&has.des[[k]]
              if(i==1){
                conmaps[["RES_mapped_xvar"]][[j]][k]<-TRUE
                Xvar[[j]][[k]]<-do.call(substitute,
                                        list(expr=Xvar[[j]][[k]],
                                             env=setNames(lapply(paste0(tmp.conrefs,"[",tmp.nm,"]"),str2lang),tmp.conrefs)))
              }else{
                conmaps[["RES_mapped_mu"]][[j]][k]<-TRUE
                Mu[[j]][[k]]<-do.call(substitute,
                                      list(expr=Mu[[j]][[k]],
                                           env=setNames(lapply(paste0(tmp.conrefs,"[",tmp.nm,"]"),str2lang),tmp.conrefs)))
              }
            }
          }
        }
      }
    }
    
    #and finally add that traits in...
    edge.seq<-match(as.numeric(gsub("^N","",dimnames(tree)[[1]]))+1,seq_len(nedges))
    tmp<-unclass(tree)[edge.seq,,,drop=FALSE]
    for(i in conrefs.nms){
      conmaps[[i]]<-unlist(tmp[,i,,drop=FALSE],use.names=FALSE)
    }
    par.nms<-par.nms[!conrefs]
  }else{
    conmaps.flag<-FALSE
  }
  
  ####VERY GENERIC HELPER FUNCTIONS####
  
  #multiplies each slice of m*dim*n array by each column of n-row matrix (columns recycled)
  .mult.arr.mat<-function(arr,mat,m,dim,n){
    arr[,1,]<-as.vector(mat[rep(1,m),,drop=FALSE])*arr[,1,,drop=FALSE]
    if(dim>1){
      for(i in seq_len(dim-1)){
        arr[,1,]<-arr[,1,,drop=FALSE]+as.vector(mat[rep(i+1,m),,drop=FALSE])*arr[,i+1,,drop=FALSE]
      }
    }
    matrix(arr[,1,,drop=FALSE],m,n)
  }
  #does tcrossprod for each column of mat
  .tcross.mat<-function(mat,k,dim){
    tmp<-array(mat,c(k,dim,k))
    aperm.default(tmp,c(1,3,2))*aperm.default(tmp,c(3,1,2))
  }
  #does tcrossprod for each column-slice combo of array
  .tcross.arr<-function(arr,k,m,n){
    tmp<-array(arr,c(k,m,n,k))
    aperm.default(tmp,c(1,4,3,2))*aperm.default(tmp,c(4,1,3,2))
  }
  #sums elements of a dim-length list (equivalent to Reduce("+",ls))
  .sum.ls<-function(ls,dim){
    if(dim>1){
      for(i in seq_len(dim-1)){
        ls[[1]]<-ls[[1]]+ls[[i+1]]
      }
    }
    ls[[1]]
  }
  #for vectorized calculation of determinants...
  .log.prod.diag<-function(x,k,dim,inf=NULL){
    x<-x[rep(rep(c(TRUE,FALSE),c(1,k)),length.out=k^2)]
    if(!is.null(inf)) x[inf]<-1
    .colSums(log(x),k,dim)
  }
  #convert TRUE/FALSE vectors to 0/1 sequences across matrices and arrays...
  .inds2codes.mat<-function(mat,m,n){
    out<-rep(list(rep("0",m)),n)
    for(i in seq_len(n)){
      out[[i]][mat[,i]]<-"1"
    }
    do.call(paste0,out)
  }
  .inds2codes.arr<-function(arr,k,m,n){
    out<-rep(list(matrix("0",m,n)),k)
    for(i in seq_len(k)){
      out[[i]][arr[i,,]]<-"1"
    }
    matrix(do.call(paste0,out),m,n)
  }
  .tcross.OR.mat<-function(mat,k,dim){
    tmp<-array(mat,c(k,dim,k))
    aperm.default(tmp,c(1,3,2))|aperm.default(tmp,c(3,1,2))
  }
  
  
  ####PRE-PROCESSING STUFF FOR FINAL FUNCTION####
  
  ##setting tree/traitID stuff##
  treeID<-tree.info[['treeID']]
  traitID<-trait.info[['traitID']]
  traitID.inds<-lapply(unique(traitID),'==',traitID)
  treeID.inds<-lapply(unique(treeID),'==',treeID)
  ntraitID<-length(traitID.inds)
  ntreeID<-length(treeID.inds)
  traitID.seq<-seq_len(ntraitID)
  treeID.seq<-seq_len(ntreeID)
  sims.per.traitID<-unlist(lapply(traitID.inds,sum),use.names=FALSE)
  sims.per.treeID<-unlist(lapply(treeID.inds,sum),use.names=FALSE)
  traitID.per.treeID<-lapply(treeID.inds,function(ii) unique(traitID[ii]))
  
  ##other tree stuff##
  states<-tree.info[['states']]
  prune.seq<-topo.info[['prune.seq']]
  des<-topo.info[['des']]
  ndes<-lengths(des)
  has.des<-ndes>0
  if(conmaps.flag){
    maps<-.get.maps(tree,element="dts")
    tmp<-.get.maps(tree,element="state")
    maps[]<-lapply(seq_along(maps),function(ii) setNames(maps[[ii]],tmp[[ii]]))
  }else{
    maps<-do.call(cbind,lapply(tree.info[["tree"]],'[[','maps'))
    maps[]<-lapply(maps,function(ii) tapply(ii,names(ii),sum))
  }
  maps<-rbind(list(NULL),maps)
  map.lens<-lengths(maps)
  nedges<-nrow(maps)
  #just indicators to see if a given edge/treeID combo include mapped xvars or mus
  xvar.maps<-matrix(unlist(lapply(maps,function(ii) any(unlist(conmaps[["RES_mapped_xvar"]][inds[["Xvar"]][unique(names(ii))]],use.names=FALSE))),use.names=FALSE),nedges,ntreeID)
  mu.maps<-matrix(unlist(lapply(maps,function(ii) any(unlist(conmaps[["RES_mapped_mu"]][inds[["mu"]][unique(names(ii))]],use.names=FALSE))),use.names=FALSE),nedges,ntreeID)
  
  ##other trait stuff##
  nobs<-trait.info[['nobs']]
  has.obs<-nobs>0
  parsed.mis<-trait.info[['parsed.mis']]
  parsed.obs<-trait.info[['parsed.obs']]
  parsed.obs[]<-lapply(parsed.obs,function(ii) if(!is.null(ii)) t(ii))
  #let's first convert parsed.mis to codes...
  #eventually may just want to integrate this into get.trait.info!
  #I feel like this could be simplified...but it seems to get the job done for now...
  lens<-lengths(parsed.mis)
  tmp.seq<-seq_along(parsed.mis)
  splits<-factor(rep(tmp.seq,lens/ntraits),levels=tmp.seq)
  tmp.m<-sum(lens)/ntraits
  tmp<-matrix(unlist(parsed.mis,use.names=FALSE),tmp.m,ntraits,byrow=TRUE)
  obs.codes<-.inds2codes.mat(tmp,tmp.m,ntraits)
  not.dups<-!duplicated(obs.codes)
  code.inds<-tmp[not.dups,,drop=FALSE]
  tmp.obs.codes<-obs.codes[not.dups]
  obs.codes<-match(obs.codes,tmp.obs.codes)
  obs.codes<-matrix(split(obs.codes,splits),nedges,ntraitID)
  tmp<-array(FALSE,dim=c(ntraits,nedges,ntraitID))
  has.x<-has.des|has.obs
  for(e in prune.seq){
    for(i in traitID.seq){
      if(has.x[e,i]){
        tmp[,e,i]<-apply(matrix(c(parsed.mis[[e,i]],tmp[,des[[e]],i,drop=FALSE]),ntraits,nobs[e,i]+ndes[e]),
                         1,
                         any)
      }
    }
  }
  treeID.tmp<-array(unlist(lapply(traitID.per.treeID,function(ii) apply(tmp[,,ii,drop=FALSE],c(1,2),any)),use.names=FALSE),
                    c(ntraits,nedges,ntreeID))
  des.codes<-.inds2codes.arr(tmp,ntraits,nedges,ntraitID)
  treeID.des.codes<-.inds2codes.arr(treeID.tmp,ntraits,nedges,ntreeID)
  not.dups<-!duplicated(c(des.codes,treeID.des.codes))
  new.codes<-!(c(des.codes,treeID.des.codes)[not.dups]%in%tmp.obs.codes)
  if(any(new.codes)){
    tmp.obs.codes<-c(tmp.obs.codes,c(des.codes,treeID.des.codes)[not.dups][new.codes])
    code.inds<-rbind(code.inds,t(matrix(c(tmp,treeID.tmp),ntraits,nedges*(ntraitID+ntreeID))[,not.dups,drop=FALSE][,new.codes,drop=FALSE]))
  }
  des.codes<-matrix(match(des.codes,tmp.obs.codes),nedges,ntraitID)
  treeID.des.codes<-matrix(match(treeID.des.codes,tmp.obs.codes),nedges,ntreeID)
  code.k<-rowSums(code.inds)
  all.code<-which(code.k==ntraits)
  if(!length(all.code)){
    code.inds<-rbind(code.inds,TRUE)
    code.k<-c(code.k,ntraits)
    all.code<-nrow(code.inds)
    tmp.obs.codes<-c(tmp.obs.codes,paste0(rep("1",ntraits),collapse=""))
  }
  non.code<-which(code.k==0)
  if(!length(non.code)){
    code.inds<-rbind(code.inds,FALSE)
    code.k<-c(code.k,0)
    non.code<-nrow(code.inds)
    tmp.obs.codes<-c(tmp.obs.codes,paste0(rep("0",ntraits),collapse=""))
  }
  #setting up inds with obs.code system...
  tmp.inds<-attr(Ycor,'inds')
  tmp<-length(Ycor)
  codes.per.cor<-vector("list",tmp)
  for(i in seq_len(tmp)){
    codes.per.cor[[i]]<-unique(unlist(obs.codes[tmp.inds==i,,drop=FALSE],use.names=FALSE))
  }
  #let's also reduce obs.codes a bit to save on computation later...
  foo<-function(x){
    tmp<-unique(x)
    list(tmp,match(x,tmp))
  }
  obs.codes[]<-lapply(obs.codes,foo)
  
  ##holders##
  #some reminders:
  #R = remainders (scalar stuff--likelihood is a scaled multivariate normal)
  #Z = trait value vectors
  #V = variance-covariance matrices
  #P = precision matrices (inverse of V)
  #vectors (for scalars)
  RR.holder<-numeric(nsim)
  RR<-rep(list(RR.holder),nedges)
  if(mult.traits){
    #matrices (for vectors)
    PZ.holder<-matrix(0,ntraits,nsim)
    ZZ<-rep(list(PZ.holder),nedges)
    #arrays (for covariance matrices)
    VV.holder<-array(0,c(ntraits,ntraits,nsim))
    VV<-rep(list(VV.holder),nedges)
  }else{
    VV.holder<-PZ.holder<-RR.holder
    VV<-ZZ<-RR
  }
  
  ##more specific helper functions##
  #may want to add some kind of functionality to support Ycor/var varying by state, I realized...
  #or it may be best just to allow s_RES and t_RES to be in contsimmap traits again...
  #   except that means Ycor/Yvar must be reevaluated for each possible s_RES/t_RES...
  #   maybe best in these cases to let these things be manually specified, but that would make it a lot more frustrating to let Ysig2 vary by state...
  #   yeah, not important enough at the moment, maybe just something to think about going forward...
  calc.mu<-function(formula,env){
    eval(formula,env=env)
  }
  if(mult.traits){
    calc.xvar<-function(formula,env){
      out<-eval(formula,env=env)
      if(any(out<0|is.na(out)|is.infinite(out))){
        NULL
      }else{
        out[out<vartol]<-0
        sqrt(out)
      }
    }
    calc.yvar<-function(formula,env){
      out<-eval(formula,env=env)
      if(any(out<0|is.na(out)|is.infinite(out))){
        NULL
      }else{
        out[out<vartol]<-0
        1/sqrt(out)
      }
    }
  }else{
    calc.xvar<-function(formula,env){
      out<-eval(formula,env=env)
      if(any(out<0|is.na(out)|is.infinite(out))){
        NULL
      }else{
        out[out<vartol]<-0
        out
      }
    }
    calc.yvar<-function(formula,env){
      out<-eval(formula,env=env)
      if(any(out<0|is.na(out)|is.infinite(out))){
        NULL
      }else{
        out[out<vartol]<-0
        1/out
      }
    }
  }
  calc.mus<-function(env,tree.flag,nsim){
    for(i in seq_along(Mu)){
      for(j in seq_len(ntraits)){
        if(tree.flag){
          env[[paste0("RES_mu_",i,"_",j,"_inds")]]<-
            env[[paste0("RES_mu_",i,"_",j,"_inds")]]&env[["RES_trees"]]
        }
        tmp<-calc.mu(Mu[[i]][[j]],env)
        if(is.null(tmp)){
          return(NULL)
        }else{
          if(env[["RES_mapped_mu"]][[i]][j]){
            Mu[[i]][[j]]<-env[["RES_mu_holder"]]
            if(tree.flag){
              Mu[[i]][[j]][env[[paste0("RES_mu_",i,"_",j,"_inds")]][env[["RES_trees"]]]]<-tmp
            }else{
              Mu[[i]][[j]][env[[paste0("RES_mu_",i,"_",j,"_inds")]]]<-tmp
            }
            Mu[[i]][[j]]<-matrix(split(Mu[[i]][[j]],env[["RES_splits"]]),nedges,nsim)
          }else{
            Mu[[i]][[j]]<-tmp
          }
        }
      }
    }
    Mu
  }
  calc.yvars<-function(env,tree.flag){
    for(i in seq_along(Yvar)){
      for(j in seq_len(ntraits)){
        if(tree.flag){
          env[[paste0("RES_Yvar_",i,"_",j,"_inds")]]<-
            env[[paste0("RES_Yvar_",i,"_",j,"_inds")]]&env[["RES_nodes_trees"]]
        }
        tmp<-calc.yvar(Yvar[[i]][[j]],env)
        if(is.null(tmp)){
          return(NULL)
        }else{
          if(env[["RES_mapped_yvar"]][[i]][j]){
            Yvar[[i]][[j]]<-env[["RES_yvar_holder"]]
            if(tree.flag){
              Yvar[[i]][[j]][env[[paste0("RES_Yvar_",i,"_",j,"_inds")]][env[["RES_nodes_trees"]]]]<-tmp
            }else{
              Yvar[[i]][[j]][env[[paste0("RES_Yvar_",i,"_",j,"_inds")]]]<-tmp
            }
          }else{
            Yvar[[i]][[j]]<-tmp
          }
        }
      }
    }
    Yvar
  }
  #could potentially do smarter splits based by reorganizing simulation to be in a better order...
  #but I'm not gonna worry about it for now
  calc.xvars<-function(env,tree.flag,nsim){
    for(i in seq_along(Xvar)){
      for(j in seq_len(ntraits)){
        if(tree.flag){
          env[[paste0("RES_Xvar_",i,"_",j,"_inds")]]<-
            env[[paste0("RES_Xvar_",i,"_",j,"_inds")]]&env[["RES_trees"]]
        }
        tmp<-calc.xvar(Xvar[[i]][[j]],env)
        if(is.null(tmp)){
          return(NULL)
        }else{
          if(env[["RES_mapped_xvar"]][[i]][j]){
            Xvar[[i]][[j]]<-env[["RES_xvar_holder"]]
            if(tree.flag){
              Xvar[[i]][[j]][env[[paste0("RES_Xvar_",i,"_",j,"_inds")]][env[["RES_trees"]]]]<-tmp
            }else{
              Xvar[[i]][[j]][env[[paste0("RES_Xvar_",i,"_",j,"_inds")]]]<-tmp
            }
            Xvar[[i]][[j]]<-matrix(split(Xvar[[i]][[j]],env[["RES_splits"]]),nedges,nsim)
          }else{
            Xvar[[i]][[j]]<-tmp
          }
        }
      }
    }
    Xvar
  }
  if(mult.traits){
    cor.holder<-matrix(0,nrow=ntraits,ncol=ntraits)
    diag.inds<-rep(rep(c(TRUE,FALSE),c(1,ntraits)),length.out=ntraits^2)
    odiag.inds<-!diag.inds
    calc.cor<-function(cor,env){
      if(is.matrix(cor)){
        cor.holder[odiag.inds]<-unlist(lapply(cor,eval,envir=env),use.names=FALSE)
        cor.holder[diag.inds]<-1
        cor.holder
      }else{
        eval(cor,envir=env)
      }
    }
    #technically these checks are unnecessary if cor was constructed via make.cor function...
    #but it really doesn't add to the computational burden of the function anyways in the face of all the other stuff going on
    check.cor<-function(cor){
      any(eigen(cor,symmetric=TRUE,only.values=TRUE)[['values']]<=0)
    }
    check.cors<-function(cors){
      any(unlist(lapply(cors,check.cor),use.names=FALSE))
    }
    #converting cors.holder to a kind of "cache"
    ncodes<-nrow(code.inds)
    tmp<-setNames(rep(list(cor.holder),ncodes),tmp.obs.codes)
    attr(tmp,"log.dets")<-setNames(numeric(ncodes),tmp.obs.codes)
    attr(tmp,"calc.flag")<-setNames(logical(ncodes),tmp.obs.codes)
    cors.holder<-rep(list(tmp),length(Ycor))
    invert.ycor<-function(cor){
      for(i in seq_along(cor)){
        for(j in codes.per.cor[[i]]){
          cors.holder[[i]][[j]][code.inds[j,],code.inds[j,]]<-solve.default(cor[[i]][code.inds[j,],code.inds[j,],drop=FALSE])
          attr(cors.holder[[i]],"log.dets")[j]<-determinant.matrix(cors.holder[[i]][[j]][code.inds[j,],code.inds[j,],drop=FALSE])[["modulus"]]
        }
        #can add dimensional correction directly into determinants!
        attr(cors.holder[[i]],"log.dets")<-attr(cors.holder[[i]],"log.dets")-code.k*log(2*pi)
        attr(cors.holder[[i]],"calc.flag")[c(codes.per.cor[[i]],non.code)]<-TRUE
      }
      cors.holder
    }
    calc.cors<-function(env){
      xcor<-lapply(Xcor,calc.cor,env=env)
      if(check.cors(xcor)){
        return(NULL)
      }
      ycor<-lapply(Ycor,calc.cor,env=env)
      if(check.cors(ycor)){
        return(NULL)
      }
      list(xcor,invert.ycor(ycor),ycor)
    }
    .get.edgewise.xvar.mu<-function(ls,e,m,n,k,sim.inds,trait.inds,map.inds){
      out<-array(dim=c(k,m,n))
      #could potentially pre-compute some of this stuff...
      for(i in seq_along(ls)){
        tmp<-map.inds==i
        if(any(tmp)){
          for(j in seq_len(k)){
            if(is.list(ls[[i]][[trait.inds[j]]])){
              out[j,tmp,]<-unlist(ls[[i]][[trait.inds[j]]][e,sim.inds],use.names=FALSE)[tmp]
            }else{
              out[j,tmp,]<-ls[[i]][[trait.inds[j]]]
            }
          }
        }
      }
      out
    }
    #calculates V and Z for given descendant edge
    calc.des<-function(e,i,m,n,k,sim.inds,trait.inds,xvar,mu,xcor,V,Z){
      tmp.map<-maps[[e,i]]
      tmp.states<-names(tmp.map)
      if(xvar.maps[e,i]){
        tmp.xvar<-.tcross.arr(.get.edgewise.xvar.mu(xvar,e,m,n,k,sim.inds,which(trait.inds),inds[["Xvar"]][tmp.states]),k,m,n)
      }else{
        tmp.xvar<-.tcross.arr(matrix(unlist(xvar[inds[["Xvar"]][tmp.states]],use.names=FALSE)[trait.inds],c(k,m,1)),k,m,1)
      }
      if(mu.maps[e,i]){
        tmp.mu<-.get.edgewise.xvar.mu(mu,e,m,n,k,sim.inds,which(trait.inds),inds[["mu"]][tmp.states])
      }else{
        tmp.mu<-array(unlist(mu[inds[["mu"]][tmp.states]],use.names=FALSE)[trait.inds],c(k,m,1))
      }
      tmp.xcor<-lapply(xcor,function(ii) ii[trait.inds,trait.inds,drop=FALSE])[inds[["Xcor"]][tmp.states]]
      
      tmp.xvar[,,,1]<-tmp.map[[1]]*tmp.xvar[,,,1,drop=FALSE]*as.vector(tmp.xcor[[1]])
      tmp.mu[,1,]<-tmp.map[[1]]*tmp.mu[,1,,drop=FALSE]
      if(m>1){
        for(j in seq_len(m)[-1]){
          tmp.xvar[,,,1]<-tmp.xvar[,,,1,drop=FALSE]+tmp.map[[j]]*tmp.xvar[,,,j,drop=FALSE]*as.vector(tmp.xcor[[j]])
          tmp.mu[,1,]<-tmp.mu[,1,,drop=FALSE]-tmp.map[[j]]*tmp.mu[,j,,drop=FALSE]
        }
      }
      list(V+array(tmp.xvar[,,,1,drop=FALSE],c(k,k,n)),Z-matrix(tmp.mu[,1,,drop=FALSE],k,n))
    }
    #calculates sums of P, PZ, and associated scalars for all observations of given edge
    calc.obs<-function(e,i,n,nobs,yvar,in.ycor,inf,des.inf,base.ycor){
      codes<-obs.codes[[e,i]]
      matches<-codes[[2]]
      codes<-codes[[1]]
      any.obs.inf<-!is.null(inf)
      any.des.inf<-!is.null(des.inf)
      if(any.des.inf){
        n<-ncol(des.inf)
        yvar<-array(yvar,c(ntraits,ntraits,n))
      }
      tmp<-rep(log(yvar[diag.inds]),length.out=n*ntraits)
      if(any.obs.inf|any.des.inf){
        P<-det<-vector("list",length(codes))
        if(any.obs.inf&any.des.inf){
          tmp.inf<-!(as.vector(inf)|des.inf)
        }else if(any.obs.inf){
          tmp.inf<-!inf
        }else{
          tmp.inf<-!des.inf
        }
        cached<-attr(in.ycor,"calc.flag")
        tmp.nms<-names(cached)
        counter<-0
        for(j in codes){
          counter<-counter+1
          tmp.inds<-code.inds[j,]&tmp.inf
          tmp.codes<-.inds2codes.mat(t(tmp.inds),n,ntraits)
          not.dups<-!duplicated(tmp.codes)
          new.codes<-which(not.dups)[!(tmp.codes[not.dups]%in%tmp.nms[cached])]
          #maybe an inefficient expansion here, but hopefully doesn't occur enough to warrant worrying about it...
          #I'm more worried about having to recalculate codes for each observation... :S
          if(length(new.codes)){
            for(k in new.codes){
              tmp.nm<-tmp.codes[k]
              tmp.inds2<-tmp.inds[,k]
              in.ycor[[tmp.nm]]<-in.ycor[[non.code]]
              in.ycor[[tmp.nm]][tmp.inds2,tmp.inds2]<-tmp.cor<-solve(base.ycor[tmp.inds2,tmp.inds2,drop=FALSE])
              if(tmp.nm%in%tmp.nms){
                attr(in.ycor,"log.det")[tmp.nm]<-determinant.matrix(tmp.cor)[["modulus"]]-sum(tmp.inds2)*log(2*pi)
                attr(in.ycor,"calc.flag")[tmp.nm]<-TRUE
              }else{
                attr(in.ycor,"log.det")<-c(attr(in.ycor,"log.det"),
                                           setNames(determinant.matrix(tmp.cor)[["modulus"]]-sum(tmp.inds2)*log(2*pi),tmp.nm))
                attr(in.ycor,"calc.flag")<-c(attr(in.ycor,"calc.flag"),
                                             setNames(TRUE,tmp.nm))
              }
            }
            cached<-attr(in.ycor,"calc.flag")
            tmp.nms<-names(cached)
          }
          tmp.matches<-match(tmp.codes,tmp.nms)
          #maybe not the most efficient...but you have get 1s in here...
          tmp.ycor<-unlist(in.ycor[tmp.matches],use.names=FALSE)
          if(any.obs.inf){
            tmp.ycor[diag.inds][code.inds[j,]&inf]<-1
          }
          P[[counter]]<-tmp.ycor*yvar
          det[[counter]]<-attr(in.ycor,"log.det")[tmp.matches]+.colSums(tmp[code.inds[j,]],code.k[j],n)
        }
      }else{
        det<-lapply(codes,function(ii) attr(in.ycor,"log.det")[ii]+.colSums(tmp[code.inds[ii,]],code.k[ii],n))
        ycor<-in.ycor[codes]
        P<-lapply(seq_along(ycor),function(ii) as.vector(ycor[[ii]])*yvar)
      }
      tmp.zs<-parsed.obs[[e,i]]
      tmp.n<-dim(yvar)[3]
      tmp.z<-tmp.zs[,1,drop=FALSE]
      out.P<-tmp.P<-P[[1]]
      out.PZ<-tmp.PZ<-.mult.arr.mat(tmp.P,tmp.z,ntraits,ntraits,tmp.n)
      if(any.obs.inf) tmp.PZ[inf]<-0
      out.det<-det[[1]]-.colSums(as.vector(tmp.z)*tmp.PZ,ntraits,tmp.n)
      if(nobs>1){
        for(j in seq_len(nobs)[-1]){
          tmp.z<-tmp.zs[,j,drop=FALSE]
          tmp.match<-matches[j]
          tmp.P<-P[[tmp.match]]
          out.P<-out.P+tmp.P
          tmp.PZ<-.mult.arr.mat(tmp.P,tmp.z,ntraits,ntraits,tmp.n)
          out.PZ<-out.PZ+tmp.PZ
          #need to do this such that things cancel nicely in the next step...
          if(any.obs.inf) tmp.PZ[inf]<-0
          out.det<-out.det+det[[tmp.match]]-.colSums(as.vector(tmp.z)*tmp.PZ,ntraits,tmp.n)
        }
      }
      #if tmp.n is 1, should be fine because it will get replicated when inserting into larger array
      list(out.P,out.PZ,out.det,in.ycor)
    }
    #calculates V as well as sums PZ and associated scalars for a given edge
    #(also adds determinant of V and dimensional correction to associated scalars)
    des.inf.holder<-matrix(FALSE,ntraits,nsim)
    calc.all<-function(e,xvar,yvar,cors,mu,VV,ZZ,RR,double.RR,
                       traitID.inds,treeID.inds,sims.per.traitID,sims.per.treeID,traitID.seq,treeID.seq,nsim,
                       full.has.obs,full.has.x,traitID,
                       RR.holder,PZ.holder,VV.holder,des.inf.holder){
      #need to actually start keeping track of infinite precisions (aka 0 variances) here!
      #I think this all works, but hard to tell--basically just keeps track of variance-covariance matrices with 0 along diagonal
      #Then uses this information in get.obs() below (but importantly avoids counting observations as infinite)
      if(has.des[e]){
        des.list<-rep(list(list(VV.holder,PZ.holder,des.inf.holder)),ndes[e])
        counter<-0
        for(d in des[[e]]){
          counter<-counter+1
          for(i in treeID.seq){
            tmp.code<-treeID.des.codes[d,i]
            if(tmp.code!=non.code){
              tmp.inds<-treeID.inds[[i]]
              tmp.inds2<-code.inds[tmp.code,]
              tmp<-calc.des(d,i,
                            map.lens[d,i],sims.per.treeID[i],code.k[tmp.code],
                            tmp.inds,tmp.inds2,
                            xvar,mu,cors[[1]],
                            VV[[d]][tmp.inds2,tmp.inds2,tmp.inds,drop=FALSE],
                            ZZ[[d]][tmp.inds2,tmp.inds,drop=FALSE])
              des.list[[counter]][[1]][tmp.inds2,tmp.inds2,tmp.inds]<-tmp[[1]]
              des.list[[counter]][[2]][tmp.inds2,tmp.inds]<-tmp[[2]]
              des.list[[counter]][[3]][tmp.inds2,tmp.inds]<-tmp[[1]][rep(rep(c(TRUE,FALSE),c(1,code.k[tmp.code])),length.out=code.k[tmp.code]^2)]==0
              des.inf.holder<-des.inf.holder|des.list[[counter]][[3]]
            }
          }
        }
        has.des.inf<-any(des.inf.holder)
      }else{
        has.des.inf<-FALSE
      }
      tmp.nobs<-nobs[e,]
      has.obs.inf<-FALSE
      if(any(tmp.nobs>0)){
        tmp.yvar<-do.call(rbind,lapply(yvar[[inds[["Yvar"]][e]]],function(ii) if(length(ii)>1) ii[e,,drop=FALSE] else ii))
        #taking care of infinite precisions with 2 modifications...
        #ad-hoc --> set infinite precisions to 0 to get 0 correlation with other trait dimensions
        #       --> BUT set diagonal entry to 1 afterwards (cancels out when calculating determinants that way)
        #more fundamental --> unfortunately still need to account for how infinite precisions affect observation codes
        #                 --> exact measurements essentially don't "count"
        #                 --> solution: if any infinite precisions are present, recalculate observation codes
        #                 --> then add new ycor inverses/determinants as necessary to a cache
        #                 --> can be carried over as needed in subsequent steps!
        #                 --> requires "base.ycor" which is just the uninverted correlation matrix with no dropped dimensions
        single.yvar<-ncol(tmp.yvar)==1
        obs.inf<-is.infinite(tmp.yvar)
        has.obs.inf<-any(obs.inf)
        if(has.obs.inf){
          tmp.yvar[obs.inf]<-0
        }
        if(single.yvar){
          tmp.yvar<-array(tcrossprod(tmp.yvar),c(ntraits,ntraits,1))
        }else{
          tmp.yvar<-.tcross.mat(tmp.yvar,ntraits,nsim)
        }
        if(has.obs.inf){
          for(i in seq_len(ntraits)){
            tmp.yvar[i,i,obs.inf[i,]]<-1
          }
        }
        tmp.ind<-inds[["Ycor"]][e]
        tmp.ycor<-cors[[2]][[tmp.ind]]
        base.ycor<-cors[[3]][[tmp.ind]]
      }
      for(i in traitID.seq){
        tmp.code<-des.codes[e,i]
        if(tmp.code!=non.code){
          tmp.inds<-traitID.inds[[i]]
          tmp.n<-sims.per.traitID[i]
          tmp.obs.inf<-if(has.obs.inf) {if(single.yvar) obs.inf else obs.inf[,tmp.inds,drop=FALSE]} else NULL
          tmp.has.obs.inf<-any(tmp.obs.inf)
          tmp.des.inf<-if(has.des.inf) des.inf.holder[,tmp.inds,drop=FALSE] else NULL
          tmp.has.des.inf<-any(tmp.des.inf)
          if(tmp.nobs[i]){
            tmp<-calc.obs(e,i,
                          if(single.yvar) 1 else tmp.n,
                          tmp.nobs[i],
                          if(single.yvar) tmp.yvar else tmp.yvar[,,tmp.inds,drop=FALSE],
                          tmp.ycor,
                          if(tmp.has.obs.inf) tmp.obs.inf else NULL,
                          if(tmp.has.des.inf) tmp.des.inf else NULL,
                          base.ycor)
            VV.holder[,,tmp.inds]<-tmp[[1]]
            PZ.holder[,tmp.inds]<-tmp[[2]]
            RR.holder[tmp.inds]<-tmp[[3]]
            cors[[2]][[tmp.ind]]<-tmp.ycor<-tmp[[4]]
          }
          #now the issue is to account for infinite precisions below...
          if(has.des[e]){
            counter<-0
            for(d in des[[e]]){
              counter<-counter+1
              tmp.code2<-des.codes[d,i]
              if(tmp.code2!=non.code){
                tmp.inds2<-code.inds[tmp.code2,]
                tmp.k<-code.k[tmp.code2]
                tmp.PP<-des.list[[counter]][[1]][tmp.inds2,tmp.inds2,tmp.inds,drop=FALSE]
                tmp.PZ<-tmp.ZZ<-des.list[[counter]][[2]][tmp.inds2,tmp.inds,drop=FALSE]
                tmp.det<-numeric(tmp.n)
                
                #finding infinite precisions for focal and non-focal observations...
                tmp.foc.inf<-if(tmp.has.des.inf) des.list[[counter]][[3]][tmp.inds2,tmp.inds,drop=FALSE] else NULL
                tmp.has.foc.inf<-any(tmp.foc.inf)
                if(tmp.has.des.inf|tmp.has.obs.inf){
                  if(tmp.has.des.inf&tmp.has.obs.inf){
                    tmp.oth.inf<-as.vector(tmp.obs.inf[tmp.inds2,,drop=FALSE])|tmp.des.inf[tmp.inds2,,drop=FALSE]
                  }else if(tmp.has.des.inf){
                    tmp.oth.inf<-tmp.des.inf[tmp.inds2,,drop=FALSE]
                  }else{
                    tmp.oth.inf<-tmp.obs.inf[tmp.inds2,,drop=FALSE]
                  }
                }else{
                  tmp.oth.inf<-NULL
                }
                tmp.has.oth.inf<-any(tmp.oth.inf)
                if(tmp.has.foc.inf|tmp.has.oth.inf){ #this is where we handle cases with infinite precisions
                  if(tmp.has.foc.inf&tmp.has.oth.inf){
                    tmp.inf<-!(as.vector(tmp.oth.inf)|tmp.foc.inf)
                  }else if(tmp.has.foc.inf){
                    tmp.inf<-!tmp.foc.inf
                  }else{
                    tmp.inf<-!tmp.oth.inf
                  }
                  tmp.PP[.tcross.OR.mat(!tmp.inf,tmp.k,ncol(tmp.inf))]<-0
                  if(tmp.has.foc.inf){
                    tmp.PP[rep(rep(c(TRUE,FALSE),c(1,tmp.k)),length.out=tmp.k^2)][tmp.foc.inf]<-1
                  }
                  if(tmp.k==1){
                    tmp.inds3<-as.vector(tmp.inf)
                    tmp.PP[tmp.inds3]<-1/tmp.PP[tmp.inds3]
                    tmp.PZ<-as.vector(tmp.PP)*tmp.ZZ
                    tmp.det[tmp.inds3]<-log(tmp.PP[tmp.inds3])-log(2*pi)
                  }else{
                    tmp.chol<-tmp.PP
                    tmp.ks<-.colSums(tmp.inf,tmp.k,ncol(tmp.inf))
                    for(j in seq_len(tmp.n)[tmp.ks>0]){
                      tmp.inds3<-tmp.inf[,j]
                      tmp.chol[tmp.inds3,tmp.inds3,j]<-chol(tmp.PP[tmp.inds3,tmp.inds3,j])
                      tmp.PP[tmp.inds3,tmp.inds3,j]<-chol2inv(tmp.chol[tmp.inds3,tmp.inds3,j],tmp.ks[j])
                    }
                    tmp.PZ<-.mult.arr.mat(tmp.PP,tmp.ZZ,tmp.k,tmp.k,tmp.n)
                    tmp.det<- -2*.log.prod.diag(tmp.chol,tmp.k,tmp.n,if(tmp.has.oth.inf) tmp.oth.inf else NULL)-tmp.ks*log(2*pi)
                  }
                }else{ #this is where we handle normal cases
                  if(tmp.k==1){
                    tmp.PP<-1/tmp.PP
                    tmp.PZ<-as.vector(tmp.PP)*tmp.ZZ
                    tmp.det<-log(tmp.PP)-log(2*pi)
                  }else{
                    tmp.chol<-tmp.PP
                    for(j in seq_len(tmp.n)){
                      tmp.chol[,,j]<-chol(tmp.PP[,,j])
                      tmp.PP[,,j]<-chol2inv(tmp.chol[,,j],tmp.k)
                    }
                    tmp.PZ<-.mult.arr.mat(tmp.PP,tmp.ZZ,tmp.k,tmp.k,tmp.n)
                    tmp.det<- -2*.log.prod.diag(tmp.chol,tmp.k,tmp.n)-tmp.k*log(2*pi)
                  }
                }
                VV.holder[tmp.inds2,tmp.inds2,tmp.inds]<-VV.holder[tmp.inds2,tmp.inds2,tmp.inds,drop=FALSE]+tmp.PP
                PZ.holder[tmp.inds2,tmp.inds]<-PZ.holder[tmp.inds2,tmp.inds,drop=FALSE]+tmp.PZ
                if(tmp.has.foc.inf) tmp.PZ[tmp.foc.inf]<-0
                RR.holder[tmp.inds]<-RR.holder[tmp.inds]+tmp.det-.colSums(tmp.ZZ*tmp.PZ,tmp.k,tmp.n)
              }
            }
          }
          #an annoying quirk here is that tmp.VV is fully expanded even when not necessary (since it's expanded following calc.obs)
          #would be annoying to fix, but not impossible--just would need to check for single.yvar and no descendant edges
          tmp.inds2<-code.inds[tmp.code,]
          tmp.k<-code.k[tmp.code]
          tmp.VV<-VV.holder[tmp.inds2,tmp.inds2,tmp.inds,drop=FALSE]
          tmp.det<-numeric(tmp.n)
          if(tmp.has.des.inf|tmp.has.obs.inf){
            if(tmp.has.des.inf&tmp.has.obs.inf){
              tmp.inf<-!(as.vector(tmp.obs.inf[tmp.inds2,,drop=FALSE])|tmp.des.inf[tmp.inds2,,drop=FALSE])
            }else if(tmp.has.des.inf){
              tmp.inf<-!(tmp.des.inf[tmp.inds2,,drop=FALSE])
            }else{
              tmp.inf<-!(tmp.obs.inf[tmp.inds2,,drop=FALSE])
            }
            if(tmp.k==1){
              tmp.inds3<-as.vector(!tmp.inf)
              #check for invalid precision matrices
              #if precision is greater than 1 for an exact measurement, then there are multiple exact measurements of the same node!
              #technically would be fine if measurements agree, but that's such an edge case I don't think it's worth considering at this point
              if(any(tmp.VV[tmp.inds3]>1)) return(NULL)
              tmp.det[tmp.inf]<-log(tmp.VV[tmp.inf])-log(2*pi)
              tmp.VV<-1/tmp.VV
            }else{
              #check for invalid precision matrices as above
              if(any(tmp.VV[rep(rep(c(TRUE,FALSE),c(1,tmp.k)),length.out=tmp.k^2)][!tmp.inf]>1)) return(NULL)
              tmp.chol<-tmp.VV
              tmp.ks<-.colSums(tmp.inf,tmp.k,ncol(tmp.inf))
              for(j in seq_len(tmp.n)){
                tmp.chol[,,j]<-chol.default(tmp.VV[,,j])
                tmp.VV[,,j]<-chol2inv(tmp.chol[,,j],tmp.k)
              }
              tmp.det<-2*.log.prod.diag(tmp.chol,tmp.k,tmp.n,!tmp.inf)-tmp.ks*log(2*pi)
            }
          }else{
            if(tmp.k==1){
              tmp.det<-log(tmp.VV)-log(2*pi)
              tmp.VV<-1/tmp.VV
            }else{
              tmp.chol<-tmp.VV
              for(j in seq_len(tmp.n)){
                tmp.chol[,,j]<-chol.default(tmp.VV[,,j])
                tmp.VV[,,j]<-chol2inv(tmp.chol[,,j],tmp.k)
              }
              tmp.det<-2*.log.prod.diag(tmp.chol,tmp.k,tmp.n)-tmp.k*log(2*pi)
            }
          }
          VV.holder[tmp.inds2,tmp.inds2,tmp.inds]<-tmp.VV
          RR.holder[tmp.inds]<-RR.holder[tmp.inds]-tmp.det
        }
      }
      ZZ<-.mult.arr.mat(VV.holder,PZ.holder,ntraits,ntraits,nsim)
      if(has.obs.inf|has.des.inf){
        if(has.obs.inf&has.des.inf){
          tmp.inf<-(as.vector(obs.inf)|des.inf.holder)
        }else if(has.des.inf){
          tmp.inf<-des.inf.holder
        }else{
          tmp.inf<-obs.inf
        }
        VV.holder[diag.inds][tmp.inf]<-0
        PZ.holder[tmp.inf]<-0
      }
      
      list(VV.holder,
           ZZ,
           (if(double.RR) (.colSums(ZZ*PZ.holder,ntraits,nsim)+RR.holder) else (.colSums(ZZ*PZ.holder,ntraits,nsim)+RR.holder)/2)+
             (if(has.des[e]) .sum.ls(RR[des[[e]]],ndes[e]) else 0),
           cors[[2]])
    }
    full.has.obs<-NULL
    full.has.x<-NULL
  }else{
    #make simpler univariate versions...
    sum.obs<-matrix(unlist(lapply(parsed.obs,sum),use.names=FALSE),nedges,ntraitID)
    sum.sq.obs<-matrix(unlist(lapply(parsed.obs,function(ii) sum(ii^2)),use.names=FALSE),nedges,ntraitID)
    full.has.obs<-has.obs[,traitID,drop=FALSE]
    full.has.x<-(des.codes!=non.code)[,traitID,drop=FALSE]
    .get.edgewise.xvar.mu<-function(ls,e,m,n,sim.inds,map.inds,map){
      out<-matrix(nrow=m,ncol=n)
      #could potentially pre-compute some of this stuff...
      for(i in seq_along(ls)){
        tmp<-map.inds==i
        if(any(tmp)){
          if(is.list(ls[[i]][[1]])){
            out[tmp,]<-unlist(ls[[i]][[1]][e,sim.inds],use.names=FALSE)[tmp]
          }else{
            out[tmp,]<-ls[[i]][[1]]
          }
        }
      }
      .colSums(map*out,m,n)
    }
    #could better simplify the observation code system in the case of 1 trait...
    des.inf.holder<-logical(nsim)
    calc.all<-function(e,xvar,yvar,cors,mu,VV,ZZ,RR,double.RR,
                       traitID.inds,treeID.inds,sims.per.traitID,sims.per.treeID,traitID.seq,treeID.seq,nsim,
                       full.has.obs,full.has.x,traitID,
                       RR.holder,PZ.holder,VV.holder,des.inf.holder){
      #need to do descendants first (good checkpoint for infinite precisions)
      if(has.des[e]){
        des.list<-rep(list(list(VV.holder,PZ.holder,des.inf.holder)),ndes[e])
        counter<-0
        for(d in des[[e]]){
          counter<-counter+1
          for(i in treeID.seq){
            if(treeID.des.codes[d,i]!=non.code){
              tmp.inds<-treeID.inds[[i]]
              tmp.map<-maps[[d,i]]
              tmp.states<-names(tmp.map)
              tmp.map<-as.vector(tmp.map)
              tmp.m<-map.lens[d,i]
              tmp.n<-sims.per.treeID[i]
              if(xvar.maps[d,i]){
                des.list[[counter]][[1]][tmp.inds]<-.get.edgewise.xvar.mu(xvar,d,tmp.m,tmp.n,tmp.inds,inds[["Xvar"]][tmp.states],tmp.map)
              }else{
                des.list[[counter]][[1]][tmp.inds]<-sum(tmp.map*unlist(xvar[inds[["Xvar"]][tmp.states]],use.names=FALSE))
              }
              if(mu.maps[d,i]){
                des.list[[counter]][[2]][tmp.inds]<-.get.edgewise.xvar.mu(mu,d,tmp.m,tmp.n,tmp.inds,inds[["mu"]][tmp.states],tmp.map)
              }else{
                des.list[[counter]][[2]][tmp.inds]<-sum(tmp.map*unlist(mu[inds[["mu"]][tmp.states]],use.names=FALSE))
              }
            }
          }
          tmp.inds<-full.has.x[d,]
          des.list[[counter]][[1]][tmp.inds]<-VV[[d]][tmp.inds]+des.list[[counter]][[1]][tmp.inds]
          des.list[[counter]][[2]][tmp.inds]<-ZZ[[d]][tmp.inds]-des.list[[counter]][[2]][tmp.inds]
          des.list[[counter]][[3]][tmp.inds]<-des.list[[counter]][[1]][tmp.inds]==0
          des.inf.holder[tmp.inds]<-des.inf.holder[tmp.inds]|des.list[[counter]][[3]][tmp.inds]
        }
        has.des.inf<-any(des.inf.holder)
      }else{
        has.des.inf<-FALSE
      }
      has.obs.inf<-FALSE
      tmp.inds<-has.obs[e,]
      if(any(tmp.inds)){
        tmp.nobs<-nobs[e,]
        tmp.yvar<-yvar[[inds[["Yvar"]][e]]][[1]]
        if(length(tmp.yvar)>1){
          single.yvar<-FALSE
          tmp.inds2<-full.has.obs[e,]
          tmp.inds3<-traitID[tmp.inds2]
          tmp.yvar<-tmp.yvar[e,tmp.inds2]
          tmp.nobs<-tmp.nobs[tmp.inds3]
        }else{
          single.yvar<-TRUE
        }
        log.yvar<-log(tmp.yvar)
        obs.inf<-is.infinite(tmp.yvar)
        if(has.des.inf){
          if(single.yvar){
            tmp.inds2<-full.has.obs[e,]
            tmp.inds3<-traitID[tmp.inds2]
            tmp.nobs<-tmp.nobs[tmp.inds3]
            tmp.n<-sum(tmp.inds2)
            tmp.yvar<-rep(tmp.yvar,length.out=tmp.n)
            log.yvar<-rep(log.yvar,length.out=tmp.n)
            obs.inf<-rep(obs.inf,length.out=tmp.n)
            single.yvar<-FALSE
          }
          tmp.yvar[des.inf.holder]<-0
          log.yvar[des.inf.holder]<-0
        }
        if(any(obs.inf)){
          has.obs.inf<-TRUE
          tmp.yvar[obs.inf]<-1
          log.yvar[obs.inf]<-0
        }
        if(single.yvar){
          VV.holder<-(tmp.nobs*tmp.yvar)[traitID]
          PZ.holder<-(sum.obs[e,]*tmp.yvar)[traitID]
          RR.holder<-(tmp.nobs*log.yvar-sum.sq.obs[e,]*tmp.yvar-tmp.nobs*log(2*pi))[traitID]
          if(has.obs.inf|has.des.inf){
            RR.holder[obs.inf|des.inf.holder]<-0
          }
        }else{
          VV.holder[tmp.inds2]<-tmp.nobs*tmp.yvar
          PZ.holder[tmp.inds2]<-sum.obs[e,tmp.inds3]*tmp.yvar
          RR.holder[tmp.inds2]<-tmp.nobs*log.yvar-sum.sq.obs[e,tmp.inds3]*tmp.yvar-tmp.nobs*log(2*pi)
          if(has.obs.inf|has.des.inf){
            RR.holder[obs.inf|des.inf.holder]<-0
          }
        }
      }
      if(has.des[e]){
        counter<-0
        for(d in des[[e]]){
          counter<-counter+1
          tmp.inds<-full.has.x[d,]
          if(any(tmp.inds)){
            if(has.des.inf){
              foc.inf<-des.list[[counter]][[3]][tmp.inds]
              has.foc.inf<-any(foc.inf)
            }else{
              has.foc.inf<-FALSE
            }
            if(has.des.inf|has.obs.inf){
              if(has.des.inf&has.obs.inf){
                if(single.yvar){
                  oth.inf<-des.inf.holder[tmp.inds]
                }else{
                  oth.inf<-des.inf.holder[tmp.inds]|obs.inf[tmp.inds]
                }
              }else if(has.des.inf){
                oth.inf<-des.inf.holder[tmp.inds]
              }else{
                if(single.yvar){
                  oth.inf<-TRUE
                }else{
                  oth.inf<-obs.inf[tmp.inds]
                }
              }
              has.oth.inf<-any(oth.inf)
            }else{
              has.oth.inf<-FALSE
            }
            tmp.VV<-1/des.list[[counter]][[1]][tmp.inds]
            if(has.oth.inf) tmp.VV[oth.inf]<-0
            if(has.foc.inf) tmp.VV[foc.inf]<-1
            tmp.ZZ<-des.list[[counter]][[2]][tmp.inds]
            tmp.PZ<-tmp.VV*tmp.ZZ
            VV.holder[tmp.inds]<-VV.holder[tmp.inds]+tmp.VV
            PZ.holder[tmp.inds]<-PZ.holder[tmp.inds]+tmp.PZ
            if(has.oth.inf) tmp.VV[oth.inf]<-1
            if(has.foc.inf) tmp.PZ[foc.inf]<-0
            RR.holder[tmp.inds]<-RR.holder[tmp.inds]+log(tmp.VV)-tmp.ZZ*tmp.PZ
            if(has.foc.inf|has.oth.inf){
              if(has.foc.inf&has.oth.inf){
                tmp.inf<-!(foc.inf|oth.inf)
              }else if(has.foc.inf){
                tmp.inf<-!foc.inf
              }else{
                tmp.inf<-!oth.inf
              }
              RR.holder[tmp.inds][tmp.inf]<-RR.holder[tmp.inds][tmp.inf]-log(2*pi)
            }else{
              RR.holder[tmp.inds]<-RR.holder[tmp.inds]-log(2*pi)
            }
          }
        }
      }
      tmp.inds<-full.has.x[e,]
      VV.holder[tmp.inds]<-1/VV.holder[tmp.inds]
      ZZ<-VV.holder*PZ.holder
      if(has.des.inf|has.obs.inf){
        if(has.des.inf&has.obs.inf){
          if(single.yvar){
            tmp.inf<-des.inf.holder[tmp.inds]
          }else{
            tmp.inf<-des.inf.holder[tmp.inds]|obs.inf[tmp.inds]
          }
        }else if(has.des.inf){
          tmp.inf<-des.inf.holder[tmp.inds]
        }else{
          if(single.yvar){
            tmp.inf<-TRUE
          }else{
            tmp.inf<-obs.inf[tmp.inds]
          }
        }
        #check for invalid variances less than 1 (for exact measurements)! This should work...
        if(any(VV.holder[tmp.inds][tmp.inf]<1)) return(NULL)
        log.VV<-numeric(sum(tmp.inds))
        log.VV[!tmp.inf]<-log(VV.holder[tmp.inds][!tmp.inf])+log(2*pi)
        VV.holder[tmp.inds][tmp.inf]<-0
        PZ.holder[tmp.inds][tmp.inf]<-0
      }else{
        log.VV<-log(VV.holder[tmp.inds])+log(2*pi)
      }
      RR.holder[tmp.inds]<-RR.holder[tmp.inds]+ZZ[tmp.inds]*PZ.holder[tmp.inds]+log.VV
      
      list(VV.holder,
           ZZ,
           (if(double.RR) RR.holder else RR.holder/2)+
             (if(has.des[e]) .sum.ls(RR[des[[e]]],ndes[e]) else 0))
    }
  }
  
  ####FINAL FUNCTION TO OUTPUT####
  
  int.out<-function(par,tree.inds=NULL){
    #calculate tree indices
    if(!is.null(tree.inds)&conmaps.flag){
      tree.flag<-TRUE
      
      #correct mapping of contsimmapped variables
      nsim<-length(tree.inds)
      wgts<-wgts[tree.inds]
      conmaps[["RES_trees"]]<-conmaps[["RES_trees"]]%in%tree.inds
      conmaps[["RES_nodes_trees"]]<-conmaps[["RES_nodes_trees"]]%in%tree.inds
      conmaps[["RES_splits"]]<-conmaps[["RES_splits"]][conmaps[["RES_trees"]]]
      conmaps[["RES_splits"]]<-factor(conmaps[["RES_splits"]],levels=which(tabulate(as.numeric(conmaps[["RES_splits"]]))>0))
      if(any(unlist(conmaps[["RES_mapped_mu"]],use.names=FALSE))){
        conmaps[["RES_mu_holder"]]<-conmaps[["RES_mu_holder"]][conmaps[["RES_trees"]]]
      }
      if(any(unlist(conmaps[["RES_mapped_yvar"]],use.names=FALSE))){
        conmaps[["RES_yvar_holder"]]<-matrix(conmaps[["RES_yvar_holder"]][conmaps[["RES_nodes_trees"]]],nedges,nsim)
      }
      if(any(unlist(conmaps[["RES_mapped_xvar"]],use.names=FALSE))){
        conmaps[["RES_xvar_holder"]]<-conmaps[["RES_xvar_holder"]][conmaps[["RES_trees"]]]
      }
      
      #correct indexing of tree/traitIDs
      #better idea--just change traitID.seq and treeID.seq--I think this should work, actually!
      #just need to update the indices and nsims per accordingly...
      traitID.inds<-lapply(traitID.inds,'[',tree.inds)
      treeID.inds<-lapply(treeID.inds,'[',tree.inds)
      sims.per.traitID<-unlist(lapply(treeID.inds,sum),use.names=FALSE)
      sims.per.treeID<-unlist(lapply(treeID.inds,sum),use.names=FALSE)
      traitID.seq<-which(sims.per.traitID>0)
      treeID.seq<-which(sims.per.treeID>0)
      
      if(!mult.traits){
        full.has.obs<-full.has.obs[,tree.inds,drop=FALSE]
        full.has.x<-full.has.x[,tree.inds,drop=FALSE]
        traitID<-traitID[tree.inds]
      }
      
      #correct sizes of containers
      RR.holder<-numeric(nsim)
      RR<-rep(list(RR.holder),nedges)
      if(mult.traits){
        PZ.holder<-matrix(0,ntraits,nsim)
        ZZ<-rep(list(PZ.holder),nedges)
        VV.holder<-array(0,c(ntraits,ntraits,nsim))
        VV<-rep(list(VV.holder),nedges)
        des.inf.holder<-matrix(FALSE,ntraits,nsim)
      }else{
        VV.holder<-PZ.holder<-RR.holder
        VV<-ZZ<-RR
        des.inf.holder<-logical(nsim)
      }
      
    }else{
      tree.flag<-FALSE
    }
    env<-c(as.list(setNames(par,par.nms)),conmaps)
    xvar<-calc.xvars(env,tree.flag,nsim)
    if(is.null(xvar)) return(-Inf)
    yvar<-calc.yvars(env,tree.flag)
    if(is.null(yvar)) return(-Inf)
    if(mult.traits){
      cors<-calc.cors(env)
      if(is.null(cors)) return(-Inf)
    }else{
      cors<-NULL
    }
    mu<-calc.mus(env,tree.flag,nsim)
    for(e in prune.seq){
      tmp<-calc.all(e,xvar,yvar,cors,mu,VV,ZZ,RR,root.nuisance.prior&e==1,
                    traitID.inds,treeID.inds,sims.per.traitID,sims.per.treeID,traitID.seq,treeID.seq,nsim,
                    full.has.obs,full.has.x,traitID,
                    RR.holder,PZ.holder,VV.holder,des.inf.holder)
      if(is.null(tmp)) return(-Inf)
      VV[[e]]<-tmp[[1]]
      ZZ[[e]]<-tmp[[2]]
      RR[[e]]<-tmp[[3]]
      if(mult.traits) cors[[2]]<-tmp[[4]]
    }
    LL<-wgts+RR[[1]]
    
    #is there something you can do here to make numerical gradient calculation better???
    #seems difficult with these more exotic priors...
    
    #probably could make this a bit more elegant and less hacky...
    # if(LL.cap){
    #   cutoff<-quantile(LL,prob=cap.quant)
    #   if(cutoff>min(LL)){
    #     LL.inds<-LL<cutoff
    #     LL<-LL[LL.inds]
    #     LL<-LL+log(nsim-sum(!LL.inds))
    #   }else{
    #     LL.inds<-rep(TRUE,nsim)
    #     LL<-LL+log(nsim)
    #   }
    # }else{
    #   LL<-LL+log(nsim)
    # }
    #an alt procedure would be to simply round down to cutoff...would this be better?
    #may have some advantages wrt to dummy models not being overly "stunted", yet...
    #makes the likelihood surface seem WAY more multimodal, which would be undesirable here...
    
    #now allowing both options (trimming/winsorization) in two-tailed manner
    #trimming = ignore outliers
    #winsorization = censor outliers by rounding to cutoff
    #(trimming extremely small outliers should help with dummy model competitiveness, right?)
    if(outlier.cutoff){
      cutoffs<-quantile(LL,prob=quant.cutoffs)
      if(winsorize){
        LL[LL<cutoffs[1]]<-cutoffs[1]
        LL[LL>cutoffs[2]]<-cutoffs[2]
      }else{
        LL.inds<-LL<cutoffs[1]|LL>cutoffs[2]
        n.otl<-sum(LL.inds)
        if(n.otl>0&n.otl<nsim){
          LL<-LL[LL.inds]
          nsim<-nsim-n.otl
          LL.ind.flag<-TRUE
        }
      }
    }
    LL<-LL+log(nsim)
    
    if(tree.balance.prior){
      max.LL<-max(LL)
      tmp<-log(sum(exp(LL-max.LL)))+max.LL-LL
      LL<-LL+log(tmp)-log(sum(tmp))
    }else{
      #should be numerically stable now
      tmp<-(tree.prior.exp-1)*LL
      max.LL<-max(tmp)
      LL<-tree.prior.exp*LL-log(sum(exp(tmp-max.LL)))-max.LL
    }
    
    #one last thing to try...
    #you can try "capping" the top log(nsim) weights
    #this is common in importance sampling
    #could allow you to do the nuisance prior while still "balancing" likelihoods...
    #but would that be any better than the flat prior? feels like it wouldn't...
    #but still might be worth a try...
    
    max.LL<-max(LL)
    LL<-exp(LL-max.LL)
    out.lik<-log(sum(LL))+max.LL
    
    if(LL.ind.flag){
      old.LL<-LL
      LL<-RR.holder
      LL[LL.inds]<-old.LL
    }
    
    if(is.na(out.lik)){
      -Inf
    }else{
      attr(out.lik,"LL")<-LL
      out.lik
    }
  }
  
  npar<-length(par.nms)
  gg<-numeric(npar)
  out<-function(par,
                grad=FALSE,grad.qual=0.9,grad.step=(.Machine$double.eps)^(1/3),
                invert=FALSE,
                return.LL=FALSE){
    if(length(par)==npar){
      lik<-tryCatch(int.out(par),error=function(e) -Inf)
      if(invert) lik<- -lik
      if(return.LL){
        if(is.infinite(lik)){
          rep(lik,nsim)
        }else{
          if(invert) -attr(lik,"LL") else attr(lik,"LL")
        }
      }else{
        if(grad){
          if(is.infinite(lik)){
            #probably a terrible idea...but oh well :/
            gg<-rnorm(npar)
          }else{
            if(grad.qual<1){
              LL<-attr(lik,"LL")
              LL<-LL/sum(LL)
              ord<-order(LL,decreasing=TRUE)
              #trees that contribute to 90% of likelihood!
              #might even be able to get away with 50%...
              #1.5 is just a random guess as to an appropriate sample size in this situation...
              #2 seems safer
              tmp.size<-min(sum(LL>0),round(grad.qual*nsim),2*min(which(cumsum(LL[ord])>grad.qual)))
              inds<-sort.int(sample.int(nsim,tmp.size,prob=LL))
            }else{
              inds<-NULL
            }
            cur<-tryCatch(int.out(par,inds),error=function(e) -Inf)
            #random gradients upon errors didn't seem awful...so I think I'll just do that from now on
            if(is.infinite(cur)){
              gg<-c(-0.1,0.1)[rbinom(npar,1,0.5)+1]
            }else{
              hh<-rep(grad.step,npar)
              tmp.inds<-abs(par)>1
              hh[tmp.inds]<-hh[tmp.inds]*abs(par[tmp.inds])
              pp<-par
              for(i in seq_along(par)){
                pp[i]<-par[i]+hh[i]/2
                fw<-tryCatch(int.out(pp,inds),error=function(e) -Inf)
                fw.inf<-is.infinite(fw)
                pp[i]<-par[i]-hh[i]/2
                bw<-tryCatch(int.out(pp,inds),error=function(e) -Inf)
                bw.inf<-is.infinite(bw)
                if(fw.inf&bw.inf){
                  gg[i]<-c(-0.1,0.1)[rbinom(1,1,0.5)+1]
                }else{
                  if(fw.inf){
                    gg[i]<-2*(fw-cur)/hh[i]
                  }else if(bw.inf){
                    gg[i]<-2*(cur-bw)/hh[i]
                  }else{
                    gg[i]<-(fw-bw)/hh[i]
                  }
                }
                pp[i]<-par[i]
              }
            }
            if(invert) gg<- -gg
          }
          list("objective"=lik[1],"gradient"=gg)
        }else{
          lik[1]
        }
      }
    }else{
      stop("Wrong number of parameters: double-check input!")
    }
  }
  attr(out,"par.nms")<-par.nms
  attr(out,"call")<-list("formulae"=setNames(list(Xvar,Xcor,Mu,Yvar,Ycor),
                                             c("Xvar","Xcor","mu","Yvar","Ycor")),
                         "trait.data"=trait.info[["trait.data"]],
                         "tree"=tree,
                         "treeID"=treeID,
                         "traitID"=traitID,
                         "conmaps"=conmaps)
  
  out
}

.fix.bds<-function(bds,upper,par.nms){
  if(is.list(bds)) bds<-unlist(bds)
  out<-bds[par.nms]
  nms<-names(bds)
  probs<-is.na(out)
  if(any(probs)){
    if(is.null(nms)){
      n.bds<-length(bds)
      nms<-character(n.bds)
      prob.nms<-!logical(n.bds)
    }else{
      nms[is.na(nms)]<-""
      prob.nms<-!(nms%in%par.nms)
    }
    if(any(prob.nms)){
      out[probs]<-rep(bds[prob.nms],length.out=sum(probs))
    }else if(upper){
      out[probs]<-Inf
    }else{
      out[probs]<- -Inf
    }
  }
  if(upper){
    out[is.infinite(out)&out<0]<-Inf
  }else{
    out[is.infinite(out)&out>0]<- -Inf
  }
  out
}

.get.init.fxn<-function(lb,ub,npar,init.width){
  #now forces all parameters to initialize around 0 as much as possible
  if(is.null(lb)) lb<-rep(-Inf,npar)
  if(is.null(ub)) ub<-rep(Inf,npar)
  lb.inf<-is.infinite(lb)
  ub.inf<-is.infinite(ub)
  lb[lb.inf]<-pmin(-init.width/2,ub[lb.inf]-init.width)
  ub[ub.inf]<-pmax(init.width/2,lb[ub.inf]+init.width)
  cents<-cbind((lb+ub)/2,lb+init.width/2,ub-init.width/2)
  whichs<-apply(abs(cents),1,which.min)
  cents<-cents[cbind(seq_len(npar),whichs)]
  lb<-pmax(lb,cents-init.width/2)
  ub<-pmin(ub,cents+init.width/2)
  function(inds){
    runif(length(inds),lb[inds],ub[inds])
  }
}

#allow for sequences of NLOPT calls by passing vectors to opts...
#' @export
find.mle<-function(lik.fun,init=NULL,times=1,lb=NULL,ub=NULL,...,
                   recycle.init=FALSE,init.width=10,
                   on.fail=c("random.restart","NA"),max.tries=100,
                   grad.qual=0.9,grad.step=.Machine$double.eps^(1/3),
                   verbose=FALSE){
  par.nms<-attr(lik.fun,"par.nms")
  npar<-length(par.nms)
  if(!is.null(lb)){
    lb<-.fix.bds(lb,upper=FALSE,par.nms)
  }
  if(!is.null(ub)){
    ub<-.fix.bds(ub,upper=TRUE,par.nms)
  }
  #swap any bounds needing it
  if(!is.null(lb)&!is.null(ub)){
    tmp.lb<-lb
    tmp.ub<-ub
    tmp<-tmp.lb>tmp.ub
    lb[tmp]<-tmp.ub[tmp]
    tmp<-tmp.ub<tmp.lb
    ub[tmp]<-tmp.lb[tmp]
  }
  out.init<-matrix(nrow=npar,ncol=times,
                   dimnames=list(par.nms,NULL))
  if(!is.null(init)){
    if(!is.list(init)){
      ndims<-length(dim(init))
      if(ndims<2){
        if(is.null(names(init))){
          init<-list(init)
        }else{
          init<-split(init,names(init))
        }
      }else if(ndims>2){
        stop("Please input init as either a vector, matrix (with columns corresponding to different parameters), or list")
      }else{
        init<-asplit(init,2)
      }
    }
    nms<-names(init)
    if(is.null(nms)){
      n.init<-length(init)
      nms<-character(n.init)
      prob.nms<-!logical(n.init)
    }else{
      nms[is.na(nms)]<-""
      prob.nms<-!(nms%in%par.nms)
    }
    for(i in nms[!prob.nms]){
      if(recycle.init){
        out.init[i,]<-rep(init[[i]],length.out=times)
      }else{
        tmp.seq<-seq_len(min(times,length(init[[i]])))
        out.init[i,tmp.seq]<-init[[i]][tmp.seq]
      }
    }
    probs<-which(!(par.nms%in%nms[!prob.nms]))
    if(any(probs)){
      counter<-0
      if(recycle.init){
        for(i in rep(which(prob.nms),length.out=length(probs))){
          counter<-counter+1
          out.init[probs[counter],]<-rep(init[[i]],length.out=times)
        }
      }else{
        for(i in which(prob.nms)[seq_len(min(sum(prob.nms),length(probs)))]){
          counter<-counter+1
          tmp.seq<-seq_len(min(times,length(init[[i]])))
          out.init[probs[counter],tmp.seq]<-init[[i]][tmp.seq]
        }
      }
    }
  }
  init<-out.init
  nas<-is.na(init)
  inds<-(which(nas)-1)%%npar+1
  init.fxn<-.get.init.fxn(lb,ub,npar,init.width)
  init[nas]<-init.fxn(inds)
  opts<-list(...)
  if(is.null(opts[["algorithm"]])){
    opts[["algorithm"]]<-c(if(npar>2) "NLOPT_LD_TNEWTON_PRECOND_RESTART" else NULL,
                           "NLOPT_LN_SBPLX",
                           if(npar>2) "NLOPT_LN_PRAXIS" else NULL)
  }
  if(is.null(opts[["maxeval"]])){
    opts[["maxeval"]]<-c(1e4,1e3)[as.numeric(grepl("LD",opts[["algorithm"]]))+1]
  }
  if(is.null(opts[["ftol_rel"]])){
    opts[["ftol_rel"]]<-sqrt(.Machine$double.eps)
  }
  if(is.null(opts[["xtol_res"]])){
    opts[["ftol_rel"]]<-sqrt(.Machine$double.eps)
  }
  if(length(on.fail)==1){
    if(is.na(on.fail)) on.fail<-"NA"
  }
  if(is.character(on.fail)) on.fail<-pmatch(on.fail[1],c("random.restart","NA"))
  if(is.na(on.fail)|!is.numeric(on.fail)) on.fail<-1
  if(verbose&is.null(opts[["print_level"]])){
    opts[["print_level"]]<-3
  }
  runs.per.time<-max(lengths(opts))
  opts<-lapply(opts,function(ii) rep(ii,length.out=runs.per.time))
  pars<-matrix(nrow=npar,ncol=times)
  codes<-liks<-numeric(times)
  for(i in seq_len(times)){
    if(verbose) cat("Maximizing likelihood function... (rep ",i," out of ",times,"):\n\n",sep="")
    for(j in seq_len(runs.per.time)){
      if(verbose) cat("Optimization round ",j," out of ",runs.per.time,"... (rep ",i," out of ",times,"):\n\n",sep="")
      tmp.opts<-lapply(opts,'[[',j)
      if(grepl("_LD_",tmp.opts[["algorithm"]])){
        grad<-TRUE
      }else{
        grad<-FALSE
      }
      res<-nloptr::nloptr(init[,i],lik.fun,lb=lb,ub=ub,opts=tmp.opts,
                          grad=grad,grad.qual=grad.qual,grad.step=grad.step,invert=TRUE,return.LL=FALSE)
      #check for early terminations
      if((is.infinite(res[["objective"]])|res[["status"]]<0)&grepl("LBFGS|TNEWTON",tmp.opts[["algorithm"]])&res[["iterations"]]<20){
        break
      }
      init[,i]<-res[["solution"]]
    }
    if(on.fail==1){
      counter<-1
      while(counter<=max.tries&(is.infinite(res[["objective"]])|res[["status"]]<0)){
        if(verbose) cat("Random restart ",counter," out of ",max.tries,"... (rep ",i," out of ",times,"):\n\n",sep="")
        counter<-counter+1
        init[,i]<-init.fxn(seq_len(npar))
        for(j in seq_len(runs.per.time)){
          if(verbose) cat("Optimization round ",j," out of ",runs.per.time,"... (rep ",i," out of ",times,"):\n\n",sep="")
          tmp.opts<-lapply(opts,'[[',j)
          if(grepl("_LD_",tmp.opts[["algorithm"]])){
            grad<-TRUE
          }else{
            grad<-FALSE
          }
          res<-nloptr::nloptr(init[,i],lik.fun,lb=lb,ub=ub,opts=tmp.opts,
                              grad=grad,grad.qual=grad.qual,grad.step=grad.step,invert=TRUE,return.LL=FALSE)
          if((is.infinite(res[["objective"]])|res[["status"]]<0)&grepl("LBFGS|TNEWTON",tmp.opts[["algorithm"]])&res[["iterations"]]<20){
            break
          }
          init[,i]<-res[["solution"]]
        }
      }
    }
    liks[i]<- -res[["objective"]]
    pars[,i]<-res[["solution"]]
    codes[i]<-res[["status"]]
  }
  tmp.inds<-which(!(codes<0))
  #should also check for max iteration warnings
  #should rounding errors (-4) also be reported separately?
  if(!length(tmp.inds)){
    warning("Optimization algorithm failed to properly converge")
    final.max<-which.max(liks)
  }else{
    final.max<-tmp.inds[which.max(liks[tmp.inds])]
  }
  list("lnLik"=liks[final.max],
       "estimates"=setNames(pars[,final.max],par.nms),
       "return_code"=codes[final.max])
  #eventually make a function to format output nicely
}
