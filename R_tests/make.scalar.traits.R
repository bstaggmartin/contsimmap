#more efficient alternative to traits2rates in the case that some other trait is only affecting the rate at which
#a new trait evolves
#may be nice to make this eventually store the extra parameters created
#as it stands right now, I think this might be messed up by reordering simulations?
#I think this is the way to go potentially--can do anything traits2rates does if Xcor doesn't vary as well
#' @export
make.scalar.trait<-function(contsimmap,...,
                            Xsig2=NULL,
                            X0=0){
  if(nedge(contsimmap)<Nedge(contsimmap)){
    if(length(contsimmap[['nodes']])==1){
      stop('The make.scalar.trait() function will support single subtrees in the future, but not yet')
    }else{
      stop('The make.scalar.trait() function does not work with contsimmaps consisting of multiple subtrees--did you mean to subset the edges in your contsimmaps?')
    }
  }
  
  formulae<-list(...)
  ntraits<-length(formulae)
  traits<-unlist(lapply(formulae,function(ii) as.character(ii[[2]])),use.names=FALSE)
  nsims<-nsim(contsimmap)
  states<-contsimmap:::.get.states(contsimmap[['tree']])
  nstates<-length(states)
  if(is.null(Xsig2)){
    Xsig2<-diag(ntraits)
  }
  Xsig2<-contsimmap:::.fix.mat.list(Xsig2,ntraits,traits,nstates,states)
  traits<-rownames(Xsig2[[1]])
  X0<-matrix(rep(X0,length.out=ntraits*nsims),ntraits,nsims)
  tmp.env<-as.list(setNames(paste0(traits,'_rate'),traits))
  tmp.env<-lapply(tmp.env,str2lang)
  tmp.formulae<-lapply(formulae,function(ii) do.call(substitute,list(expr=ii,env=tmp.env)))
  contsimmap<-do.call(make.traits,c(list(contsimmap),tmp.formulae))
  
  tmp.ntraits<-ntrait(contsimmap)
  treeID<-contsimmap[['perm.inds']][,1,drop=FALSE]
  treeID<-as.integer(factor(treeID,unique(treeID)))
  sims.per.tree<-tabulate(treeID)
  nts<-do.call(cbind,lapply(contsimmap[['x']],function(ii) lengths(ii)/tmp.ntraits/sims.per.tree))
  seeds.per.tree.edge<-nts*ntraits
  seeds.per.edge<-colSums(seeds.per.tree.edge)
  seed<-contsimmap:::.get.seed(seeds.per.edge,ntraits)
  anc<-anc.edges(contsimmap[['tree']][[1]])
  anc[lengths(anc)==0]<-0
  anc<-as.character(unlist(anc))
  names(anc)<-seq_along(anc)
  cladewise.ord<-contsimmap:::.get.cladewise.ord(anc)
  contsimmap<-subset(contsimmap,traits=c(seq_len(tmp.ntraits),rep(NA,ntraits)))
  contsimmap[['traits']][seq_len(ntraits)+tmp.ntraits]<-traits
  rate.inds<-match(as.character(tmp.env),contsimmap[['traits']])
  scalar.inds<-match(traits,contsimmap[['traits']])
  contsimmap[['x']]<-contsimmap:::.accumulate.seed(seed,
                                                   anc,cladewise.ord,contsimmap$maps,
                                                   Xsig2,X0,
                                                   ntraits,traits,nstates,states,
                                                   treeID,nts,sims.per.tree,seeds.per.tree.edge,
                                                   scalar.trait=TRUE,original=contsimmap[['x']],rate.inds=rate.inds,scalar.inds=scalar.inds)
  contsimmap[['nodes']][['0']][scalar.inds,]<-X0
  contsimmap
}