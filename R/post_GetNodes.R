#notes: internal.only only matters if nodes is NULL, in which case all available (internal if TRUE) nodes are selected
#       partial only affects matching of character selections to tips, since node labels are numeric in contsimmap as of 2/28/23

#' @export
get.nodes<-function(contsimmap,nodes=NULL,traits=NULL,partial=FALSE,internal.only=TRUE){
  if(!is.null(traits)){
    contsimmap<-contsimmap[,traits,]
  }
  edges<-dimnames(contsimmap)[[1]]
  roots<-edges=='N0'
  tmp<-as.numeric(gsub('N','',edges))
  tmp[roots]<-NA
  avail.nodes<-attr(contsimmap,'tree')[[1]][['edge']][tmp,2]
  avail.nodes[roots]<-Ntip(contsimmap)+1
  w.flag<-TRUE
  if(is.null(nodes)){
    nodes<-(Ntip(contsimmap)+1):(Ntip(contsimmap)+Nnode(contsimmap))
    if(!internal.only) nodes<-c(seq_len(Ntip(contsimmap)),nodes)
    w.flag<-FALSE
  }
  if(is.numeric(nodes)){
    matches<-match(nodes,avail.nodes)
    tips<-nodes<=Ntip(contsimmap)
    nodes[tips]<-attr(contsimmap,'tree')[[1]][['tip.label']][nodes[tips]]
  }else if(is.character(nodes)){
    if(partial){
      strings<-is.na(as.numeric(nodes))
      nodes<-as.list(nodes)
      tmp<-attr(contsimmap,'tree')[[1]][['tip.label']]
      nodes[strings]<-lapply(nodes[strings],grep,x=tmp,value=TRUE)
      nodes<-unlist(nodes,use.names=FALSE)
    }
    tips<-avail.nodes<=Ntip(contsimmap)
    avail.nodes[tips]<-attr(contsimmap,'tree')[[1]][['tip.label']][avail.nodes[tips]]
    matches<-match(nodes,avail.nodes)
  }
  probs<-is.na(matches)
  if(any(probs)&w.flag){
    warning("placeholder warning about nodes specified that either don't exist in tree or contsimmap")
  }
  matches<-matches[!probs]
  nodes<-nodes[!probs]
  lens<-rep(.get.ns(contsimmap,'nts',uncompress=TRUE)[matches,,drop=FALSE],dim(contsimmap)[2])
  tmp<-aperm(unclass(contsimmap)[matches,,,drop=FALSE],c(1,3,2))
  aperm(array(unlist(lapply(seq_along(tmp),function(ii) tmp[[ii]][lens[ii]]),use.names=FALSE),
              c(length(matches),dim(tmp)[-1]),c(node=list(nodes),dimnames(tmp)[-1])),
        c(1,3,2))
}

#' @export
get.tips<-function(contsimmap,tips=NULL,traits=NULL,partial=TRUE){
  if(is.null(tips)){
    tips<-tip.inds(contsimmap,unique.only=TRUE)
    tips<-sort(tips[!is.na(tips)])
  }
  get.nodes(contsimmap,nodes=tips,traits=traits,partial=partial)
}