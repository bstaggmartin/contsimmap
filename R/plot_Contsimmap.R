.capitalize<-function(x){
  tmp<-match(substr(x,0,1),letters)
  if(!is.na(tmp)){
    paste0(LETTERS[tmp],substr(x,2,nchar(x)))
  }else{
    x
  }
}

.parse.bys<-function(choices,
                     Col.by,Layer.by,Alpha.by,Mix.by,Wgt.by,Lty.by,Lwd.by,
                     one.trait){
  #determine bys
  bys<-list(col=Col.by,layer=Layer.by,alpha=Alpha.by,mix=Mix.by,wgt=Wgt.by,lty=Lty.by,lwd=Lwd.by)
  trait.flag<-setNames(rep(FALSE,7),names(bys))
  foo<-function(i){
    switch(i,
           col=match('simulation',choices),
           layer=if(one.trait) match('simulation',choices) else match('time',choices),
           alpha=bys[['col']],
           mix=bys[['col']],
           wgt=bys[['mix']],
           lty=bys[['col']],
           lwd=bys[['lty']])
  }
  for(i in names(bys)){
    if(is.null(bys[[i]])){
      bys[[i]]<-foo(i)
    }else{
      bys[[i]]<-pmatch(bys[[i]][1],choices)
      if(is.na(bys[[i]])){
        bys[[i]]<-foo(i)
        warning(paste0(.capitalize(i),'.by'),' should be one of the following:',.report.names(nms=choices,combine='or'),
                ': it was set to ',choices[bys[[i]]],' by default')
      }
    }
    if(bys[[i]]>4){
      trait.flag[i]<-TRUE
    }
  }
  bys<-unlist(bys)
  bys<-setNames(choices[bys],names(bys))
  attr(bys,'trait.flag')<-trait.flag
  bys
}

# .get.priority<-function(x,tmp.by){
#   
# }


.parse.args<-function(args.ls,extractions,bys,traits,one.trait,sims,states,edges){
  defaults<-list('col'=palette(),
                 'layer'=integer(0),
                 'alpha'=NA,
                 'mix'=NA,
                 'wgt'=0.5,
                 'lty'=par('lty'),
                 'lwd'=par('lwd'),
                 'xlab'=if(one.trait) 'time' else traits[1],
                 'ylab'=traits[length(traits)])
  for(i in names(defaults)){
    if(is.null(args.ls[[i]])){
      args.ls[[i]]<-defaults[[i]]
    }
  }
  #should I put nodes together in a single array first?
  if(is.null(args.ls[['xlim']])){
    args.ls[['xlim']]<-range(extractions[[if(one.trait) 'time' else traits[1]]],na.rm=TRUE)
  }
  if(is.null(args.ls[['ylim']])){
    args.ls[['ylim']]<-range(extractions[[traits[length(traits)]]],na.rm=TRUE)
  }
  #form final argument list/recycling
  inds<-list()
  len<-length(extractions[[1]])
  #need to figure out how I want layering to work though...
  resolve.layering<-function(x,nms,breaks=NULL){
    if(is.logical(x)){
      x<-nms[x]
    }
    if(is.character(x)|is.factor(x)){
      x<-match(x,nms)
    }else if(!is.null(breaks)&is.numeric(x)&!is.integer(x)){
      x<-findInterval(x,breaks)+1
    }
    x<-unique(x[!is.na(x)])
    tmp<-vector('integer',length(nms))
    tmp[x]<-seq_along(x)
    #lower numbers/picked items = higher priority, therefore items later in default ordering (hence plotted last) should receive lower numbers!
    tmp[tmp==0]<-length(tmp):(length(x)+1)
    tmp
  }
  foo<-function(i,tmp.by){
    x<-args.ls[[i]]
    nms2match<-switch(tmp.by,
                      'simulation'=sims,
                      'state'=states,
                      'edge'=edges,
                      NULL)
    if(is.null(nms2match)){
      #prepping args to match up to continuous scale
      #NA values are ignored for now, though perhaps a better system could be envisioned?
      breaks<-paste0(tmp.by,'.breaks')
      nbreaks<-length(args.ls[[breaks]])+1
      if(i=='layer'){
        x<-resolve.layering(x,seq_len(nbreaks),args.ls[[breaks]])
      }else if(i=='col'){
        x<-colorRampPalette(x,alpha=TRUE)(nbreaks)
      }else{
        if(is.numeric(x)){
          x<-approx(x,n=nbreaks)[['y']]
        }else{
          x<-x[approx(seq_along(x),n=nbreaks)[['y']]]
        }
      }
    }else{
      #prepping args to match up to discrete factor
      if(i=='layer'){
        x<-resolve.layering(x,nms2match)
      }else{
        nms<-names(x)
        if(is.null(nms)){
          x<-setNames(rep(x,length.out=length(nms2match)),nms2match)
        }else{
          tmp<-x[nms2match]
          probs<-!(nms%in%nms2match)
          nas<-is.na(tmp)
          tmp[nas]<-if(any(probs)) x[probs] else NA
          x<-tmp
        }
      }
    }
    #check to see if we have states from summarized contsimmap!
    if(!is.null(attr(inds[[tmp.by]],'summarized_states'))){
      #color averaging
      if(i=="col"|i=="mix"){
        tmp<-col2rgb(x,alpha=TRUE)
        tmp<-lapply(seq_len(4),
                    function(ii) 
                      .rowSums(sweep(extractions[[tmp.by]][attr(inds[[tmp.by]],'notnas'),,drop=FALSE],2,tmp[ii,],'*'),attr(inds[[tmp.by]],'nnotnas'),length(states)))
        x<-rep(as.character(NA),len)
        x[attr(inds[[tmp.by]],'notnas')]<-do.call(rgb,c(tmp,maxColorValue=255))
      #"normal" averaging
      }else if(is.numeric(x)){
        tmp<-x
        x<-rep(as.numeric(NA),len)
        x[attr(inds[[tmp.by]],'notnas')]<-
          .rowSums(sweep(extractions[[tmp.by]][attr(inds[[tmp.by]],'notnas'),,drop=FALSE],2,tmp,'*'),attr(inds[[tmp.by]],'nnotnas'),length(states))
        if(i=="lty") x<-round(x)
      #just take maximum probs state for discrete, character-style arguments
      }else{
        x<-x[inds[[tmp.by]]]
      }
    }else{
      x<-x[inds[[tmp.by]]]
    }
    x
  }
  #potentially totally screwed up by NA indexing!!!
  #it will be difficult to generalize this to cases with NAs, but let's just assume for now this isn't a problem...
  nas<-which(is.na(extractions[[1]]))
  nasm1<-nas-1
  has.alpha<-TRUE
  has.mix<-TRUE
  for(i in names(bys)){
    if(i=='alpha'){
      if(all(is.na(args.ls[[i]]))){
        args.ls[[i]]<-NA
        has.alpha<-FALSE
      }
    }else if(i=='mix'){
      if(all(is.na(args.ls[[i]]))){
        args.ls[[i]]<-NA
        has.mix<-FALSE
      }
    }else if(i=='wgt'){
      if(all(!args.ls[[i]])){
        args.ls[[i]]<-0
        has.mix<-FALSE
      }
    }
    if(length(args.ls[[i]])>1|i=='layer'){
      tmp.by<-bys[[i]]
      #set up indices
      if(is.null(inds[[tmp.by]])){
        #check to see if we have states from a summarized contsimmap!
        if(length(dim(extractions[[tmp.by]]))>1){
          tmp<-apply(extractions[[tmp.by]],1,which.max)
          notnas<-lengths(tmp)>0
          tmp[!notnas]<-list(NA)
          tmp<-unlist(tmp,use.names=FALSE)
          attr(tmp,"summarized_states")<-TRUE
          attr(tmp,"notnas")<-notnas
          attr(tmp,"nnotnas")<-sum(notnas)
          inds[[tmp.by]]<-tmp
        }else if(any(tmp.by==c('simulation','state','edge'))){
          inds[[tmp.by]]<-match(extractions[[tmp.by]],
                                switch(tmp.by,simulation=sims,state=states,edge=edges))
        }else{
          breaks<-paste0(tmp.by,'.breaks')
          if(is.null(args.ls[[breaks]])){
            args.ls[[breaks]]<-100L
          }
          if(is.integer(args.ls[[breaks]])){
            lims<-range(extractions[[tmp.by]],na.rm=TRUE)
            args.ls[[breaks]]<-seq(lims[1],lims[2],length.out=args.ls[[breaks]]+2)[-c(1,args.ls[[breaks]]+2)]
          }
          args.ls[[breaks]]<-sort(args.ls[[breaks]])
          args.ls[[breaks]]<-args.ls[[breaks]][!is.na(args.ls[[breaks]])]
          tmp<-extractions[[tmp.by]]
          tmp[nas]<-tmp[nasm1]
          tmp[-len]<-(tmp[-len]+tmp[-1])/2
          tmp[nas]<-NA
          inds[[tmp.by]]<-findInterval(tmp,args.ls[[breaks]])+1
        }
      }
      #prep arguments
      args.ls[[i]]<-foo(i,tmp.by)
    }
  }
  if(has.alpha|has.mix){
    args.ls[['col']]<-alter.cols(args.ls[['col']],args.ls[['alpha']],args.ls[['mix']],args.ls[['wgt']])
  }
  args.ls[c('alpha','mix','wgt')]<-NULL
  foo<-function(x){
    if(length(x)>1){
      x[nas]<-x[nasm1]
    }
    x
  }
  args.ls[c('col','layer','lty','lwd')]<-lapply(args.ls[c('col','layer','lty','lwd')],foo)
  args.ls
}

#as it stands, will break if any traits are named simulation, state, edge, or time!
#NA indexing could severely break things
#breaks arguments are parsed by what's being broken (i.e., time/trait.breaks), which is a pretty good system, but then multiple args
##end up depending on the same breaks, which may not be desirable
#would be good to allow Layer.by to be set to multiple things so layering may have sort-within-sort style rules!
#might also be good to eventually add Pch.by and Bg.by for times when folks do type='o', etc.
##although the above would ideally require making a separate col, etc. for the points
##also a little weird since all splitting is based on edge parameters, rather than point parameters...
##thinking about it, allowing point parameters to reflect actual pointwise quantities would be an annoying amount of work...
##ideally would mean handling points as a separate thing manually, which seems annoying
##probably best to just make something akin to nodelabels()/tiplabels() that can be called afterwards
#also should probably add tip labelling capability eventually...
#automated legend making might also be a good idea down the line
#Decided to drop keys for now and see how it goes...but this function will surely break when you try to plot character-valued traits
#Should states be made traits in this case? Nah, because they vary per tree, not necessarily per sim...right?

#' Plot continuous stochastic character maps
#'
#' This function runs the default plotting method for continuous stochastic character maps (class "\code{contsimmap}").
#'
#'
#'
#' @param contsimmap An object of class "\code{contsimmap}".
#' @param traits A character/numeric vector specifying which traits within \code{contsimmap} to plot. Specifying a single trait results in a
#' "phenogram"-style plot with the x-axis corresponding to time; specifying multiple traits results in a "phylomorphospace"-style plot with the
#' first two traits corresponding to the x and y-axes, respectively. Additional traits are ignored for now, though a \code{pairs()} or 3D plotting
#' method may implemented in the future. Also, \code{NA} traits are \emph{not yet supported and could break the function spectacularly}.
#' 
#' Alternatively, phylogenies may be plotted in more traditional styles akin to 
#' the \code{plot.phylo()} function in \bold{ape} by setting 
#' \code{trait = "phylogram"} or \code{trait = "cladogram"}. In conjunction with the 
#' \code{Col.by}, \code{Mix.by}, etc. arguments, this allows one to 
#' annotate phylogenies with colors according to mapped continuous and/or 
#' discrete states, much like the \code{contMap()} or \code{plotSimmap()} functions 
#' from \bold{phytools}. Fan-style phylogenies may also be plotted by 
#' setting the \code{polarize} argument to \code{TRUE}.
#' @param sims A character/numeric vector specifying which simulations within \code{contsimmap} to plot. Set to \code{NULL} to plot all
#' simulations. \code{NA} simulations are not yet supported and could break the function spectacularly.
#' @param edges A character/numeric vector specifying which edges within \code{contsimmap} to plot. Set to \code{NULL} to plot all edges.
#' \code{NA} simulations are not yet supported and could break the function spectacularly.
#' @param Col.by,Layer.by,Alpha.by,Mix.by,Wgt.by,Lty.by,Lwd.by Any unambiguous abbreviation of "simulation", "state", "edge", "time", or trait
#' names within \code{contsimmap}. These arguments specify how their respective graphical arguments are allocated across the plot. For example,
#' setting \code{Col.by = "simulation"} results in plotting each simulation in a different color, while setting it to
#' \code{"edge"} or \code{"time"} instead causes edges or time slices to be colored differently. \code{Layer.by} controls how parts of the
#' plot are layered on top of one another; for example, setting \code{Layer.by = "time"} results in plotting younger parts of the phylogeny last.
#' \code{Alpha.by}, \code{Mix.by}, and \code{Wgt.by} control modifications of the base colors using \code{alter.cols()}; for example, setting
#' \code{Col.by = "simulation"} and \code{Alpha.by = "time"} would result in plotting each simulation in a different color but adjusting the
#' transparency of these colors according to time slice. Setting any of these to \code{NULL} or \code{NA} results in defaults: "simulation"
#' for \code{Col.by}, "simulation" for \code{Layer.by} if producing a "phenogram"-style plot and "time" otherwise, whatever \code{Col.by}
#' is for \code{Alpha.by}/\code{Mix.by}/\code{Lty.by}, whatever \code{Mix.by} is for \code{Wgt.by}, and whatever \code{Lty.by} is for
#' \code{Lwd.by}.
#' @param add \code{TRUE} or \code{FALSE}: should the plot be added to the existing plot window or initiate a new plotting window?
#' @param reverse.layers \code{TRUE} or \code{FALSE}: should the layering order specified by \code{Layer.by} (and the optional \code{layer}
#' argument passed to \code{...}) be reversed? For example, if \code{Layer.by = "time"} and \code{reverse.layering = TRUE}, the
#' \emph{older}, rather than younger, parts of the phylogeny will be plotted last.
#' @param polarize \code{TRUE} or \code{FALSE}: should the x and y coordinates be 
#' transformed to polar coordinates r and \eqn{\theta}, respectively? This is mainly 
#' for plotting phylogenies in conventional fan style by specifying \code{traits = "phylogram"} and 
#' \code{polarize = TRUE}, though maybe you want to "polarize" a phenogram/phylomorphospace for some reason? Seems kinda like a
#' neat idea, and I don't judge!
#' @param ang.min If \code{polarize = TRUE}: specifies the angle (in radians) that the minimum y coordinate gets converted to.
#' @param ang.max If \code{polarize = TRUE}: specifies the angle (in radians) that the maximum y coordinate gets converted to.
#' @param curviness If \code{traits = "cladogram"}: a positive number specifying the degree to which
#' edge lines are "biased" towards the y coordinate of ancestral nodes (if <1) versus 
#' descendant nodes (if >1). The default is 1, resulting in perfectly straight lines that don't curve.
#' @param ... Just about any base R graphical argument you can think of--it should theoretically all work! Notably, entries of \code{col},
#' \code{lty}, and \code{lwd} are recycled according to their respective \code{<xxx>.by} arguments:
#' \itemize{
#' \item{If "simulation", entries are assigned to simulations in the same order as given by \code{sims}. Entries may also be named to assign
#' them to specific simulations.}
#' \item{If "state", entries are assigned to discrete states in the \code{multiSimmap} object \code{contsimmap} is based on, in alphanumeric
#' order. Entries may also be named to assign them to specific states.}
#' \item{If "edge", entries are assigned to edges available in \code{contsimmap} in order of their numeric indices. Entries may also be named
#' to assign them to specific edges.}
#' \item{If "time" or a trait name, entries are interpolated to create a continuous gradient (as best as R can manage).
#' As a result, any \code{NA} entries are generally ignored. You can use \code{<xxx>.breaks} arguments to control how this gradient is
#' constructed (see below).}
#' }
#' There are also a few additional and potentially important arguments:
#' \itemize{
#' \item{\code{layer}, which is used to alter the layering order by specifying "high priority" elements. For example, setting \code{layer = 
#' "tree1_sim1"} ensures that simulation "tree1_sim1" is plotted last, assuming that \code{Layer.by = "simulation"}. \code{layer} can also be a
#' numeric or logical vector. If \code{Layer.by} is set to "time" or a trait name, one can even use this argument to specify particular time
#' slices or trait values that should be plotted last. Note that that this instead specifies "low priority" elements that are plotted \emph{first} if
#' \code{reverse.layering = TRUE}.}
#' \item{\code{alpha}, \code{mix}, and \code{wgt}, which are recycled as described above and correspond to the same arguments in
#' \code{alter.cols()}.}
#' \item{\code{<xxx>.breaks} (where \code{<xxx>} is either "time" or a trait name), which controls how continuous gradients
#' are constructed. This is similar to the \code{breaks} argument in the \code{image()} function, and basically ends up getting passed to the
#' \code{vec} argument in a call to \code{findInterval()}. Essentially, you provide thresholds where corresponding graphical arguments change.
#' This can also be set to an integer, in which case \code{<xxx>.breaks + 2} equally-spaced intervals are constructed across an appropriate range.
#' If not specified, these arguments are set to \code{100L} by default.}
#' }
#' 
#' 
#' 
#' @return Invisibly assigns plot information to "\code{last_plot.phylo}" in 
#' \code{.PlotPhyloEnv} environment, just like the \code{plot.phylo()} 
#' function from \bold{ape}. Most importantly, this allows the use of the 
#' \bold{ape} functions \code{nodelabels()}, \code{tiplabels()}, and 
#' \code{edgelabels()} for further plot annotation. Note that tip/node 
#' coordinates are averaged if multiple simulations are plotted at once.
#' 
#' 
#' 
#' @export
plot.contsimmap<-function(contsimmap,
                          traits=1,sims=sample(dim(contsimmap)[3],min(20,dim(contsimmap)[3])),edges=NULL,
                          Col.by=c('simulation','state','edge','time',dimnames(contsimmap)[[2]]),
                          Layer.by=NULL,Alpha.by=NULL,Mix.by=NULL,Wgt.by=NULL,Lty.by=NULL,Lwd.by=NULL,
                          add=FALSE,reverse.layers=FALSE,
                          polarize=FALSE,ang.min=0,ang.max=2*pi,
                          curviness=1,
                          ...){
  
  #be nice to eventually treat time similarly to this
  if(is.character(traits)){
    if(any(traits=="phylogram"|traits=="cladogram")){
      tmp.dev<-dev.cur()
      pdf(file=NULL)
      plot.phylo(attr(contsimmap,"tree")[[1]],show.tip.label=FALSE)
      dev.off()
      dev.set(tmp.dev)
      tmp<-get("last_plot.phylo",envir=.PlotPhyloEnv)[c("edge","yy")]
      nts<-.get.ns(contsimmap,uncompress=TRUE)
      tmp.edges<-as.numeric(gsub("N","",dimnames(contsimmap)[[1]]))
      tmp.edges[!tmp.edges]<-NA
      nodes<-tmp$edge[tmp.edges,2]
      nodes[is.na(nodes)]<-Ntip(contsimmap)+1
      #below could probs be made more efficient...
      if(any(traits=="phylogram")){
        contsimmap<-contsimmap[,c(seq_len(dim(contsimmap)[2]),NA)]
        dimnames(contsimmap)[[2]][dim(contsimmap)[2]]<-"phylogram"
        params<-attr(contsimmap,"params")
        params<-cbind(params,rep(list(NULL),nrow(params)))
        params[["call_info",ncol(params)]]<-list("fxn"="phylogram","formula"=NA)
        attr(contsimmap,"params")<-params
        attr(contsimmap,"traits")[length(attr(contsimmap,"traits"))]<-ncol(params)
        names(attr(contsimmap,"traits"))[length(attr(contsimmap,"traits"))]<-"phylogram"
        tmp.class<-class(contsimmap)
        contsimmap<-unclass(contsimmap)
        contsimmap[,"phylogram",]<-do.call(rbind,lapply(seq_along(tmp$yy),function(ii) lapply(nts[ii,],rep,x=tmp$yy[nodes[ii]])))
        class(contsimmap)<-tmp.class
      }
      if(any(traits=="cladogram")){
        anc.nodes<-tmp$edge[tmp.edges,1]
        anc.nodes[is.na(anc.nodes)]<-Ntip(contsimmap)+1
        diffs<-tmp$yy[nodes]-tmp$yy[anc.nodes]
        dts<-contsimmap:::.get.maps(contsimmap,"dts",uncompress=TRUE)
        dts[]<-lapply(dts,function(ii) if(sum(ii)==0) 0 else cumsum(ii)/sum(ii))
        contsimmap<-contsimmap[,c(seq_len(dim(contsimmap)[2]),NA)]
        dimnames(contsimmap)[[2]][dim(contsimmap)[2]]<-"cladogram"
        params<-attr(contsimmap,"params")
        params<-cbind(params,rep(list(NULL),nrow(params)))
        params[["call_info",ncol(params)]]<-list("fxn"="cladogram","formula"=NA)
        attr(contsimmap,"params")<-params
        attr(contsimmap,"traits")[length(attr(contsimmap,"traits"))]<-ncol(params)
        names(attr(contsimmap,"traits"))[length(attr(contsimmap,"traits"))]<-"cladogram"
        tmp.class<-class(contsimmap)
        contsimmap<-unclass(contsimmap)
        contsimmap[,"cladogram",]<-do.call(rbind,lapply(seq_along(tmp$yy),function(ii) lapply(dts[ii,],function(jj) diffs[ii]*jj^(1/curviness)+tmp$yy[anc.nodes[ii]])))
        class(contsimmap)<-tmp.class
      }
    }
    traits<-pmatch(traits,dimnames(contsimmap)[[2]])
  }
  traits<-dimnames(contsimmap)[[2]][traits]
  traits<-traits[!is.na(traits)]
  if(!length(traits)) traits<-dimnames(contsimmap)[[2]][1] #need to account for possibility of NA trait dimension eventually...
  one.trait<-length(traits)==1
  
  #parse the by arguments...
  bys<-.parse.bys(c('simulation','state','edge','time',dimnames(contsimmap)[[2]]),
                  Col.by,Layer.by,Alpha.by,Mix.by,Wgt.by,Lty.by,Lwd.by,one.trait)
  traits.to.extract<-unique(c(traits,bys[attr(bys,'trait.flag')]))
  if(is.null(traits.to.extract)) traits.to.extract<-0
  contsimmap<-contsimmap[edges,traits.to.extract,sims]
  foo<-function(x){
    tmp<-contsimmap[,x,]
    tmp[seq_along(tmp)]
  }
  parsed.contsimmap<-setNames(lapply(traits.to.extract,foo),traits.to.extract)
  len<-length(parsed.contsimmap[[1]])
  tmp.seq<-seq_len(len)
  lens<-lengths(lapply(parsed.contsimmap[[1]],'[[','values'))
  things.to.extract<-unique(c('edge',if(one.trait) 'time',bys,traits.to.extract)) #always need to extract to edges to keep track of separation!
  #interleaved NAs are intended for ensuring separation of unjoined lines...still need to some processing when splitting things up
  ##by various variables...
  #they may not be necessary ultimately
  extractions<-setNames(rep(list(rep(list(NA),2*len)),length(things.to.extract)),things.to.extract)
  if(any(traits=="phylogram")){
    if(polarize){
      state.foo<-if(inherits(contsimmap,"summarized_contsimmap"))
        lapply(parsed.contsimmap[[1]],function(ii) ii[["states"]][c(rep(1,100),seq_len(nrow(ii[["states"]]))),,drop=FALSE])
      else
        lapply(parsed.contsimmap[[1]],function(ii) ii[["states"]][c(rep(1,100),seq_along(ii[["states"]]))])
      lens<-lens+100
      extractions[[1]][c(TRUE,FALSE)]<-lapply(tmp.seq,function(ii) rep(attr(parsed.contsimmap[[1]][[ii]],'info')[1],lens[ii]))
      for(i in things.to.extract[-1]){
        extractions[[i]][c(TRUE,FALSE)]<-switch(i,
                                                simulation=lapply(tmp.seq,function(ii) rep(attr(parsed.contsimmap[[1]][[ii]],'info')[3],lens[ii])),
                                                state=state.foo,
                                                time=lapply(parsed.contsimmap[[1]],function(ii) ii[["ts"]][c(rep(1,100),seq_along(ii[["ts"]]))]),
                                                phylogram=lapply(parsed.contsimmap[[i]],function(ii) c(ii[["values"]][1],
                                                                                                       (ii[["values"]][2]-ii[["values"]][1])*seq(0.01,1,0.01)+ii[["values"]][1],
                                                                                                       ii[["values"]][-1])),
                                                lapply(parsed.contsimmap[[i]],function(ii) ii[["values"]][c(rep(1,100),seq_along(ii[["values"]]))]))
      }
    }else{
      state.foo<-if(inherits(contsimmap,"summarized_contsimmap"))
        lapply(parsed.contsimmap[[1]],function(ii) ii[["states"]][c(1,seq_len(nrow(ii[["states"]]))),,drop=FALSE])
      else
        lapply(parsed.contsimmap[[1]],function(ii) ii[["states"]][c(1,seq_along(ii[["states"]]))])
      lens<-lens+1
      extractions[[1]][c(TRUE,FALSE)]<-lapply(tmp.seq,function(ii) rep(attr(parsed.contsimmap[[1]][[ii]],'info')[1],lens[ii]))
      for(i in things.to.extract[-1]){
        extractions[[i]][c(TRUE,FALSE)]<-switch(i,
                                                simulation=lapply(tmp.seq,function(ii) rep(attr(parsed.contsimmap[[1]][[ii]],'info')[3],lens[ii])),
                                                state=state.foo,
                                                time=lapply(parsed.contsimmap[[1]],function(ii) ii[["ts"]][c(1,seq_along(ii[["ts"]]))]),
                                                phylogram=lapply(parsed.contsimmap[[i]],function(ii) ii[["values"]][c(1,2,seq_along(ii[["values"]])[-1])]),
                                                lapply(parsed.contsimmap[[i]],function(ii) ii[["values"]][c(1,seq_along(ii[["values"]]))]))
      }
    }
  }else{
    extractions[[1]][c(TRUE,FALSE)]<-lapply(tmp.seq,function(ii) rep(attr(parsed.contsimmap[[1]][[ii]],'info')[1],lens[ii]))
    for(i in things.to.extract[-1]){
      extractions[[i]][c(TRUE,FALSE)]<-switch(i,
                                              simulation=lapply(tmp.seq,function(ii) rep(attr(parsed.contsimmap[[1]][[ii]],'info')[3],lens[ii])),
                                              state=lapply(parsed.contsimmap[[1]],'[[','states'),
                                              time=lapply(parsed.contsimmap[[1]],'[[','ts'),
                                              lapply(parsed.contsimmap[[i]],'[[','values'))
    }
  }
  state.ind<-which(things.to.extract=="state")
  if(inherits(contsimmap,"summarized_contsimmap")&length(state.ind)){
    extractions[[state.ind]]<-do.call(rbind,extractions[[state.ind]])
    extractions[-state.ind]<-lapply(extractions[-state.ind],unlist,use.names=FALSE)
    states<-colnames(attr(contsimmap,"maps")[[1,1,"state"]])
  }else{
    extractions<-lapply(extractions,unlist,use.names=FALSE)
    states<-colnames(attr(contsimmap,'tree')[[1]][["mapped.edge"]])
  }
  sims<-dimnames(contsimmap)[[3]] #Definitely don't sort sims
  edges<-sort(edge.inds(contsimmap)) #I feel like edges should be sorted for now, but may want to revisit in future
  args.ls<-.parse.args(list(...),extractions,bys,traits,one.trait,sims,states,edges)
  
  nedges<-Nedge(contsimmap)
  nnodes<-Nnode(contsimmap,internal.only=FALSE)
  tmp.inds<-cumsum(lens+1)
  tmp.inds1<-tmp.inds-1
  tmp.inds0<-c(1,tmp.inds[-len]+1)
  tmp.edges<-as.numeric(extractions[["edge"]][tmp.inds-1])
  edge.matches<-lapply(seq_len(nedges),function(ii) which(tmp.edges==ii))
  endpt.inds<-vector("list",nnodes)
  tmp.edge.mat<-attr(contsimmap,"tree")[[1]][["edge"]]
  has.matches<-lengths(edge.matches)>0
  for(i in which(has.matches)){
    endpt.inds[[tmp.edge.mat[i,2]]]<-tmp.inds1[edge.matches[[i]]]
  }
  tmp.inds<-match(which(lengths(endpt.inds)==0),tmp.edge.mat[,1])
  tmp.inds<-tmp.inds[!is.na(tmp.inds)]
  tmp.inds<-tmp.inds[has.matches[tmp.inds]]
  for(i in tmp.inds){
    endpt.inds[[tmp.edge.mat[i,1]]]<-tmp.inds0[edge.matches[[i]]]
  }
  # edge.matches<-match(seq_len(nedges),tmp.edges)
  # endpt.inds<-rep(as.integer(NA),nnodes)
  # endpt.inds[attr(contsimmap,"tree")[[1]][["edge"]][,2]]<-
  #   (tmp.inds-1)[edge.matches]
  # tmp.nas<-is.na(endpt.inds[attr(contsimmap,"tree")[[1]][["edge"]][,1]])&
  #   !is.na(c(1,tmp.inds[-len]+2)[edge.matches])
  # endpt.inds[attr(contsimmap,"tree")[[1]][["edge"]][tmp.nas,1]]<-
  #   c(1,tmp.inds[-len]+1)[edge.matches][tmp.nas]
  
  #split into unique argument combos and plot!
  par.fac<-do.call(paste,c(args.ls[c('col','layer','lty','lwd')],sep='@'))
  splits<-rle(par.fac)
  len<-length(splits[['values']])
  tmp.seq<-seq_len(len)
  lens<-splits[['lengths']]
  levs<-unique(splits[['values']])
  inds<-do.call(cbind,lapply(levs,function(ii) splits[['values']]==ii))
  splits<-rep(tmp.seq,lens)
  foo<-function(x){
    lapply(tmp.seq,function(ii) c(if(ii>1&!is.na(x[[ii]][1])) c(NA,x[[ii-1]][lens[ii-1]]),
                                  x[[ii]]))#,
                                  #if(ii<len&!is.na(x[[ii]][lens[ii]])) c(x[[ii+1]][1],NA)))
  }
  if(polarize){
    tmp.xx<-extractions[[if(one.trait) 'time' else traits[1]]]
    tmp.yy<-extractions[[traits[length(traits)]]]
    tmp.yy<-(tmp.yy-min(tmp.yy,na.rm=TRUE))/diff(range(tmp.yy,na.rm=TRUE))*
      (ang.max-ang.min)+ang.min
    xx<-tmp.xx*cos(tmp.yy)
    yy<-tmp.xx*sin(tmp.yy)
  }else{
    xx<-extractions[[if(one.trait) 'time' else traits[1]]]
    yy<-extractions[[traits[length(traits)]]]
  }
  node.xx<-unlist(lapply(endpt.inds,function(ii) mean(xx[ii],na.rm=TRUE)),
                  use.names=FALSE)
  node.yy<-unlist(lapply(endpt.inds,function(ii) mean(yy[ii],na.rm=TRUE)),
                  use.names=FALSE)
  # node.xx<-xx[endpt.inds]
  # node.yy<-yy[endpt.inds]
  xx<-foo(split(xx,splits))
  yy<-foo(split(yy,splits))
  tmp<-strsplit(levs,'@')
  col<-unlist(lapply(tmp,'[[',1),use.names=FALSE)
  layer<-as.numeric(unlist(lapply(tmp,'[[',2),use.names=FALSE))
  ord<-order(layer,decreasing=!reverse.layers)
  lty<-unlist(lapply(tmp,'[[',3),use.names=FALSE)
  if(all(grepl('^\\d+$',lty))){
    lty<-as.numeric(lty)
  }
  lwd<-as.numeric(unlist(lapply(tmp,'[[',4),use.names=FALSE))
  final.args<-args.ls[!(grepl('\\.breaks$',names(args.ls))|names(args.ls)%in%c('col','layer','lty','lwd','x','y'))]
  if(!add){
    if(polarize){
      final.args[["xlim"]]<-c(-1,1)*max(final.args[["xlim"]])
      final.args[["ylim"]]<-final.args[["xlim"]]
    }
    do.call(plot,c(final.args,x=NA))
  }
  for(i in ord){
    final.args[c('x','y','col','lty','lwd')]<-list(unlist(xx[inds[,i]],use.names=FALSE),
                                                   unlist(yy[inds[,i]],use.names=FALSE),
                                                   col[i],lty[i],lwd[i])
    do.call(lines,final.args)
  }
  #figure out some way of grabbing node coordinates
  PP<-list(type=if(polarize) "fan"
             else if(any(traits=="phylogram")) "phylogram"
             else if(any(traits=="cladogram")) "cladogram"
             else if(one.trait) "cladogram"
             else "unrooted",
           use.edge.length=TRUE,
           node.pos=if(any(traits=="phylogram"|traits=="cladogram")) 1 else NULL,
           node.depth=1,
           show.tip.label=FALSE,
           show.node.label=FALSE,
           font=3,
           cex=par("cex"),
           adj=0,
           srt=0,
           no.margin=FALSE,
           label.offset=0,
           x.lim=args.ls[["xlim"]],
           y.lim=args.ls[["ylim"]],
           direction="rightwards", #only allow rightwards at moment, but may change at some point...
           tip.color=rep(par("col"),length.out=Ntip(contsimmap)),
           Ntip=Ntip(contsimmap),
           Nnode=Nnode(contsimmap),
           root.time=NULL,
           align.tip.label=FALSE,
           edge=attr(contsimmap,"tree")[[1]][["edge"]],
           xx=node.xx,
           yy=node.yy)
  assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
}
