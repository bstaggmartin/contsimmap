#helper function to merge lists in contsimmap by trait
#then automatically makes a function that resplits and integrates those merged lists back into contsimmap object
.extract.traits<-function(contsimmap,foc.traits,lhs,traits2rates=FALSE){
  flags<-setNames(rep(TRUE,4),paste0(c('n','t','dt','s'),'_PH'))
  inds<-rep(TRUE,length(foc.traits))
  for(i in names(flags)){
    tmp<-foc.traits==i
    inds<-inds&!tmp
    if(i!='n_PH'){
      flags[i]<-any(tmp)
    }
  }
  foc.traits<-foc.traits[inds]
  lhs<-lhs[inds]
  ntraits<-ntrait(contsimmap)
  nts<-lapply(contsimmap[['x']],function(ii) lengths(ii)/ntraits)
  edge.seq<-seq_along(nts)
  edge.split.inds<-rep(edge.seq,unlist(lapply(nts,sum),use.names=FALSE))
  tree.seq<-seq_along(contsimmap[['tree']])
  tree.split.inds<-lapply(nts,function(ii) rep(tree.seq,ii))
  trait.ind<-match(foc.traits,contsimmap[['traits']])
  #use &lhs to exclude non-existent traits that aren't being created in call
  #(likely represent variables stored in global environment!)
  prob.traits<-is.na(trait.ind)&lhs
  if(any(prob.traits)){
    new.traits<-unique(foc.traits[prob.traits])
    n.new.traits<-length(new.traits)
    #use subset to extend data to include new traits
    contsimmap<-subset(contsimmap,traits=c(seq_len(ntraits),rep(NA,n.new.traits)))
    contsimmap[['traits']][seq_len(n.new.traits)+ntraits]<-new.traits
    if(!is.null(contsimmap[['trait.data']])){
      dimnms<-dimnames(contsimmap[['trait.data']])
      dimnms[[2]][seq_len(n.new.traits)+ntraits]<-new.traits
    }
    trait.ind<-match(foc.traits,contsimmap[['traits']])
  }
  inds<-!is.na(trait.ind)
  foc.traits<-foc.traits[inds]
  trait.ind<-trait.ind[inds]
  lhs<-lhs[inds]
  trait.seq<-seq_along(foc.traits)
  if(traits2rates){
    #exclude nodes (they don't affect regimes) and make resplit.and.integrate "fxn" just a contsimmap with indices for constructing simmaps
    out<-setNames(rep(contsimmap['x'],length(foc.traits)),foc.traits)
    resplit.and.integrate<-unname(contsimmap['x'])
    counter<-0
    for(i in seq_along(foc.traits)){
      for(j in edge.seq){
        for(k in tree.seq){
          tmp<-out[[i]][[j]][[k]][,trait.ind[i],,drop=FALSE]
          out[[i]][[j]][[k]]<-tmp
          if(i==1){
            tmp.nt<-nts[[j]][k]
            tmp[]<-seq_len(tmp.nt)+counter
            counter<-counter+tmp.nt
            resplit.and.integrate[[1]][[j]][[k]]<-tmp
          }
        }
      }
      out[[i]]<-unlist(out[[i]],use.names=FALSE)
    }
  }else{
    nodes<-contsimmap[['nodes']]
    nnodes<-length(nodes)
    node.array<-array(unlist(nodes,use.names=FALSE),
                      c(dim(nodes[[1]]),nnodes),
                      list(NULL,NULL,names(nodes)))
    tmp.nodes<-asplit(node.array[trait.ind,,,drop=FALSE],1)
    out<-setNames(rep(contsimmap['x'],length(foc.traits)),foc.traits)
    for(i in seq_along(foc.traits)){
      for(j in edge.seq){
        for(k in tree.seq){
          out[[i]][[j]][[k]]<-out[[i]][[j]][[k]][,trait.ind[i],,drop=FALSE]
        }
      }
      out[[i]]<-c(tmp.nodes[[i]],unlist(out[[i]],use.names=FALSE))
    }
    init.split<-rep.int(2,length(out[[1]]))
    init.split[seq_along(tmp.nodes[[1]])]<-1
    resplit.and.integrate<-function(foc.traits,nms,contsimmap){
      tmp.foc.inds<-match(nms,names(foc.traits))
      tmp.trait.inds<-match(nms,contsimmap[['traits']])
      nodes<-contsimmap[['nodes']]
      node.array<-array(unlist(nodes,use.names=FALSE),
                        c(dim(nodes[[1]]),nnodes),
                        list(NULL,NULL,names(nodes)))
      for(i in seq_along(nms)){
        init<-split(foc.traits[[tmp.foc.inds[i]]],init.split)
        edges<-split(init[[2]],edge.split.inds)
        nodes<-init[[1]]
        for(j in edge.seq){
          tmp<-split(edges[[j]],tree.split.inds[[j]])
          for(k in tree.seq){
            contsimmap[['x']][[j]][[k]][,tmp.trait.inds[i],]<-tmp[[k]]
          }
        }
        node.array[tmp.trait.inds[i],,]<-nodes
      }
      contsimmap[['nodes']]<-asplit(node.array,3)
      contsimmap
    }
  }
  if(flags['n_PH']){
    out[['n_PH']]<-length(out[[1]])
  }
  focs<-names(flags[-1])[flags[-1]]
  if(length(focs)){
    foc.seq<-seq_along(focs)
    grabs<-c('ts','dts','state')[flags[-1]]
    maps<-contsimmap[['maps']]
    out[focs]<-list(rep(vector('list',length(contsimmap[['tree']])),length(nts)))
    treeID<-contsimmap[['perm.inds']][,1]
    zeros<-rep(0,length(treeID)*nnodes)
    treeID<-as.integer(factor(treeID,unique(treeID)))
    per.tree<-tabulate(treeID)
    inds<-rep(tree.seq,per.tree)
    for(i in edge.seq){
      for(j in tree.seq){
        for(k in foc.seq){
          out[[focs[k]]][[i]][[j]]<-maps[[i]][[j]][[grabs[k]]]
        }
      }
      for(k in foc.seq){
        out[[focs[k]]][[i]]<-out[[focs[k]]][[i]][inds]
      }
    }
    for(i in focs){
      if(i=='t_PH'|i=='dt_PH'){
        out[[i]]<-c(zeros,unlist(out[[i]],use.names=FALSE))
      }else if(i=='s_PH'){
        node.states<-as.character(zeros)
        node.labs<-names(contsimmap[['nodes']])
        ancs<-contsimmap[['anc']]
        #implicitly picks first edge and sticks with that--might cause some issues in weird edge cases
        #for example, if descending edges from node have different initial states
        #that's impossible under a normal markov model, but in the case of regime mapping things could get weird...
        node.des<-match(node.labs,ancs)
        out[[i]]<-unlist(c(lapply(node.des,function(ii) lapply(out[[i]][[ii]],'[',1)),out[[i]]),use.names=FALSE)
      }
    }
  }
  out<-c(out,contsimmap=list(contsimmap),fxn=resplit.and.integrate)
  out
}

.fix.formulae<-function(formulae,traits){
  lens<-lengths(formulae)
  probs<-lens<3
  if(any(probs)){
    all.traits<-c(unique(unlist(lapply(formulae,all.vars),use.names=FALSE)),
                  traits)
    def.nms<-grepl('^trait_\\d+$',all.traits)
    if(any(def.nms)){
      def.nms.offset<-max(as.numeric(gsub('^trait_(\\d+)','\\1',all.traits[def.nms])))
    }else{
      def.nms.offset<-0
    }
    new.def.nms<-paste0('trait_',seq_len(sum(probs))+def.nms.offset)
    deparsed.formulae<-unlist(lapply(formulae[probs],deparse),use.names=FALSE)
    formulae[probs]<-lapply(paste0(new.def.nms,deparsed.formulae),str2lang)
  }
  formulae
}

#for making new traits given arbitrary formula of existing traits
#' @export
make.traits<-function(contsimmap,...){
  #old formula deparser:
  # deparsed.formulae<-unlist(lapply(formulae,deparse),use.names=FALSE)
  # allowed.prefixes<-c(' ','^',',','\\(')
  # lhs.patterns<-lapply(foc.traits,function(ii) paste0(allowed.prefixes,ii))
  # lhs.patterns<-lapply(lhs.patterns,function(ii) paste0(ii,rep(c('~',' ~'),each=length(allowed.prefixes))))
  # lhs.patterns<-unlist(lapply(lhs.patterns,paste0,collapse='|'),use.names=FALSE)
  # lhs<-unlist(lapply(seq_along(lhs.patterns),function(ii) any(grepl(lhs.patterns[ii],deparsed.formulae))),use.names=FALSE)
  formulae<-.fix.formulae(list(...),contsimmap[['traits']])
  resps<-lapply(formulae,function(ii) all.vars(ii[[2]]))
  calls<-lapply(formulae,'[[',3)
  all.traits<-unique(unlist(lapply(formulae,all.vars),use.names=FALSE))
  lhs<-all.traits%in%unlist(resps,use.names=FALSE)
  call.envir<-.extract.traits(contsimmap,all.traits,lhs)
  resplit.and.integrate<-call.envir[['fxn']]
  contsimmap<-call.envir[['contsimmap']]
  call.envir<-call.envir[seq_len(length(call.envir)-2)]
  check.length<-function(x){
    len<-length(x)
    if(len==1){
      keys<-attr(x,'keys')
      x<-rep(x,call.envir[['n_PH']])
      if(length(keys)>0){
        attr(x,'keys')<-keys
      }
    }else if(len!=call.envir[['n_PH']]){
      stop("New trait is of wrong length: double-check provided formulae!")
    }
    x
  }
  for(i in seq_along(formulae)){
    tmp.resp<-resps[[i]]
    tmp.res<-eval(calls[[i]],envir=call.envir)
    if(length(tmp.resp)==1){
      if(is.list(tmp.res)){
        tmp.res<-tmp.res[[1]]
      }
      call.envir[[tmp.resp]]<-check.length(tmp.res)
    }else{
      if(!is.list(tmp.res)){
        tmp.res<-list(tmp.res)
      }
      tmp.res<-lapply(tmp.res,check.length)
      #maybe should make this check stricter--tmp.res is recycled to match length of tmp.resp
      call.envir[tmp.resp]<-tmp.res
    }
    #this function's unfortunately a little sluggish...but I don't know what else to do
    contsimmap<-resplit.and.integrate(call.envir,tmp.resp,contsimmap)
  }
  
  #check if any foc.traits have a 'keys' attribute and add to contsimmap
  keys<-lapply(call.envir,function(ii) attr(ii,'keys'))
  keys<-keys[lengths(keys)>0]
  if(length(keys)){
    matches<-match(names(keys),names(contsimmap[['keys']]))
    new<-is.na(matches)
    contsimmap[['keys']][matches[!new]]<-keys[!new]
    contsimmap[['keys']]<-c(contsimmap[['keys']],keys[new])
  }
  contsimmap
}