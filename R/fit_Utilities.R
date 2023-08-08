.proc.trait.data<-function(trait.data,tips,nodes,edges,ntrees,treeID,nsim){
  trait.data<-.fix.trait.data(trait.data,tips,nodes)
  ntraits<-ncol(trait.data[[1]])
  nms<-c(nodes[1],c(tips,nodes)[edges[,2]])
  mis<-lapply(trait.data,is.na)
  foo<-function(x){
    trait.data[[x]][mis[[x]]]<-0
    split(trait.data[[x]],rownames(trait.data[[x]]))[nms]
  }
  parsed.obs<-do.call(cbind,lapply(seq_along(trait.data),foo))
  nobs<-lengths(parsed.obs)/ntraits
  tmp<-which(nobs>0)
  parsed.obs[tmp]<-lapply(tmp,function(ii) matrix(parsed.obs[[ii]],nobs[ii],ntraits))
  parsed.mis<-do.call(cbind,lapply(mis,function(ii) split(!ii,rownames(ii))[nms]))
  parsed.mis[tmp]<-lapply(tmp,function(ii) matrix(parsed.mis[[ii]],ntraits,nobs[ii],byrow=TRUE))
  ntrait.data<-length(trait.data)
  trait.data.treeID<-(seq_len(ntrait.data)-1)%%ntrees+1
  traitID<-vector('integer',nsim)
  tree.inds<-vector('list',ntrees)
  for(i in seq_len(ntrees)){
    tree.inds[[i]]<-treeID==i
    tmp.n<-sum(tree.inds[[i]])
    traitID[tree.inds[[i]]]<-rep(which(trait.data.treeID%%ntrait.data==i%%ntrait.data),length.out=tmp.n)
  }
  list(traits=colnames(trait.data[[1]]),
       trait.data=trait.data,
       nobs=nobs,
       parsed.obs=parsed.obs,
       parsed.mis=parsed.mis,
       traitID=traitID)
}

#for fixing Xvar, Yvar, and mu formulae
.fix.var.mu<-function(var,nice.name,base,Yvar.flag,contsimmap.traits){
  if(!is.list(var)){
    if(inherits(var,"formula")){
      var<-list(var)
    }else{
      var<-as.list(var)
    }
  }
  #parse out formula response variables
  formulae<-unlist(lapply(var,function(ii) inherits(ii,"formula")),use.names=FALSE)
  tmp.nms<-NULL
  resps<-formulae&lengths(var)==3
  if(any(resps)){
    tmp.nms<-rep("",length(var))
    tmp.nms[resps]<-unlist(lapply(var[resps],function(ii) all.vars(ii[[2]])[1]),use.names=FALSE)
  }
  var[formulae]<-lapply(var[formulae],function(ii) as.character(ii[length(ii)]))
  var[]<-lapply(var,function(ii) if(length(ii)) str2lang(as.character(ii)) else NULL)
  #cut out NULL formulae!
  var<-var[lengths(var)>0]
  nms<-names(var)
  if(!is.null(tmp.nms)){
    if(!length(nms)){
      nms<-tmp.nms
    }else{
      tmp<-!nzchar(nms)&nzchar(tmp.nms)
      nms[tmp]<-tmp.nms[tmp]
    }
  }
  n<-length(var)
  out.nms<-names(base)
  if(!length(nms)){
    nms<-names(var)<-rep('',n)
    prob.nms<-rep(TRUE,n)
    has.prob.nms<-TRUE
    quiet.flag<-TRUE
  }else{
    prob.nms<-!(nms%in%out.nms)
    has.prob.nms<-any(prob.nms)
    quiet.flag<-FALSE
  }
  if(any(!prob.nms)){
    base[nms[!prob.nms]]<-var[!prob.nms]
  }
  prob.traits<-!(out.nms%in%nms)
  if(any(prob.traits)){
    message<-" and set to default"
    if(has.prob.nms){
      nrem<-sum(prob.traits)
      tmp.n<-sum(prob.nms)
      if(nrem>=tmp.n){
        tmp.prob.traits<-which(prob.traits)[seq_len(tmp.n)]
        tmp.prob.nms<-prob.nms
      }else{
        tmp.prob.traits<-prob.traits
        tmp.prob.nms<-which(prob.nms)[seq_len(nrem)]
      }
      base[tmp.prob.traits]<-var[tmp.prob.nms]
      if(!quiet.flag){
        message<-paste0(if(identical(out.nms[prob.traits],out.nms[tmp.prob.traits]))
          " and "
          else 
            .report.names(c('; missing formula for trait','; missing formulae for traits'),out.nms[tmp.prob.traits],),
          "filled with unmatched formulae from ",
          nice.name)
        if(nrem>tmp.n){
          message<-paste0(message," with the rest set to default")
        }
      }
    }
    if(!quiet.flag){
      warning(.report.names(c('Formula for trait','Formulae for traits'),out.nms[prob.traits]),
              'missing in ',
              nice.name,
              message,
              immediate.=TRUE)
    }
  }
  base
}

#for fixing Xcor, Ycor formulae
.fix.cor<-function(cor,nice.name,base,holder,diag.inds,contsimmap.traits){
  if(!is.list(cor)){
    if(inherits(cor,"formula")){
      cor<-list(cor)
    }else{
      cor<-as.list(cor)
    }
  }
  #parse out formula response variables
  formulae<-unlist(lapply(cor,function(ii) inherits(ii,"formula")),use.names=FALSE)
  tmp.nms<-NULL
  resps<-formulae&lengths(cor)==3
  if(any(resps)){
    tmp.nms<-rep("",length(cor))
    tmp.nms[resps]<-unlist(lapply(cor[resps],function(ii) paste(all.vars(ii[[2]],unique=FALSE),collapse="_BY_")),use.names=FALSE)
  }
  cor[formulae]<-lapply(cor[formulae],function(ii) as.character(ii[length(ii)]))
  cor[]<-lapply(cor,function(ii) if(length(ii)) str2lang(as.character(ii)) else NULL)
  ndims<-length(dim(cor))
  if(ndims>2){
    stop(nice.name,' has more than two dimensions; multiple correlation matrices should be stored in lists, not arrays')
  }else if(ndims<2){
    nms<-names(cor)
    if(is.null(nms)) nms<-rep("",length(cor))
    if(!is.null(tmp.nms)){
      tmp<-!nzchar(nms)&nzchar(tmp.nms)
      nms[tmp]<-tmp.nms[tmp]
    }
    #alright, I THINK this is the best way to do things...
    #truly assumes column-major order for lower triangle unless otherwise specified, I think
    nms<-strsplit(nms,"_BY_")
    nms<-do.call(rbind,lapply(nms,function(ii) ii[1:2]))
    nas<-is.na(nms)
    tmp.vars<-unique(nms[!nas])
    tmp.nms<-unique(t(apply(nms,1,sort,na.last=TRUE)))
    tmp.nms<-tmp.nms[!is.na(tmp.nms[,1]),,drop=FALSE]
    #better calculation of implied number of correlations
    #ignores specified diagonal entries and avoids double-counting specified symmetric pairs
    tmp.n<-nrow(tmp.nms)+sum(nas[,1])-sum(tmp.nms[,1]==tmp.nms[,2],na.rm=TRUE)
    n.tmp.vars<-length(tmp.vars)
    if(n.tmp.vars){
      max.tmp.vars<-max(unlist(lapply(tmp.vars,function(ii) sum(apply(tmp.nms,1,function(jj) any(jj==ii)),na.rm=TRUE)),use.names=FALSE))
    }else{
      max.tmp.vars<-0
    }
    n.tmp.vars<-length(tmp.vars)
    #if n choose k = x, then n is between sqrt(k!*x) and sqrt(k!*x)+k
    #so check sqrt(2*x)+0:2, then take the minimum integer x
    checks<-ceiling(sqrt(2*tmp.n)+0:2)
    tmp.n<-max(checks[choose(checks,2)>=tmp.n][1],n.tmp.vars,max.tmp.vars)
    nrem<-tmp.n-n.tmp.vars
    if(nrem>0){
      out.nms<-unique(as.vector(nms),incomparables=NA)
      tmp.nas<-is.na(out.nms)
      n.nas<-sum(tmp.nas)
      tmp.seq<-seq_len(min(n.nas,nrem))
      out.nms[tmp.nas][tmp.seq]<-paste0("TEMPORARY_",tmp.seq)
      if(nrem>n.nas)  out.nms<-c(out.nms,paste0("TEMPORARY_",(n.nas+1):nrem))
      out.nms<-unique(out.nms[!is.na(out.nms)])
    }else{
      out.nms<-tmp.vars
    }
    #now just have to figure out how to inset lower triangle column major sequence into existing nms matrix...
    cols<-.col(c(tmp.n,tmp.n))
    rows<-t(cols)
    def.nms<-cbind(out.nms[rows[rows>cols]],
                   out.nms[cols[cols<rows]])
    n.nas<-rowSums(nas)
    completes<-n.nas==0
    if(any(completes)){
      tmp.nms<-apply(nms[completes,,drop=FALSE],1,function(ii) paste(ii[order(match(ii,out.nms),decreasing=TRUE)],collapse="_"))
      inds<-paste(def.nms[,1],def.nms[,2],sep="_")%in%tmp.nms
      def.nms<-def.nms[!inds,,drop=FALSE]
    }
    semis<-n.nas==1
    if(any(semis)){
      tmp.nms<-nms[semis,1]
      for(i in seq_along(tmp.nms)){
        tmp<-match(tmp.nms[i],def.nms[,1])
        tmp.ind<-2
        #I think, based on things above, that this  should always work...
        if(is.na(tmp)){
          tmp<-match(tmp.nms[i],def.nms[,2])
          tmp.ind<-1
        }
        #just in case...
        if(is.na(tmp)) stop("dang it my correlation matrix formatting function didn't work")
        nms[semis,2][i]<-def.nms[tmp,tmp.ind]
        def.nms<-def.nms[-tmp,,drop=FALSE]
      }
    }
    #again, I THINK this should always work!
    incompletes<-!(completes|semis)
    n.incompletes<-sum(incompletes)
    nms[incompletes,]<-def.nms[seq_len(n.incompletes),,drop=FALSE]
    tmp<-matrix(list(),length(out.nms),length(out.nms),
                dimnames=list(out.nms,out.nms))
    for(i in seq_along(cor)){
      tmp[[nms[i,1],nms[i,2]]]<-cor[[i]]
    }
    cor<-tmp
    rownames(cor)[grepl("^TEMPORARY_\\d+$",rownames(cor))]<-""
    colnames(cor)[grepl("^TEMPORARY_\\d+$",colnames(cor))]<-""
  }else if(!is.null(tmp.nms)){
    dims<-dim(cor)
    nms<-dimnames(cor)
    if(!length(nms)) nms<-rep(list(NULL),2)
    for(i in c(1,2)){
      if(is.null(nms[[i]])) nms[[i]]<-rep('',dims[i])
    }
    tmp.nms<-strsplit(tmp.nms,"_BY_")
    tmp.nms[!lengths(tmp.nms)]<-list("")
    tmp.nms<-lapply(tmp.nms,function(ii) if(length(ii)==1) c(ii,"") else ii[c(1,2)])
    tmp.rownms<-matrix(unlist(lapply(tmp.nms,'[[',1),use.names=FALSE),dims[1],dims[2])
    tmp.rownms<-apply(tmp.rownms,1,function(ii) c(ii[nzchar(ii)],"")[1])
    tmp.colnms<-matrix(unlist(lapply(tmp.nms,'[[',2),use.names=FALSE),dims[1],dims[2])
    tmp.colnms<-apply(tmp.colnms,1,function(ii) c(ii[nzchar(ii)],"")[1])
    tmp<-!nzchar(nms[[1]])&nzchar(tmp.rownms)
    nms[[1]][tmp]<-tmp.rownms[tmp]
    tmp<-!nzchar(nms[[2]])&nzchar(tmp.colnms)
    nms[[2]][tmp]<-tmp.colnms[tmp]
  }
  dims<-dim(cor)
  nms<-dimnames(cor)
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
    nms[[i]][minseq][inds[minseq]&other.inds[minseq]]<-other[inds[minseq]&other.inds[minseq]]
    inds<-!nzchar(nms[[i]])
    nms[[i]][inds]<-paste0('TEMPORARY_',seq_len(sum(inds))+offset)
    offset<-sum(inds)
  }
  dimnames(cor)<-nms
  #append rows corresponding to cols
  missing.rows<-!(nms[[2]]%in%nms[[1]])
  cor<-do.call(rbind,c(list(cor),rep(list(NULL),sum(missing.rows))))
  rownames(cor)<-c(nms[[1]],nms[[2]][missing.rows])
  nms<-dimnames(cor)
  #append cols corresponding to rows
  missing.cols<-!(nms[[1]]%in%nms[[2]])
  cor<-do.call(cbind,c(list(cor),rep(list(NULL),sum(missing.cols))))
  colnames(cor)<-c(nms[[2]],nms[[1]][missing.cols])
  #sort
  nms<-sort(nms[[1]])
  cor<-cor[nms,nms,drop=FALSE]
  #check for invalid formulae including contsimmap traits
  #barebones check for now--could make more specific for users to pinpoint issue...
  probs<-unlist(lapply(cor,function(ii) any(all.vars(ii)%in%contsimmap.traits)),use.names=FALSE)
  if(any(probs)){
    cor[probs]<-list(NULL)
    warning("One or more formulae in ",
            nice.name,
            " were ignored because they referenced a mapped continuous trait; ",
            "due to certain mathematical difficulties, correlations are not (yet) allowed to vary as ",
            "functions of continuous variables (though you can make them vary according to state)")
  }
  #reflect
  lo<-lower.tri(cor)
  up<-upper.tri(cor)
  tcor<-t(cor)
  inds<-lengths(cor)==0
  inds<-inds&t(!inds)
  cor[lo&inds]<-tcor[lo&inds]
  cor[up&inds]<-tcor[up&inds]
  prob.nms<-grepl('^TEMPORARY_\\d+$',nms)
  out.nms<-dimnames(base)[[1]]
  if(all(prob.nms)){
    has.prob.nms<-TRUE
    quiet.flag<-TRUE
  }else{
    prob.nms<-prob.nms|!(nms%in%out.nms)
    has.prob.nms<-any(prob.nms)
    quiet.flag<-FALSE
  }
  if(any(!prob.nms)){
    holder[nms[!prob.nms],nms[!prob.nms]]<-cor[!prob.nms,!prob.nms]
  }
  prob.traits<-!(out.nms%in%nms)
  #I think no replication makes the most sense here
  if(any(prob.traits)&has.prob.nms){
    nrem<-sum(prob.traits)
    tmp.n<-sum(prob.nms)
    if(nrem>=tmp.n){
      tmp.prob.traits<-which(prob.traits)[seq_len(tmp.n)]
      tmp.prob.nms<-prob.nms
    }else{
      tmp.prob.traits<-prob.traits
      tmp.prob.nms<-which(prob.nms)[seq_len(nrem)]
    }
    holder[tmp.prob.traits,tmp.prob.traits]<-cor[tmp.prob.nms,tmp.prob.nms] #this is fine
    #but now have to think about how to insert other blocks...
    if(any(!prob.nms)){
      holder[tmp.prob.traits,nms[!prob.nms]]<-cor[tmp.prob.nms,!prob.nms]
      holder[nms[!prob.nms],tmp.prob.traits]<-cor[!prob.nms,tmp.prob.nms]
    }
    if(!quiet.flag){
      warning(.report.names(c('Trait','Traits'),out.nms[prob.traits],c('has','have')),
              ' no associated formulae in ',
              nice.name,
              "; missing formulae ",
              if(!identical(out.nms[tmp.prob.traits],out.nms[prob.traits])) .report.names(c('for trait','for traits'),out.nms[tmp.prob.traits]),
              "filled with unmatched formulae from ",
              nice.name,
              immediate.=TRUE)
    }
  }
  prob.entries<-lengths(holder)==0
  prob.diags<-!prob.entries[diag.inds]
  prob.entries[diag.inds]<-FALSE
  if(any(prob.diags)){
    holder[diag.inds]<-list(NULL)
    warning('Diagonal ',
            .report.names(c('formula for trait','formulae for traits'),traits[prob.diags]),
            'in ',
            nice.name,
            " ignored because ",
            nice.name,
            " specifies a correlation matrix which MUST have 1s along the diagonal",
            immediate.=TRUE)
  }
  if(any(prob.entries)){
    tmp<-unique(t(apply(which(prob.entries,arr.ind=TRUE),1,sort,decreasing=TRUE)))
    tmp[]<-out.nms[tmp]
    tmp<-paste0(tmp[,1]," and ",tmp[,2])
    holder[prob.entries]<-base[prob.entries]
    warning(.report.names(c("Formula for correlation between traits","Formulae for correlations between traits"),tmp),
            "missing in ",
            nice.name,
            " and set to default",
            immediate.=TRUE)
  }
  if(!isSymmetric(base)){
    stop('Failed to create symmetric correlation matrix formula from ',nice.name,
         immediate.=TRUE)
  }
  holder
}

#to fix calls to the make.cor formula!
.fix.make.cor<-function(cor,nice.name,base,contsimmap.traits){
  cor<-cor[[length(cor)]]
  args<-as.list(cor[-1])
  ntraits<-dim(base)[1]
  lt<-lower.tri(base)
  out.len<-sum(lt)
  #not bothering with name matching for now--but will want to add in future
  #will look similar to stuff in the .fix.cor function--args will have names corresponding to user inputs
  #(no arguments should be formulae, but may want to add some checks for this in the future...)
  #Oh, an also need to check for any references to contsimmap traits in the formulae! Cause that will break things for sure
  if(length(args)>out.len){
    args<-args[seq_len(out.len)]
  }else if(length(args)<out.len){
    args<-c(args,base[lt][(length(args)+1):out.len])
  }
  args<-c(args,list(ntraits))
  cor[2:(length(args)+1)]<-args
  cor
}

# cor<-list("g",a_BY_b~cool,c(c,c)~sweet)
# holder<-matrix(list(),3,3,dimnames=list(letters[1:3],letters[1:3]))
# nice.name<-"Xcor[[1]]"
# par<-list(.fix.cor(cor,nice.name,tmp.base,holder,!holder.inds))
# 
# par<-list(list(0.5,c~g,b~3*f),list(a~z))
# traits<-letters[1:3]
# states<-c("0","1")
.fix.par.formulae<-function(par,traits,states,contsimmap.traits){
  ntraits<-length(traits)
  nstates<-length(states)
  type<-deparse(substitute(par))
  cor.flag<-type=="Xcor"|type=="Ycor"
  #wrap par in list if it's a formula (otherwise as.list() below might break things)
  if(inherits(par,"formula")) par<-list(par)
  list.flag<-is.list(par)
  ndims<-length(dim(par))
  #I think this should work, but it's a confusing list of lists type deal!
  #if this is for a cor mat, wrap as list only if input is an array of lists, a list including formulae with c() response variables and/or with names including "_BY_", OR is not a list
  #don't have to worry about checking for ntraits being 1 or anything like that since Xcor/Ycor don't matter in such cases and will be ignored by parent function
  if(cor.flag){
    probs.single.state.flag<-FALSE
    if(list.flag){
      if(ndims<2){
        tmp<-unlist(lapply(par,function(ii) inherits(ii,"formula")),use.names=FALSE)
        tmp<-tmp&lengths(par)==3
        if(any(tmp)){
          tmp<-lapply(par[tmp],function(ii) all.vars(ii[[2]],unique=FALSE))
        }else{
          tmp<-NULL
        }
        tmp<-c(tmp,if(!is.null(names(par))) strsplit(names(par),"_BY_"))
        if(!is.null(tmp)){
          probs.single.state.flag<-any(lengths(tmp)>1)
        }
      }else{
        probs.single.state.flag<-TRUE
      }
    }
    if(probs.single.state.flag|!list.flag){
      par<-list(par)
    }
    tmp.holder<-holder<-matrix(list(),ntraits,ntraits,dimnames=list(traits,traits))
    rows<-.row(c(ntraits,ntraits))
    diag(rows)<-0
    cols<-t(rows)
    ut<-rows<cols
    rows[ut]<-t(rows)[ut]
    cols[ut]<-t(cols)[ut]
    base<-paste0("tanh(PLACEHOLDER_",traits[rows],'_BY_',traits[cols],')')
    holder.inds<-ut|t(ut)
  }else{
    #this is a bit more confusing than anticipated...
    #alright, so if it is a list input, check if it's a matrix of lists and split by columns if so (always assume columns correspond to states here--may want to revisit at some point?)
    #otherwise, if it's not a list, check number of dimensions
    #if there less than two, coerce to list if there's only 1 trait and multiple states (assume multiple entries correspond to different states), otherwise just wrap as list (assume entries correspond to different traits)
    #if there is two, split by columns
    #if there is more than two return error
    #assume multiple entries correspond to traits if there's only 1 state
    if(list.flag){
      if(ndims==2){
        par<-asplit(par,2)
      }
    }else{
      if(ndims<2){
        if(ntraits==1&nstates>1&length(par)>0){
          par[]<-as.list(par)
        }else{
          par<-list(par)
        }
      }else if(ndims==2){
        par<-asplit(par,2)
      }else{
        stop("PLACEHOLDER: unrecognized par format for ",type)
      }
    }
    if(nstates==1&length(par)>1){
      par<-list(par)
    }
    tmp.holder<-holder<-setNames(vector('list',ntraits),traits)
    if(type=="mu"){
      base<-setNames(paste0("PLACEHOLDER_",traits),traits)
    }else{
      base<-setNames(paste0("exp(PLACEHOLDER_",traits,")"),traits)
    }
    holder.inds<-rep(TRUE,ntraits)
  }
  nms<-names(par)
  n<-length(par)
  if(!length(nms)){
    nms<-names(par)<-rep('',n)
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
      nice.name<-paste0(type,"[['",nm,"']]")
      tmp.holder[holder.inds]<-lapply(gsub("PLACEHOLDER",paste0(type,"_",nm),base),str2lang)
      tmp.base<-tmp.holder
    }else{
      nice.name<-paste0(type,"[[",i,"]]")
      tmp.holder[holder.inds]<-lapply(gsub("PLACEHOLDER",paste0(type,"_",i),base),str2lang)
      tmp.base<-tmp.holder
    }
    if(length(par[[i]])){
      if(cor.flag){
        make.cor.flag<-FALSE
        #check for new, special, make.cor() function reference
        if(inherits(par[[i]],"formula")){
          char.form<-as.character(par[[i]])
          make.cor.flag<-grepl("^make.cor\\(.*\\)$",char.form[length(char.form)])
        }
        if(make.cor.flag){
          par[[i]]<-.fix.make.cor(par[[i]],nice.name,tmp.base,contsimmap.traits)
        }else{
          par[[i]]<-.fix.cor(par[[i]],nice.name,tmp.base,holder,!holder.inds,contsimmap.traits)
        }
      }else{
        par[[i]]<-.fix.var.mu(par[[i]],nice.name,tmp.base)
      }
    }else{
      par[[i]]<-tmp.base
    }
  }
  out<-setNames(numeric(nstates),states)
  if(any(!prob.nms)){
    out[nms[!prob.nms]]<-which(!prob.nms)
  }
  prob.states<-!(states%in%nms)
  if(any(prob.states)){
    message<-"set to default"
    if(has.prob.nms){
      out[prob.states]<-rep(which(prob.nms),length.out=sum(prob.states))
      message<-paste0("filled with unmatched ",
                      if(cor.flag) 'matrices' else 'vectors',
                      ' from ',
                      type)
    }else{
      tmp.holder[holder.inds]<-lapply(gsub("PLACEHOLDER",paste0(type,"_DEF"),base),str2lang)
      par<-c(par,list(tmp.holder))
      out[prob.states]<-rep(length(par),length.out=sum(prob.states))
    }
    if(!quiet.flag){
      tmp1<-if(cor.flag) c("matrix","matrices") else c("vector","vectors")
      tmp2<-if(type=='Ysig2'|type=='Ycor') c("tip/node","tips/nodes") else c("state","states")
      prefix<-paste0("Formula ",tmp1," for ",tmp2)
      warning(.report.names(prefix,states[prob.states]),
              'missing from ',
              type,
              ' and ',
              message,
              immediate.=TRUE)
    }
  }
  attr(par,"inds")<-out
  par
}

#' @export
make.cor<-function(...){
  tmp<-unlist(list(...),use.names=FALSE)
  len<-length(tmp)
  ntraits<-tmp[len]
  tmp<-tmp[-len]
  out<-input<-matrix(0,ntraits,ntraits)
  lt<-lower.tri(input)
  input[lt]<-tmp
  tally<-numeric(ntraits)
  for(i in seq_len(ntraits)){
    out[i,i]<-sqrt(1-tally[i])
    out[lt[,i],i]<-input[lt[,i],i]*sqrt(1-tally[lt[,i]])
    tally<-tally+out[,i]^2
  }
  out%*%t(out)
}