#' @title Derives an \code{FLBRP} from life history parameters
#' 
#' @description 
#' Takes an \code{FLPar} object with life history and selectivity parameters
#' and generates a corresponding \code{FLBRP} object. Can use a range of functional forms.
#'
#' @param params an \code{FLPar} object with life history parameters
#' @param growth function that takes an \code{FLPar} object with parameters, by default \code{vonB}
#' @param m \code{character} takes the natural mortality model name, by default gislason 
#' @param mat function that takes an \code{FLPar} object with parameters, by default \code{logistic}
#' @param sel function that takes an \code{FLPar} object with parameters, by default \code{dnormal}
#' @param sr \code{character} value, "bevholt" by default
#' @param range \code{numeric} with age range by default from 0 to 40
#' @param spwn \code{numeric} give propotion of year when spawning occurs, by default is params["a50"]-floor(params["a50"])
#' @param fish \code{numeric} give propotion of year when fishing occurs, by default 0.5        
#' @param units \code{character} for vectors in \code{FLBRP} returned by method
#' @param midyear when growth measured, default 0.5
#' @param ... any other arguments 
#' 
#' @aliases lhEql lhEql-method lhEql,FLPar-method
#' 
#' @return \code{FLBRP} object
#'
#' @seealso \code{\link{lhPar}}, \code{\link{lhRef}}
#'
#' @export lhEql
#' @docType methods
#' @rdname lhEql
#'
#' @seealso  \code{\link{vonB}} \code{\link{lorenzen}} \code{\link{sigmoid}}
#' 
#' @examples
#' \dontrun{
#' data(teleost)
#' alb=teleost[,"Thunnus alalunga"]
#' eql=lhEql(lhPar(alb))
#' }
#' 
#' 
#' 

setMethod("lhEql", signature(params='FLPar'),
  function(params, 
    growth =FLifeE::vonB,
    m      ="gislason",
    sr     ="bevholt",
    mat    =logistic,
    sel    =dnormal,
    range  =c(min=0,max=40,minfbar=1,maxfbar=40,plusgroup=40),
    spwn   =0, #c(params["a50"]-floor(params["a50"])),
    fish   =0.5, # proportion of year when fishing happens
    midyear=0.5, ...){
  
  pNms=dimnames(params)$params
  if ("sl"%in%pNms&!("sel1"%in%pNms))
    dimnames(params)$params["sel1"==pNms]="sel1"
  if ("sr"%in%pNms&!("sel2"%in%pNms))
    dimnames(params)$params["sel2"==pNms]="sel2"
  if ("a1"%in%pNms&!("sel3"%in%pNms))
    dimnames(params)$params["sel3"==pNms]="sel3"
     
  # Check that spwn and fish are [0, 1]
  if (any(spwn > 1) | any(spwn < 0) | any(fish > 1) | any(fish < 0))
    stop("spwn and fish must be in the range 0 to 1\n")

  args<-list(...)

  if (("m.spwn" %in% names(args)))
    m.spwn =args[["m.spwn"]]
  else
    m.spwn=FLQuant(rep(spwn,each=length(range["min"]:range["max"])), dimnames=list(age=range["min"]:range["max"],iter=seq(length(spwn))), units="")

  if (("harvest.spwn" %in% names(args)))
    harvest.spwn =args[["harvest.spwn"]]
  else
    harvest.spwn=FLQuant(rep(spwn,each=length(range["min"]:range["max"])), dimnames=list(age=range["min"]:range["max"],iter=seq(length(spwn))), units="")

  age=FLQuant(range["min"]:range["max"],
              dimnames=list(age =range["min"]:range["max"],
                            iter=dimnames(params)$iter), units="")
  # Get the lengths through different times of the year
  slen   <- growth(age+m.spwn,params) # slen is length at spawning time
  clen   <- growth(age+fish,  params) # clen is length when fishing happens
  midyearlen <- growth(age+midyear,params) # midyear length used for natural mortality
  
  # Corresponding weights
  # bug warning cos of log(NA)
  cwt=FLifeE::len2wt(clen,params)
  swt=len2wt(slen,params)

  mat. =mat(age + m.spwn,params) # maturity is biological therefore + m.spwn
  if (dims(mat.)["min"]==0) mat.[1]=0

  sel. =sel(age + fish,  params) # selectivty is fishery  based therefore + fish

  ## create a FLBRP object to   calculate expected equilibrium values and ref pts
  dms=dimnames(swt)
  
  res=FLBRP(stock.wt       =swt,
            landings.wt    =cwt,
            discards.wt    =cwt,
            bycatch.wt     =cwt,
            mat            =FLQuant(mat.,         dimnames=dms, units=""),
            landings.sel   =FLQuant(sel.,         dimnames=dms, units=""),
            discards.sel   =FLQuant(0,            dimnames=dms, units=""),
            bycatch.harvest=FLQuant(0,            dimnames=dms, units="f"),
            harvest.spwn   =FLQuant(harvest.spwn, dimnames=dms, units=""),
            m.spwn         =FLQuant(m.spwn,       dimnames=dms, units=""),
            availability   =FLQuant(1,            dimnames=dms, units=""),
            range          =range)

  #warning("FLPar%*%FLQuant operator sets 1st dim name to quant regardless")
  #if ("numeric" %in% is(m)) m.=FLQuant(m,dimnames=dimnames(age)) else{
  #  if ("length" %in% names(formals(m)))
  #    m.   =m(length=midyearlen,params=params) # natural mortality is always based on mid year length
  #  else if ("age" %in% names(formals(m))){
  #    m.   =m(age=age+midyear,params=params) # natural mortality is always based on mid year length
  #  }else if ("wt" %in% names(formals(m)))
  #    {
  #    m.   =m(swt,params=params[c("m1","m2")])
  #    }
  #
#res    <<-res
#m.     <<-m
#params.<<-params

  #names(dimnames(m.))[1]="age"}
  if ("character"%in%is(m))
    m(res)=setM(res,m,params)
  else if ("numeric"%in%is(m))
    m(res)[]=m
  else if ("function"%in%is(m))
    m(res)=m(res,params)
  else if ("FLQuant"%in%is(m))
    m(res)=m
  
  ## FApex
  #if (!("range" %in% names(args))) range(res,c("minfbar","maxfbar"))[]<-as.numeric(dimnames(landings.sel(res)[landings.sel(res)==max(landings.sel(res))][1])$age)

  ## replace any slot passed in as an arg
  for (slt in names(args)[names(args) %in% names(getSlots("FLBRP"))[names(getSlots("FLBRP"))!="fbar"]])
    slot(res, slt)<-args[[slt]]
  
  params(res)=propagate(params(res),dims(res)$iter)
    
  ## Stock recruitment relationship
  model(res) =do.call(sr,list())$model
  #return(res)
  
  if (sr=="shepherd" & !("c" %in% dimnames(params)[[1]])){

    dmns=dimnames(params)

    dmns$params=c(dmns$params,"c")

    par.=FLPar(NA,dimnames=dmns)
    par.[dimnames(params)$params]=params
    par.["c"]=1
    params=par.}
  
  if (dims(params)$iter>1) {

    #warning("Scarab, iters dont work for SRR:sv/ab etc")
    #warning("Should be no need to specify mode of FLPar element")
  
    if (sr=="shepherd")
      params(res)=FLPar(c(a=as.numeric(NA),b=as.numeric(NA),c=as.numeric(NA)),iter=dims(params)$iter)
    else
      params(res)=FLPar(c(a=as.numeric(NA),b=as.numeric(NA)),iter=dims(params)$iter)
   
    for (i in seq(dims(params)$iter))
      if (sr=="shepherd")
        params(res)[,i][]=unlist(c(FLCore::ab(params[c("s","v","c"),i],sr,spr0=iter(spr0(res),i))[c("a","b","c")]))
      else {
        params(res)[,i][]=unlist(c(FLCore::ab(params[c("s","v"),i],sr,spr0=iter(spr0(res),i))[c("a","b")]))  
      }
    
    #warning("iter(params(res),i)=ab(params[c(s,v),i],sr,spr0=iter(spr0(res),i))[c(a,b)] assignment doesnt work")
    #warning("iter(FLBRP,i) doesn't work")
  } else{
    
    if (sr=="shepherd"){
      params(res)=FLCore::ab(params[c("s","v","c")],sr,spr0=spr0(res))[c("a","b","c")]
      }else{
      params(res)=FLCore::ab(params[c("s","v")],sr,spr0=spr0(res))[c("a","b")]
      }
    }
  
  refpts(res)=propagate(refpts(res)[c("virgin","msy","crash","f0.1","fmax")],dims(params)$iter)
  res=brp(res)
  
  if ("fbar" %in% names(args))
    fbar(res)<-args[["fbar"]] else
      if (any((!is.nan(refpts(res)["crash","harvest"]))))
        fbar(res)<-FLQuant(seq(0,1,length.out=101),quant="age")%*%refpts(res)["crash","harvest"]
  
  names(dimnames(fbar(res)))[1]="age"
  res=brp(res)
  
  if (!("units" %in% names(attributes(params))))  return(res)
  if (all(is.na(attributes(params)$units)))  return(res)

  try(res <- setUnits(res, params),silent=TRUE)

  return(res)})

#
# setMethod('ab', signature(x='FLparams', model='character'),
#   function(x, model, spr0=NA){
#
#    s=x["a"]
#    v=x["b"]
#    a=FLPar(a=1,dimnames=dimnames(s))
#    b=FLPar(b=1,dimnames=dimnames(v))
#
#    if ("spr0" %in% dimnames(x)$params)
#       spr0=x["spr0"]  else
#       spr0=FLPar(spr0,dimnames=dimnames(a))
#
#    c=FLPar(c=1,dimnames=dimnames(a))
#    d=FLPar(d=1,dimnames=dimnames(a))
#    if (("c" %in% dimnames(x)$params))  c=x["c"]
#    if (("d" %in% dimnames(x)$params))  d=x["d"]
#
#    v <- v*spr2v(model, spr0, a, b, c, d)
#    s <- s*srr2s(model, ssb=v*.2, a=a, b=b, c=c, d=d) / srr2s(model, ssb=v, a=a, b=b, c=c, d=d)
#
#    res=rbind(s, v, spr0)
#
#    if ("c" %in% dimnames(x)$params)
#      res=rbind(res, c)
#
#    if ("d" %in% dimnames(x)$params)
#      res=rbind(res, d)
#
#    res=rbind(res, spr0)
#
#    return(res)})

matFn=function(age,params) {
  a50=FLQuant(ceiling(rep(c(params["a50"]),each=dim(age)[1])),
              dimnames=dimnames(age))
  res=FLQuant(0.5,dimnames=dimnames(age))
  res[age> a50]=1
  res[age< a50]=0
  res}
