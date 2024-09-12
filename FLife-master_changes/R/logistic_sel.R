setGeneric("logistic_Sel",function(age,params,...)
  standardGeneric("logistic_Sel"))

setMethod("logistic_Sel", signature(age='FLQuant',params='FLPar'),
          function(age,params,...) { 
            res=logisticFn_Sel(age,params)
            res@units='proportion'
            res})

logisticFn_Sel <- function(age,params){
  
  pNms=dimnames(params)$params
  
  if ("sel1"%in%pNms){
    params=params[dimnames(params)$params!="a1"]
    dimnames(params)$params["sel1"==pNms]="a1"
  }
  
  a1=FLQuant(1,dimnames=dimnames(age))%*%params["a1"]
  
  res =params["asym"]%/%(1.0+pow(19.0,(params["a1"]%-%age)%/%params["ato95"]))
  res[is.na(res)]=0
  asym=FLQuant(1,dimnames=dimnames(age))%*%params["asym"]
  grt =(params["a1"]%-%age)%/%params["ato95"] >  5
  lss =(params["a1"]%-%age)%/%params["ato95"] < -5
  
  res[grt]=0
  res@.Data[lss]=asym@.Data[lss]
  
  dmns          =dimnames(res)
  names(dmns)[1]="age"
  dimnames(res) =dmns
  
  return(res)
}