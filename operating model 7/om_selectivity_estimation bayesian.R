Catch_proportion <- function(data){
  catch_prop <- apply(data,2,function(x){res <- x/sum(x)})
  return(catch_prop)
}

om <- function(params,sel_f,fhist,inp_file){
  
  n_iter <- 500
  yrs_hist <- 100
  yrs_proj <- 50
  
  ######################
  
  if (identical(fhist, "random")) {
    start <- rep(0, n_iter)
    middle <- runif(n = n_iter, min = 0, max = 1)
    end <- runif(n = n_iter, min = 0, max = 1)
    df <- t(sapply(seq(n_iter), 
                   function(x) {
                     c(approx(x = c(1, yrs_hist/2), 
                              y = c(start[x], middle[x]), 
                              n = yrs_hist/2)$y,
                       approx(x = c(yrs_hist/2, yrs_hist + 1), 
                              y = c(middle[x], end[x]), 
                              n = (yrs_hist/2) + 1)$y[-1])
                   }))
    df2 <- as.data.frame(df)
    rownames(df2) <- seq(n_iter)
    colnames(df2) <- seq(yrs_hist)
    df2$iter <- seq(n_iter)
    #df2 %>% 
    #  gather(key = "year", value = "value", 1:100) %>%
    #  mutate(year = as.numeric(as.character(year))) %>%
    #  ggplot(aes(x = year, y = value, group = as.factor(iter))) +
    #  geom_line(alpha = 0.5) +
    #  theme_bw()
    
    f_array <- array(dim = c(yrs_hist, 3, n_iter),
                     dimnames = list(seq(yrs_hist), c("min","val","max"),
                                     iter = seq(n_iter)))
    f_array[, "val", ] <- c(t(df))
  }
  
  ######################
  
  input_params <- readRDS(inp_file)
  lh_pars1 <- input_params$lh_pars
  
  lh_pars1["sel1"] <- params[1]
  lh_pars1["sel2"] <- params[2]
  lh_pars1["sel3"] <- params[3]
  #lh_pars1["a50"] <- params[1]

  ### initiate the biological functions 
  params <- lh_pars1
  growth = FLifeE::vonB
  m      = "gislason"
  sr     = "bevholt"
  mat    = logistic
  sel    = sel_f
  range  = c(min=1,max=input_params$Amax,minfbar=2,maxfbar=10,plusgroup=input_params$Amax)
  spwn   = 0 
  fish   = 0.5 # proportion of year when fishing happens
  midyear= 0.5
  
  pNms=dimnames(params)$params
  if ("sl"%in%pNms&!("sel1"%in%pNms)) {
    dimnames(params)$params["sel1"==pNms]="sel1"
  } 
  if ("sr"%in%pNms&!("sel2"%in%pNms)) {
    dimnames(params)$params["sel2"==pNms]="sel2"
  }
  if ("a1"%in%pNms&!("sel3"%in%pNms)) {
    dimnames(params)$params["sel3"==pNms]="sel3"
  }
  
  if (("m.spwn" %in% names(args))) {
    m.spwn =args[["m.spwn"]]
  } else {
    m.spwn=FLQuant(rep(spwn,each=length(range["min"]:range["max"])), dimnames=list(age=range["min"]:range["max"],iter=seq(length(spwn))), units="")
  }
  if (("harvest.spwn" %in% names(args))) {
    harvest.spwn =args[["harvest.spwn"]]
  } else {
    harvest.spwn=FLQuant(rep(spwn,each=length(range["min"]:range["max"])), dimnames=list(age=range["min"]:range["max"],iter=seq(length(spwn))), units="")
  }
  
  age=FLQuant(range["min"]:range["max"],
              dimnames=list(age =range["min"]:range["max"],
                            iter=dimnames(params)$iter), units="")
  
  # Get the lengths through different times of the year
  slen   <- growth(age+m.spwn,params) # slen is length at spawning time
  clen   <- growth(age+fish,  params) # clen is length when fishing happens
  midyearlen <- growth(age+midyear,params)
  
  # Corresponding weights
  # bug warning cos of log(NA)
  cwt=FLifeE::len2wt(clen,params)
  swt=len2wt(slen,params)
  
  # maturity
  mat. =mat(age,params) 
  # selectivty
  sel. =sel(age,params)
  
  ## create a FLBRP object to   calculate expected equilibrium values and ref pts
  dms=dimnames(swt)
  
  FLBRP_object =  FLBRP(stock.wt       =swt,
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
  
  # gislason model
  m(FLBRP_object) =  exp(0.55 - 1.6*log(slen) + 1.44*log(lh_pars1["linf"]) + log(lh_pars1["k"])) 

  ## replace any slot passed in as an arg
  for (slt in names(args)[names(args) %in% names(getSlots("FLBRP"))[names(getSlots("FLBRP"))!="fbar"]]){
    slot(FLBRP_object, slt)<-args[[slt]]
  }
  
  params(FLBRP_object)=propagate(params(FLBRP_object),dims(FLBRP_object)$iter)
  model(FLBRP_object) =do.call(sr,list())$model
  
  params(FLBRP_object)=FLPar(c(a=as.numeric(NA),b=as.numeric(NA)),iter=dims(params)$iter)
  for (i in seq(dims(params)$iter)) {
    params(FLBRP_object)[,i][]=unlist(c(FLCore::ab(params[c("s","v"),i],sr,spr0=iter(spr0(FLBRP_object),i))[c("a","b")]))  
  }
  
  refpts(FLBRP_object)=propagate(refpts(FLBRP_object)[c("virgin","msy","crash","f0.1")],dims(params)$iter)
  FLBRP_object=brp(FLBRP_object)
  
  ### stk ###
  
  stk <- as(FLBRP_object,"FLStock")
  refpts <- refpts(FLBRP_object)
  stk <- qapply(stk, function(x) {#browser()
    dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
  })
  stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
  stk <- propagate(stk, 500)
  
  ### create stock recruitment model
  stk_sr <- FLSR(params = FLBRP_object@params, model = FLBRP_object@model)
  ### create residuals for (historical) projection
  set.seed(0)
  residuals(stk_sr) <- rlnoise(dim(stk)[6], rec(stk) %=% 0, 
                            sd = 0.6, b = 0)
  ### replicate residuals from catch rule paper for historical period
  set.seed(0)
  residuals <- rlnoise(dim(stk)[6], (rec(stk) %=% 0)[, ac(1:100)], 
                       sd = 0.6, b = 0)
  residuals(stk_sr)[, ac(1:100)] <- residuals[, ac(1:100)]

  #### fishing history
  if (isTRUE(fhist == "one-way")) {
    
    ### 0.5Fmsy until year 75, then increase to 0.8Fcrash
    fs <- sapply(c(refpts["msy","harvest"]),function(x){rep(x,74)*0.5})
    f0 <- sapply(c(refpts["msy","harvest"]),function(x){x*0.5})
    fmax <- sapply(c(refpts["crash","harvest"]),function(x){x*0.8})
    rate <- exp((log(fmax) - log(f0)) / (25))
    
    fs_array <- array(NA,dim = c(99,500))
    for (i in 1:500){
      fs_array[,i] <- c(fs[,i],rate[i]^(1:25)*f0[i])
    }
    
    arr <- array(NA,dim = c(99,3,500),dimnames = list(1:99,c("min","val","max"),iter=1:500))
    arr[,"val",] <- fs_array
    
    ctrl <- fwdControl(data.frame(year = c(2:100), 
                                  quantity = c("f"), value = NA))
    ### control object
    ctrl@trgtArray <- arr 
    
  } else if (isTRUE(fhist == "random")) {
    
    ### control object template
    ctrl <- fwdControl(data.frame(year = seq(yrs_hist), 
                                  quantity = c("f"), value = NA))
    ### add iterations
    ctrl@trgtArray <- f_array
    ### target * Fcrash
    ctrl@trgtArray [,"val",] <- ctrl@trgtArray [,"val",] * 
      c(refpts["crash", "harvest"]) * 1
  
  }  
  
  ### project fishing history
  stk_stf <- fwd(stk, ctrl, sr = stk_sr, sr.residuals = residuals(stk_sr),
                 sr.residuals.mult = TRUE, maxF = 5) 
  
  invisible(gc())
  
  return(list(stk = stk_stf, sr = stk_sr,sel = sel))
}


fitness_catch <- function(params,sel_f,catch_ipma_prop,fhist,path,inp_file){
  
  ### check for files 
  params[1] <- round(params[1],1)
  params[2:3] <- round(params[2:3])
  run_i <- paste0(params, collapse = "_")

  ### check if path exists
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  ### check if run already exists
  if (isTRUE(file.exists(paste0(path, run_i, ".rds")))) {
    ### load stats
    new_error <- readRDS(paste0(path, run_i, ".rds"))
    ### set flag for running MP
    run_mp <- FALSE
    
  } else {
    run_mp <- TRUE
  }
  
  if (isTRUE(run_mp)){
    
    new_om <- om(params,sel_f,fhist,inp_file)
    catch_means <- iterMeans(new_om$stk@catch.n)
    
    catch_prop_om <- Catch_proportion(catch_means[,86:100])
    new_error <- sum((catch_ipma_prop - catch_prop_om)^2)
    
    saveRDS(new_error,paste0(path,run_i,".rds"))
    
    if (getDoParWorkers() > 1){
      . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
    } else {
      invisible(gc())
    }
  
    return(new_error)
            
  }
  
  return(new_error)
  
}

