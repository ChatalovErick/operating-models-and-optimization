
cumulative_catch <- function(data){
  catch_prop <- apply(data,2,function(x){res <- x/sum(x)})
  catch_cumulative <- apply(catch_prop,2,function(x){res <- cumsum(x)})
  return(catch_cumulative)
  
}

Catch_proportion <- function(data){
  catch_prop <- apply(data,2,function(x){res <- x/sum(x)})
  return(catch_prop)
}

om <- function(params,sel_f,fhist){
  
  n_iter <- 500
  yrs_hist <- 100
  yrs_proj <- 50
  
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

  ##############################################################
  ############### parameters for the specie ##################
  ### extended list
  stocks_lh1 <- read.csv("input/stock_list_full_pt.csv")
  ### subset to non-NA stocks
  stocks_lh1 <- stocks_lh1[!is.na(stocks_lh1$a), ]
  names(stocks_lh1)[1] <- "name"
  
  ### set fbar range
  ### only used for default OMs
  stocks_lh1$minfbar <- 2
  stocks_lh1$maxfbar <- 10
  
  ### OM scenarios (M, etc.)
  ### stored in csv file, load here
  OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
  OM_scn <- OM_scns[38,];OM_scn
  stock <- stocks_lh1[1,]

  ### create brp
  ### get lh params
  lh_res <- c(dimnames(lhPar(FLPar(linf = 1)))$params, "l50");lh_res
  lh_avail <- intersect(lh_res, names(stock))
  lh_pars1 <- stock[, lh_avail]
  lh_pars1 <- lh_pars1[, !is.na(lh_pars1)]
  lh_pars1 <- as(lh_pars1, "FLPar")
  
  l50linf <- lh_pars1["l50"] / lh_pars1["linf"]
  ### calculate steepness according to Wiff et al. 2018
  x <- l50linf
  y <- 2.706 - 3.698*x
  invLogit <- function(y) (0.2 + exp(y))/(1 + exp(y))
  lh_pars1["s"] <- invLogit(y)
  # spr0 #
  #lh_pars1$spr0 <- 20
  ### create missing pars ###
  lh_pars1 <- lhPar(lh_pars1); lh_pars1
  # a1  selectivity #
  #lh_pars1["a1"] <- -(log(1-(pop$l50/linf))/0.117)-0.617
  
  lh_pars1["sel1"] <- params[1]
  lh_pars1["sel2"] <- params[2]
  lh_pars1["sel3"] <- params[3]
  #lh_pars1["a50"] <- params[1]
  ###########################
  
  ############# brp fucntion to initialize the population ######
  # for fixed parameters #
  
  ##############################################################
  #####   stock recruitment model   #####
  ####### S0 = 1000, R0 = 0.3276 #########
  #beta <- 1000*0.2*(0.7568-1)/(0.2-0.7568);beta
  #alpha <- -0.8*0.7568*0.3276/(0.2-0.7568);alpha
  
  brp <- lhEql(lh_pars1, range = c(min = 1, max = 25, 
                                   minfbar = 1, maxfbar = 10, plusgroup = 25), mat = logistic,sel = sel_f,
               growth=vonB,m="gislason",sr="bevholt",fish=0.5,midyear=0.5)
  
  
  sel <- do.call(sel_f,list(FLQuant(1:25),dimnames= list("age"=25),lh_pars1))
  
  stk <- as(brp,"FLStock")
  refpts <- refpts(brp)
  stk <- qapply(stk, function(x) {#browser()
    dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
  })
  stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
  stk <- propagate(stk, 500)
  
  ### create stock recruitment model
  stk_sr <- FLSR(params = brp@params, model = brp@model)
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
    fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 74)
    f0 <- c(refpts["msy", "harvest"]) * 0.5
    fmax <- c(refpts["crash", "harvest"]) * 0.8
    rate <- exp((log(fmax) - log(f0)) / (25))
    fs <- c(fs, rate ^ (1:25) * f0)
    
    ### control object
    ctrl <- fwdControl(data.frame(year = 2:100, quantity = "f", val = fs))  
    
  } else {
    
    ### control object template
    ctrl <- fwdControl(data.frame(year = seq(yrs_hist), 
                                  quantity = c("f"), val = NA))
    ### add iterations
    ctrl@trgtArray <- f_array
    ### target * Fcrash
    ctrl@trgtArray[,"val",] <- ctrl@trgtArray[,"val",] * 
      c(refpts["crash", "harvest"]) * 1
    
  }
  
  ### project fishing history
  stk_stf <- fwd(stk, ctrl, sr = stk_sr, sr.residuals = residuals(stk_sr),
                 sr.residuals.mult = TRUE, maxF = 5) 
  
  invisible(gc())
  
  return(list(stk = stk_stf, sr = stk_sr,sel = sel))
}


fitness_catch <- function(params,sel_f,cumsum,ipma_catch_age,fhist,path){
  
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
    
    new_om <- om(params,sel_f,fhist)
    catch_means <- matrix(iterMedians(new_om$stk@catch.n)[,1:100],nrow=25,ncol=100)
    print(c(catch_means[,99]))
    
    if (isTRUE(cumsum)){
      catch_cumsum_prop_om <- cumulative_catch(catch_means[,86:100])
      catch_ipma_cumsum_prop <- cumulative_catch(ipma_catch_age)
      new_error <- sum((catch_ipma_cumsum_prop - catch_cumsum_prop_om)^2)

    } else { 
      catch_prop_om <- Catch_proportion(catch_means[,86:100])
      catch_ipma_prop <- Catch_proportion(ipma_catch_age)
      new_error <- sum((catch_ipma_prop - catch_prop_om)^2)
    }
    
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

