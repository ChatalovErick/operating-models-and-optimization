# installing the parallelpso ##############################################

# parallel pso #
#install.packages("devtools")
#library("devtools")
#remotes::install_github("https://github.com/ChatalovErick/parallelPSO.git")

# Flife with changes 
#install.packages("FLife-master_changes", repos = NULL, type = "source")

###########################################################################
###########################################################################

scenario <- "OM 5"
fhist <- "one-way"

#install.packages("C:/Users/erick/OneDrive/Ambiente de Trabalho/github-thesis V2 - Copy/FLBRP-master", repos = NULL, type="source",
#                 dependencies = TRUE, INSTALL_opts = '--no-lock')


#install.packages("C:/Users/erick/OneDrive/Ambiente de Trabalho/github-thesis V2 - Copy/TropFishR-master", repos = NULL, type="source",
#                 dependencies = TRUE, INSTALL_opts = '--no-lock')

#remotes::install_github("flr/mse")
#install.packages("FLash", repos="http://flr-project.org/R")
#################################################################
#### r file to use the Thompson and bell model to make the 
#### reference points

library(FLifeE)
library(parallelPSO)
library(parallel)
library(doParallel)
library(dplyr)
library(ggplot2)
library(mse)
library(foreach)

library(TropFishR)
library(fishdynr)
source("funs.R")

# ------------------------------------------------------------------------ #
# data 
# ------------------------------------------------------------------------ #

req_pckgs <- c("FLCore", "FLash", "mse", "GA","FLifeE", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)
library(FLBRP)

params <- read.csv("input/stock_list_full_pt.csv")[1,]
stock_lh <- read.csv("input/stock_list_full_pt.csv")[1,]
data_ipma <- read.csv("input/RJC_numbers_IC_all.csv",header=TRUE, sep = ";")
data_ipma <- data_ipma[,c(3:17)]

### ------------------------------------------------------------------------ ###
### fishing history dimensions ####
### ------------------------------------------------------------------------ ###

n_iter <- 500
yrs_hist <- 100
yrs_proj <- 50

### ------------------------------------------------------------------------ ###
### with uniform distribution and random F trajectories ####
### ------------------------------------------------------------------------ ###

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


# ------------------------------------------------------------------------ #
# Aggregate Catch data for length data using different bin_sizes  
# ------------------------------------------------------------------------ #

# add more lengths to the data with zeros until linf
m <- matrix(0,nrow=params$linf-nrow(data_ipma),ncol=ncol(data_ipma))
colnames(m) <- names(data_ipma)
data_ipma <- rbind(data_ipma,m)

# aggregate catch data for different bin sizes in catch data

aggregate_data <- function(data,bin_size,linf){
  C <- matrix(NA,nrow = linf/bin_size,ncol = ncol(data))
  interval_val <- c(seq(1,linf,bin_size),linf)
  mid_lengths <- c()
  for (len in 1:(length(interval_val)-1)){
    mid_lengths[len] <- (interval_val[len]+interval_val[len+1])/2
    C[len,] <- colSums(data[c(interval_val[len]:(interval_val[len+1]-1)),])
  }
  
  rownames(C) <- mid_lengths
  return(C)
}

data_ipma_bin_size_4 <- aggregate_data(data_ipma,4,128)

# ------------------------------------------------------------------------ #
# calculate some parameters 
# ------------------------------------------------------------------------ #

#------------------------#
# steepness 
#------------------------#

### OM scenarios (M, etc.)
### stored in csv file, load here
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
# only using the scenario 38
OM_scn <- OM_scns[38,]
OM_scn$steepness <- "l50linf"

OM_scn$selectivity <- "dnormal"

stock <- params[1,]

### handle steepness
if (!is.na(OM_scn$steepness)) {
  if (!is.na(as.numeric(OM_scn$steepness))) {
    steepness <- as.numeric(OM_scn$steepness)
  } else if (OM_scn$steepness == "l50linf") {
    l50linf <- params$l50 / params$linf
    ### calculate steepness according to Wiff et al. 2018
    x <- l50linf
    y <- 2.706 - 3.698*x
    invLogit <- function(y) (0.2 + exp(y))/(1 + exp(y))
    steepness <- invLogit(y)
  }
} else {
  steepness <- 0.75
}

### Max age: age at l = 0.95 * linf
max_age <- ceiling(log(0.05)/(-c(stock_lh$k)) + c(stock_lh$t0))

# ------------------------------------------------------------------------ #
# prepare the om, calculate the reference points
# ------------------------------------------------------------------------ #

# create the input for the models in tropFish

data <- list()
data$Linf <- params$linf
data$K <- params$k
data$t0 <- params$t0
data$Lmat <- params$l50
data$wmat <- params$l50*0.2
data$M <- params$M
data$s <- steepness

data$sample.no <- 1:nrow(data_ipma_bin_size_4)
data$midLengths <- as.numeric(rownames(data_ipma_bin_size_4))
data$catch <- data_ipma_bin_size_4

# estimate the selectivity ogive and the total mortality
res_cc <- catchCurve(data,catch_columns = c(1:15),calc_ogive=TRUE,plot=TRUE,auto = TRUE)
res_cc$wqs <- (res_cc$L75 - res_cc$L50)*2

data$a <- params$a
data$b <- params$b
data$L50 <- res_cc$L50
data$FM <- res_cc$Z - data$M
data$E <- data$FM/res_cc$Z

#--------------------------#
# selectivity 
#--------------------------#

# mortality by length group, selectivity by length
sel_CC <-  logisticSelect(data$midLengths,res_cc$L50,res_cc$wqs)

### ------------------------------------------------------------------------ ###
### create final brps ####
### ------------------------------------------------------------------------ ###

### get params

lh_res <- c(dimnames(lhPar(FLPar(linf = 1)))$params, "l50")
lh_avail <- intersect(lh_res, names(stock))
lh_pars <- stock[, lh_avail]
lh_pars <- lh_pars[, !is.na(lh_pars)]
lh_pars <- as(lh_pars, "FLPar")
lh_pars["s"] <- steepness

### create missing pars
lh_pars <- lhPar(lh_pars)

### ------------------------------------------ #
### set-up selectivity for FLBRP (if specified)
### ------------------------------------------ #
lh_pars[c("sel1","sel2","sel3")] <- c(res_cc$t50,1,5000)

### initiate the biological functions 
params <- lh_pars
growth = FLifeE::vonB
m      = "gislason"
sr     = "bevholt"
mat    = logistic
sel    = FLifeE::logistic_Sel
range  = c(min=1,max=25,minfbar=2,maxfbar=10,plusgroup=25)
spwn   = 0 
fish   =0.5 # proportion of year when fishing happens
midyear=0.5

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
sel. =sel(age,  params) 

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

m(FLBRP_object) = data$M 

## replace any slot passed in as an arg
for (slt in names(args)[names(args) %in% names(getSlots("FLBRP"))[names(getSlots("FLBRP"))!="fbar"]]){
  slot(FLBRP_object, slt)<-args[[slt]]
}

params(FLBRP_object)=propagate(params(FLBRP_object),dims(FLBRP_object)$iter)
model(FLBRP_object) =do.call(sr,list())$model
params(FLBRP_object)=FLCore::ab(params[c("s","v")],sr,spr0=spr0(FLBRP_object))[c("a","b")]

## Stock recruitment relationship
#params(FLBRP_object)[1] <- population_values$virgin_population$alpha
#params(FLBRP_object)[2] <- population_values$virgin_population$beta

refpts(FLBRP_object)=propagate(refpts(FLBRP_object)[c("virgin","msy","crash","f0.1")],dims(params)$iter)
FLBRP_object=brp(FLBRP_object)

refpts(FLBRP_object)
##############################################################################
##############################################################################

### save life-history parameters in FLBRP
attr(FLBRP_object, "lhpar") <- lh_pars
attr(FLBRP_object, "OM_scn") <- OM_scn
attr(FLBRP_object, "lhpar_calc") <- FLPar(
  M = mean(m(FLBRP_object)),
  M_mat = weighted.mean(x = m(FLBRP_object), w = FLBRP_object@mat),
  MK = weighted.mean(x = m(FLBRP_object), w = FLBRP_object@mat) / c(lh_pars["k"]),
  amax = dims(FLBRP_object)$plusgroup
)

name(FLBRP_object) <- stock$stock
desc(FLBRP_object) <- stock$stock

### ------------------------------------------------------------------------ ###
### create OM ####
### ------------------------------------------------------------------------ ###

# Blim value 
attr(FLBRP_object, "Blim") <- 1000*0.163 

# convert brp into FLstock() and initialize a population
stk <- as(FLBRP_object, "FLStock")
refpts <- refpts(FLBRP_object)

stk <- qapply(stk, function(x) {#browser()
  dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
})

stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
stk <- propagate(stk, n_iter)

### create stock recruitment model
stk_sr <- FLSR(params = params(FLBRP_object), model = model(FLBRP_object))

### create residuals for (historical) projection
set.seed(0)
residuals(stk_sr) <- rlnoise(dim(stk)[6], rec(stk) %=% 0, 
                             sd = 0.6, b = 0)
### replicate residuals from catch rule paper for historical period
set.seed(0)
residuals <- rlnoise(dim(stk)[6], (rec(stk) %=% 0)[, ac(1:100)], 
                     sd = 0.6, b = 0)
residuals(stk_sr)[, ac(1:100)] <- residuals[, ac(1:100)]

### fishing history
if (isTRUE(fhist == "one-way")) {
  
  ### 0.5Fmsy until year 75, then increase to 0.8Fcrash
  fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 74)
  f0 <- c(refpts["msy", "harvest"]) * 0.5
  fmax <- c(refpts["crash", "harvest"]) * 0.8
  rate <- exp((log(fmax) - log(f0)) / (25))
  fs <- c(fs, rate ^ (1:25) * f0)
  
  ### control object
  ctrl <- FLash::fwdControl(data.frame(year = 2:100,val = fs,quantity = c("f")))
  
  ### random F trajectories
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


path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", fhist, "/")
dir.create(path, recursive = TRUE)
saveRDS(list(stk = stk_stf, sr = stk_sr), file = paste0(path,stock$stock,"_",scenario,".rds"))

### ------------------------------------------------------------------------ ###
### prepare OMs for flr/mse MP ####
### ------------------------------------------------------------------------ ###

stock <- "rjc_1"
stocks <- stock_lh

### load stock
tmp <- readRDS(paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", fhist,
                      "/", stock,"_",scenario,".rds"))
stk_fwd <- tmp$stk
stk_sr <- tmp$sr
### life-history data
lhist <- stocks[stocks$stock == stock, ]
#range(stk_stf)
### cut of history
stk_fwd <- window(stk_fwd, start = 50)
stk_sr@residuals <- window(stk_sr@residuals, start = 50)
### length data
pars_l <- FLPar(a = lhist$a,
                b = lhist$b,
                Lc = calc_lc(stk = stk_fwd[, ac(75:100)], 
                             a = lhist$a, b = lhist$b))

### indices
#q <- 1/(1 + exp(-1*(an(dimnames(stk_fwd)$age) - dims(stk_fwd)$max/10)))
sel_dnormal <- c(FLifeE::dnormal(FLQuant(1:max_age,dimnames=list("age"=1:max_age)),lh_pars))

idx <- FLQuants(
  sel = stk_fwd@mat %=% sel_dnormal,
  idxB = quantSums(stk_fwd@stock.n * stk_fwd@stock.wt * (stk_fwd@mat)),
  idxL = lmean(stk = stk_fwd, params = pars_l),
  PA_status = ssb(stk_fwd) %=% NA_integer_)

### index deviation
PA_status_dev <- FLQuant(NA, dimnames = list(age = c("positive", "negative"), 
                                             year = dimnames(stk_fwd)$year, 
                                             iter = dimnames(stk_fwd)$iter))
set.seed(1)
PA_status_dev["positive"] <- rbinom(n = PA_status_dev["positive"], 
                                    size = 1, prob = 0.9886215)
set.seed(2)
PA_status_dev["negative"] <- rbinom(n = PA_status_dev["negative"], 
                                    size = 1, prob = 1 - 0.4216946)
set.seed(696)
idx_dev <- FLQuants(sel = stk_fwd@mat %=% 1,
                    idxB = rlnoise(n = dims(idx$idxB)$iter, idx$idxB %=% 0, 
                                   sd = 0.2, b = 0),
                    idxL = rlnoise(n = dims(idx$idxL)$iter, idx$idxL %=% 0, 
                                   sd = 0.2, b = 0),
                    PA_status = PA_status_dev)
### iem deviation
set.seed(205)
iem_dev <- FLQuant(rlnoise(n = dims(stk_fwd)$iter,  catch(stk_fwd) %=% 0,
                           sd = 0.1, b = 0))
### lowest observed index in last 50 years
I_loss <- list()
I_loss$SSB_idx <- apply(ssb(stk_fwd)[, ac(50:100)], 6, min)
I_loss$SSB_idx_dev <- apply((ssb(stk_fwd) * idx_dev$idxB)[, ac(50:100)], 
                            6, min)
I_loss$idx <- apply(idx$idxB[, ac(50:100)], 6, min)
I_loss$idx_dev <- apply((idx$idxB * idx_dev$idxB)[, ac(50:100)], 6, min)

### parameters for components
pars_est <- list(
  comp_r = TRUE, comp_f = TRUE, comp_b = TRUE,
  comp_c = TRUE, comp_m = 1,
  idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3, idxB_range_3 = 1,
  catch_lag = 1, catch_range = 1,
  interval = 2,
  idxL_lag = 1, idxL_range = 1,
  exp_r = 1, exp_f = 1, exp_b = 1,
  Lref = rep((lhist$linf + 2*1.5*c(pars_l["Lc"])) / (1 + 2*1.5), n_iter),
  B_lim = rep(FLBRP_object@Blim, n_iter),
  I_trigger = c(I_loss$idx_dev * 1.4), ### default, can be overwritten later
  pa_buffer = FALSE, pa_size = 0.8, pa_duration = 3,
  upper_constraint = Inf,
  lower_constraint = 0
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = stk_sr, ### stock recruitment and precompiled residuals
           fleetBehaviour = mseCtrl(),
           projection = mseCtrl(method = fwd_attr,
                                args = list(dupl_trgt = TRUE)))
tracking = c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
             "multiplier", "exp_r", "exp_f", "exp_b")
oem <- FLoem(method = obs_generic,
             observations = list(stk = stk_fwd, idx = idx), 
             deviances = list(stk = FLQuant(), idx = idx_dev),
             args = list(idx_dev = TRUE, ssb = FALSE,
                         lngth = TRUE, lngth_dev = TRUE,
                         lngth_par = pars_l,
                         PA_status = FALSE, PA_status_dev = FALSE,
                         PA_Bmsy = c(refpts(FLBRP_object)["msy", "ssb"]), 
                         PA_Fmsy = c(refpts(FLBRP_object)["msy", "harvest"])))

ctrl <- mpCtrl(list(
  est = mseCtrl(method = est_comps,
                args = pars_est),
  phcr = mseCtrl(method = phcr_comps,
                 args = pars_est),
  hcr = mseCtrl(method = hcr_comps,
                args = pars_est),
  isys = mseCtrl(method = is_comps,
                 args = pars_est)
))
iem <- FLiem(method = iem_comps,
             args = list(use_dev = FALSE, iem_dev = iem_dev))
### args
args <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
             y0 = dims(stk_fwd)$minyear, ### first data year
             iy = 100, ### first simulation (intermediate) year
             nsqy = 3, ### not used, but has to provided
             nblocks = 1, ### block for parallel processing
             seed = 1, ### random number seed before starting MSE
             seed_part = FALSE
)
### get reference points
refpts <- refpts(FLBRP_object)
Blim <- attr(FLBRP_object, "Blim")

### list with input to mp()
input <- list(om = om, oem = oem, iem = iem, ctrl = ctrl, 
              args = args,
              scenario = "GA", tracking = tracking, 
              verbose = TRUE,
              refpts = refpts, Blim = Blim, I_loss = I_loss)

### save OM
path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_2_mp_input/", fhist, "/")
dir.create(path, recursive = TRUE)
saveRDS(object = input, file = paste0(path, stock,"_",scenario,".rds"))
