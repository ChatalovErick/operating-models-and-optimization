# installing the parallelpso ##############################################

# parallel pso #
#install.packages("devtools")
#library("devtools")
#remotes::install_github("https://github.com/ChatalovErick/parallelPSO.git")

# Flife with changes 
#install.packages("FLife-master_changes", repos = NULL, type = "source")

###########################################################################
###########################################################################

scenario <- "OM 6"
fhist <- "one-way"

library(FLifeE)
library(parallelPSO)
library(parallel)
library(doParallel)
library(dplyr)
library(ggplot2)
library(mse)
library(foreach)
library(FLa4a)
source("funs.R")

################################################################
################################################################
req_pckgs <- c("FLCore", "FLash", "mse", "GA","FLifeE", "doParallel", "doRNG", "FLBRP","FLa4a")
for (i in req_pckgs) library(package = i, character.only = TRUE)

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

########################################################
########################################################
### specie and list of parameters for some of the variables.

stocks_lh <- read.csv("input/stock_list_full_pt.csv")
stocks_lh[1,]
stocks_lh <- stocks_lh[!is.na(stocks_lh$a), ]
names(stocks_lh)[1] <- "name"

### use lmax as proxy for linf, if linf not provided but lmax exits
pos_lmax <- which(!is.na(stocks_lh$lmax) & is.na(stocks_lh$linf))
stocks_lh$linf[pos_lmax] <- stocks_lh$lmax[pos_lmax]

### change name for stocks that appear twice
stocks_lh$stock <- as.character(stocks_lh$stock)
if (nrow(stocks_lh) > 15) {
  stocks_lh$stock[29] <- paste0(stocks_lh$stock[29], "_2")
}

### set fbar range
### only used for default OMs
stocks_lh$minfbar <- 2
stocks_lh$maxfbar <- 10

### OM scenarios (M, etc.)
### stored in csv file, load here
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
# only using the scenario 38
OM_scn <- OM_scns[38,]
OM_scn$selectivity <- "dnormal"

stock <- stocks_lh[1,]

########################
#### steepness ####
OM_scn$steepness <- "l50linf"
#########################

# The conversion of length data to age is performed through a growth model.
# The implementation is done through the a4aGr class.

vbObj <- a4aGr(
  grMod=~linf*(1-exp(-k*(t-t0))),      
  grInvMod=~t0-1/k*log(1-len/linf),      
  params= FLPar(linf=128, k=0.1170, t0=-0.6170))

# -------------------------------------------------------------------------------- #
# Adding uncertainty to growth parameters with a multivariate normal distribution
# -------------------------------------------------------------------------------- #

# Make an empty cor matrix
cm <- diag(c(1,1,1))
# k and linf are negatively correlated while t0 is independent
cm[1,2] <- cm[2,1] <- -0.5
# scale cor to var using CV=0.05
cv <- 0.05
p <- c(linf=128, k=0.1170, t0=-0.6170)
vc <- matrix(1, ncol=3, nrow=3)
l <- vc
l[1,] <- l[,1] <- p[1]*cv
k <- vc
k[,2] <- k[2,] <- p[2]*cv
t <- vc
t[3,] <- t[,3] <- p[3]*cv
mm <- t*k*l
diag(mm) <- diag(mm)^2
mm <- mm*cm
# check that we have the intended correlation
all.equal(cm, cov2cor(mm))

# Create the a4aGr object as before, but now we also include the
# vcov argument for the variance-covariance matrix.

vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), 
               params=FLPar(linf=p['linf'], k=p['k'], t0=p['t0']), vcov=mm)


# We simulate 1000 iterations from the a4aGr object by calling mvrnorm() using the 
# variance-covariance matrix we created earlier.
set.seed(0)
vbNorm <- mvrnorm(500,vbObj)

# Now we have 1000 iterations of each parameter, randomly sampled from the 
# multivariate normal distribution

vbNorm@params$l50 <- array(stock$l50,dim = c(1,500))
vbNorm@params$a50 <- array(stock$a50,dim = c(1,500))
vbNorm@params$a <- array(stock$a,dim = c(1,500))
vbNorm@params$b <- array(stock$b,dim = c(1,500))

lh_pars <- lhPar(FLPar(vbNorm@params,iter=500))
lh_pars["sel1"] <- array(7.5,dim = c(1,500))

### get params
### handle steepness
if (!is.na(OM_scn$steepness)) {
  if (!is.na(as.numeric(OM_scn$steepness))) {
    lh_pars["s"] <- as.numeric(OM_scn$steepness)
  } else if (OM_scn$steepness == "l50linf") {
    l50linf <- lh_pars["l50"] / lh_pars["linf"]
    ### calculate steepness according to Wiff et al. 2018
    x <- l50linf
    y <- 2.706 - 3.698*x
    invLogit <- function(y) (0.2 + exp(y))/(1 + exp(y))
    lh_pars["s"] <- invLogit(y)
  }
} else {
  lh_pars["s"] <- 0.75
}


### Max age: age at l = 0.95 * linf
max_age <- ceiling(log(0.05)/(-c(stock$k)) + c(stock$t0))

lh_pars1 <- list(lh_pars=lh_pars,Amax=max_age)

### initiate the biological functions 
params <- lh_pars
growth = FLifeE::vonB
m      = "gislason"
sr     = "bevholt"
mat    = logistic
sel    = dnormal
range  = c(min=1,max=25,minfbar=2,maxfbar=10,plusgroup=25)
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

# Get the lengths through different times of the year for the 500 iter
length_at_age   <- growth(age+m.spwn,params) # slen is length at spawning time

# lenght-weight 
swt <- len2wt(length_at_age,params)

params <- lh_pars

##############################################################################
##############################################################################

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
m(FLBRP_object) =  exp(0.55 - 1.6*log(slen) + 1.44*log(params["linf"]) + log(params["k"])) 

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

### plots from brp ###
par(mfrow=c(3,2))
hist(c(refpts(FLBRP_object)["msy","harvest"]),main = "Fmsy")
hist(c(refpts(FLBRP_object)["crash","harvest"]),main = "Fcrash")
hist(c(refpts(FLBRP_object)["msy","yield"]),main = "msy yield")
hist(c(refpts(FLBRP_object)["msy","rec"]),main = "msy rec")
hist(c(refpts(FLBRP_object)["msy","ssb"]),main="msy ssb")
hist(c(refpts(FLBRP_object)["msy","biomass"]),main="msy biomass")
######################

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
stk_stf <- fwd(stk, ctrl, sr = stk_sr,sr.residuals = residuals(stk_sr),
               sr.residuals.mult = TRUE, maxF = 5)

path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", fhist, "/")
dir.create(path, recursive = TRUE)
saveRDS(list(stk = stk_stf, sr = stk_sr), file = paste0(path,stock$stock,"_",scenario,".rds"))

### ------------------------------------------------------------------------ ###
### prepare OMs for flr/mse MP ####
### ------------------------------------------------------------------------ ###

stock <- "rjc_1"
stocks <- stocks_lh

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

