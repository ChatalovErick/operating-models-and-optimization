##########################################################################################

# Flife with changes 
#install.packages("FLife-master_changes", repos = NULL, type = "source")

##########################################################################################
##########################################################################################

scenario <- "OM 2"
fhist <- "one-way"

library(FLifeE)
library(parallelPSO)
library(parallel)
library(doParallel)
library(dplyr)
library(ggplot2)
library(mse)
library(foreach)
source("funs.R")

################################################################
################################################################
req_pckgs <- c("FLCore", "FLash", "mse", "GA","FLifeE", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)
library(FLBRP)

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

### ------------------------------------------------------------------------ ###
### create final brps ####
### ------------------------------------------------------------------------ ###

### create brp
### get params

lh_res <- c(dimnames(lhPar(FLPar(linf = 1)))$params, "l50")
lh_avail <- intersect(lh_res, names(stock))
lh_pars <- stock[, lh_avail]
lh_pars <- lh_pars[, !is.na(lh_pars)]
lh_pars <- as(lh_pars, "FLPar")

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
 
### create missing pars
lh_pars <- lhPar(lh_pars)
lh_pars["sel1"] <- 7.5

##############################################################################
##############################################################################

### Max age: age at l = 0.95 * linf
max_age <- ceiling(log(0.05)/(-c(lh_pars["k"])) + c(lh_pars["t0"]))

sel_def <- FLifeE::dnormal
### create brp
brp <- lhEql(lh_pars, range = c(min = 1, max = max_age, 
                                minfbar = 2, maxfbar = 10, 
                                plusgroup = max_age), 
             m = "gislason", sel = sel_def)

### save life-history parameters in FLBRP
attr(brp, "lhpar") <- lh_pars
attr(brp, "OM_scn") <- OM_scn
attr(brp, "lhpar_calc") <- FLPar(
  M = mean(m(brp)),
  M_mat = weighted.mean(x = m(brp), w = mat(brp)),
  MK = weighted.mean(x = m(brp), w = mat(brp)) / c(lh_pars["k"]),
  amax = dims(brp)$plusgroup
)

name(brp) <- stock$stock
desc(brp) <- stock$stock

### ------------------------------------------------------------------------ ###
### extract the reference points ####
### ------------------------------------------------------------------------ ###

### extract and save refpts
refpts <- refpts(brp)
saveRDS(refpts, file = paste0("input/refpts_",scenario,".rds"))

### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

# Blim value 
attr(brp, "Blim") <- 1000*0.163 

# convert brp into FLstock() and initialize a population
stk <- as(brp, "FLStock")
refpts <- refpts(brp)

stk <- qapply(stk, function(x) {#browser()
  dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
})

stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
stk <- propagate(stk, n_iter)

### create stock recruitment model
stk_sr <- FLSR(params = params(brp), model = model(brp))


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
  ctrl <- fwdControl(data.frame(year = 2:100, quantity = "f", val = fs))
  
  ### random F trajectories
} else if (isTRUE(fhist == "random")) {
  
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

path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", fhist, "/")
dir.create(path, recursive = TRUE)
saveRDS(list(stk = stk_stf, sr = stk_sr), file = paste0(path, stock,"_",scenario, ".rds"))

### ------------------------------------------------------------------------ ###
### prepare OMs for flr/mse MP ####
### ------------------------------------------------------------------------ ###
stock <- "rjc_1"
stocks <- stocks_lh[1,]

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
  B_lim = rep(brp@Blim, n_iter),
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
                         PA_Bmsy = c(refpts(brp)["msy", "ssb"]), 
                         PA_Fmsy = c(refpts(brp)["msy", "harvest"])))
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
refpts <- refpts(brp)
Blim <- attr(brp, "Blim")

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
