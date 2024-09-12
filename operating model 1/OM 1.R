####################################################################################
####################################################################################
# R file used to create all the necessary data, input for the optimization process.
####################################################################################
####################################################################################

###############################################################################
########### installing packages                              ##################
###############################################################################

devtools::install_github(repo = "flr/FLCore", ref = "d55bc6570c0134c6bea6c3fc44be20378691e042")
devtools::install_github(repo = "flr/FLash", ref = "7c47560cf57627068259404bb553f2b644682726")
devtools::install_github(repo = "flr/FLBRP", ref = "142d5e14137c5ceb4526afd6718c26269ad81e7c")
devtools::install_github(repo = "flr/ggplotFL", ref = "9b502a1aa01524637f4f269a3353a92c7d452db0")
#devtools::install_github(repo = "flr/FLife", ref = "d0cca5e574a77fb52ec607a25c244969b9c8dd38")
devtools::install_github(repo = "shfischer/mse", ref = "80b5cf18dc9611f7307f599564ccdfbad433948d", INSTALL_opts = "--no-multiarch")
install.packages("FLasher", repos="http://flr-project.org/R")
install.packages("FLAssess", repos="http://flr-project.org/R")
install.packages("doParallel")
remotes::install_github("shfischer/mse", ref = "mseDL2.0")
install.packages("foreach")
remotes::install_github(repo = "shfischer/GA")

# Flife with changes package #
install.packages("FLife-master_changes", repos = NULL, type = "source")

###############################################################################
###############################################################################
###############################################################################

library(FLife)
library(FLBRP)
#library(FLasher)
library(foreach)
library(doParallel)
#library(FLAssess)
library(mse)
library(foreach)
library(FLCore)
#library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(FLash)

source("funs.R")
source("funs_GA.R")
### set up cluster for parallel computing
cl <- makeCluster(1)
registerDoParallel(cl)
### ------------------------------------------------------------------------ ###
### fishing history dimensions ####
### ------------------------------------------------------------------------ ###
n_iter <- 500
yrs_hist <- 100
yrs_proj <- 50
scenario <- "scenario_2.1"
fhist <- "one-way" #"one-way"
### ------------------------------------------------------------------------ ###
### with uniform distribution and random F trajectories ####
### ------------------------------------------------------------------------ ###
fhist <- "random"
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
  #df2$iter <- seq(n_iter)
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

f_array

### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###
?FLife


### specie and list of parameters for some of the variables.
stocks_lh <- read.csv("input/stock_list_full_pt.csv")
stocks_lh[1,]
stocks_lh <- stocks_lh[!is.na(stocks_lh$a), ]
names(stocks_lh)[1] <- "name"

### use lmax as proxy for linf, if linf not provided but lmax exits
pos_lmax <- which(!is.na(stocks_lh$lmax) & is.na(stocks_lh$linf))
stocks_lh$linf[pos_lmax] <- stocks_lh$lmax[pos_lmax]

### set steepness to 0.75
stocks_lh$s <- 0.75

### change name for stocks that appear twice
stocks_lh$stock <- as.character(stocks_lh$stock)
if (nrow(stocks_lh) > 15) {
  stocks_lh$stock[29] <- paste0(stocks_lh$stock[29], "_2")
}

### OM scenarios (M, etc.)
### stored in csv file, load here
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
# only using the scenario 38
OM_scns <- OM_scns[38,]

### ------------------------------------------------------------------------ ###
### create final brps ####
### ------------------------------------------------------------------------ ###

i <-  stocks_lh[1,]
OM_scn <- OM_scns[38,]
  
brps <- foreach(OM_scn = split(OM_scns, 1:nrow(OM_scns)),
                .final = function(x) {
                  names(x) <- OM_scns$id
                  return(x)
                }) %:%
  foreach(i = split(stocks_lh, 1:nrow(stocks_lh)),
          .errorhandling = "pass", .packages = c("FLife", "FLBRP"),
          .final = function(x) {
            names(x) <- stocks_lh$stock
            return(x)
          }) %dopar% {
            
            ### create brp
            ### get params
            lh_res <- c(dimnames(lhPar(FLPar(linf = 1)))$params, "l50")
            lh_avail <- intersect(lh_res, names(i))
            lh_pars <- i[, lh_avail]
            lh_pars <- lh_pars[, !is.na(lh_pars)]
            lh_pars <- as(lh_pars, "FLPar")
            ### set default t0, if missing
            if (!"t0" %in% dimnames(lh_pars)$params) {
              lh_pars <- rbind(lh_pars, FLPar(t0 = -0.1))
            }
            ### calculate l50, if missing and a50 exists
            if (!"l50" %in% dimnames(lh_pars)$params & 
                all(c("a50", "k", "linf", "t0") %in% dimnames(lh_pars)$params)) {
              lh_pars <- rbind(lh_pars, 
                               FLPar(l50 = vonB(age = c(lh_pars["a50"]), 
                                                params = lh_pars)))
            }
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
              } else if (OM_scn$steepness == "k") {
                ### h dependent on k...
                dat <- data.frame(ks = c(min(stocks_lh$k), max(stocks_lh$k)),
                                  hs = c(0.5, 0.9))
                model <- lm(hs ~ ks, data = dat)
                new_s <- predict(object = model, new = data.frame(ks = c(lh_pars["k"])))
                lh_pars["s"] <- new_s
              } else if (OM_scn$steepness == "Myers") {
                ### h from Myers et al. 1999
                lh_pars["s"] <- ifelse(!is.na(i$h_Myers), i$h_Myers, c(lh_pars["s"]))
              }
            }
            
            ### create missing pars
            lh_pars <- lhPar(lh_pars)
            
            ### Max age: age at l = 0.95 * linf
            max_age <- ceiling(log(0.05)/(-c(lh_pars["k"])) + c(lh_pars["t0"]))
            ### max from F
            ### alternative for Lorenzen mortality
            if (OM_scn$m == "lorenzen") {
              ### age with 5% survival
              max_ageP <- ceiling(-log(0.05)/(OM_scn$M1))
              if (max_age < max_ageP) max_age <- max_ageP
            }
            
            ### set up M
            ### either return function or vector with values
            m_def <- if (OM_scn$m == "gislason") {
              ### hardcode Gislason natural mortality function here
              ### used to be default in FLife
              function(length,params) {
                exp(0.55)*(length^-1.61) %*% (params["linf"]^1.44) %*% params["k"]
              }
            } else if (OM_scn$m == "lorenzen") {
              
              ages <- 1:max_age
              M1 <- OM_scn$M1
              M2 <- OM_scn$M2
              age_M1 <- 1
              age_M2 <- 20
              
              a <- exp(log((vonB(age_M2, lh_pars))^log(M2) / 
                             (vonB(age_M1, lh_pars))^log(M1)) / 
                         log((vonB(age_M2, lh_pars)) / (vonB(age_M1, lh_pars))))
              b <- log(M1 / M2) /
                log( (vonB(age_M2, lh_pars)) / (vonB(age_M1, lh_pars)))
              a * vonB(ages, lh_pars)^b
              
            } else if (is.numeric(as.numeric(as.character(OM_scn$m)))) {
              as.numeric(as.character(OM_scn$m))
            }
            
            ### set-up selectivity (if specified)
            if (!is.na(OM_scn$selectivity)) {
              if (isTRUE(OM_scn$selectivity == "before")) {
                ### move maturity curve to the left 
                ### 2x the difference between a50 and a95
                sel_def <- function(age, params) {
                  ### mimic age definition as used in maturity function
                  new_age <- age - 0.5 + c(params["a50"] - floor(params["a50"]) +
                                             2 * params["ato95"])
                  sel. <- FLife::logistic(new_age, params)
                  return(sel.)
                }
              } else if (isTRUE(OM_scn$selectivity == "maturity")) {
                sel_def <- function(age, params) {
                  ### mimic age definition as used in maturity function
                  new_age <- age - 0.5 + c(params["a50"] - floor(params["a50"]))
                  sel. <- FLife::logistic(new_age, params)
                  return(sel.)
                }
              } else if (isTRUE(OM_scn$selectivity == "after")) {
                ### move maturity curve to the right 
                ### 2x the difference between a50 and a95
                sel_def <- function(age, params) {
                  ### mimic age definition as used in maturity function
                  new_age <- age - 0.5 + c(params["a50"] - floor(params["a50"]) -
                                             2 * params["ato95"])
                  sel. <- FLife::logistic(new_age, params)
                  return(sel.)
                }
              } else if (isTRUE(grepl(x = OM_scn$selectivity, pattern = "shift_"))) {
                ### shift selectivity curve by ages:
                shift <- an(strsplit(x = OM_scn$selectivity, split = "_")[[1]][2])
                sel_def <- function(age, params) {
                  age_new <- age - shift
                  sel. <- FLife::dnormal(age = age_new, params = params)
                  return(sel.)
                }
              }
            } else {
              sel_def <- FLife::dnormal
            }
            
            ### fbar range: 
            ### for default OMs: use manual definition 
            ### for modified OMs: use age at full selection
            if (isTRUE(!OM_scn$idSEQ %in% c(2:16))) {
              fbar_age <- c(i$minfbar, i$maxfbar)
            } else {
              fbar_age <- rep(c(round(lh_pars["a1"])), 2)
            }
            
            ### create brp
            brp <- lhEql(lh_pars, range = c(min = 1, max = max_age, 
                                            minfbar = fbar_age[1], maxfbar = fbar_age[2], 
                                            plusgroup = max_age), 
                         sel = sel_def)

          
            ### save life-history parameters in FLBRP
            attr(brp, "lhpar") <- lh_pars
            attr(brp, "OM_scn") <- OM_scn
            attr(brp, "lhpar_calc") <- FLPar(
              M = mean(m(brp)),
              M_mat = weighted.mean(x = m(brp), w = mat(brp)),
              MK = weighted.mean(x = m(brp), w = mat(brp)) / c(lh_pars["k"]),
              amax = dims(brp)$plusgroup
            )
            
            name(brp) <- i$stock
            desc(brp) <- i$stock
            
            return(brp)
            
          }

### ------------------------------------------------------------------------ ###
### extract the reference points ####
### ------------------------------------------------------------------------ ###

### extract and save refpts
refpts <- lapply(brps, function(x) {
  lapply(x, refpts)
})
  
refpts_tmp <- lapply(brps, function(x) lapply(x, refpts))
saveRDS(refpts_tmp, file = paste0("input/refpts_",scenario,".rds"))

### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

# Blim value 
attr(brps$new_baseline$rjc_1, "Blim") <- 1000*0.163 

stocks_subset <- stocks_lh$stock[1]
stock <- "rjc_1"

stk <- as(brp,"FLStock")
refpts <- refpts(brp)
stk_sr <- FLSR(params = params(brp), model = model(brp))

stks_hist <- foreach(stock = stocks_subset, .errorhandling = "pass", 
                     .packages = c("FLCore", "FLash", "FLBRP")) %dopar% {
  
   # convert brp into FLstock() and initialize a population
   stk <- as(brps$new_baseline[[stock]], "FLStock")
   refpts <- refpts(brps$new_baseline[[stock]])
   
   stk <- qapply(stk, function(x) {#browser()
     dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
   })
   
   stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
   stk <- propagate(stk, n_iter)
   
   ### create stock recruitment model
   stk_sr <- FLSR(params = params(brps$new_baseline[[stock]]), model = model(brps$new_baseline[[stock]]))
   
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
   

   #plot(stk_stf, iter = 1:50)
   #plot(ssb(stk_stf), iter = 1:50)
   ### run a few times to get closer to target
   # for (i in 1:5) {
   #   stk_stf <- fwd(stk_stf, ctrl, sr = stk_sr,
   #                  sr.residuals.mult = TRUE, maxF = 4)
   # }
   
   name(stk_stf) <- stock
   path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", fhist, "/")
   dir.create(path, recursive = TRUE)
   saveRDS(list(stk = stk_stf, sr = stk_sr), file = paste0(path, stock,"_",scenario, ".rds"))
   
   #return(NULL)
   return(list(stk = stk_stf, sr = stk_sr))
}


### ------------------------------------------------------------------------ ###
### prepare OMs for flr/mse MP ####
### ------------------------------------------------------------------------ ###

stock <- "rjc_1"
brps <- brps$new_baseline
stocks <- stocks_lh

stks_mp <- foreach(stock = stocks_subset, .errorhandling = "pass", 
                   .packages = c("FLCore", "mse")) %do% {
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
   q <- 1/(1 + exp(-1*(an(dimnames(stk_fwd)$age) - dims(stk_fwd)$max/10)))
   idx <- FLQuants(
     sel = stk_fwd@mat %=% q,
     idxB = quantSums(stk_fwd@stock.n * stk_fwd@stock.wt * (stk_fwd@mat %=% q)),
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
     B_lim = rep(brps[[stock]]@Blim, n_iter),
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
                            PA_Bmsy = c(refpts(brps[[stock]])["msy", "ssb"]), 
                            PA_Fmsy = c(refpts(brps[[stock]])["msy", "harvest"])))
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
   refpts <- refpts(brps[[stock]])
   Blim <- attr(brps[[stock]], "Blim")
   
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
   return(input)
}

