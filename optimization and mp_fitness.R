args = c("use_MPI=FALSE", "n_workers=3", "n_blocks=3", "maxiter=10", "popSize=20", "run=10", "stock_id=c(1)", "n_iter=500", "n_yrs=50", "multiplier=TRUE")

setwd("C:/Users/erick/OneDrive/Ambiente de Trabalho/operating models and optimization")
# cl1 = NULL
### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### arguments ####
### ------------------------------------------------------------------------ ###

#args = c()
#args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  ### set default arguments
  ### parallelization
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("n_blocks")) n_blocks <- 3
  if (!exists("n_workers")) n_workers <- 3
  ### scenario definition
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 50
  if (!exists("fhist")) fhist <- "one-way" # fishing history 
  if (!exists("scenario")) scenario1 <- "OM 1" # operating model
  if (!exists("catch_rule")) catch_rule <- "catch_rule"
  if (!exists("comp_r")) comp_r <- TRUE
  if (!exists("comp_f")) comp_f <- TRUE
  if (!exists("comp_b")) comp_b <- TRUE
  if (!exists("scenario")) scenario <- "PA"
  if (!exists("cap_below_b")) cap_below_b <- TRUE
  ### GA search
  if (!exists("ga_search")) ga_search <- TRUE
  if (isTRUE(ga_search)) {
    if (!exists("popSize")) stop("popSize missing")
    if (!exists("maxiter")) stop("maxiter missing")
    if (!exists("stock_id")) stop("stock_id missing")
    if (!exists("run")) run <- maxiter
    if (!exists("collate")) collate <- TRUE
    ### objective function elements
    if (!exists("obj_SSB")) obj_SSB <- TRUE
    if (!exists("obj_F")) obj_F <- FALSE
    if (!exists("obj_C")) obj_C <- TRUE
    if (!exists("obj_risk")) obj_risk <- FALSE
    if (!exists("obj_ICV")) obj_ICV <- TRUE 
    if (!exists("obj_ICES_PA")) obj_ICES_PA <- FALSE
    if (!exists("obj_ICES_PA2")) obj_ICES_PA2 <- FALSE
    if (!exists("obj_ICES_MSYPA")) obj_ICES_MSYPA <- TRUE
    if (!exists("risk_threshold")) risk_threshold <- 0.05
    ### GA
    if (!exists("add_suggestions")) add_suggestions <- TRUE
    if (!exists("stat_yrs")) stat_yrs <- "all"
  }
  
} else {
  
  stop("no argument passed to R")
  
}

## lower and upper limits for the pso ##
ga_lower <- c(0, 1, 1, 1, 0, 0, 0, 1, 0.5, 1, 0.7)
ga_upper <- c(1, 5, 5, 1, 2, 2, 2, 5, 1.2, 1.2, 1)

np <- 20  #number of particles
mi <- 20  #max number of iterations

### seed ###
set.seed(4366)

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###
req_pckgs <- c("FLCore", "FLash", "mse", "GA", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)
library(FLBRP)

source("funs.R")
source("funs_ga.R")
#####################################################################################################

if (isTRUE(use_MPI)) {
  ### 1st: doMPI cluster with 1 worker per node
  message("starting doMPI")
  library(doMPI)
  cl1 <- startMPIcluster()
  message("startMPIcluster() succeeded")
  print(cl1)
  registerDoMPI(cl1)
  cl_length_1 <- cl1$workerCount
  print(cl_length_1)
  
  ### 2nd: doParallel workers inside doMPI workers
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    ### load packages and functions into MPI workers
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
  }
  message("MPI package loading succeeded")
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    source("funs.R", echo = FALSE)
    source("funs_ga.R", echo = FALSE)
  }
  message("MPI script loading succeeded")
  ### start doParallel inside MPI processes
  if (isTRUE(n_workers > 1)) {
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      cl2 <- makeCluster(n_workers)
      registerDoParallel(cl2)
      cl_length_2 <- length(cl2)
      ### load packages and functions into parallel workers
      . <- foreach(i = seq(cl_length_2)) %dopar% {
        for (i in req_pckgs) library(package = i, character.only = TRUE,
                                     warn.conflicts = FALSE, verbose = FALSE,
                                     quietly = TRUE)
        source("funs.R", echo = FALSE)
        source("funs_ga.R", echo = FALSE)
      }
    }
  }
  message("setting up doParallel inside MPI succeeded")
} else {
  if (isTRUE(n_workers > 1)) {
    ### start doParallel cluster
    cl1 <- makeCluster(n_workers)
    registerDoParallel(cl1)
    cl_length_1 <- length(cl1)
    ### load packages and functions into parallel workers
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      for (i in req_pckgs) library(package = i, character.only = TRUE,
                                   warn.conflicts = FALSE, verbose = FALSE,
                                   quietly = TRUE)
      source("funs.R", echo = FALSE)
      source("funs_GA.R", echo = FALSE)
    }
  } else {
    cl1 <- FALSE
  }
}
###########################################################################

### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###

stocks <- read.csv("input/stock_list_full_pt.csv", stringsAsFactors = FALSE)
stock <- stocks$stock[c(1)]
names(stock) <- stock

input <- lapply(stock, function(x) {
  readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/", x,
                 "_",scenario1,".rds"))
})

### ------------------------------------------------------------------------ ###
### specify scenario ####
### ------------------------------------------------------------------------ ###

### default catch rule
input <- lapply(input, function(x) {
  ### OEM: activate uncertainty
  x$oem@args$idx_dev <- TRUE
  x$oem@args$ssb_idx <- FALSE
  x$oem@args$tsb_idx <- FALSE
  x$oem@args$lngth <- TRUE
  x$oem@args$lngth_dev <- TRUE
  ### IEM: do not activate uncertainty
  x$iem@args$use_dev <- FALSE
  ### catch rule components
  x$ctrl$est@args$comp_r <- comp_r
  x$ctrl$est@args$comp_f <- comp_f
  x$ctrl$est@args$comp_b <- comp_b
  ### catch lag fixed
  x$ctrl$est@args$catch_lag <- 1
  ### turn of uncertainty cap when index below Itrigger?
  if (isFALSE(cap_below_b)) {
    x$ctrl$isys@args$cap_below_b <- cap_below_b
    #x$ctrl$isys@method <- is_comps
  }
  
  return(x)
})

### default ICES rule: 2 over 3
if (isTRUE(catch_rule == "2over3")) {
  input <- lapply(input, function(x) {
    ### OEM: turn of length index
    x$oem@args$lngth <- FALSE
    x$oem@args$lngth_dev <- FALSE
    ### add PA buffer stock status and deviation
    x$oem@args$PA_status <- TRUE
    x$oem@args$PA_status_dev <- TRUE
    ### catch rule components: turn of f & b, 2 over 3 rule
    x$ctrl$est@args$comp_f <- FALSE
    x$ctrl$est@args$comp_b <- FALSE
    x$ctrl$est@args$idxB_lag <- 1
    x$ctrl$est@args$idxB_range_1 <- 2
    x$ctrl$est@args$idxB_range_2 <- 3
    ### PA buffer
    x$ctrl$est@args$pa_buffer <- TRUE
    ###
    x$ctrl$phcr@args$exp_r <- 1
    x$ctrl$phcr@args$exp_f <- 0
    x$ctrl$phcr@args$exp_b <- 1 ### PA buffer
    ### biennial
    #x$ctrl$hcr@args$interval <- 2
    ### uncertainty cap
    x$ctrl$isys@args$upper_constraint <- 1.2
    x$ctrl$isys@args$lower_constraint <- 0.8
    return(x)
  })
  # input$pol$oem@method <- wklife_3.2.1_obs
  # input$pol$ctrl$est@method <- wklife_3.2.1_est
  # input$pol$ctrl$isys@method <- is_r
  #debugonce(goFishDL)
  #res <- do.call(mp, c(input$pol, cut_hist = FALSE))
  
} else if (isTRUE(catch_rule == "hr")) {
  input <- lapply(input, function(x) {
    ### OEM
    x$oem@method <- obs_generic
    x$oem@args$lngth <- FALSE
    x$oem@args$lngth_dev <- FALSE
    x$oem@args$ssb_idx <- FALSE
    x$oem@args$tsb_idx <- TRUE
    ### est
    x$ctrl$est@method <- est_hr
    x$ctrl$est@args$idxB_lag <- 1
    x$ctrl$est@args$idxB_range <- 1
    ### phcr
    x$ctrl$phcr@method <- phcr_hr
    x$ctrl$phcr@args$rate <- hr_rate
    ### hcr
    x$ctrl$hcr@method <- hcr_hr
    x$ctrl$hcr@args$interval <- 1
    ### is
    x$ctrl$isys@method <- is_r
    x$ctrl$isys@args$interval <- 1
    x$ctrl$isys@args$upper_constraint <- Inf
    x$ctrl$isys@args$lower_constraint <- 0
    return(x)
  })
}

### within scenario parallelisation?
if (isTRUE(n_workers > 1) & isTRUE(n_blocks > 1)) {
  ### use Iloss
  input <- lapply(input, function(x) {
    x$args$nblocks <- n_blocks
    return(x)
  })
}

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#remotes::install_github("https://github.com/ChatalovErick/parallelPSO.git")
library("parallelPSO")
library(parallel)

inp_file <- tempfile()
saveRDS(object = input, file = inp_file, compress = FALSE)
rm(input)
gc(reset = TRUE)

path_out <- paste0("output_pso/", n_iter, "_", n_yrs, "/",
                   fhist, "/",
                   paste0(stock,"_",scenario1), "/")

dir.create(path_out, recursive = TRUE)
################################################################################
################################################################################
ga_names <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
              "exp_r", "exp_f", "exp_b", "interval", "multiplier",
              "upper_constraint", "lower_constraint")

pso_res <- pso(fitness_function=mp_fitness,number_parameters = 11,number_of_partiples = np,W.1=0.9,W.2=0.6,
               max_number_iterations = mi,parameters_bounds = cbind(ga_lower,ga_upper),parallel = cl1,path=path_out,inp_file=inp_file)

pso_res$par
pso_res$value
saveRDS(object = pso_res, file = paste0(path_out,"sol.rds"))

################################################################################
################################################################################

mp_fitness1 <- function(par, inp_file, path, check_file = TRUE,scenario,
                       return_res = FALSE,collapse_correction = TRUE,
                       obj_SSB = TRUE, obj_C = TRUE, obj_F = FALSE,
                       obj_risk = FALSE, obj_ICV = TRUE, obj_ICES_PA = FALSE,
                       obj_ICES_PA2 = FALSE, obj_ICES_MSYPA = TRUE,
                       stat_yrs = "all",
                       risk_threshold = 0.05,
                       ...) {
  params <- par
  ### housekeeping
  invisible(gc())
  if (exists("res_mp")) {
    rm(res_mp)
    invisible(gc())
  }
  
  ### rounding of arguments
  params[1:4] <- round(params[1:4])
  params[5:7] <- round(params[5:7], 1)
  params[8] <- round(params[8])
  params[9] <- round(params[9], 2)
  params[10:11] <- round(params[10:11], 2)
  ### fix NaN for upper_constraint
  if (is.nan(params[10])) params[10] <- Inf
  
  ### check for files?
  if (isTRUE(check_file)) {
    ### current run
    run_i <- paste0(params, collapse = "_")
    ### get current stock(s)
    stock_i <- strsplit(x = tail(strsplit(x = path, split = "/")[[1]], 1), 
                        split = "_")[[1]]
    base_path <- paste0(paste0(head(strsplit(x = path, split = "/")[[1]], -1), 
                               collapse = "/"), "/")
    ### check if path exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    ### check if run already exists
    if (isTRUE(file.exists(paste0(path, run_i, ".rds")))) {
      ### load stats
      stats <- readRDS(paste0(path, run_i, ".rds"))
      ### set flag for running MP
      run_mp <- FALSE
      ### use different period to calculate stats?
      if (!any(stat_yrs %in% c("all", "more"))) {
        if (!any(grepl(x = rownames(stats), pattern = stat_yrs))) run_mp <- TRUE
      }
    } else {
      ### check if run exist in larger group
      dir_i <- paste0(stock_i, collapse = "_")
      dirs_i <- setdiff(x = dir(path = base_path, pattern = dir_i),
                        y = dir_i)
      if (isTRUE(length(dirs_i) > 0)) {
        dirs_i <- dirs_i[which(sapply(dirs_i, function(x) {
          tmp <- strsplit(x = x, split = "_")[[1]]
          ifelse(isFALSE(dir_i %in% tmp), FALSE, TRUE)
        }))]
        files_tmp <- lapply(dirs_i, function(x) {
          #browser()
          path_tmp <- paste0(base_path, x, "/", run_i, ".rds")
          if (isTRUE(file.exists(path = path_tmp))) {
            return(path_tmp)
          } else {
            return(NA)
          }
        })
        files_tmp[is.na(files_tmp)] <- NULL
        if (isTRUE(length(files_tmp) > 0)) {
          ### load stats from larger group
          stats <- readRDS(files_tmp[[1]])
          ### subset to current group
          stats <- stats[, stock_i]
          ### do not run MP
          run_mp <- FALSE
        } else {
          run_mp <- TRUE
        }
      } else {
        run_mp <- TRUE
      }
    }
  } else {
    run_mp <- TRUE
  }
  
  if (isTRUE(run_mp)) {
    
    ### load input file from disk
    input <- readRDS(inp_file)
    
    ### insert arguments into input object for mp
    input <- lapply(input, function(x) {
      x$ctrl$est@args$idxB_lag     <- params[1]
      x$ctrl$est@args$idxB_range_1 <- params[2]
      x$ctrl$est@args$idxB_range_2 <- params[3]
      x$ctrl$est@args$catch_range  <- params[4]
      x$ctrl$est@args$comp_m <- params[9]
      x$ctrl$phcr@args$exp_r <- params[5]
      x$ctrl$phcr@args$exp_f <- params[6]
      x$ctrl$phcr@args$exp_b <- params[7]
      x$ctrl$hcr@args$interval <- params[8]
      x$ctrl$isys@args$interval <- params[8]
      x$ctrl$isys@args$upper_constraint <- params[10]
      x$ctrl$isys@args$lower_constraint <- params[11]
      
      return(x)
    })
    
    ### if group of stocks, check if results for individual stocks exist
    group <- ifelse(isTRUE(length(input) > 1) & isTRUE(check_file), TRUE, FALSE)
    if (group) {
      ### get paths
      group_stocks <- names(input)
      path_base <- gsub(x = path, 
                        pattern = paste0(paste0(group_stocks, collapse = "_"), 
                                         "/"),
                        replacement = "")
      path_stocks <- paste0(path_base, group_stocks, "/")
      ### check for files
      run_exists <- file.exists(paste0(path_stocks, run_i, ".rds"))
      group <- ifelse(any(run_exists), TRUE, FALSE)
      
      ### do some results exist?
      if (group) {
        ### load results
        files_exist <- paste0(path_stocks, run_i, ".rds")[run_exists]
        stats_group <- lapply(files_exist, readRDS)
        names(stats_group) <- group_stocks[run_exists]
        ### get stocks which require simulation
        run_stocks <- group_stocks[!run_exists]
        ### subset input
        input <- input[run_stocks]
        
      }
      
    }
    
    res_mp <- lapply(input, function(x) { do.call(mp,x) }) 
    
    return(res_mp)
  }
  
}

################################################################################
################################################################################
################################################################################
params <- pso_res$par

path_out_sol <- paste0("best_solution from each scenario/", n_iter, "_", n_yrs, "/",
                       fhist, "/",
                       paste0(stock, collapse = "_"), "/")
dir.create(path_out_sol, recursive = TRUE)


res <- mp_fitness1(params,inp_file, path=path_out_sol)

saveRDS(res,paste0(path_out_sol,stock,"_",scenario1,".rds"))
################################################################################
################################################################################
