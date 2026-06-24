################################################################################
# TITLE: case_study_simulations.R
#
# PURPOSE: run simulations based on case study data to give optimal MA method
#          and corresponding weight cutoffs
#
# SECTIONS: Section 0 - load packages 
#           Section 1 - helper functions (risk distr., model templates)
#           Section 2 - model averaging functions
#           Section 3 - calibration function (simulations in parallel)
#           Section 4 - run case study simulations
#
# AUTHOR: Shannon Thomas
# NOTES: The progress bar feature is more experimental than it is useful because
#        the parallelization of the simulations makes the values jump back and
#        forth rather than steadily increasing/decreasing as the simulation runs.
################################################################################

##########################################################
################ SECTION 0 : LOAD PACKAGES ############### 
##########################################################

library(brms)
library(progressr)

##########################################################
############### SECTION 1: HELPER FUNCTIONS ############## 
##########################################################

# risk distribution functions
lambda_exponential <- function(p_event, t) {
  -log(1 - p_event) / t
}

scale_weibull <- function(p_event, t, shape) {
  t / (-log(1 - p_event))^(1 / shape)
}

mu_lognormal <- function(p_event, t, sigma) {
  log(t) - qnorm(p_event) * sigma
}

make_templates <- function(N) {
  
  # minimal dummy data (just to compile model)
  dummy_dat_count <- data.frame(
    trt = factor(rep(c(0,1), each = N/2)),
    event_ind = rbinom(N, 1, 0.5)
  )
  
  dummy_dat_surv <- data.frame(
    trt = factor(rep(c(0,1), each = N/2)),
    time_surv = rexp(N, 1),
    event_ind = rbinom(N, 1, 0.5)
  )
  
  #COUNT MODELS
  fit_count_int <- brm(
    event_ind ~ 1,
    data = dummy_dat_count,
    family = bernoulli(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_count_trt <- brm(
    event_ind ~ trt,
    data = dummy_dat_count,
    family = bernoulli(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  #SURVIVAL MODELS
  fit_survint_exp <- brm(
    time_surv | cens(1 - event_ind) ~ 1,
    data = dummy_dat_surv,
    family = exponential(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_survint_weib <- brm(
    time_surv | cens(1 - event_ind) ~ 1,
    data = dummy_dat_surv,
    family = weibull(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_exp_trt <- brm(
    time_surv | cens(1 - event_ind) ~ trt,
    data = dummy_dat_surv,
    family = exponential(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_weib_trt <- brm(
    time_surv | cens(1 - event_ind) ~ trt,
    data = dummy_dat_surv,
    family = weibull(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  #output
  list(
    count_int = fit_count_int,
    count_trt = fit_count_trt,
    survint_exp = fit_survint_exp,
    survint_weib = fit_survint_weib,
    exp_trt = fit_exp_trt,
    weib_trt = fit_weib_trt
  )
}

##########################################################
########## SECTION 2: MODEL AVERAGING FUNCTIONS ########## 
##########################################################

#count model fucntion
runmodels_count <- function(N, trt_effect_event, control_event_rate, templates) {
  
  X_trt <- rep(c(0,1), each = N/2)
  
  p_ctrl <- control_event_rate
  p_trt  <- plogis(qlogis(control_event_rate) + trt_effect_event)
  
  y_ctrl <- rbinom(N/2, 1, p_ctrl)
  y_trt  <- rbinom(N/2, 1, p_trt)
  
  dat <- data.frame(
    trt = factor(X_trt),
    event_ind = c(y_ctrl, y_trt)
  )
  
  if(sum(dat$event_ind) <= 1) return(NULL)
  
  fit_int <- update(templates$count_int,
                    newdata = dat,
                    recompile = FALSE)
  
  fit_trt <- update(templates$count_trt,
                    newdata = dat,
                    recompile = FALSE)
  
  rbind(
    stacking = model_weights(fit_int, fit_trt, weights = "stacking"),
    pBMA     = model_weights(fit_int, fit_trt, weights = "pseudobma"),
    BMA      = model_weights(fit_int, fit_trt, weights = "bma")
  )
}



#survival model function
runmodels_survival <- function(N, trt_effect_event, distr,
                               control_event_rate, censoring_time, templates) {
  
  X_trt <- rep(c(0,1), each = N/2)
  hr <- exp(trt_effect_event)
  
  if(distr == "exponential"){
    lambda_c <- lambda_exponential(control_event_rate, censoring_time)
    lambda_t <- lambda_c * hr
    time_c <- rexp(N/2, lambda_c)
    time_t <- rexp(N/2, lambda_t)
  }
  
  if(distr == "weibull"){
    shape <- 2
    scale_c <- scale_weibull(control_event_rate, censoring_time, shape)
    scale_t <- scale_c / hr^(1/shape)
    time_c <- rweibull(N/2, shape, scale_c)
    time_t <- rweibull(N/2, shape, scale_t)
  }
  
  if(distr == "lognormal"){
    sigma <- 0.7
    mu_c <- mu_lognormal(control_event_rate, censoring_time, sigma)
    mu_t <- mu_c - trt_effect_event
    time_c <- rlnorm(N/2, mu_c, sigma)
    time_t <- rlnorm(N/2, mu_t, sigma)
  }
  
  event_times <- c(time_c, time_t)
  censor_times <- runif(N, 0, censoring_time)
  
  dat <- data.frame(trt = factor(X_trt),
                    time_e = event_times,
                    time_c = censor_times)
  
  dat$time_surv <- pmin(dat$time_e, dat$time_c)
  dat$event_ind <- ifelse(dat$time_e <= dat$time_c, 1, 0)
  
  if(sum(dat$event_ind) <= 1) return(NULL)
  
  fit_survint <- update(templates$survint_exp,
                        newdata = dat,
                        recompile = FALSE)
  
  fit_weibint <- update(templates$survint_weib,
                        newdata = dat,
                        recompile = FALSE)
  
  fit_exp <- update(templates$exp_trt,
                    newdata = dat,
                    recompile = FALSE)
  
  fit_weibull <- update(templates$weib_trt,
                        newdata = dat,
                        recompile = FALSE)
  
  rbind(
    stacking = model_weights(fit_survint, fit_weibint, fit_exp, fit_weibull, weights = "stacking"),
    pBMA     = model_weights(fit_survint, fit_weibint, fit_exp, fit_weibull, weights = "pseudobma"),
    BMA      = model_weights(fit_survint, fit_weibint, fit_exp, fit_weibull, weights = "bma")
  )
}


#function to extract signal from count or survival models
extract_signal <- function(weights, data_type){
  
  if(is.null(weights)) return(c(NA,NA,NA))
  
  if(data_type == "count"){
    return(c(
      stacking = weights["stacking","fit_trt"],
      pBMA     = weights["pBMA","fit_trt"],
      BMA      = weights["BMA","fit_trt"]
    ))
  }
  
  if(data_type == "survival"){
    return(c(
      stacking = sum(weights["stacking", c("fit_exp","fit_weibull")]),
      pBMA     = sum(weights["pBMA", c("fit_exp","fit_weibull")]),
      BMA      = sum(weights["BMA", c("fit_exp","fit_weibull")])
    ))
  }
}

##########################################################
############# SECTION 3: CALIBRATION FUNCTION ############ 
##########################################################

calibrate_scenario <- function(N,control_event_rate,data_type = c("count","survival"),distr = "exponential",
                               censoring_time = 3, trt_effect = NULL, nsim_null = 500, nsim_power = 500, seed = 123,
                               parallel = FALSE, workers = 4, templates){
  
  set.seed(seed)

  sim_fun <- if(data_type == "count") runmodels_count else runmodels_survival
  
  #require packages
  if (!requireNamespace("progressr", quietly = TRUE)) stop("install.packages('progressr')")
  if (parallel && !requireNamespace("future.apply", quietly = TRUE)) stop("install.packages('future.apply')")
  
  options(progressr.enable = TRUE)
  progressr::handlers(global = TRUE)
  progressr::handlers("progress")
  
  #parallel set-up
  old_plan <- NULL
  if (parallel) {
    old_plan <- future::plan()
    future::plan(future::multisession, workers = workers)
  }
  
  
  #define a function to convert the NULL outputs to NA
  safe_extract <- function(w) {
    if (is.null(w)) {
      return(c(NA, NA, NA))
    } else {
      return(extract_signal(w, data_type))
    }
  }
  
  #define function to run simulation with progress bars
  run_with_progress <- function(nsim, sim_expr, label) {
    
    progressr::with_progress({
      
      p <- progressr::progressor(steps = nsim)
      
      start_time <- Sys.time()
      results <- vector("list", nsim)
      
      run_one <- function(i, templates) {
        
        val <- sim_expr(templates = templates)
        results[[i]] <<- val
        
        # ---- KEEP NA, summarize conditional ----
        mat <- do.call(rbind, results)
        
        means <- if (!is.null(mat) && nrow(mat) > 5) {
          colMeans(mat, na.rm = TRUE)
        } else {
          rep(NA, 3)
        }
        
        # ---- ETA ----
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        avg_time <- elapsed / i
        eta <- avg_time * (nsim - i)
        eta_str <- sprintf("%02d:%02d", floor(eta/60), round(eta %% 60))
        
        stat_str <- paste0(
          "mean weights: ",
          paste(names(means), round(means,3), collapse = ", ")
        )
        
        p(sprintf("%s %d/%d | ETA: %s | %s",
                  label, i, nsim, eta_str, stat_str))
        
        return(val)
      }
      
      if (parallel) {
        out <- future.apply::future_lapply(
          seq_len(nsim),
          function(i) run_one(i, templates = templates),
          future.seed = TRUE
        )
      } else {
        out <- lapply(seq_len(nsim), run_one)
      }
      
      do.call(rbind, out)
    })
  }
  
  #run null scenario simulation
  null_signals <- run_with_progress(
    nsim_null,
    sim_expr = function(templates) {
      
      w <- if(data_type == "count"){
        sim_fun(N, log(1.0), control_event_rate, templates = templates)
      } else {
        sim_fun(N, log(1.0), distr, control_event_rate, censoring_time, templates = templates)
      }
      
      safe_extract(w)
    },
    label = "Null sims"
  )
  
  #set NA to 0 weight (non-reject)
  null_signals_filled <- null_signals
  null_signals_filled[is.na(null_signals_filled)] <- 0
  
  thresholds <- apply(null_signals_filled, 2, quantile, 0.95)
  
  #power calculation
  power <- NULL
  
  if(!is.null(trt_effect)){
    
    alt_signals <- run_with_progress(
      nsim_power,
      sim_expr = function(templates) {
        
        w <- if(data_type == "count"){
          sim_fun(N, trt_effect, control_event_rate, templates = templates)
        } else {
          sim_fun(N, trt_effect, distr, control_event_rate, censoring_time, templates = templates)
        }
        
        safe_extract(w)
      },
      label = "Power sims"
    )
    
    #define decision and set NA to FALSE (non-reject)
    decision_matrix <- alt_signals > matrix(thresholds,
                                            nrow = nrow(alt_signals),
                                            ncol = length(thresholds),
                                            byrow = TRUE)
    
    decision_matrix[is.na(decision_matrix)] <- FALSE
    
    power <- colMeans(decision_matrix)
  }
  
  #NA diagnostics
  prop_na_null <- colMeans(is.na(null_signals))
  
  prop_na_alt <- NULL
  if(!is.null(trt_effect)){
    prop_na_alt <- colMeans(is.na(alt_signals))
  }
  

  #reset parallel  
  if (parallel && !is.null(old_plan)) {
    future::plan(old_plan)
  }
  

  #output
  list(
    thresholds     = thresholds,
    power          = power,
    prop_na_null   = prop_na_null,
    prop_na_alt    = prop_na_alt
  )
}

##########################################################
########## SECTION 4: RUN CASE STUDY SIMULATIONS ######### 
##########################################################

#### MAKE MODEL TEMPLATES
templates <- make_templates(N = 478)

cal_allae <- calibrate_scenario(N = 478, control_event_rate = 0.97,
                                data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                                parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_allae, file='casestudy_sim_results/cal_allae.rds')

cal_arrhyANDcvaANDmi <- calibrate_scenario(N = 478, control_event_rate = 0.01,
                                data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                                parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_arrhyANDcvaANDmi, file='casestudy_sim_results/cal_arrhyANDcvaANDmi.rds')

cal_cv <- calibrate_scenario(N = 478, control_event_rate = 0.24,
                             data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                             parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_cv, file='casestudy_sim_results/cal_cv.rds')

cal_emb <- calibrate_scenario(N = 478, control_event_rate = 0.18,
                              data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                              parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_emb, file='casestudy_sim_results/cal_emb.rds')

cal_hyp <- calibrate_scenario(N = 478, control_event_rate = 0.05,
                              data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                              parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_hyp, file='casestudy_sim_results/cal_hyp.rds')

cal_ser <- calibrate_scenario(N = 478, control_event_rate = 0.40,
                              data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                              parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_ser, file='casestudy_sim_results/cal_ser.rds')

cal_sev <- calibrate_scenario(N = 478, control_event_rate = 0.56,
                              data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                              parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_sev, file='casestudy_sim_results/cal_sev.rds')

cal_trtrel <- calibrate_scenario(N = 478, control_event_rate = 0.06,
                                 data_type = "survival", distr = "weibull", trt_effect = log(1.1), 
                                 parallel = TRUE, workers = 5, templates = templates)
saveRDS(cal_trtrel, file='casestudy_sim_results/cal_trtrel.rds')



