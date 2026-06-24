#################################################################################
# TITLE: simulation.R
#
# PURPOSE: simulate AE data, run BMA and stacking, and plot model results
#
# OUTPUT: .rds files with all simulation output
#
# SECTIONS: Section 0 - load packages
#           Section 1 - risk distributions
#           Section 2 - simulation functions
#           Section 3 - run simulation
#           Section 4 - plot results
#
# AUTHOR: Shannon Thomas
# DATE CREATED: MAR 27, 2026
# NOTES: takes days to run
#################################################################################

##########################################################
################ SECTION 0: LOAD PACKAGES ################
##########################################################

library(brms)
library(survival)
library(splines)
library(loo)
library(tidyverse)
library(patchwork)
library(future)
library(future.apply)

##########################################################
############# SECTION 1: Risk Distributions ##############
##########################################################

# #consider lognormal distribution for weird bump in prob that alex is interested in 
# #SURVIVAL DATA SIMULATION
# N <- 20
# trt_effect_event <- log(1.0)
# shape_par <- 2 
# hr <- exp(trt_effect_event)
# control_event_rate <- 0.99
# censoring_time <- 3

#test baseline risks 
  
  # Exponential: P(T <= t) = 1 - exp(-lambda * t)
  lambda_exponential <- function(p_event, t) {
    -log(1 - p_event) / t
  }
  
  # Weibull: P(T <= t) = 1 - exp(-(t / scale)^shape)
  scale_weibull <- function(p_event, t, shape) {
    t / (-log(1 - p_event))^(1 / shape)
  }
  
  # Log-normal: P(T <= t) = Phi((log(t) - mu) / sigma)
  mu_lognormal <- function(p_event, t, sigma) {
    log(t) - qnorm(p_event) * sigma
  }
  
  # #########################
  # # Exponential
  # #########################
  # lambda_c <- lambda_exponential(control_event_rate, censoring_time)
  # lambda_t <- lambda_c * hr
  # 
  # time_c <- rexp(n=N/2, rate = lambda_c)
  # time_t <- rexp(n=N/2, rate = lambda_t)
  # 
  # dist_c <- dexp(x=seq(0,5,0.1), rate = lambda_c)
  # dist_t <- dexp(x=seq(0,5,0.1), rate = lambda_t)
  # 
  # plot(x = seq(0,5,0.1), y = dist_c, ylim = c(0,2)) +
  #   lines(x = seq(0,5,0.1), y = dist_t)
  # 
  # hist(time_c)
  # hist(time_t)
  # 
  # #########################
  # # Weibull (shape = 2)
  # #########################
  # shape_w <- 2
  # 
  # scale_c <- scale_weibull(control_event_rate, censoring_time, shape_w)
  # scale_t <- scale_c / hr^(1 / shape_w)
  # 
  # time_c <- rweibull(n=N/2, shape = shape_w, scale = scale_c)
  # time_t <- rweibull(n=N/2, shape = shape_w, scale = scale_t)
  # 
  # dist_c <- dweibull(x=seq(0,5,0.1), shape = shape_w, scale = scale_c)
  # dist_t <- dweibull(x=seq(0,5,0.1), shape = shape_w, scale = scale_t)
  # 
  # plot(x = seq(0,5,0.1), y = dist_c, ylim = c(0,2)) +
  #   lines(x = seq(0,5,0.1), y = dist_t)
  # 
  # hist(time_c)
  # hist(time_t)
  # 
  # #########################
  # # Log-normal
  # #########################
  # sigma_ln <- 0.7
  # 
  # mu_c <- mu_lognormal(control_event_rate, censoring_time, sigma_ln)
  # mu_t <- mu_c - trt_effect_event
  # 
  # time_c <- rlnorm(n=N/2, meanlog = mu_c, sdlog = sigma_ln)
  # time_t <- rlnorm(n=N/2, meanlog = mu_t, sdlog = sigma_ln)
  # 
  # dist_c <- dlnorm(x=seq(0,5,0.1), meanlog = mu_c, sdlog = sigma_ln)
  # dist_t <- dlnorm(x=seq(0,5,0.1), meanlog = mu_t, sdlog = sigma_ln)
  # 
  # plot(x = seq(0,5,0.1), y = dist_c, ylim = c(0,2)) +
  #   lines(x = seq(0,5,0.1), y = dist_t)
  # 
  # hist(time_c)
  # hist(time_t)
  # 

##########################################################
############ SECTION 2: Simulation Functions #############
##########################################################


########## Predicted value function used in simulation function

get_averaged_pred <- function(type = "logistic", weights, new_data,
                              fit_int = NULL, fit_trt = NULL,
                              fit_survint = NULL, fit_exp = NULL, 
                              fit_weibint = NULL, fit_weibull = NULL) {

  if(type == "logistic") {
    # 1. Extract Point-in-Time Probabilities
    lp1 <- plogis(posterior_linpred(fit_int, newdata = new_data))
    lp2 <- plogis(posterior_linpred(fit_trt,  newdata = new_data))
    
    # 2. Calculate Weighted Average
    avg <- (lp1 * weights["fit_int"]) + (lp2 * weights["fit_trt"])   
    
  } else {
    # 1. Extract parameters for Survival Models
    mu_int     <- posterior_linpred(fit_survint, newdata = new_data)
    mu_exp     <- posterior_linpred(fit_exp,     newdata = new_data)
    mu_weibint <- posterior_linpred(fit_weibint, newdata = new_data)
    mu_wei     <- posterior_linpred(fit_weibull, newdata = new_data)
    
    k_weibint  <- as.matrix(fit_weibint, variable = "shape")
    k_wei  <- as.matrix(fit_weibull, variable = "shape")
    
    # 2. Calculate Hazards h(t)
    h_int <- 1 / exp(mu_int) # Constant hazard for intercept-only
    h_exp <- 1 / exp(mu_exp) # Constant hazard for exp + trt
    
    t_mat <- matrix(new_data$time_surv, nrow = nrow(mu_wei), ncol = nrow(new_data), byrow = TRUE)
    
    lambda_weibint <- exp(mu_weibint) / gamma(1 + 1/as.vector(k_weibint))
    h_weibint      <- (as.vector(k_weibint) / lambda_weibint) * (t_mat / lambda_weibint)^(as.vector(k_weibint) - 1)
    
    lambda_wei <- exp(mu_wei) / gamma(1 + 1/as.vector(k_wei))
    h_wei      <- (as.vector(k_wei) / lambda_wei) * (t_mat / lambda_wei)^(as.vector(k_wei) - 1)
    
    # 3. Calculate Weighted Average
    avg <- (h_int * weights["fit_survint"]) + 
      (h_exp * weights["fit_exp"]) + 
      (h_weibint * weights["fit_weibint"]) +
      (h_wei * weights["fit_weibull"])
    
  }
  
  # Summarize 
  res <- data.frame(prob_event = colMeans(avg),
                    low = apply(avg, 2, quantile, 0.025),
                    high = apply(avg, 2, quantile, 0.975))
  
  return(res)
}





########## Main Simulation Function
runmodels <- function(i, N, trt_effect_event, distr, control_event_rate, censoring_time, brms_models){
  
  ################# SIMULATE DATA ################
  #set.seed(i+406) #don't need this when future.seed = TRUE
  
  X_trt <- rep(c(0,1), each = N/2) #treatment assignment (placebo for first half and treatment for second half)
  
  if(distr == "exponential"){
    hr <- exp(trt_effect_event)
    lambda_c <- lambda_exponential(control_event_rate, censoring_time)
    lambda_t <- lambda_c * hr
    
    time_c <- rexp(n=N/2, rate = lambda_c)
    time_t <- rexp(n=N/2, rate = lambda_t)
    
  }
  if(distr == "weibull"){
    hr <- exp(trt_effect_event)
    shape_par <- 2
    
    scale_c <- scale_weibull(control_event_rate, censoring_time, shape_par)
    scale_t <- scale_c / hr^(1 / shape_par)
    
    time_c <- rweibull(n=N/2, shape = shape_par, scale = scale_c)
    time_t <- rweibull(n=N/2, shape = shape_par, scale = scale_t)
  
  }
  if(distr == "lognormal"){
    sigma_ln <- 0.7
    
    mu_c <- mu_lognormal(control_event_rate, censoring_time, sigma_ln)
    mu_t <- mu_c - trt_effect_event
    
    time_c <- rlnorm(n=N/2, meanlog = mu_c, sdlog = sigma_ln)
    time_t <- rlnorm(n=N/2, meanlog = mu_t, sdlog = sigma_ln)
  }
  
  event_times <- c(time_c, time_t)
  
  
  censor_times <- runif(n = N, min = 0, max = censoring_time)  #random censoring (i.e. drop out)
  censor_times[censor_times >= censoring_time] <- censoring_time #administrative censoring
  
  #combine simulated data into data frame
  dat <- data.frame(trt = factor(X_trt, levels = c(0,1)), 
                    time_e = event_times, 
                    time_c = censor_times)
  
  #define event time based on censoring and event times
  dat$time_surv <- pmin(dat$time_e, dat$time_c)
  
  #define event indicator
  dat$event_ind <- ifelse(dat$time_surv == dat$time_e, 1, 0)

  #calculate event counts
  trt_eventcount  <- sum(dat$event_ind[dat$trt == "1"])
  ctrl_eventcount <- sum(dat$event_ind[dat$trt == "0"])
  
  # print(paste("N =", N))
  # print(paste("treatment event count =", trt_eventcount))
  # print(paste("control event count =", ctrl_eventcount))
  
  if(trt_eventcount + ctrl_eventcount > 1){
    ############ FIT STANDARD MODELS ##############
    lmfit <- glm(event_ind ~ trt, family = binomial(link = "logit"), data = dat) 
    
    mainfit <- glm(event_ind~1, family = binomial(link = "logit"), data = dat)
    lrtresult <- lmtest::lrtest(mainfit, lmfit)
    
    lrtpval <- lrtresult$`Pr(>Chisq)`[2] #p-value from LRT
    
    
    chisqpval <- tryCatch({  chisq.test(matrix(c(ctrl_eventcount, trt_eventcount, 
                                                 N/2 - ctrl_eventcount, N/2 - trt_eventcount), 
                                               ncol = 2))$p.value}, 
                          warning = function(w){
                            fisher.test(matrix(c(ctrl_eventcount, trt_eventcount, 
                                                 N/2 - ctrl_eventcount, N/2 - trt_eventcount),
                                               ncol = 2))$p.value})
    
    ############ FIT LOGISTIC MODELS ##############
    fit_int <- update(brms_models$fit_ind, 
                   newdata = dat, 
                   chains = 2, cores = 1, iter = 2000,
                   silent = 2, refresh = 0)
    
    fit_trt <- update(brms_models$fit_trt, 
                      newdata = dat, 
                      chains = 2, cores = 1, iter = 2000,
                      silent = 2, refresh = 0)
    
    #compute weights 
    logistic_w_stack <-  model_weights(fit_int, fit_trt,  
                                       weights = "stacking")
    logistic_w_pBMA  <-  model_weights(fit_int, fit_trt,  
                                       weights = "pseudobma")
    logistic_w_BMA   <-  model_weights(fit_int, fit_trt,  
                                       weights = "bma")
    logistic_weights <- rbind(logistic_w_stack, logistic_w_pBMA, logistic_w_BMA)
    
    #get posterior model summaries
    logistic_summary_stack <- posterior_summary(posterior_average(fit_int, fit_trt,  
                                                                  weights = "stacking",  missing = 0))
    logistic_summary_pBMA  <- posterior_summary(posterior_average(fit_int, fit_trt,  
                                                                  weights = "pseudobma", missing = 0))
    logistic_summary_BMA   <- posterior_summary(posterior_average(fit_int, fit_trt,  
                                                                  weights = "bma",       missing = 0))
    
    # ######### GENERATE PREDICTED VALUES ###########   
    # time_seq <- seq(0, censoring_time, length.out = 200)
    # new_data <- expand.grid(time_surv = time_seq, trt = unique(dat$trt))
    new_data <- NULL
    # 
    # logistic_pred_stack <- get_averaged_pred(type = "logistic", weights = logistic_weights["logistic_w_stack",], 
    #                                          new_data = new_data, fit_int = fit_int, fit_trt = fit_trt)
    # logistic_pred_pBMA  <- get_averaged_pred(type = "logistic", weights = logistic_weights["logistic_w_pBMA",], 
    #                                          new_data = new_data, fit_int = fit_int, fit_trt = fit_trt)
    # logistic_pred_BMA   <- get_averaged_pred(type = "logistic", weights = logistic_weights["logistic_w_BMA",], 
    #                                          new_data = new_data, fit_int = fit_int, fit_trt = fit_trt)
    
    
    ############ FIT SURVIVAL MODELS ##############
    fit_survint <- update(brms_models$fit_survint, 
                          newdata = dat, 
                          chains = 2, cores = 1, iter = 2000,
                          silent = 2, refresh = 0)
    
    fit_weibint <- update(brms_models$fit_weibint, 
                          newdata = dat, 
                          chains = 2, cores = 1, iter = 2000,
                          silent = 2, refresh = 0)
    
    fit_exp <- update(brms_models$fit_exp, 
                      newdata = dat, 
                      chains = 2, cores = 1, iter = 2000,
                      silent = 2, refresh = 0)
    
    fit_weibull <- update(brms_models$fit_weibull, 
                          newdata = dat, 
                          chains = 2, cores = 1, iter = 2000,
                          silent = 2, refresh = 0)
    
      #compute weights 
      w_stack <-  model_weights(fit_survint, fit_weibint, fit_exp, fit_weibull, 
                                weights = "stacking")
      w_pBMA  <-  model_weights(fit_survint, fit_weibint, fit_exp, fit_weibull, 
                                weights = "pseudobma")
      w_BMA   <-  model_weights(fit_survint, fit_weibint, fit_exp, fit_weibull, 
                                weights = "bma")
      survival_weights <- rbind(w_stack, w_pBMA, w_BMA)
      
      #get posterior model summaries
      summary_stack <- posterior_summary(posterior_average(fit_survint, fit_weibint, fit_exp, fit_weibull, 
                                                           weights = "stacking",  missing = 0))
      summary_pBMA  <- posterior_summary(posterior_average(fit_survint, fit_weibint, fit_exp, fit_weibull, 
                                                           weights = "pseudobma", missing = 0))
      summary_BMA   <- posterior_summary(posterior_average(fit_survint, fit_weibint, fit_exp, fit_weibull, 
                                                           weights = "bma",       missing = 0))
  
      # ######### GENERATE PREDICTED VALUES ###########
      # survival_pred_stack <- get_averaged_pred(type = "survival", weights = survival_weights["w_stack",], 
      #                                          new_data = new_data, fit_survint = fit_survint, fit_exp = fit_exp,
      #                                          fit_weibint = fit_weibint, fit_weibull = fit_weibull)
      # survival_pred_pBMA  <- get_averaged_pred(type = "survival", weights = survival_weights["w_pBMA",], 
      #                                          new_data = new_data, fit_survint = fit_survint, fit_exp = fit_exp,
      #                                          fit_weibint = fit_weibint, fit_weibull = fit_weibull)
      # survival_pred_BMA   <- get_averaged_pred(type = "survival", weights = survival_weights["w_BMA",], 
      #                                          new_data = new_data, fit_survint = fit_survint, fit_exp = fit_exp,
      #                                          fit_weibint = fit_weibint, fit_weibull = fit_weibull)
  } else{
    time_seq <- seq(0, censoring_time, length.out = 200)
    new_data <- expand.grid(time_surv = time_seq, trt = unique(dat$trt))
    
    chisqpval = NA
    lrtpval   = NA
    logistic_weights       = data.frame(fit_int = rep(NA,3), fit_trt = rep(NA,3), 
                                        row.names = c("logistic_w_stack","logistic_w_pBMA", "logistic_w_BMA"))
    logistic_summary_stack = data.frame(Estimate = rep(NA,4), Est.Error = rep(NA,4), Q2.5 = rep(NA,4), Q97.5 = rep(NA,4),
                                        row.names = c("b_Intercept","Intercept", "lprior", "b_trt1"))
    logistic_summary_pBMA  = data.frame(Estimate = rep(NA,4), Est.Error = rep(NA,4), Q2.5 = rep(NA,4), Q97.5 = rep(NA,4),
                                        row.names = c("b_Intercept","Intercept", "lprior", "b_trt1"))
    logistic_summary_BMA   = data.frame(Estimate = rep(NA,4), Est.Error = rep(NA,4), Q2.5 = rep(NA,4), Q97.5 = rep(NA,4),
                                        row.names = c("b_Intercept","Intercept", "lprior", "b_trt1"))
    # logistic_pred_stack    = data.frame(prob_event = rep(NA, nrow(new_data)),
    #                                     low = rep(NA, nrow(new_data)),
    #                                     high = rep(NA, nrow(new_data)))
    # logistic_pred_pBMA     = data.frame(prob_event = rep(NA, nrow(new_data)),
    #                                     low = rep(NA, nrow(new_data)),
    #                                     high = rep(NA, nrow(new_data)))
    # logistic_pred_BMA      = data.frame(prob_event = rep(NA, nrow(new_data)),
    #                                     low = rep(NA, nrow(new_data)),
    #                                     high = rep(NA, nrow(new_data)))
    survival_weights       = data.frame(fit_survint = rep(NA,3), fit_weibint = rep(NA,3), fit_exp = rep(NA,3), fit_weibull = rep(NA,3),
                                        row.names = c("w_stack","w_pBMA", "w_BMA"))
    summary_stack = data.frame(Estimate = rep(NA,5), Est.Error = rep(NA,5), Q2.5 = rep(NA,5), Q97.5 = rep(NA,5),
                                        row.names = c("b_Intercept","Intercept", "lprior", "shape","b_trt1"))
    summary_pBMA  = data.frame(Estimate = rep(NA,5), Est.Error = rep(NA,5), Q2.5 = rep(NA,5), Q97.5 = rep(NA,5),
                                        row.names = c("b_Intercept","Intercept", "lprior", "shape","b_trt1"))
    summary_BMA   = data.frame(Estimate = rep(NA,5), Est.Error = rep(NA,5), Q2.5 = rep(NA,5), Q97.5 = rep(NA,5),
                                        row.names = c("b_Intercept","Intercept", "lprior", "shape","b_trt1"))
    # survival_pred_stack    = data.frame(prob_event = rep(NA, nrow(new_data)),
    #                                     low = rep(NA, nrow(new_data)),
    #                                     high = rep(NA, nrow(new_data)))
    # survival_pred_pBMA     = data.frame(prob_event = rep(NA, nrow(new_data)),
    #                                     low = rep(NA, nrow(new_data)),
    #                                     high = rep(NA, nrow(new_data)))
    # survival_pred_BMA      = data.frame(prob_event = rep(NA, nrow(new_data)),
    #                                     low = rep(NA, nrow(new_data)),
    #                                     high = rep(NA, nrow(new_data)))
  }
    
  ############# CREATE OUTPUT LIST ##############
  out <- list(samp_num           = i, 
              N                  = N, 
              distr              = distr,
              trt_effect         = trt_effect_event,
              control_event_rate = control_event_rate,
              censoring_time     = censoring_time, 
              trt_eventcount  = trt_eventcount,
              ctrl_eventcount = ctrl_eventcount,
              chisqpval = chisqpval,
              lrtpval   = lrtpval,
              logistic_weights       = logistic_weights,
              logistic_summary_stack = logistic_summary_stack,
              logistic_summary_pBMA  = logistic_summary_pBMA,
              logistic_summary_BMA   = logistic_summary_BMA,
              # logistic_pred_stack    = logistic_pred_stack,
              # logistic_pred_pBMA     = logistic_pred_pBMA,
              # logistic_pred_BMA      = logistic_pred_BMA,
              survival_weights       = survival_weights,
              survival_summary_stack = summary_stack,
              survival_summary_pBMA  = summary_pBMA,
              survival_summary_BMA   = summary_BMA,
              # survival_pred_stack    = survival_pred_stack,
              # survival_pred_pBMA     = survival_pred_pBMA,
              # survival_pred_BMA      = survival_pred_BMA,
              new_data = new_data)
  
  return(out)
}



##########################################################
################ SECTION 3: RUN FUNCTIONS ################ 
##########################################################

nsamp <- 1
nsim <- 500
sampnums <- seq(1:nsamp)
trt_effs <- c(log(1.0),log(1.1),log(2.0))
sampsizes <- c(20,50,100,500,1000,5000)
distr <- c("exponential", "weibull", "lognormal")
control_event_rate <- c(0.1,0.99)
censoring_times <- c(3)

allparamcombos <- expand.grid(samp_num = sampnums, 
                              trt_effect = trt_effs, 
                              N = sampsizes, 
                              distr = distr, 
                              control_event_rate = control_event_rate,
                              censoring_time = censoring_times)
allparamcombos <- apply(allparamcombos,1,as.list)

numeric_fields <- c("samp_num", "trt_effect", "N", "control_event_rate","censoring_time")

allparamcombos <- lapply(allparamcombos, function(x) {
  x[numeric_fields] <- lapply(x[numeric_fields], as.numeric)
  x
})

# Run parallelization
plan(multisession, workers = 5)

results_all <- future_lapply(
  c(88,107,108),
  function(scen) {
    
    library(brms)
    
    # Compile brms models in advance
    formulas <- list(
      fit_ind = bf(event_ind ~ 1, family=bernoulli()),
      fit_trt = bf(event_ind ~ trt, family=bernoulli()),
      fit_survint = bf(time_surv | cens(1 - event_ind) ~ 1, family = exponential()),
      fit_weibint = bf(time_surv | cens(1 - event_ind) ~ 1, family = weibull()),
      fit_exp = bf(time_surv | cens(1 - event_ind) ~ trt, family = exponential()),
      fit_weibull = bf(time_surv | cens(1 - event_ind) ~ trt, family = weibull())
    )
    
    # Small data to compile
    small_dummy_data <- data.frame(
      event_ind = c(0, 1),
      trt = factor(c(0, 1), levels = c(0, 1)),
      time_surv = c(1, 2)
    )
    
    # Precompile models
    brms_models <- lapply(formulas, function(f) {
      brm(
        formula = f,
        data = small_dummy_data,
        chains = 0,
        backend = "rstan"
      )
    })
    
    vals <- allparamcombos[[scen]]
    
    sim_results <- lapply(
      1:nsim,
      function(iter) {
        runmodels(
          i = iter,
          N = vals$N,
          trt_effect_event = vals$trt_effect,
          distr = vals$distr,
          control_event_rate = vals$control_event_rate,
          censoring_time = vals$censoring_time,
          brms_models = brms_models
        )
      }
    )
    
    saveRDS(
      sim_results,
      file = paste0("rareoutcome_bothmodels_sim_results_scenario", scen, ".rds")
    )
    
    NULL
  },
  future.seed = TRUE,
  future.packages = 'brms'
)

plan(sequential)


# 1. Group the 10,000 flat lists into 10 settings (chunks of 1,000)
# We use rep() to create 10 groups of 1000 indices each
settings_grouped <- split(sim_results, rep(1:(length(allparamcombos)/nsamp), each = nsamp))

# 2. Process all scenarios
mean_with_na <- function(x_list) {
  # x_list: list of numeric data frames / matrices
  n_non_na <- Reduce(
    '+',
    lapply(x_list, function(x) !is.na(x))
  )
  
  summed <- Reduce(
    '+',
    lapply(x_list, function(x) replace(x, is.na(x), 0))
  )
  
  out <- summed / n_non_na
  out[n_non_na == 0] <- NA   # keep NA if no replication had data
  out
}

mean_predictions <- map(settings_grouped, function(one_setting) {
  
  # A. Use the first iteration as a template for scenario variables (1-9, 13-20)
  # This preserves all non-averaged elements exactly as they are
  setting_template <- one_setting[[1]][c(2,3,4,5,6,25)]
  
  # B. Calculate averages for the 3 target data frames (10, 11, 12)
  target_indices <- c(chisqpval = 9,
                      lrtpval = 10,
                      logistic_weights = 11,
                      logistic_stack = 15, logistic_pBMA = 16, logistic_BMA = 17,
                      survival_weights = 18,
                      survival_stack = 22, survival_pBMA = 23, survival_BMA = 24)
  
  averages <- map(target_indices, function(idx) {
    values <- map(one_setting, ~ .x[[idx]])
    mean_with_na(values)
  })
  
  
  return(c(setting_template, averages))
})

saveRDS(mean_predictions, file = "rareoutcome_bothmodels_sim_results_avgpred.RData")


# ##########################################################
# ################# SECTION 4: PLOT RESULTS ################
# ##########################################################
# 
# # Example: Compare Survival Stacking vs Logistic Stacking
# plot_comparison <- function(type, predictions, weights, title) {
#   #df_plot <- get_averaged_hazard(type, weights, new_data)
#   ggplot(predictions, aes(x = time_surv, y = prob_event, color = factor(trt), fill = factor(trt))) +
#     geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.3, color = NA) +
#     geom_line(linewidth = 1) +
#     {if(type == "survival")labs(title = title, 
#                                 subtitle = paste("Model Weights:", 
#                                                  "\nW_Int = ",     signif(weights["fit_survint"],3),
#                                                  ", W_WeibInt = ", signif(weights["fit_weibint"],3),
#                                                  "\nW_Exp = ",     signif(weights["fit_exp"],    3),
#                                                  ", W_Weib = ",    signif(weights['fit_weibull'],3), 
#                                                  sep = ""), 
#                                 y = "Hazard", x = "Time")} +
#     {if(type == "logistic")labs(title = title, 
#                                 subtitle = paste("Model Weights:", 
#                                                  "\nW_Int = ",      signif(weights["fit_int"],3),
#                                                  ", W_Trt = ",      signif(weights["fit_trt"],3), 
#                                                  sep = ""), 
#                                 y = "Event Probability", x = "Time")} +
#     {if(type == "survival")ylim(c(0,6.5))} +
#     {if(type == "logistic")ylim(c(0,1))} +
#     scale_color_discrete(name = "Group", labels = c("Pla.","Trt")) +
#     scale_fill_discrete(name = "Group", guide = FALSE) +
#     theme_minimal()
# }
# 
# 
# for(i in 1:length(mean_predictions)){
# 
#   df <- mean_predictions[[i]]
#   outplot <- (plot_comparison("survival", cbind(df$new_data,df$survival_stack),
#                               df$survival_weights["w_stack",], "Surv: Stacking") | 
#                 plot_comparison("survival", cbind(df$new_data,df$survival_pBMA),
#                                 df$survival_weights["w_pBMA",],  "Surv: Pseudo-BMA") |
#                 plot_comparison("survival", cbind(df$new_data,df$survival_BMA),
#                                 df$survival_weights["w_BMA",],   "Surv: BMA")) /
#     (plot_comparison("logistic", cbind(df$new_data,df$logistic_stack),
#                      df$logistic_weights["logistic_w_stack",], "Logit: Stacking") |
#        plot_comparison("logistic",  cbind(df$new_data,df$logistic_pBMA),
#                        df$logistic_weights["logistic_w_pBMA",],  "Logit: Pseudo-BMA") |
#        plot_comparison("logistic",  cbind(df$new_data,df$logistic_BMA),
#                        df$logistic_weights["logistic_w_BMA",],   "Logit: BMA"))
#   
#   
#   jpeg(paste("rareoutcome_bothmodels_AVGResults_","trteff", round(exp(df$trt_effect),2), 
#              "_N", df$N, 
#              "_", df$distr,
#              "_controleventrate", df$control_event_rate,
#              "_nsamp", nsamp,
#              ".jpg",sep = ""),
#        height = 800, width = 1000)
#   print(outplot +
#     plot_layout(guides = "collect") & 
#     theme(legend.position = "bottom") &
#     plot_annotation(title = as.character(paste("Average Predicted Values and Weights for Sim Trt Effect = ", round(exp(df$trt_effect),2), 
#                                                ", N = ", df$N, 
#                                                ", ", df$distr,
#                                                ", control_event_rate = ", df$control_event_rate,
#                                                ", nsamp = ", nsamp,
#                                                sep = ''))))
#   dev.off()
# 
# }
# 



