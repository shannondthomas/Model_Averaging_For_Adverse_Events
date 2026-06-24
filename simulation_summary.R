################################################################################
# TITLE: simulation_summary.R
#
# PURPOSE: process simulation results, compute performance metrics, compare 
#          simulated methods, run CoxPH power calculations, and generate
#          plots/tables for publication.
#
# INPUT: all .RData files in sim_results folder
#        (rareoutcome_bothmodels_sum_results_scenario(1-108).RData or .rds)
#
# OUTPUT: Simulation_Hazard_Functions.jpeg 
#         Power_CalibratedThresholds_withCoxPH_highlightbest_nsamp500.jpeg
#         Pareto_FrontCurve_nsamp500.jpeg
#         calibrated_results_best_table_wCoxPH.html
#
# SECTIONS: Section 0 - load packages and data
#           Section 1 - create simulation results data frame
#           Section 2 - define calibrated cutoffs and corresponding power
#           Section 3 - perform CoxPH power simulations
#           Section 4 - create dfs for tables and figures
#           Section 5 - create risk distribution figure
#           Section 6 - create power figure and Pareto front curve figure
#           Section 7 - create table of Pareto-optimal methods
#
# AUTHOR: Shannon Thomas
#
# NOTES: The first 9 simulation scenarios were output as .RData files rather 
#        than .rds files and have additional variables. 
################################################################################

############################################################
############ SECTION 0: LOAD PACKAGES AND DATA ##############
############################################################

library(tidyverse)
library(purrr)
library(stringr)
library(patchwork)
library(kableExtra)
library(powerSurvEpi)
library(survival)
library(future)
library(furrr)

##### Load simulation files

path <- "C:/Users/mushanno/OneDrive - The University of Colorado Denver/Desktop/Work/Dissertation3/Model_Averaging_For_Adverse_Events/sim_results"
files <- list.files(path, full.names = TRUE)

read_any <- function(f) {
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.null(obj)) return(obj)
  
  obj <- tryCatch({
    e <- new.env()
    load(f, envir = e)
    as.list(e)
  }, error = function(e) NULL)
  
  if (!is.null(obj)) return(obj)
  stop(paste("Cannot read file:", f))
}

all_data <- map(files, read_any)
names(all_data) <- str_extract(basename(files), "\\d+")
all_data <- all_data[order(as.integer(names(all_data)))]

nsamp <- length(all_data[[1]])

############################################################
########### SECTION 1: CREATE RESULTS DATA FRAME ########### 
############################################################

allweights_df <- map_dfr(seq_along(all_data), function(i) {
  
  scenario <- all_data[[i]]
  
  map_dfr(seq_along(scenario), function(j) {
    
    sim <- scenario[[j]]
    
    base_info <- tibble(
      scenario = i,
      sim_draw = j,
      N = sim$N,
      distr = sim$distr,
      trt_effect = paste("Trt Effect =", round(exp(sim$trt_effect), 3)),
      control_event_rate = sim$control_event_rate
    )
    
    # Standard tests
    chisq <- sim$chisqpval %||% NA_real_
    lrt   <- sim$lrtpval  %||% NA_real_
    
    # Logistic MA
    logit <- sim$logistic_weights
    logit_stack <- tryCatch(logit["logistic_w_stack","fit_trt"], error=function(e) NA)
    logit_pBMA  <- tryCatch(logit["logistic_w_pBMA","fit_trt"], error=function(e) NA)
    logit_BMA   <- tryCatch(logit["logistic_w_BMA","fit_trt"], error=function(e) NA)
    
    # Survival MA
    surv <- sim$survival_weights
    surv_stack <- tryCatch(surv["w_stack","fit_exp"] + surv["w_stack","fit_weibull"], error=function(e) NA)
    surv_pBMA  <- tryCatch(surv["w_pBMA","fit_exp"] + surv["w_pBMA","fit_weibull"], error=function(e) NA)
    surv_BMA   <- tryCatch(surv["w_BMA","fit_exp"] + surv["w_BMA","fit_weibull"], error=function(e) NA)
    
    bind_rows(
      base_info %>% mutate(model="standard", method="chisq", trtmodel_weight=chisq),
      base_info %>% mutate(model="standard", method="LRT",   trtmodel_weight=lrt),
      
      base_info %>% mutate(model="logit", method="stacking",   trtmodel_weight=logit_stack),
      base_info %>% mutate(model="logit", method="pseudo-BMA", trtmodel_weight=logit_pBMA),
      base_info %>% mutate(model="logit", method="BMA",        trtmodel_weight=logit_BMA),
      
      base_info %>% mutate(model="survival", method="stacking",   trtmodel_weight=surv_stack),
      base_info %>% mutate(model="survival", method="pseudo-BMA", trtmodel_weight=surv_pBMA),
      base_info %>% mutate(model="survival", method="BMA",        trtmodel_weight=surv_BMA)
    )
    
  })
})


##### convert NA scenarios to "non-reject" weight/p values

allweights_df <- allweights_df %>%
  mutate(
    trt_effect_n = as.numeric(str_extract(trt_effect, "\\d+\\.?\\d*")),
    
    # NA handling
    trtmodel_weight_noNA =
      case_when(
        model == "standard" ~ replace_na(trtmodel_weight, 1),
        TRUE                ~ replace_na(trtmodel_weight, 0)
      )
  )

############################################################
############ SECTION 2: PERFORMANCE METRICS ################
############################################################

scenario_vars <- c("N","distr","control_event_rate","model","method")

##### Calibrated Thresholds for 0.05 Type 1 Error Rate
thresholds_df <- allweights_df %>%
  filter(trt_effect_n == 1) %>%
  group_by(across(all_of(scenario_vars))) %>%
  summarize(
    threshold = if (
      first(model) == "standard"
    ) {
      0.05
    } else {
      quantile(trtmodel_weight_noNA, 0.95)
    },
    .groups = "drop"
  )

##### Compute Power
perf_calibrated <- allweights_df %>%
  left_join(thresholds_df, by = scenario_vars) %>%
  mutate(
    decision = if_else(
      model == "standard",
      trtmodel_weight_noNA < threshold,
      trtmodel_weight_noNA > threshold
    )
  ) %>%
  group_by(across(all_of(scenario_vars)), trt_effect_n) %>%
  summarize(prob_reject = mean(decision), .groups="drop") %>%
  mutate(threshold_type = "calibrated")

############################################################
################## SECTION 3: COX PH POWER ################# 
############################################################

set.seed(406)
alpha <- 0.05
nsim_cox <- nsamp   # match simulation count

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

simulate_cox_pval <- function(N, logHR, distr, control_event_rate, censoring_time){
  
  X_trt <- rep(c(0,1), each = N/2)
  
  if(distr == "exponential"){
    hr <- exp(logHR)
    
    lambda_c <- lambda_exponential(control_event_rate, censoring_time)
    lambda_t <- lambda_c * hr
    
    time_c <- rexp(n = N/2, rate = lambda_c)
    time_t <- rexp(n = N/2, rate = lambda_t)
  }
  
  if(distr == "weibull"){
    hr <- exp(logHR)
    shape_par <- 2
    
    scale_c <- scale_weibull(control_event_rate, censoring_time, shape_par)
    scale_t <- scale_c / hr^(1 / shape_par)
    
    time_c <- rweibull(n = N/2, shape = shape_par, scale = scale_c)
    time_t <- rweibull(n = N/2, shape = shape_par, scale = scale_t)
  }
  
  if(distr == "lognormal"){
    sigma_ln <- 0.7
    
    mu_c <- mu_lognormal(control_event_rate, censoring_time, sigma_ln)
    mu_t <- mu_c - logHR   # matches your DGP
    
    time_c <- rlnorm(n = N/2, meanlog = mu_c, sdlog = sigma_ln)
    time_t <- rlnorm(n = N/2, meanlog = mu_t, sdlog = sigma_ln)
  }
  
  event_times <- c(time_c, time_t)
  
  censor_times <- runif(n = N, min = 0, max = censoring_time)
  censor_times[censor_times >= censoring_time] <- censoring_time
  
  dat <- data.frame(
    trt = factor(X_trt, levels = c(0,1)),
    time_e = event_times,
    time_c = censor_times
  )
  
  dat$time_surv <- pmin(dat$time_e, dat$time_c)
  dat$event_ind <- ifelse(dat$time_surv == dat$time_e, 1, 0)
  
  # Fit Cox PH
  fit <- try(coxph(Surv(time_surv, event_ind) ~ trt, data = dat), silent = TRUE)
  
  if(inherits(fit, "try-error")) return(NA_real_)
  
  summary(fit)$coef[1,5]
}


options(future.apply.debug=TRUE)
options(future.globals.onReference = "error")
plan(multisession, workers = 5)

cox_results <- allweights_df %>%
  distinct(N, distr, control_event_rate, trt_effect_n) %>%
  mutate(
    logHR = log(trt_effect_n),
    
    prob_reject = future_pmap_dbl(
      list(N, logHR, distr, control_event_rate),
      function(N, logHR, distr, control_event_rate){
        
        pvals <- replicate(
          nsim_cox,
          simulate_cox_pval(N, logHR, distr, control_event_rate, 3)
        )
        
        pvals[is.na(pvals)] <- 1
        mean(pvals < alpha)
        
      },
      .options = furrr_options(seed = TRUE)
      
    )
  ) %>%
  mutate(
    model = "standard",
    method = "CoxPH",
    threshold_type = "calibrated"
  ) %>%
  select(N, distr, control_event_rate, trt_effect_n,
         model, method, threshold_type, prob_reject)

plan(sequential)

############################################################
###### SECTION 4: MERGE + BUILD TRADEOFF AND POWER DFS #####
############################################################

results_all <- bind_rows(perf_calibrated, cox_results)

# Build tradeoff df
tradeoff_df <- results_all %>%
  group_by(N,distr,control_event_rate,model,method,threshold_type,trt_effect_n) %>%
  summarize(prob_reject=mean(prob_reject), .groups="drop") %>%
  pivot_wider(names_from=trt_effect_n,
              values_from=prob_reject,
              names_prefix="effect_")

find_pareto_3d <- function(df){
  df %>%
    rowwise() %>%
    mutate(
      dominated = any(
        (df$effect_1 <= effect_1 &
           df$effect_1.1 >= effect_1.1 &
           df$effect_2 >= effect_2) &
          (df$effect_1 < effect_1 |
             df$effect_1.1 > effect_1.1 |
             df$effect_2 > effect_2)
      )
    ) %>%
    ungroup() %>%
    filter(!dominated)
}

pareto_results <- tradeoff_df %>%
  filter(!is.na(effect_1), !is.na(effect_1.1)) %>%
  group_by(N,distr,control_event_rate,threshold_type) %>%
  group_modify(~find_pareto_3d(.x)) %>%
  ungroup()

pareto_results <- pareto_results %>%
  arrange(N, distr, control_event_rate, effect_1) %>%
  mutate(distr = factor(distr, levels = c("exponential","weibull","lognormal")))

# Build power df
power_plot_df <- results_all %>%
  
  filter(trt_effect_n != 1,
         threshold_type == "calibrated") %>%
  
  mutate(
    model_method = paste(model, ":", method),
    
    model_method = factor(
      model_method,
      levels = c(
        "standard : chisq",
        "standard : LRT",
        "logit : BMA",
        "logit : pseudo-BMA",
        "logit : stacking",
        "standard : CoxPH",
        "survival : BMA",
        "survival : pseudo-BMA",
        "survival : stacking"
      )
    ),
    distr = factor(distr, levels = c("exponential","weibull","lognormal"))
  )


power_plot_df_best <- power_plot_df %>%
  group_by(distr, N, control_event_rate, trt_effect_n) %>%
  mutate(is_best = (prob_reject >= max(prob_reject) - 1e-10)) %>%
  ungroup()


############################################################
############ SECTION 5: HAZARD FUNCTION FIGURE #############
############################################################

set.seed(805)

# Lognormal hazard = f(t) / S(t)
lognormal_hazard <- function(t, mu, sigma) {
  dlnorm(t, meanlog = mu, sdlog = sigma) /
    (1 - plnorm(t, meanlog = mu, sdlog = sigma))
}

censoring_time <- 3
design_vals <- expand_grid(
  distr = c("exponential", "weibull", "lognormal"),
  control_event_rate = c(0.1, 0.99),
  trt_effect = c(1, 1.1, 2)
)
t_grid <- seq(0.01, censoring_time, length.out = 200)

# ---- Build hazard data ----
hazard_df <- design_vals %>%
  mutate(distr = factor(distr, levels = c("exponential","weibull","lognormal"))) %>%
  crossing(
    t = t_grid,
    group = c("Control","Treatment")
  ) %>%
  mutate(
    logHR = log(trt_effect),
    HR = trt_effect,
    
    hazard = case_when(
      
      # ---- EXPONENTIAL ----
      distr == "exponential" & group == "Control" ~ {
        lambda_exponential(control_event_rate, censoring_time)
      },
      
      distr == "exponential" & group == "Treatment" ~ {
        lambda_exponential(control_event_rate, censoring_time) * HR
      },
      
      # ---- WEIBULL ----
      distr == "weibull" & group == "Control" ~ {
        shape_par <- 2
        scale_c <- scale_weibull(control_event_rate, censoring_time, shape_par)
        (shape_par / scale_c) * (t / scale_c)^(shape_par - 1)
      },
      
      distr == "weibull" & group == "Treatment" ~ {
        shape_par <- 2
        scale_c <- scale_weibull(control_event_rate, censoring_time, shape_par)
        scale_t <- scale_c / HR^(1 / shape_par)
        (shape_par / scale_t) * (t / scale_t)^(shape_par - 1)
      },
      
      # ---- LOGNORMAL ----
      distr == "lognormal" & group == "Control" ~ {
        sigma_ln <- 0.7
        mu_c <- mu_lognormal(control_event_rate, censoring_time, sigma_ln)
        lognormal_hazard(t, mu_c, sigma_ln)
      },
      
      distr == "lognormal" & group == "Treatment" ~ {
        sigma_ln <- 0.7
        mu_c <- mu_lognormal(control_event_rate, censoring_time, sigma_ln)
        mu_t <- mu_c - logHR
        lognormal_hazard(t, mu_t, sigma_ln)
      }
    )
  )


hazard_df_plot <- hazard_df %>%
  mutate(
    curve_label = case_when(
      group == "Control" ~ "Control",
      TRUE ~ paste0("Treatment (Trt Effect=", round(HR, 2), ")")
    ),
    curve_label = factor(curve_label,
                         levels = c("Control",
                                    "Treatment (Trt Effect=1)",
                                    "Treatment (Trt Effect=1.1)",
                                    "Treatment (Trt Effect=2)"))
  )

# ---- Plot ----
jpeg("Figures/Simulation_Hazard_Functions.jpg",
     height = 700, width = 1000)

ggplot(hazard_df_plot,
       aes(x = t, y = hazard,
           color = curve_label,
           linetype = curve_label)) +
  
  geom_line(linewidth = 1) +
  
  facet_grid(
    rows = vars(control_event_rate),
    cols = vars(distr),
    scales = "free_y",
    labeller = labeller(
      control_event_rate = function(x) paste("ER =", x)
    )
  ) +
  
  labs(
    title = "Hazard Functions for Each Data Generating Mechanism",
    x = "Time",
    y = "Hazard",
    color = "Curve",
    linetype = "Curve"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )


dev.off()


############################################################
############ SECTION 6: SUMMARY PLOTS ######################
############################################################

jpeg(paste("Figures/Power_CalibratedThresholds_withCoxPH_highlightbest",
           "_nsamp", nsamp,
           ".jpg",sep = ""),
     height = 800, width = 1200)
ggplot(power_plot_df_best,aes(x = model_method,y = prob_reject,
                         fill = model_method)) +
  # Base layer: all bars with NO outline
  geom_col(position = "dodge", colour = NA, width = 0.8) +
  # # Overlay layer: ONLY standard bars, with black outline
  # geom_col(
  #   data = subset(power_plot_df_best, model == "standard"),
  #   aes(x = model_method, y = prob_reject),
  #   position = position_dodge(width = 0.6),
  #   fill = NA,          # don't redraw fill
  #   colour = "black",
  #   linewidth = 0.4) +
  geom_col(
    data = power_plot_df_best %>% dplyr::filter(is_best),
    aes(x = model_method, y = prob_reject),
    position = position_dodge(width = 0.8),
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  facet_grid(
    rows = vars(control_event_rate,trt_effect_n),
    cols = vars(distr, N),
    labeller = labeller(
      control_event_rate = function(x) paste("ER =", x),
      trt_effect_n = function(x) paste("HR =", x),
      N = function(x) paste("N =", x))) +
  labs(title = "Power Comparison Across Methods and Scenarios",
       x = "",y = "Power",fill = "Model Type") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()

jpeg(paste("Figures/Pareto_FrontCurve",
           "_nsamp", nsamp,
           ".jpg",sep = ""),
     height = 500, width = 800)
ggplot(pareto_results,
       aes(x = effect_1, y = effect_1.1)) +
  
  geom_point(aes(color = method, shape = model), size = 3) +
  
  geom_line(
    data = pareto_results %>%
      group_by(N, distr, control_event_rate) %>%
      filter(n() > 1),   # ✅ only draw when >1 point
    aes(group = 1),
    linewidth = 0.8
  ) +
  
  facet_grid(rows = vars(N),
             cols = vars(control_event_rate, distr),
             scales = "free_x") +
  
  labs(
    x = "Type I Error",
    y = "Power",
    title = "Pareto Front Including Cox PH"
  )
dev.off()

############################################################
############ SECTION 7: TABLE ##############################
############################################################

pareto_table <- pareto_results %>%
  filter(threshold_type=="calibrated") %>%
  select(N,distr,control_event_rate,model,method,
         effect_1,effect_1.1,effect_2) %>%
  arrange(N,distr,control_event_rate,
          desc(effect_1.1),desc(effect_2))

pareto_kable <- pareto_table %>%
  kable(caption="Pareto-Optimal Methods Including Cox PH",
        digits=3) %>%
  kable_styling(full_width=FALSE) %>%
  collapse_rows(columns=c(1,2,3))

save_kable(pareto_kable, "calibrated_results_best_table_wCoxPH.html")

