#################################################################################
# TITLE: case_study_runMA1.R
#
# PURPOSE: create TTE data, run all models, and output final table with AE 
#          counts and analysis
#
# OUTPUT: case_study_summary.csv - case study table with all AEs, counts, and 
#                                  p-value/weight
#
# SECTIONS: Section 0 - load packages and read in data
#           Section 1 - process data
#           Section 2 - define tte data set function and make tte data sets
#           Section 3 - define functions for standard tests and model avg. with
#                       count or survival data
#           Section 4 - initialize table
#           Section 5 - run models on all tte data sets
#           Section 6 - create and export table
#           Section 7 - optional ggplot comparing weights from all methods 
#
# AUTHOR: Shannon Thomas
# NOTES: 
#################################################################################

##########################################################
################ SECTION 0 : LOAD PACKAGES ############### 
##########################################################

library(haven)
library(dplyr)
library(tidyr)
library(survival)
library(brms)
library(loo)
library(lmtest)
library(tidyverse)

set.seed(12345)

#read in data (available on Project Data Sphere)
data_dir <- "C:/Users/mushanno/Downloads/AllProvidedFiles_266/SAS dataset - 20010145/SAS dataset"

c_keyvar <- read_sas(file.path(data_dir, "C_KEYVAR.sas7bdat"))
c_ae     <- read_sas(file.path(data_dir, "C_AE.sas7bdat"))

############################################################
## SECTION 1: PROCESS DATA
############################################################

pop_safety <- c_keyvar %>%
  filter(!is.na(SAFGROUP)) %>%
  select(SUBJID, SAFGROUP, B_ECOGN, B_LDHN)

ae <- c_ae %>%
  inner_join(pop_safety, by = c("SUBJID","SAFGROUP")) %>%
  filter(AEYN == 1, TEAEYN == 1) %>%
  mutate(
    SERIOUS = SRIOUSYN == 1,
    SEVERE  = SEVRCD %in% c(3,4),
    FATAL   = SEVRCD == 5,
    TRTREL  = RELATEYN == 1
  )

safe_min <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) min(x) else NA_real_
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) max(x) else NA_real_
}


cv_terms <- c("EMBOL","THROMB","MYOCARD","STROKE",
              "CEREBROVASC","ARRHYTHM","HEART FAILURE","HYPERTENS")

ae_subject_flags <- ae %>%
  group_by(SUBJID, SAFGROUP) %>%
  summarise(
    any_ae      = TRUE,
    any_serious = any(SERIOUS),
    any_severe  = any(SEVERE | FATAL),
    any_trtrel  = any(TRTREL),
    
    # ---- ADD THESE ----
    any_cv   = any(grepl(paste(cv_terms, collapse="|"), AEPTERM, ignore.case=TRUE)),
    any_arr  = any(grepl("ARRHYTHM", AEPTERM, ignore.case=TRUE)),
    any_cva  = any(grepl("STROKE|CEREBROVASC", AEPTERM, ignore.case=TRUE)),
    any_mi   = any(grepl("MYOCARD|CORONARY", AEPTERM, ignore.case=TRUE)),
    any_emb  = any(grepl("EMBOL|THROMB", AEPTERM, ignore.case=TRUE)),
    any_htn  = any(grepl("HYPERTENS", AEPTERM, ignore.case=TRUE)),
    
    .groups = "drop"
  )

############################################################
## SECTION 2: GENERIC TTE BUILDER
############################################################

build_tte <- function(flag_var = NULL, terms = NULL) {
  
  ## STEP 1: BUILD EVENTS 
  if (!is.null(flag_var)) {
    
    # subject-level definition
    events <- ae_subject_flags %>%
      filter(.data[[flag_var]]) %>%
      select(SUBJID) %>%   # ✅ keep ONLY SUBJID
      left_join(
        ae %>%
          group_by(SUBJID) %>%
          summarise(DAY = safe_min(STUDYDAY), .groups = "drop"),
        by = "SUBJID"
      ) %>%
      mutate(EVENT = 1L)
    
  } else {
    
    # term-based definition
    events <- ae %>%
      filter(grepl(paste(terms, collapse = "|"),
                   AEPTERM, ignore.case = TRUE)) %>%
      group_by(SUBJID) %>%
      summarise(
        DAY = safe_min(STUDYDAY),
        EVENT = 1L,
        .groups = "drop"
      )
  }
  
  ## STEP 2: CENSORING 
  last_day <- ae %>%
    group_by(SUBJID) %>%
    summarise(
      LAST_DAY = safe_max(STUDYDAY),
      .groups = "drop"
    )
  
  global_censor <- max(ae$STUDYDAY, na.rm = TRUE)
  
  ## STEP 3: ATTACH TO POPULATION 
  tte <- pop_safety %>%
    select(SUBJID, SAFGROUP, B_ECOGN, B_LDHN) %>%  # ✅ ONLY here
    left_join(events, by = "SUBJID") %>%
    left_join(last_day, by = "SUBJID") %>%
    mutate(
      EVENT = ifelse(is.na(EVENT), 0L, EVENT),
      time  = ifelse(EVENT == 1,
                     DAY,
                     ifelse(!is.na(LAST_DAY), LAST_DAY, global_censor))
    ) %>%
    mutate(SAFGROUP = factor(SAFGROUP))   # ✅ ensure proper modeling
  
  return(tte)
}


## BUILD INDIVIDUAL TTE DATASETS
tte_all    <- build_tte(flag_var="any_ae")
tte_serious<- build_tte(flag_var="any_serious")
tte_severe <- build_tte(flag_var="any_severe")
tte_trtrel <- build_tte(flag_var="any_trtrel")

tte_cv   <- build_tte(terms=cv_terms)
tte_arr  <- build_tte(terms=c("ARRHYTHM"))
tte_cva  <- build_tte(terms=c("STROKE","CEREBROVASC"))
tte_mi   <- build_tte(terms=c("MYOCARD","CORONARY"))
tte_emb  <- build_tte(terms=c("EMBOL","THROMB"))
tte_htn  <- build_tte(terms=c("HYPERTENS"))

############################################################
## SECTION 3: TEST FUNCTIONS
############################################################

make_templates <- function(N) {
  
  # minimal dummy data (just to compile model)
  dummy_dat_count <- data.frame(
    SAFGROUP = factor(rep(c(0,1), each = N/2)),
    EVENT = rbinom(N, 1, 0.5)
  )
  
  dummy_dat_surv <- data.frame(
    SAFGROUP = factor(rep(c(0,1), each = N/2)),
    time = rexp(N, 1),
    EVENT = rbinom(N, 1, 0.5)
  )
  
  #COUNT MODELS
  fit_count_int <- brm(
    EVENT ~ 1,
    data = dummy_dat_count,
    family = bernoulli(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_count_trt <- brm(
    EVENT ~ SAFGROUP,
    data = dummy_dat_count,
    family = bernoulli(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  #SURVIVAL MODELS
  fit_survint_exp <- brm(
    time | cens(1 - EVENT) ~ 1,
    data = dummy_dat_surv,
    family = exponential(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_survint_weib <- brm(
    time | cens(1 - EVENT) ~ 1,
    data = dummy_dat_surv,
    family = weibull(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_exp_trt <- brm(
    time | cens(1 - EVENT) ~ SAFGROUP,
    data = dummy_dat_surv,
    family = exponential(),
    chains = 1, iter = 300, refresh = 0, silent = 2
  )
  
  fit_weib_trt <- brm(
    time | cens(1 - EVENT) ~ SAFGROUP,
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

#standard tests
logistic_tests <- function(df){
  
  m1 <- glm(EVENT ~ SAFGROUP, df, family=binomial)
  m0 <- glm(EVENT ~ 1, df, family=binomial)
  
  lrt <- lrtest(m0,m1)$`Pr(>Chisq)`[2]
  
  tab <- table(df$SAFGROUP,df$EVENT)
  
  chi <- tryCatch(chisq.test(tab)$p.value,
                  warning=function(w) fisher.test(tab)$p.value)
  
  c(lrt_p=lrt, chi_p=chi)
}

#model averaging for count data
get_logistic_weights_unadj <- function(df, templates) {
  
  fit0 <- update(templates$count_int,
                    newdata = df,
                    recompile = FALSE)
  
  fit1 <- update(templates$count_trt,
                    newdata = df,
                    recompile = FALSE)
  
  list(
    stack = model_weights(fit0, fit1, weights = "stacking"),
    pBMA  = model_weights(fit0, fit1, weights = "pseudobma"),
    BMA   = model_weights(fit0, fit1, weights = "bma")
  )
}

#model averaging for survival data
get_survival_weights_unadj <- function(df, templates) {

  
  fit_exp0 <- update(templates$survint_exp,
                     newdata = df,
                     recompile = FALSE)
  
  fit_weib0 <- update(templates$survint_weib,
                      newdata = df,
                      recompile = FALSE)
  
  fit_exp1 <- update(templates$exp_trt,
                     newdata = df,
                     recompile = FALSE)
  
  fit_weib1 <- update(templates$weib_trt,
                      newdata = df,
                      recompile = FALSE)
  
  list(
    stack = model_weights(fit_exp0, fit_weib0, fit_exp1, fit_weib1, weights = "stacking"),
    pBMA  = model_weights(fit_exp0, fit_weib0, fit_exp1, fit_weib1, weights = "pseudobma"),
    BMA   = model_weights(fit_exp0, fit_weib0, fit_exp1, fit_weib1, weights = "bma")
  )
}


############################################################
## SECTION 4 : INITIALIZE TABLE
############################################################

denoms <- pop_safety %>% count(SAFGROUP,name="N")

t3_counts <- ae_subject_flags %>%
  pivot_longer(starts_with("any_"),
               names_to="type",values_to="flag") %>%
  filter(flag) %>%
  count(SAFGROUP,type,name="events") %>%
  left_join(denoms,by="SAFGROUP") %>%
  mutate(
    percent = 100*events/N,
    value   = sprintf("%d (%.0f%%)",events,percent),
    Category = recode(type,
                      any_ae     = "All adverse events",
                      any_serious= "Serious",
                      any_severe = "Severe, life-threatening, or fatal",
                      any_trtrel = "Treatment related",
                      any_cv     = "CV",
                      any_arr    = "Arrhythmia",
                      any_cva    = "CVA",
                      any_mi     = "MI",
                      any_emb    = "Embolism",
                      any_htn    = "Hypertension"
    )
  ) %>%
  select(Category,SAFGROUP,value) %>%
  pivot_wider(names_from=SAFGROUP,values_from=value)

############################################################
## SECTION 5: RUN MODELS
############################################################

run_all <- function(df, label){
  
  ## standard tests 
  lt <- logistic_tests(df)
  
  #generate MA templates
  templates <- make_templates(N = 100)
  
  ## WEIGHTS
  lw_unadj <- get_logistic_weights_unadj(df, templates)   # length 2
  sw_unadj <- get_survival_weights_unadj(df, templates)   # length 4

  #output
  data.frame(
    Category = label,
    
    LRT_p   = as.numeric(lt["lrt_p"]),
    ChiSq_p = as.numeric(lt["chi_p"]),
    
    ###### LOGISTIC UNADJUSTED ######
    logit_unadj_stack_null = lw_unadj$stack[1],
    logit_unadj_stack_trt  = lw_unadj$stack[2],
    
    logit_unadj_pBMA_null = lw_unadj$pBMA[1],
    logit_unadj_pBMA_trt  = lw_unadj$pBMA[2],
    
    logit_unadj_BMA_null = lw_unadj$BMA[1],
    logit_unadj_BMA_trt  = lw_unadj$BMA[2],
    
    ###### SURVIVAL UNADJUSTED (4 MODELS) ######
    surv_unadj_stack_exp_null = sw_unadj$stack[1],
    surv_unadj_stack_weib_null = sw_unadj$stack[2],
    surv_unadj_stack_exp_trt  = sw_unadj$stack[3],
    surv_unadj_stack_weib_trt = sw_unadj$stack[4],
    
    surv_unadj_pBMA_exp_null = sw_unadj$pBMA[1],
    surv_unadj_pBMA_weib_null = sw_unadj$pBMA[2],
    surv_unadj_pBMA_exp_trt  = sw_unadj$pBMA[3],
    surv_unadj_pBMA_weib_trt = sw_unadj$pBMA[4],
    
    surv_unadj_BMA_exp_null = sw_unadj$BMA[1],
    surv_unadj_BMA_weib_null = sw_unadj$BMA[2],
    surv_unadj_BMA_exp_trt  = sw_unadj$BMA[3],
    surv_unadj_BMA_weib_trt = sw_unadj$BMA[4]
  )
}

#run models
set.seed(303)
allae <- run_all(tte_all,"All adverse events")
serae <- run_all(tte_serious,"Serious")
sevae <- run_all(tte_severe,"Severe")
trtrelae <- run_all(tte_trtrel,"Treatment related")
cvae <- run_all(tte_cv,"CV")
arrythae <- run_all(tte_arr,"Arrhythmia")
cvaae <- run_all(tte_cva,"CVA")
miae <- run_all(tte_mi,"MI")
embae <- run_all(tte_emb,"Embolism")
hypae <-   run_all(tte_htn,"Hypertension")

#create results df
results <- rbind(allae, serae, sevae, trtrelae, cvae, arrythae, cvaae, miae, embae, hypae)


############################################################
## SECTION 6: CREATE TABLE WITH ALL WEIGHTS OF INTEREST
############################################################

table3_full <- t3_counts %>%
  left_join(results,by="Category")

table3_full


## COLLAPSE EXISTING WEIGHTS (NO MODEL RERUN NEEDED)
results_collapsed <- results %>%
  mutate(
    
    #### SURVIVAL (SUM EXP + WEIBULL TREATMENT MODELS)
    surv_unadj_trt_weight_stack =
      surv_unadj_stack_exp_trt + surv_unadj_stack_weib_trt,
    
    surv_unadj_trt_weight_pBMA =
      surv_unadj_pBMA_exp_trt + surv_unadj_pBMA_weib_trt,
    
    surv_unadj_trt_weight_BMA =
      surv_unadj_BMA_exp_trt + surv_unadj_BMA_weib_trt,
  )

results_compact <- results_collapsed %>%
  select(
    -contains("_null"),
    -contains("_exp_"),
    -contains("_weib_")
  )

results_compact$Category[results_compact$Category == "Severe"] <- "Severe, life-threatening, or fatal"

table3_final_collapsed <- t3_counts %>%
  left_join(results_compact, by = "Category")

#View(table3_final_collapsed %>% select(-contains("unadj")))
#View(table3_final_collapsed %>% select("Category", "NESP","PLACEBO", "LRT_p", "ChiSq_p",contains("unadj")))

write.csv(table3_final_collapsed, "case_study_summary.csv")


############################################################
## SECTION 7: (OPTIONAL) PLOT RESULTS TO COMPARE WEIGHTS
############################################################

table3_long <- pivot_longer(data = table3_final_collapsed, 
                            cols = contains(c("logit","surv")),
                            )

table3_long$endpoint <- ifelse(str_detect(table3_long$name, "surv"), "survival", "logistic")
table3_long$model_adj <- ifelse(str_detect(table3_long$name, "_unadj"), "unadjusted", "adjusted")
table3_long$model_avg_type <- ifelse(str_detect(table3_long$name, "stack"), "stacking",
                                     ifelse(str_detect(table3_long$name, "pBMA"), 
                                            "pseudo-BMA","BMA"))
  
  
ggplot(data = table3_long %>%
         filter(str_detect(name, "_unadj")), 
       aes(x = endpoint, y = value, shape = endpoint, color = model_avg_type)) +
  geom_point() + facet_wrap(~paste(Category,
                                   "\nNESP:",NESP,
                                   "\nPLACEBO:",PLACEBO)
                            +paste("ChiSq p-value =",signif(ChiSq_p,2)),
                            ncol = 5) +
  ggtitle("Model Results") + ylab("Sum of Weights for Treament Models")







