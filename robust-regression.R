####=====================  Setup  ======================####
library(haven)
library(lfe)
library(L1pack)
library(MASS)
library(systemfit)
library(tidyverse)

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"
results_path <- "original-project"

decisions_df <- read_dta(
  file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta")
)

####===================== Functions =======================####
get_t_stat <- function(coefficient, standard_error, h0) {
  t_stat <- (coefficient - h0) / standard_error
  return(t_stat)
}

get_coefficient <- function(model){
  model_summary  = model %>% summary() %>% coef() %>% head(1)
  crowd_out_rate = model_summary[1]
  return(crowd_out_rate)
}

get_standard_error <- function(model){
  model_summary  = model %>% summary() %>% coef() %>% head(1)
  standard_error = model_summary[2]
  return(standard_error)
}

get_p_value_of_full_crowd_out_rejection <- function(model) {
  coefficient    <- model %>% get_coefficient()
  standard_error <- model %>% get_standard_error()
  h0 <- -1
  t_stat <- get_t_stat(coefficient, standard_error, h0)
  
  # The number of degrees of freedom is the number of observations minus 1.
  # We count 85 observations because we have participant fixed effect.
  degrees_of_freedom <- 85 - 1
  
  p_value_one_tailed <- 1 - pt(t_stat, degrees_of_freedom)
  p_value_two_tailed <- p_value_one_tailed / 2
  
  return(p_value_two_tailed)
}


####=====================  Main code  ======================####
#### --------- Reproduction
# Column (1) and (2) of table 2
# Low govt provision:  $4 -> $10 (income $46 -> $40)
Low_df <- subset(decisions_df, Budget == 2 | Budget == 5)
low_reg <- felm(h ~ Ggov | newid, data=Low_df)
summary(low_reg)
low_reg %>% get_p_value_of_full_crowd_out_rejection()

texreg(low_reg, file = file.path(results_path,"Low_rob_reg_results.tex"))

# High govt provision:  $28 -> $34 (income $46 -> $40)
High_df <- subset(decisions_df, Budget == 4 | Budget == 6)
high_reg <- felm(h ~ Ggov | newid, data=High_df)
summary(high_reg)
high_reg %>% get_p_value_of_full_crowd_out_rejection()

texreg(high_reg, file = file.path(results_path,"High_rob_reg_results.tex"))


#### --------- Robust regressions
# Least Absolute Deviation (LAD) Regression
# BR = Barrodale & Roberts method (default)
LAD_low_reg <- lad(h ~ 0 + Ggov + factor(newid), data = Low_df, method = "BR") 
summary(LAD_low_reg)
LAD_low_reg %>% get_p_value_of_full_crowd_out_rejection()

LAD_high_reg <- lad(h ~ 0 + Ggov+ factor(newid), data = High_df, method = "BR")
summary(LAD_high_reg)


# M-estimation -> kappa standard 1.345
M_low1345_reg <- rlm(
  h ~ 0 + Ggov + factor(newid),
  data = Low_df,
  init = "ls",
  psi = psi.huber,
  scale.est = "Huber",
  k2 = 1.345,
  method = "M",
  maxit = 20,
  acc = 1e-4,
  test.vec = "resid",
  lqs.control = NULL
)

summary(M_low1345_reg)
M_low1345_reg %>% get_p_value_of_full_crowd_out_rejection()

texreg(M_low1345_reg, file = file.path(results_path,
                                       "M_low1345_reg_results.tex"))

M_high1345_reg <- rlm(
  h ~ 0 + Ggov + factor(newid),
  data = High_df,
  init = "ls",
  psi = psi.huber,
  scale.est = "Huber",
  k2 = 1.345,
  method = "M",
  maxit = 20,
  acc = 1e-4,
  test.vec = "resid",
  lqs.control = NULL
)

summary(M_high1345_reg)
M_high1345_reg %>% get_p_value_of_full_crowd_out_rejection()

texreg(M_high1345_reg, file = file.path(results_path, 
                                        "M_high1345_reg_results.tex"))

#### --------- Hausman test
# Hausman test: LAD vs. OLS
hausman.systemfit(LAD_high_reg,high_reg)
hausman.systemfit(LAD_low_reg,low_reg)

# Hausman test: LAD vs. M
hausman.systemfit(LAD_high_reg,M_high_reg)
hausman.systemfit(LAD_low_reg,M_low_reg)

# Hausman test: M vs. OLS
hausman.systemfit(M_high_reg,high_reg)
hausman.systemfit(M_low_reg,low_reg)
