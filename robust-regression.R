####=====================  1. Setup  ======================####
library(haven)
library(tidyverse)
library(DescTools)
library(robustbase)
library(MASS)
library(labelled)
library(texreg)

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"
results_path <- "original-project"

decisions_df <- read_dta(file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta"))


####=====================  2. Functions  ======================####


####=====================  3. Main code  ======================####
### ---- Measures of location and spread of donations ---- ###
donation_stats_df <- decisions_df %>% 
  summarise(
    count = n(),
    mean = mean(h),
    median = median(h),
    trimmed_mean_10 = mean(h, trim = 0.10),
    hodges_lehman = HodgesLehmann(h),
    standard_deviation = sd(h),
    interquartile_range = IQR(h),
    median_absolute_deviation = mad(h),
    q_n = Qn(h)
  ) %>%
  pivot_longer(everything())


# Boxplots

# Hausman test
#   -> to decide whether to use robust methods or not

####=====================  4. Donation location (alternatives to the mean)  ======================####

# Alpha-trimmed mean


# Median

# Hodges-Lehmann estimator

####=====================  5. Donation spread (alternatives to the standard deviation)  ======================####

# Interquartile range

# Media absolute deviation

# The Q(n) coefficient

####=====================  6. Robust regression  ======================####

# Linear Fixed Effects Robust s.e. corresponds to assuming different variances at Low and High.

# Low govt provision:  $4 -> $10 (income $46 -> $40)
Low <- subset(decisions_df, Budget == 2 | Budget == 5)
reg_low <- rlm(h ~ Ggov, data=Low)
# Ho: crowd-out is -1 or a bigger negative number.
summary(reg_low)
texreg(reg_low, file = file.path(results_path,"Low_rob_reg_results.tex"))

# High govt provision:  $28 -> $34 (income $46 -> $40)
High <- subset(decisions_df, Budget == 4 | Budget == 6)
reg_high <- rlm(h ~ Ggov, data=High)
# Ho: crowd-out is -1 or a bigger negative number.
summary(reg_high)
texreg(reg_high, file = file.path(results_path,"High_rob_reg_results.tex"))

# STILL NEED TO INCLUDE FIXED-EFFECTS and ROBUST STANDARD ERRORS


# LS-estimator

# M-estimation

# LTS-estimation

# LMS-estimation

# S-estimation

# MM-estimation

# MS-estimation
