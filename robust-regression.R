####=====================  1. Setup  ======================####
library(haven)
library(tidyverse)
library(DescTools)
library(robustbase)
library(MASS)
library(labelled)
library(texreg)
library(depthTools)
library(xtable)
library(L1pack)
#install.packages("")



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
# R shows no error but export does not work
print(xtable(donation_stats_df, type = "latex", file = file.path(results_path,"Donation_stats.tex")))



# Boxplots

# Hausman test
#   -> to decide whether to use robust methods or not

####=====================  4. Donation location (alternatives to the mean)  ======================####

# Alpha-trimmed mean (Look if Remi's version not preferable)
alpha <- tmean(decisions_df,alpha=0.2,plotting=FALSE,new=TRUE,cols=c(1,4,8))

# Median
median <- decisions_df %>% summarise(median = median(h)) %>%pivot_longer(everything())

# Hodges-Lehmann estimator
Hod_Leh_est <- decisions_df %>% summarise(hodges_lehman = HodgesLehmann(h)) %>% pivot_longer(everything())

####=====================  5. Donation spread (alternatives to the standard deviation)  ======================####

# Interquartile range
int_quart_range <- decisions_df %>% summarise(interquartile_range = IQR(h)) %>% pivot_longer(everything())

# Median absolute deviation
med_abs_devia <- decisions_df %>% summarise(median_absolute_deviation = mad(h)) %>% pivot_longer(everything())

# The Q(n) coefficient
Q_n <- decisions_df %>% summarise(q_n = Qn(h)) %>% pivot_longer(everything())

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


# LS-estimator -> Least Square Regression
LS_reg <- lm(h ~ Ggov, data = decisions_df)
summary(LS_reg)
texreg(LS_reg, file = file.path(results_path,"LS_reg_results.tex"))

# M-estimation
M_reg <- lad(h~ Ggov, data = decisions_df, method = "BR") # BR = Barrodale & Roberts method (default)
summary(M_reg)
texreg(M_reg, file = file.path(results_path,"M_reg_results.tex")) # EXPORT DOESN'T WORK

# LTS-estimation -> Least Trimmed Sum of Squares
LTS_reg <- ltsReg(h ~ Ggov, data = decisions_df)
summary(LTS_reg)
texreg(LTS_reg, file = file.path(results_path,"LTS_reg_results.tex")) # EXPORT DOESN'T WORK


# LMS-estimation -> Least Median Squares Regression
LMS_reg <- lmsreg(h ~ Ggov, data = decisions_df)
summary(LMS_reg)
texreg(LMS_reg, file = file.path(results_path,"LMS_reg_results.tex")) # EXPORT DOESN'T WORK

# S-estimation
S_reg <- lmrob.S(x=decisions_df$Ggov, y=decisions_df$h, control = lmrob.control(nRes = 20), trace.lev=1)

# MM-estimation


# MS-estimation

