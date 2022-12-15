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
library(systemfit)
library(lfe)
library(base)
library(lme4)
library(officer)
library(flextable)
library(huxtable)
library(lmtest)
library(sandwich)
library(fixest)
#install.packages("lmtest")

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"
results_path <- "original-project"

decisions_df <- read_dta(file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta"))


####=====================  2. Functions  ======================####
export_unconventional_regression <- function(regression_result, output_path){
  regression_result %>%
    summary() %>%
    capture.output() %>%
    paste0(collapse = "\n") %>%
    cat(file = output_path)
}

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

####=====================  5. Linear Fixed Effects Regression  ======================####
# Column (1) and (2) of table 2

# Linear Fixed Effects Robust s.e. corresponds to assuming different variances at Low and High.

# Low govt provision:  $4 -> $10 (income $46 -> $40)
Low_df <- subset(decisions_df, Budget == 2 | Budget == 5)
low_reg <- felm(h ~ Ggov, data=Low_df)
summary(low_reg, robust=T)
texreg(low_reg, file = file.path(results_path,"Low_rob_reg_results.tex"))


#low_reg1 <- feols(h ~ Ggov, se="standard", data=Low_df)
#summary(low_reg1)
#coeftest(low_reg1, vcov = vcovHC(low_reg1, type = "HC1"))
# Ho: crowd-out is -1 or a bigger negative number.



# High govt provision:  $28 -> $34 (income $46 -> $40)
High_df <- subset(decisions_df, Budget == 4 | Budget == 6)
high_reg <- felm(h ~ Ggov, data=High_df)
# Ho: crowd-out is -1 or a bigger negative number.
summary(high_reg, robust=T)
texreg(high_reg, file = file.path(results_path,"High_rob_reg_results.tex"))


####=====================  6. Robust regression  ======================####


# LAD-estimation -> Least Absolute Deviation Regression --------------------------------------------
LAD_low_reg <- lad(h~ Ggov, data = Low_df, method = "BR") # BR = Barrodale & Roberts method (default)
summary(LAD_low_reg)
LAD_low_reg %>% export_unconventional_regression("LAD_low_reg.txt")

LAD_high_reg <- lad(h~ Ggov, data = High_df, method = "BR") # BR = Barrodale & Roberts method (default)
summary(LAD_high_reg)
print(LAD_high_reg, file = file.path(results_path,"LAD_high_reg_results.tex")) # EXPORT DOESN'T WORK



# M-estimation -> kappa standard 1.345 -----------------------------------------
M_low1345_reg <- rlm(Low_df$Ggov,Low_df$h,
                 init = "ls", psi = psi.huber,
                 scale.est = "Huber", k2 = 1.345,
                 method = "M",
                 maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_low1345_reg)
texreg(M_low1345_reg, file = file.path(results_path,"M_low1345_reg_results.tex"))

M_high1345_reg <- rlm(High_df$Ggov,High_df$h,
                 init = "ls", psi = psi.huber,
                 scale.est = "Huber", k2 = 1.345,
                 method = "M",
                 maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_high1345_reg)
texreg(M_high1345_reg, file = file.path(results_path,"M_high1345_reg_results.tex"))



# M-estimation -> kappa low 0.345 ------------------------------------------------
M_low0345_reg <- rlm(Low_df$Ggov,Low_df$h,
                 init = "ls", psi = psi.huber,
                 scale.est = "Huber", k2 = 0.345,
                 method = "M",
                 maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_low0345_reg)
texreg(M_low0345_reg, file = file.path(results_path,"M_low0345_reg_results.tex"))

M_high0345_reg <- rlm(High_df$Ggov,High_df$h,
                  init = "ls", psi = psi.huber,
                  scale.est = "Huber", k2 = 0.345,
                  method = "M",
                  maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_high0345_reg)
texreg(M_high0345_reg, file = file.path(results_path,"M_high0345_reg_results.tex"))



# M-estimation -> kappa high 2.345 -----------------------------------------
M_low2345_reg <- rlm(Low_df$Ggov,Low_df$h,
                 init = "ls", psi = psi.huber,
                 scale.est = "Huber", k2 = 2.345,
                 method = "M",
                 maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_low2345_reg)
texreg(M_low2345_reg, file = file.path(results_path,"M_low2345_reg_results.tex"))

M_high2345_reg <- rlm(High_df$Ggov,High_df$h,
                  init = "ls", psi = psi.huber,
                  scale.est = "Huber", k2 = 2.345,
                  method = "M",
                  maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_high2345_reg)
texreg(M_high2345_reg, file = file.path(results_path,"M_high2345_reg_results.tex"))





# Hausman test: LAD vs. OLS
hausman.systemfit(LAD_high_reg,high_reg)
hausman.systemfit(LAD_low_reg,low_reg)

# Hausman test: LAD vs. M
hausman.systemfit(LAD_high_reg,M_high_reg)
hausman.systemfit(LAD_low_reg,M_low_reg)

# Hausman test: M vs. OLS
hausman.systemfit(M_high_reg,high_reg)
hausman.systemfit(M_low_reg,low_reg)







# LTS-estimation -> Least Trimmed Sum of Squares -------------------------------
LTS_low_reg <- ltsReg(h ~ Ggov, data = Low_df)
summary(LTS_low_reg)
print(LTS_low_reg, file = file.path(results_path,"LTS_low_reg_results.tex")) # EXPORT DOESN'T WORK

LTS_high_reg <- ltsReg(h ~ Ggov, data = High_df)
summary(LTS_high_reg)
print(LTS_high_reg, file = file.path(results_path,"LTS_high_reg_results.tex")) # EXPORT DOESN'T WORK



# LMS-estimation -> Least Median Squares Regression ----------------------------
LMS_low_reg <- lmsreg(h ~ Ggov, data = Low_df)
summary(LMS_low_reg)
print(LMS_low_reg, file = file.path(results_path,"LMS_low_reg_results.tex")) # EXPORT DOESN'T WORK

LMS_high_reg <- lmsreg(h ~ Ggov, data = High_df)
summary(LMS_high_reg)
print(LMS_high_reg, file = file.path(results_path,"LMS_high_reg_results.tex")) # EXPORT DOESN'T WORK



# S-estimation -----------------------------------------------------------------
S_low_reg <- lmrob.S(x=Low_df$Ggov, y=Low_df$h, control = lmrob.control(nRes = 20), trace.lev=1)
summary(S_low_reg)


S_high_reg <- lmrob.S(x=High_df$Ggov, y=High_df$h, control = lmrob.control(nRes = 20), trace.lev=1)



# MM-estimation ----------------------------------------------------------------
MM_low_reg <- lmrob(h ~ Ggov, data = Low_df, method = "MM",
      model = TRUE, y = FALSE,
      singular.ok = TRUE, contrasts = NULL, offset = NULL, init = NULL)

MM_high_reg <- lmrob(h ~ Ggov, data = High_df, method = "MM",
                model = TRUE, y = FALSE,
                singular.ok = TRUE, contrasts = NULL, offset = NULL, init = NULL)
# MS-estimation
# Gets implemented in case there are dummies in rob



