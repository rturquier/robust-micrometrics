####=====================  Setup  ======================####
library(haven)
library(lfe)
library(L1pack)
library(MASS)
library(systemfit)

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"
results_path <- "original-project"

decisions_df <- read_dta(
  file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta")
)

####=====================  Main code  ======================####
#### --------- Reproduction
# Column (1) and (2) of table 2
# Low govt provision:  $4 -> $10 (income $46 -> $40)
Low_df <- subset(decisions_df, Budget == 2 | Budget == 5)
low_reg <- felm(h ~ Ggov, data=Low_df)
summary(low_reg, robust=T)
texreg(low_reg, file = file.path(results_path,"Low_rob_reg_results.tex"))

# High govt provision:  $28 -> $34 (income $46 -> $40)
High_df <- subset(decisions_df, Budget == 4 | Budget == 6)
high_reg <- felm(h ~ Ggov, data=High_df)
summary(high_reg, robust=T)
texreg(high_reg, file = file.path(results_path,"High_rob_reg_results.tex"))


#### --------- Robust regressions
# Least Absolute Deviation (LAD) Regression
# BR = Barrodale & Roberts method (default)
LAD_low_reg <- lad(h ~ Ggov, data = Low_df, method = "BR") 
summary(LAD_low_reg)

LAD_high_reg <- lad(h ~ Ggov, data = High_df, method = "BR")
summary(LAD_high_reg)

# M-estimation -> kappa standard 1.345
M_low1345_reg <- rlm(Low_df$Ggov,Low_df$h,
                 init = "ls", psi = psi.huber,
                 scale.est = "Huber", k2 = 1.345,
                 method = "M",
                 maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_low1345_reg)
texreg(M_low1345_reg, file = file.path(results_path,
                                       "M_low1345_reg_results.tex"))

M_high1345_reg <- rlm(High_df$Ggov,High_df$h,
                 init = "ls", psi = psi.huber,
                 scale.est = "Huber", k2 = 1.345,
                 method = "M",
                 maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)
summary(M_high1345_reg)
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
