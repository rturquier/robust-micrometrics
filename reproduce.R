####=====================  1. Setup  ======================####
library(haven)
library(tidyverse)
library(DescTools)
library(robustbase)

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"

decisions_df <- read_dta(
  file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta")
)


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
