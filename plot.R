library(haven)
library(tidyverse)

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"

decisions_df <- read_dta(
  file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta")
)