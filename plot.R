library(haven)
library(tidyverse)

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"

decisions_df <- read_dta(
  file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta")
)

decisions_wide_df <- read_dta(
  file.path(datasets_path,"voxA-cr001A_v01a-DecisionsWIDE.dta")
)

crowd_out_df <- decisions_df %>%
  filter(z %in% c(50, 74)) %>%  
  group_by(newid, GovtHigh) %>%
  arrange(newid, GovtHigh, y) %>%
  summarise(crowd_out = h - lag(h)) %>%
  drop_na() %>%
  mutate(relative_crowd_out = crowd_out / 6) %>%
  mutate(GovtHigh = as.factor(GovtHigh))

temp_df <- decisions_wide_df %>%
  select(newid, h04_y46, h10, h28_y46, h34)  %>%
  mutate(crowd_out_low_external  = h04_y46 - h10,
         crowd_out_high_external =  h28_y46 - h34,
         relative_crowd_out_low_external  =  crowd_out_low_external / 6,
         relative_crowd_out_high_external =  crowd_out_high_external / 6)


decisions_df %>%
  filter(Budget %in% c(2, 5)) %>%
  mutate(Budget = factor(Budget),
         Budget = fct_rev(Budget)) %>%
  ggplot(aes(x = Budget, y = h)) +
  geom_dotplot(binwidth = 0.5,
               binaxis = "y",
               stackdir = "center",
               dotsize = 0.6,
               stackratio = 1.7) +
  geom_line(aes(group = newid), alpha=.1) +
  theme_minimal() 


ggplot(crowd_out_df,aes(GovtHigh, relative_crowd_out, color=GovtHigh)) + 
  geom_boxplot(outlier.shape = NA, color="gray") +
  geom_jitter(alpha = 0.4, width = 0.1) + 
  theme_minimal()


ggplot(temp_df, aes(x = relative_crowd_out_low_external,
                    y = relative_crowd_out_high_external)) +
  geom_point(alpha=.2) +
  theme_light()
