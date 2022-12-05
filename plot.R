library(haven)
library(tidyverse)

datasets_path  <- "original-project/voxA/Work/Datasets/Derived/"

decisions_df <- read_dta(
  file.path(datasets_path,"voxA-cr021A_v01a-DecisionsLONG.dta")
)

# Prepare data and compute crowd-out
crowd_out_df <- decisions_df %>%
  filter(z %in% c(50, 74)) %>%  
  group_by(newid, GovtHigh) %>%
  arrange(newid, GovtHigh, y) %>%
  summarise(crowd_out = h - lag(h)) %>%
  drop_na() %>%
  mutate(relative_crowd_out = crowd_out / 6) %>%
  mutate(GovtHigh = as.factor(GovtHigh))

# Connected pairs plot
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

# Dot plot with boxplot
ggplot(crowd_out_df, aes(x = GovtHigh,
                         y = relative_crowd_out)) +
  geom_dotplot(aes(color = GovtHigh, fill = GovtHigh),
               binwidth = 0.05,
               binaxis = "y",
               stackdir = "center",
               dotsize = 1.5,
               stackratio = 1.7) +
  geom_boxplot(color="#00000060", fill = "#00000000", outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", alpha=0.3, width=0.25) +
  theme_minimal()
  
