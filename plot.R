library(haven)
library(tidyverse)
library(DescTools)

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
  mutate(crowd_out_rate = crowd_out / 6) %>%
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
                         y = crowd_out_rate)) +
  geom_dotplot(aes(color = GovtHigh, fill = GovtHigh),
               binwidth = 0.05,
               binaxis = "y",
               stackdir = "center",
               dotsize = 2,
               stackratio = 1.7) +
  geom_boxplot(color="#00000055",
               fill = "#00000000",
               outlier.shape = NA,
               width = 0.37) +
  stat_boxplot(geom = "errorbar", alpha=0.25, width=0.23) +
  labs(x = "External contribution", y = "") +
  ggtitle("Crowd-out rate") + 
  scale_x_discrete(labels = c("Low\n(4$ → 10$)", "High\n(28$ → 34$)")) +
  scale_y_continuous(breaks = -3:4, 
                     labels = c("-300%","", "", "0%", "100%", "", "", "400%"),
                     minor_breaks = c(0.5)) +
  theme_minimal() +
  theme(legend.position="none", 
        plot.title = element_text(size = 11, hjust = -0.1, vjust = -5))
  
# Estimations of location
location_statistics <- crowd_out_df %>%
  group_by(GovtHigh) %>%
  summarise(mean = mean(crowd_out_rate),
            trimmed_mean_05 = mean(crowd_out_rate, trim = 0.05),
            median = median(crowd_out_rate),
            hodges_lehman = DescTools::HodgesLehmann(crowd_out_rate))

# Number and proportion of people of each type
count_types <- crowd_out_df %>%
  group_by(GovtHigh) %>%
  mutate(pure_altruism   = crowd_out_rate == 1,
         pure_warm_glow  = crowd_out_rate == 0,
         impure_altruism = 0 < crowd_out_rate & crowd_out_rate < 1,
         over_crowd_out  = crowd_out_rate > 1,
         crowd_in        = crowd_out_rate < 0,
         misbehaving     = over_crowd_out | crowd_in) %>%
  summarise(across(pure_altruism:misbehaving,
                   list(count = ~ sum(., na.rm = TRUE), proportion = mean)))
