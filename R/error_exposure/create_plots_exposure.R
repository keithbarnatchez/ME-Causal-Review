# create_plots_tables.R
#
# This file allows users to create <...>
#
#
library(ggplot2)
library(tidyr)

rm(list=ls())
source('~/Github/ME-Causal-Review/R/error_exposure/plotting_functions.R')

# Merge most recent round with previous
opcharspath <- "~/Github/ME-Causal-Review/output/sim_results/exposure_results.csv"
df <- read.csv(opcharspath)

df_long_con <- df %>% pivot_longer(cols = c(bias, rmse, ci_cov), names_to = "outcome", 
                                   values_to = "value") %>%
  mutate(outcome=replace(outcome, outcome=='bias','% Bias'),
         outcome=replace(outcome, outcome=='ci_cov','C.I. Coverage'),
         outcome=replace(outcome, outcome=='rmse','RMSE'),
         ba=as.character(ba)) %>%
  mutate(ba=replace(ba,ba=='-0.5','Small Treatment Effect'),
         ba=replace(ba,ba=='1','Large Treatment Effect'),
         mis=replace(mis,mis=="ps-mis", "Misspecified Propensity Score"),
         mis=replace(mis,mis=="out-mis", "Misspecified Outcome Model"),
         mis=replace(mis,mis=="base", "Both Models Correct"))

df_long_con <- subset(df_long_con, !(method %in% c("SIMEX", "Naive") & outcome %in% c("C.I. Coverage", "RMSE")))

grid_plot_bin <- df_long_con %>% filter(n == 1000 & ba == "Large Treatment Effect") %>%
  ggplot(aes(x=sig_u,y=value,color=method))  + geom_point() +
  geom_line() + 
  facet_grid(outcome ~ as.factor(mis), scales='free') +
  theme_bw() + labs(x='Measurement error variance',
                    y='',
                    color='Method',
                    title='Simulation Results: Exposure Error') +
  theme(legend.position='bottom') ; grid_plot_bin

