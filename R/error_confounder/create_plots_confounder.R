# create_plots_tables.R
#
# This file allows users to create <...>
library(ggplot2)
library(tidyr)

rm(list=ls())
source('~/Github/ME-Causal-Review/R/error_confounder/plotting_functions.R')

# Merge most recent round with previous
# df1 <- read.csv('../../output/sim_results/res_merged.csv')
# df1 <- df1 %>% select(-X.1)
# df2 <- read.csv('../../output/sim_results/sim_results_2022-11-30_23-04-07.csv')
# df2 <- df2 %>% filter(method=='SIMEX dir.')
# df <- rbind(df1,df2)
# write.csv(df,file='../../output/sim_results/res_merged_1202.csv')
opcharspath <- "~/Github/ME-Causal-Review/output/sim_results/confounder_results.csv"
df <- read.csv(opcharspath)

# ------------------------------------------------------------------------------
df_long_con <- df %>% pivot_longer(cols = c(bias, rmse, ci_cov), names_to = "outcome", 
                        values_to = "value") %>%
          mutate(outcome=replace(outcome, outcome=='bias','% Bias'),
                 outcome=replace(outcome, outcome=='ci_cov','C.I. Coverage'),
                 outcome=replace(outcome, outcome=='rmse','RMSE'),
                 aw=as.character(aw),
                 rho=as.character(rho)) %>%
  mutate(aw=replace(aw,aw=='0.25','Low confounding'),
         aw=replace(aw,aw=='0.75','High confounding'),
         rho=replace(rho,rho=='0.15','corr(W,X) = 0.15'),
         rho=replace(rho,rho=='0.85','corr(W,X) = 0.85'),
         mis=replace(mis,mis=="ps-mis", "Misspecified Propensity Score"),
         mis=replace(mis,mis=="out-mis", "Misspecified Outcome Model"),
         mis=replace(mis,mis=="base", "Both Models Correct"))

df_long_con <- subset(df_long_con, !(method %in% c("SIMEX", "Naive") & outcome %in% c("C.I. Coverage", "RMSE")))
 
grid_plot_bin <- df_long_con %>% filter(ba == 1 & n == 1000 & mis == "Both Models Correct") %>%
  ggplot(aes(x=sig_u,y=value,color=method))  + geom_point() +
  geom_line() + 
  facet_grid(outcome ~ as.factor(rho), scales='free') +
  theme_bw() + labs(x='Measurement error variance',
                    y='',
                    color='Method',
                    title='Simulation Results: Confounder Error') +
  theme(legend.position='bottom') ; grid_plot_bin

