# create_plots_tables.R
#
# This file allows users to create <...>
#
#
library(ggplot2)

rm(list=ls())
source('exposure_plotting_functions.R')

# Merge most recent round with previous
opcharspath <- "../../output/sim_results/exposure/sim_results_2023-03-28_20-19-25.csv"
df <- read.csv(opcharspath)

df %>% ggplot(aes(x=u, y=bias, group=method, color=method)) + geom_point() +
  geom_line() + theme_bw() + 
  labs(title = 'Percent bias, varying M.E. variance',
       x='M.E. variance',
       y='Percent Bias') +
  theme(legend.position = 'bottom')

df %>% filter(method!='Naive') %>% ggplot(aes(x=u, y=bias, group=method, color=method)) + geom_point() +
  geom_line() + theme_bw() + 
  labs(title = 'Percent bias, varying M.E. variance',
       x='M.E. variance',
       y='Percent Bias') +
  theme(legend.position = 'bottom')

df %>% filter(method!='Naive') %>% ggplot(aes(x=u, y=sqrt(mse), group=method, color=method)) + geom_point() +
  geom_line() + theme_bw() + 
  labs(title = 'MSE, varying M.E. variance',
       x='M.E. variance',
       y='MSE') +
  theme(legend.position = 'bottom')