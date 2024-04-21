# create_plots_tables.R
#
# This file allows users to create <...>
#
#
rm(list=ls())
source('~/Github/ME-Causal-Review/R/error_confounder/plotting_functions.R')
# merge_op_chars <- TRUE
# if (merge_op_chars) {
#   files <- list.files(path="path/to/dir", pattern="*.txt", full.names=TRUE, recursive=FALSE)
# }

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
         rho=replace(rho,rho=='0.85','corr(W,X) = 0.85'))


grid_plot_bin <- df_long_con %>% filter(bw == 0.25 & n == 1000) %>%
  ggplot(aes(x=sig_u,y=value,color=method))  + geom_point() +
  geom_line() + 
  facet_grid(outcome ~ as.factor(aw) + as.factor(rho), scales='free') +
  theme_bw() + labs(x='Measurement error variance',
                    y='',
                    color='Method',
                    title='Simulation results: binary outcome') +
  theme(legend.position='bottom') ; grid_plot_bin

ggsave('../../output/figures/grid_plot_binY.pdf',
       grid_plot_bin,
       width=7,height = 5,units='in')

