# create_plots_tables.R
#
# This file allows users to create <...>
#
#
rm(list=ls())
source('plotting_functions.R')
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
opcharspath <- "../../output/sim_results/sim_results_2023-06-22_17-24-12.csv"
df <- read.csv(opcharspath)

df <- df %>% mutate(mse=sqrt(mse))
# df <- df %>% mutate(mse = (ATE-))
# ------------------------------------------------------------------------------
df_long_bin <- df %>% filter(b==1) %>% pivot_longer(cols = c(bias, mse, ci_cov), names_to = "outcome", 
                        values_to = "value") %>%
          mutate(outcome=replace(outcome, outcome=='bias','% Bias'),
                 outcome=replace(outcome, outcome=='ci_cov','C.I. Coverage'),
                 outcome=replace(outcome, outcome=='mse','RMSE'),
                 ax=as.character(ax),
                 rho=as.character(rho)) %>%
  mutate(ax=replace(ax,ax=='0.25','Low confounding'),
         ax=replace(ax,ax=='0.75','High confounding'),
         rho=replace(rho,rho=='0.1','corr(X,Z)=0.1'),
         rho=replace(rho,rho=='0.8','corr(X,Z)=0.8'))

df_long_con <- df %>% filter(b==0) %>% pivot_longer(cols = c(bias, mse, ci_cov), names_to = "outcome", 
                                                    values_to = "value") %>%
  mutate(outcome=replace(outcome, outcome=='bias','% Bias'),
         outcome=replace(outcome, outcome=='ci_cov','C.I. Coverage'),
         outcome=replace(outcome, outcome=='mse','RMSE'),
         ax=as.character(ax),
         rho=as.character(rho)) %>%
  mutate(ax=replace(ax,ax=='0.25','Low confounding'),
         ax=replace(ax,ax=='0.75','High confounding'),
         rho=replace(rho,rho=='0.1','corr(X,Z)=0.1'),
         rho=replace(rho,rho=='0.8','corr(X,Z)=0.8'))


grid_plot_bin <- df_long_bin %>% ggplot(aes(x=u,y=value,color=method))  + geom_point() +
  geom_line() + 
  facet_grid(outcome ~ as.factor(ax) + as.factor(rho), scales='free') +
  theme_bw() + labs(x='Measurement error variance',
                    y='',
                    color='Method',
                    title='Simulation results: binary outcome') +
  theme(legend.position='bottom') ; grid_plot_bin

ggsave('../../output/figures/grid_plot_binY.pdf',
       grid_plot_bin,
       width=7,height = 5,units='in')

grid_plot_con <- df_long_con %>% ggplot(aes(x=u,y=value,color=method))  + geom_point() +
  geom_line() + 
  facet_grid(outcome ~ as.factor(ax) + as.factor(rho), scales='free') +
  theme_bw() + labs(x='Measurement error variance',
                    y='',
                    color='Method',
                    title='Simulation results: continuous outcome') +
  theme(legend.position='bottom') ; grid_plot_con

ggsave('../../output/figures/grid_plot_contY.pdf',
       grid_plot_con,
       width=7,height = 5,units='in')


  

# Make various line plots
# c('u','a','n','b','rho','psi','ax')
defaults <- c(0.25,1,4500,0,0.25,0.5,0.25)
# ------------------------------------------------------------------------------
# Vary u
line_plot(df,'u','bias') # using default values
line_plot(df,'u','mse')
line_plot(df,'u','bias',fixed_defaults=defaults) # using default values
line_plot(df,'u','ci_cov',fixed_defaults=defaults) # using default values
line_plot(df,'u','mse',fixed_defaults=defaults) # using default values

# ------------------------------------------------------------------------------
# Vary rho
line_plot(df,'rho','bias',fixed_defaults=defaults)
line_plot(df,'rho','bias',fixed_defaults=c(0.9,1,5000,0,0.8,0.5,0.25))
line_plot(df,'rho','ci_cov',fixed_defaults=c(0.9,1,5000,0,0.8,0.5,0.25))
line_plot(df,'rho','mse',fixed_defaults=defaults) # using default values


# ------------------------------------------------------------------------------
# Vary psi (this plot makes more sense to alter with just IV)



# ------------------------------------------------------------------------------