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
df1 <- read.csv('../../output/sim_results/res_merged.csv')
df1 <- df1 %>% select(-X.1)
df2 <- read.csv('../../output/sim_results/sim_results_2022-11-30_23-04-07.csv')
df2 <- df2 %>% filter(method=='SIMEX dir.')
df <- rbind(df1,df2)
write.csv(df,file='../../output/sim_results/res_merged_1202.csv')
# opcharspath <- "../../output/sim_results/sim_results_2022-10-03_20-28-42.csv"
# df <- read.csv(opcharspath)


# Make various line plots
defaults <- c(0.5,1,5000,0,0.8,0.5,0.25)
# ------------------------------------------------------------------------------
# Vary u
line_plot(df,'u','bias') # using default values
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