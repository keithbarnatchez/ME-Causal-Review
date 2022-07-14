# create_plots_tables.R
#
# This file allows users to create <...>
#
#
source('plotting_functions.R')

opcharspath <- "../../output/sim_results/sim_results_2022-07-12_10-05-06.csv"
df <- read.csv(opcharspath)

# Make various line plots

# ------------------------------------------------------------------------------
# Vary u
line_plot(df,'u','bias')

# ------------------------------------------------------------------------------
# Vary rho
line_plot(df,'u','bias')


# ------------------------------------------------------------------------------
# Vary psi (this plot makes more sense to alter with just IV)



# ------------------------------------------------------------------------------