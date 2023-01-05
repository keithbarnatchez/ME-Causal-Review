# exposure_sim_main.R
#
# Code for implementing exposure ME simulation exercise

rm(list = ls())

## Preliminaries
library(SuperLearner)
library(parallel)
library(abind)
library(geex)

# Code for generating and fitting data
# Assuming directory is set to this code file's location
source("erf.R")
source("gen_data.R")
# source("pred.R")
source("rc.R")
source("csme.R")
source("correction_functions.R")
source("output_functions.R")

sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")

# 
data <- gen_data(sigU = 0.5)
A <- data$A ; A.star <- data$Astar ; Y <- data$Y
X1 <- data$X1 ; X2 <- data$X2 ; v_idx <- data$v_idx

# reg cal
rcal_res <- rcal(data$Y,data$A,data$X1,data$X2,data$Astar,
                 data$v_idx,sl.lib)

# simex
X <- cbind(data$X1,data$X2)
sig2_me <- var(data$Astar[v_idx==1] - data$A[v_idx==1])
simex_res <- simex_direct(A.star,Y,X,1,0,
                          family=gaussian(),
                          n.boot=100,
                          degree=2,mc.cores=1,
                          tau2=sig2_me)
                          

# csme
sig2_me <- var(data$Astar[v_idx==1] - data$A[v_idx==1])
csme_res <- csme_linear(Y,A,X1,X2,A.star,
                        sig2_me=sig2_me)


################################################################################
methods <- c('rcal')
sig_u_grid <- c(0.3,0.5)
ba_grid <- c(1)
bax1_grid <- c(-0.25) ; bax2_grid <- c(0.5)
n_grid <- c(1000) ; nsim <- 100
thecheck <- get_results(methods,
                        sig_u_grid,
                        ba_grid,bax1_grid,bax2_grid,
                        n_grid,
                        nsim)

