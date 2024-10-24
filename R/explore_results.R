# Use shinystan to explore model runs
library(shinystan)

# Load example sim data
sim_data <- read.csv("Data/S3_weeklydata.csv")

# GRW run
grw_fit <- readRDS("Results/SimData/grw_fit.rds")
launch_shinystan(grw_fit)

# IBM Run
ibm_fit <- readRDS("Results/SimData/ibm_fit_clamped.rds")
launch_shinystan(ibm_fit)

# OU Run
ou_fit <- readRDS("Results/SimData/ou_fit.rds")
launch_shinystan(ou_fit)