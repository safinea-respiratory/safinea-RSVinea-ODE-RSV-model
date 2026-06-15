########################################################## #
# LAUNCH
#
# Main launch function for RSVinea, a dynamic compartmental ODE
# model of RSV transmission and disease.
#
# Usage:
#  Via RStudio: Open RSVinea.Rproj and source this file (Ctrl+Shift+S)
#
########################################################## #

# Clear global environment
rm(list = ls())

# Load all required packages and functions
source("R/dependencies.R")

# Tidy up
if (interactive()) clf()  # Close figures
if (interactive()) clc()  # Clear console

# Set options
o = set_options(do_step = c(1:3), analysis_name = "IE")

message("Running RSVinea v2.0 (", o$analysis_name, ")")

# Step 1) Calibrate model ----
run_calibration(o) # See R/calibration.R

# Step 2) Run all scenarios ----
run_scenarios(o)   # See R/scenarios.R

# Step 3) Plot results ----
run_results(o)     # See R/results.R
      

# Finish up
message("* Finished!")

