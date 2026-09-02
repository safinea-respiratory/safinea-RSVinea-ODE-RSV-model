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
for (country in c("AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR")){ # "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LV", "MT", "NL", "NO", "PT", "RO", "SE", "SI", "SK"
  o = set_options(do_step = c(1:3), analysis_name = country)
  
  message("Running RSVinea v2.0 (", o$analysis_name, ")")
  
  # Step 1) Calibrate model ----
  # print(Sys.time())
  # run_calibration(o) # See R/calibration.R
  # print(Sys.time())
  
  # Step 2) Run all scenarios ----
  print(Sys.time())
  run_scenarios(o)   # See R/scenarios.R
  print(Sys.time())
  
  # Step 3) Plot results ----
  run_results(o)     # See R/results.R
  print(Sys.time())
  
}
# Finish up
message("* Finished!")

