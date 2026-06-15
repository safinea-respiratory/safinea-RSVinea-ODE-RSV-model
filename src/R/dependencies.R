########################################################## #
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
########################################################## #

# ---- R version check ----

# R versions for which this project has been tested and is stable
stable_versions = "4.3.0"

# R versions for which this project is stable (as a string)
stable_str = paste(stable_versions, collapse = ", ")

# Get details of R version currently running
version_info = R.Version()

# Construct version number from list details
version_num = paste0(version_info$major, ".", version_info$minor)

# Warn if the current R version has not been explicitly tested
if (!version_num %in% stable_versions)
  warning("RSVinea has been tested on R ", stable_str,
          " but you are running R ", version_num,
          ". Proceed with caution.")

# ---- Source files ----
# Load R files and functions
source("R/auxiliary.R")
source("R/options.R")
source("R/directories.R")
source("R/parse_input.R")
source("R/model.R")
source("R/calibration.R")
source("R/load_data.R")
source("R/scenarios.R")
source("R/uncertainty.R")
source("R/results.R")
source("R/plotting.R")
source("R/unit_tests.R")

# ---- Define packages ----

# Complete list of all R packages required for this project
packages = c("tidyverse",      # Includes ggplot2, dplyr, tidyr (www.tidyverse.org/packages/)
             "data.table",     # Next generation dataframes
             "useful",         # General helper functions (eg compare.list)
             "rlist",          # List-related helper functions (eg list.remove)
             "tgp",            # Latin hypercube sampler
             "deSolve",        # Numerical solutions for ODE model
             "tictoc",         # Code timer
             "yaml",           # Data loading functionality
             "lubridate",      # Data formatting functionality
             "gsubfn",         # Named list functionality
             "pals",           # Colour palettes
             "viridis",        # Colour palettes
             "RColorBrewer",   # Colour palettes
             "parallel",       # Parallelisation
             "foreach",        # To combine model output efficiently
             "pbapply",        # Progress bar
             "nanoparquet")    # Reading/writing Parquet files (used in results_evaluation)

# ---- Install and/or load R packages with pacman ----

message("* Installing required R packages")

# Check whether pacman itself has been installed
pacman_installed = "pacman" %in% rownames(installed.packages())

# If not, install it
if (!pacman_installed) 
  install.packages("pacman")

# Load pacman
library(pacman) 

# Load all required packages, installing them if required
pacman::p_load(char = packages)

# ---- Redefine or unmask particular functions ----
# Unmask certain functions to ensure they're not overwritten
select  = dplyr::select
filter  = dplyr::filter
rename  = dplyr::rename
recode  = dplyr::recode
count   = dplyr::count
union   = dplyr::union


