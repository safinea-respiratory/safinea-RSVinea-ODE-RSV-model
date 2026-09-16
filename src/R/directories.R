########################################################## #
# SET DIRECTORIES
#
# Set and get directories in one place for consistency and ease.
# Creates any directories that do not currently exist.
#
# OUTPUTS:
#	- A list of relevant directories (within o$pth) which can be 
#   referenced elsewhere.
#
########################################################## #

# -------------------------------------------------------- -
# Define paths for project inputs and outputs ----
# -------------------------------------------------------- -
set_dirs = function(o) {
  
  # Initiate file path lists
  pth = out = list()
  
  # We've already moved to code directory
  pth$code = getwd()

  # ---- Input and configuration files ----
  
  # Parent path of all input files
  pth$input  = file.path(pth$code, "input")
  pth$config = file.path(pth$code, "config")
  pth$data   = file.path(pth$code, "data")

  # ---- Set analysis name ----
  
  # Throw error if illegal characters used
  if (grepl("\\.", o$analysis_name))
    stop("Analysis name should not contain any period characters")
  
  # ---- Model parameter files ----
  
  # Default model parameters
  pth$params_default = file.path(pth$config, "default.yaml")
  
  # User-specified model parameters for this analysis
  pth$params_user = file.path(pth$input, paste0(o$analysis_name, ".yaml"))
  
  # ---- Data files ----
  # Population data provided by RespiCompass
  pth$data_pop = file.path(pth$data, "population")
  
  # Epi data provided by RespiCompass
  pth$data_epi = file.path(pth$data, "epidemiological")

  # Vaccination coverage data (user-specified)
  pth$data_vaccination = file.path(pth$data, "vaccination")
  
  # Synthetic contact matrices (Prem et al. 2021)
  pth$data_contact = file.path(pth$data, "contact_matrices")
  
  # ---- Output directories ----
  
  # Parent path of all output files, NAMESPACED BY GIT BRANCH.
  #
  # src/output/ is gitignored, and git deliberately leaves ignored files alone
  # when switching branches. Without this namespacing every branch writes the
  # SAME paths for the same analysis_name (output/1_calibration/<name>/...) and
  # silently overwrites the other branch's results. That is especially easy to
  # miss because nothing here ever deletes: make_out_dirs() only creates.
  #
  # Set o$output_tag in options.R to override (e.g. to share results between
  # branches deliberately, or when running outside a git checkout).
  tag = o$output_tag
  if (is.null(tag)) {
    tag = tryCatch(
      system2("git", c("-C", shQuote(pth$code), "rev-parse", "--abbrev-ref", "HEAD"),
              stdout = TRUE, stderr = FALSE),
      error   = function(e) NULL,
      warning = function(w) NULL)
  }

  # Fall back to a fixed name if git is unavailable or HEAD is detached
  if (length(tag) != 1 || is.na(tag) || !nzchar(tag)) tag = "default"

  # Branch names may contain characters that are illegal in paths (e.g. "/")
  tag = gsub("[^A-Za-z0-9._-]", "-", tag)

  pth_output = file.path(pth$code, "output", tag)

  # Expose it so other scripts (e.g. results_evaluation.R) can resolve paths
  pth$output = pth_output
  
  # Path to test run files
  out$testing = file.path(pth_output, "0_testing")
  
  # Path to calibration files
  out$fitting     = file.path(pth_output, "1_calibration", o$analysis_name)
  out$fit_samples = file.path(out$fitting, "fit_samples")
  
  # Paths to scenario files
  pth_scenarios   = file.path(pth_output, "2_scenarios", o$analysis_name)
  out$scenarios   = file.path(pth_scenarios, "scenarios")
  out$simulations = file.path(pth_scenarios, "simulations")
  out$uncertainty = file.path(pth_scenarios, "uncertainty")
  
  # Path to figures and other output results
  out$results = file.path(pth_output, "3_results", o$analysis_name)
  
  # ---- Create directory structure ----
  
  # Make all output directories
  make_out_dirs(out)
  
  # Append paths to o list
  o = append_dirs(o, pth, out)
  
  return(o)
}

# -------------------------------------------------------- -
# Make all output directories if they do not already exist  ----
# -------------------------------------------------------- -
make_out_dirs = function(out) {
  
  # Extract all path names in list
  pth_names = names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]
    
    # If it does not already exist, create it
    if (!dir.exists(this_pth) & !grepl("\\*", this_pth))
      dir.create(this_pth, recursive = TRUE)
  }
}

# -------------------------------------------------------- -
# Concatenate separators and append directories to o list  ----
# -------------------------------------------------------- -
append_dirs = function(o, pth, out) {
  
  # Extract all path names in list
  pth_names = names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]
    
    # We use * to denote a partial path
    if (grepl("\\*", this_pth)) {
      out[[pth_name]] = substr(this_pth, 1, nchar(this_pth) - 1)
      
    } else {  # Otherwise add a file separator to end of output paths
      out[[pth_name]] = paste0(this_pth, .Platform$file.sep)
    }
  }
  
  # Concatenate lists
  o$pth = c(pth, out)
  
  return(o)
}

