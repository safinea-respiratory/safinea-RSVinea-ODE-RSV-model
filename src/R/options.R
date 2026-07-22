########################################################## #
# OPTIONS
#
# Set key options for all things model related.
#
########################################################## #

# -------------------------------------------------------- -
# Set model options and assumptions ----
# -------------------------------------------------------- -
set_options = function(do_step = NA, quiet = FALSE, analysis_name = NA) {

  if (!quiet) message("* Setting options")

  # Reset some default R options
  options(stringsAsFactors = FALSE, scipen = 999, dplyr.summarise.inform = FALSE)

  # Initiate options list
  o = list(do_step = do_step)

  # Name of analysis to run (cannot contain period symbol)
  if (!is.na(analysis_name)) 
    o$analysis_name = analysis_name
  
  # Number of processes for parallelisation
  o$parallel = detectCores(all.tests = FALSE, logical = TRUE)
  
  # Set analysis name and create output directory system
  o = set_dirs(o)  # See directories.R

  # Age group label mapping from RespiCompass population data to RSVinea reporting bands
  o$respicompass_age_map = yaml::read_yaml(o$pth$params_default)$respicompass_age_map

  # ---- Data references ----
  # RespiCompass country list: maps full country names to ISO-2 codes.
  # Used throughout the model to look up country names from ISO codes
  # without hard-coding country names.
  o$countries_url = "https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/refs/heads/main/supporting-files/countries.csv"
  o$countries_df  = read.csv(o$countries_url)

  # ---- RespiCompass 2026/2027 RSV round-1 data ----
  # This round targets novel RSV immunisation strategies for older adults
  # NB: population/births/mortality use ISO-2 country codes; the target hospital
  # files use full country names. See load_data.R / model.R for the keying.
  respicompass_raw = "https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/refs/heads/main/"

  # Population by age band and country (ISO-2 country column)
  o$pop_url = paste0(respicompass_raw, "auxiliary-data/population/population_estimates.csv")

  # Monthly live births by country (ISO-2), covering the modelling period
  o$births_url = paste0(respicompass_raw, "auxiliary-data/births/births_by_month.csv")

  # All-cause mortality: annual DEATH COUNTS by age band and country (ISO-2).
  # Converted to a per-capita rate at model setup (deaths / population); see
  # compute_background_mortality() in auxiliary.R.
  o$mortality_url = paste0(respicompass_raw, "auxiliary-data/mortality/mortality_agegroups.csv")

  o$vaccine_url = NULL

  # RespiCompass target (observed) data links (full country-name column)
  o$respicompass =
    list(hospital_admissions       = paste0(respicompass_raw, "target-data/hospitaladmissions.csv"),
         hospital_burden_agegroups = paste0(respicompass_raw, "target-data/hospitalburden_agegroups.csv"))

  o$contact_matrices = paste0(o$pth$data_contact, "/contact_all.rdata")

  # ---- General data ----
  # Monthly births (real counts). RespiCompass supplies only the reference
  # season (2026-09 to 2027-08), but the modelling period spans two seasons, so
  # repeat the same monthly births forward one year to populate newborns in the
  # 2027/28 season (assumes births are stable year-to-year). Country is ISO-2.
  births_ref = read.csv(o$births_url, fileEncoding = "UTF-8-BOM") %>%
    mutate(date = ymd(date)) %>%
    select(country, date, births) %>%
    filter(!is.na(date)) %>%
    setDT()
  o$births = bind_rows(births_ref,
                       births_ref %>% mutate(date = date %m+% years(1))) %>%
    arrange(country, date) %>%
    setDT()

  # Raw all-cause death counts by age band (used to derive mortality rates).
  o$mortality = read.csv(o$mortality_url, fileEncoding = "UTF-8-BOM")

  # Age-stratified RSV hospital burden (seasonal totals per age group), kept
  # separately here as the source for age_relativity()'s p_hosp-by-age ratios.
  o$burden = read.csv(o$respicompass$hospital_burden_agegroups, fileEncoding = "UTF-8-BOM")
  
  
  # ---- Numerics ----
  # ODE solver method passed to deSolve::ode(). The dual-vaccination model has
  # ~3,000 states (3 tiers x W adult waning stages dominate), and for a STIFF
  # solver the Jacobian dominates runtime: a dense n x n Jacobian costs ~n extra
  # derivative calls to assemble plus an O(n^3) factorisation, repeated every few
  # steps — and that cost is per-STEP, so it does NOT shrink when you reduce
  # n_days. This is the single biggest runtime driver.
  #
  # This SEIRS system is only mildly stiff (fastest rates ~1/4-1/9 per day), so a
  # non-stiff / auto-switching solver usually wins by avoiding the Jacobian:
  #   "lsoda"  - DEFAULT. Auto-switches non-stiff <-> stiff, starting non-stiff.
  #              Best general choice for this model.
  #   "vode"   - the previous default; pure stiff BDF with a dense Jacobian.
  #              Was fine at ~775 states, is the bottleneck at ~3,000.
  #   "lsodes" - stiff but sparse; its sparsity auto-detection at t0 misses the
  #              V -> I couplings (V-stages are 0 until the first campaign), so it
  #              tends to thrash here. Avoid unless you supply the sparsity.
  #   "adams"  - explicit non-stiff; fastest IF the system is truly non-stiff.
  # Benchmark on your machine (see the timing snippet) and set the winner.
  o$ode_method = "adams"

  # ---- Calibration settings ----

  # Over-dispersion parameter for calculation of likelihood
  # (See Endo et al. 2020 Estimating the overdispersion in COVID-19
  # transmission using outbreak sizes outside China)
  o$k = 0.1

  # Re-run fitting, overwrite if TRUE, otherwise will use previous fit
  o$overwrite_samples = TRUE

  # Plot fits only, without re-running fitting
  o$plot_only = FALSE

  # ---- Scenario settings ----
  # Number of sampled parameter sets from calibration to use in scenarios
  o$n_best_samples = 10
  
  # Number of uncertainty parameter sets to sample 
  o$n_parameter_sets = 10 # Best to set to 1 if not simulating parameter uncertainty

  # Quantiles for credible intervals along the mean
  o$quantiles = c(0.05, 0.95)
  
  # ---- Plotting settings ----
  # Start date for plotting (RespiCompass 2026/2027 round modelling period)
  o$plot_start_date = "2026-09-01"

  # Zoom in start date for plotting (start of the 2027/28 season)
  o$plot_zoom_date = "2027-09-01"
  
  # Days of model output to include in fitting
  o$plot_from = 1
  o$plot_to = Inf
  
  # Saved figure size
  o$save_width  = 20
  o$save_height = 20
  
  # Units of figures sizes
  o$save_units = "cm"

  # Plotting resolution (in dpi)
  o$save_resolution = 300
  
  # Image format for saving multi-panel figures
  # NOTE: Use a character vector to save with multiple formats at once
  o$figure_format = "png" # Options: "png", "pdf", or "svg"
  # By default, the output is also saved as a multi-page PDF.
  
  # Turn scenario figures on or off
  o$plot_scenarios   = TRUE  # Plot scenarios

  return(o)
 }
