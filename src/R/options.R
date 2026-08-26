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
  
  # Number of processes for parallelisation.
  # Leave one core free so the machine stays responsive during long sweeps.
  o$parallel = max(1, detectCores(all.tests = FALSE, logical = TRUE) - 1)
  
  # Set analysis name and create output directory system
  o = set_dirs(o)  # See directories.R

  # Age group label mapping from RespiCompass population data to RSVinea reporting bands
  o$respicompass_age_map = yaml::read_yaml(o$pth$params_default)$respicompass_age_map

  # ---- Data references ----
  #
  # LOCAL CACHE. model() re-reads the population file on EVERY call, and model()
  # is called once per parameter sample, in parallel across workers - thousands
  # of times per country sweep. Fetching from raw.githubusercontent.com that
  # often gets the run rate-limited and then everything aborts with
  # "cannot open the connection". So each remote file is downloaded ONCE into a
  # local cache and the o$*_url paths point at the local copies; all the
  # downstream read.csv() calls are unchanged and simply read from disk.
  #
  # Set o$refresh_cache = TRUE (or delete the cache directory) to re-download,
  # e.g. when RespiCompass updates a round's data.
  respicompass_raw = "https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/refs/heads/main/"
  cache_dir = file.path("data", "respicompass_cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  refresh = isTRUE(o$refresh_cache)

  cache_file = function(url) {
    dest = file.path(cache_dir, basename(url))
    if (refresh || !file.exists(dest)) {
      ok = tryCatch({
        utils::download.file(url, dest, quiet = TRUE, mode = "wb")
        file.exists(dest) && file.size(dest) > 0
      }, error = function(e) FALSE, warning = function(w) file.exists(dest) && file.size(dest) > 0)
      if (!ok) {
        if (file.exists(dest)) unlink(dest)
        stop("Could not download ", url, "\n  and no usable cached copy exists at ", dest,
             "\n  (check your connection; RespiCompass may also be rate-limiting)")
      }
    }
    return(dest)
  }

  # RespiCompass country list: maps full country names to ISO-2 codes.
  # Used throughout the model to look up country names from ISO codes
  # without hard-coding country names.
  o$countries_url = cache_file(paste0(respicompass_raw, "supporting-files/countries.csv"))
  o$countries_df  = read.csv(o$countries_url)

  # ---- RespiCompass 2026/2027 RSV round-1 data ----
  # This round targets novel RSV immunisation strategies for older adults
  # NB: population/births/mortality use ISO-2 country codes; the target hospital
  # files use full country names. See load_data.R / model.R for the keying.

  # Population by age band and country (ISO-2 country column)
  o$pop_url = cache_file(paste0(respicompass_raw, "auxiliary-data/population/population_estimates.csv"))

  # Monthly live births by country (ISO-2), covering the modelling period
  o$births_url = cache_file(paste0(respicompass_raw, "auxiliary-data/births/births_by_month.csv"))

  # All-cause mortality: annual DEATH COUNTS by age band and country (ISO-2).
  # Converted to a per-capita rate at model setup (deaths / population); see
  # compute_background_mortality() in auxiliary.R.
  o$mortality_url = cache_file(paste0(respicompass_raw, "auxiliary-data/mortality/mortality_agegroups.csv"))

  # Adult vaccine waning curves: VE against infection (VE_inf) and against
  # severe disease (VE_sev) by MONTHS SINCE VACCINATION, supplied as 500
  # replicate curves ('rep') that carry the uncertainty in the waning model.
  # Applies to the ADULT product only - infant VE still comes from the yaml
  # (infant_vaccine_rel_protection / infant_vacc_IE).
  o$waning_url = cache_file(paste0(respicompass_raw, "auxiliary-data/waning-immunity/waning_curves.csv"))

  o$vaccine_url = NULL

  # RespiCompass target (observed) data links (full country-name column)
  o$respicompass =
    list(hospital_admissions       = cache_file(paste0(respicompass_raw, "target-data/hospitaladmissions.csv")),
         hospital_burden_agegroups = cache_file(paste0(respicompass_raw, "target-data/hospitalburden_agegroups.csv")))

  o$contact_matrices = paste0(o$pth$data_contact, "/contact_all.rdata")

  # ---- General data ----
  # Monthly births (real counts). RespiCompass supplies only the reference
  # season (2026-09 to 2027-08), but the modelling period spans two seasons, so
  # repeat the same monthly births forward one year to populate newborns in the
  # 2027/28 season (assumes births are stable year-to-year). Country is ISO-2.
  # NB: normalise_iso2() maps Eurostat's 'EL' to ISO-2 'GR' - without it Greece
  # silently matches zero rows here (see auxiliary.R).
  births_ref = read.csv(o$births_url, fileEncoding = "UTF-8-BOM") %>%
    normalise_iso2() %>%
    mutate(date = ymd(date)) %>%
    select(country, date, births) %>%
    filter(!is.na(date)) %>%
    setDT()
  o$births = bind_rows(births_ref,
                       births_ref %>% mutate(date = date %m+% years(1))) %>%
    arrange(country, date) %>%
    setDT()

  # Raw all-cause death counts by age band (used to derive mortality rates).
  # This file already uses 'GR'/'Czechia'; normalised defensively for consistency.
  o$mortality = read.csv(o$mortality_url, fileEncoding = "UTF-8-BOM") %>%
    normalise_iso2(col = "iso2_code") %>%
    normalise_country_name()

  # Adult vaccine waning curves, one row per (rep, month). Read once here and
  # resolved to a single curve per simulation by get_waning_curve() in
  # auxiliary.R - see there for how a replicate is chosen.
  o$waning = read.csv(o$waning_url, fileEncoding = "UTF-8-BOM") %>% setDT()

  # Age-stratified RSV hospital burden (seasonal totals per age group), kept
  # separately here as the source for age_relativity()'s p_hosp-by-age ratios.
  # Target files name Czechia "Czech Republic"; normalise to match countries.csv.
  o$burden = read.csv(o$respicompass$hospital_burden_agegroups, fileEncoding = "UTF-8-BOM") %>%
    normalise_country_name()
  
  
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

  # NB: the observation-model over-dispersion used by the likelihood is `k` in
  # config/default.yaml, NOT an option here. It behaves like any other model
  # parameter: list it under calibration_parameters to FIT it, otherwise the
  # fixed yaml value is used. It was previously duplicated as o$k, which
  # silently shadowed the yaml value - see the note on `k` in default.yaml.

  # ---- Adaptive perturbation kernel (see sample_parameters in calibration.R) ----
  # Between rounds, resampled parameter sets are jittered to create new candidates.
  # The step size for each parameter is derived from the SPREAD of the resampled
  # particles for that parameter, rather than being a fixed percentage: the spread
  # is the current estimate of how uncertain that parameter is, so it is the right
  # scale to explore at. It also anneals automatically - wide while the particles
  # are scattered, narrow once they concentrate.
  #
  # kernel_scale : multiplier on the measured (log-scale) particle spread.
  #   Values < 1 refine, > 1 explore more aggressively. ABC-SMC (Beaumont 2009)
  #   uses a kernel variance of 2x the particle variance, but that sits inside an
  #   importance-sampling scheme whose weights CORRECT for the proposal; this
  #   resample-and-perturb scheme has no such correction, so we start smaller.
  # kernel_floor : minimum log-scale sd, so the kernel can never collapse to zero
  #   and freeze the search if the particles concentrate prematurely.
  o$kernel_scale = 0.5
  o$kernel_floor = 0.01

  # Re-run fitting, overwrite if TRUE, otherwise will use previous fit
  o$overwrite_samples = TRUE

  # Plot fits only, without re-running fitting
  o$plot_only = FALSE

  # ---- Scenario settings ----
  # Number of sampled parameter sets from calibration to use in scenarios
  o$n_best_samples = 10
  
  # Number of uncertainty parameter sets to sample 
  o$n_parameter_sets = 1 # Best to set to 1 if not simulating parameter uncertainty

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
