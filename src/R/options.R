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

  # Population data link (RespiCompass per-country age-band estimates).
  # NB: RespiCompass archived the 2025/2026 round-1 auxiliary data under
  # Previous_Rounds/ when the round closed (2026-06-12), so this points there.
  # Update the round path if adapting to a newer RespiCompass round.
  o$pop_url = "https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/refs/heads/main/Previous_Rounds/2025-2026_round_1/auxiliary-data/population/population_estimates.csv"
  
  # Data on births
  o$births_url = "data/population/monthly_births.csv"
  
  o$vaccine_url = NULL
  
  # RespiCompass data links
  o$respicompass =
    list(hospital_admissions  = "data/epidemiological/RSV_weekly_counts.csv",
         hospital_burden_agegroups = "data/epidemiological/RSV_monthly_prop_age.csv")

  o$contact_matrices = paste0(o$pth$data_contact, "/contact_all.rdata")

  # ---- General data ----
  births_df = read.csv(o$births_url, fileEncoding = "UTF-8-BOM") %>%
    mutate(date = make_date(
      year  = as.integer(TIME_PERIOD),
      month = match(month, month.name),
      day   = 1),
      country = geo,
      births = OBS_VALUE
    ) %>%
    select(country, date, births) %>%
    filter(!is.na(date), date < ymd("2024-01-01")) %>%
    setDT() 
  df_2024 = births_df %>% filter(year(date) == "2023") %>% mutate(date = date + years(1))
  df_2025 = births_df %>% filter(year(date) == "2023") %>% mutate(date = date + years(2))
  df_2026 = births_df %>% filter(year(date) == "2023") %>% mutate(date = date + years(3))
  
  o$births = bind_rows(births_df, df_2024, df_2025, df_2026) %>% arrange(country, date)
  o$burden = read.csv(o$respicompass$hospital_burden_agegroups, fileEncoding = "UTF-8-BOM")
  
  
  # ---- Calibration settings ----
  
  # NB: the observation-model over-dispersion used by the likelihood is `k` in
  # config/default.yaml, NOT an option here. It behaves like any other model
  # parameter: list it under calibration_parameters to FIT it, otherwise the
  # fixed yaml value is used. It was previously duplicated as o$k, which
  # silently shadowed the yaml value - see the note on `k` in default.yaml.

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
  # Start date for plotting
  o$plot_start_date = "2023-09-01"
  
  # Zoom in start date for plotting
  o$plot_zoom_date = "2026-07-01"
  
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
