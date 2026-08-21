##############################o############################ #
# LOAD DATA
#
# Load epi data required for calibration 
#
########################################################## #

# -------------------------------------------------------- -
# Parent function for extracting relevant fitting data ----
# -------------------------------------------------------- -
load_data = function(o, fit) {

  # Extract subset of inputs relevant for loading appropriate data
  opts = c(type = fit$input$calibration_type, fit$input$calibration_options)
  
  if(!is.null(opts) && opts$data_source$epi != "user-defined"){
    # Look up the full country name from the RespiCompass country list
    opts$country_name = o$countries_df %>%
      filter(iso2_code == opts$country) %>%
      pull(country)
  }
  
  # Identify data source regardless of capitalisation
  opts$data_source = lapply(opts$data_source, toupper)
  
  # Load epidemiological data
  fit = load_epi(o, opts, fit)

  return(fit)
}

# -------------------------------------------------------- -
# Load empirical epidemiological data from source ----
# -------------------------------------------------------- -
load_epi = function(o, opts, fit) {

  # Otherwise load empirical epidemiological data...
  message(" - Loading epi data: ", opts$data_source$epi)
  
  # Get population size by age group
  if(!is.null(opts$country_name)){
    pop_data = read.csv(o$pop_url, fileEncoding = "UTF-8-BOM") %>%
      normalise_iso2() %>%                   # Eurostat 'EL' -> ISO-2 'GR' (Greece)
      filter(country == opts$country) %>%    # RespiCompass population uses ISO-2 codes
      remap_age_groups(o$respicompass_age_map)
  } else {message("No population data needed for user-defined analysis")}
  
  # All dates we're interested in
  dates_df = get_data_dates(opts) %>%
    mutate(date = format_date(date))
  fit$dates_model = dates_df
  
  # ---- Source: RespiCompass ----
  # See https://github.com/european-modelling-hubs/RespiCompass for details
  # Pull data from RespiCompass
  
  if (opts$data_source$epi == "RESPICOMPASS") {

    # --- Weekly hospital admissions (total, all ages) --- #
    # RespiCompass target-data uses full country names and reports the ISO-week
    # Sunday directly in `target_end_date`, so no +6 shift or ISO-2 remap needed.
    # normalise_country_name(): target files say "Czech Republic", countries.csv
    # says "Czechia" - without this Czechia matches zero rows (no fitting target).
    raw_data = read.csv(o$respicompass$hospital_admissions, fileEncoding = "UTF-8-BOM") %>%
      normalise_country_name()

    data_hosp_admissions = raw_data %>%
      filter(country == opts$country_name) %>%
      transmute(date  = as.Date(target_end_date),
                value = weekly_rsv_hospitalisations) %>%
      mutate(date = format_date(date),
             metric = "hospital_admissions",
             age_group = "total",
             data_freq = "weekly") %>%
      # Keep only dates of interest
      right_join(y  = dates_df,
                 by = "date") %>%
      arrange(date) %>%
      filter(!is.na(value)) %>%
      # Summarise time period if desired...
      change_time_period(dates_df, opts$data_period) %>%
      setDT()

    # Throw warning if no hospital data for this country
    if (nrow(data_hosp_admissions) == 0) {
      warning("No RespiCompass hospital admissions found for country ", opts$country)
    }


    # --- Age-stratified hospital burden (seasonal totals per age group) --- #
    # This round provides a single seasonal total per age band (not a 4-weekly
    # proportion series). We tag it data_freq = "total" with date = NA so it
    # matches the model's per-age 'total' aggregation in aggregate_model_output().
    # Age bands are remapped to the model's reporting-band labels so they align
    # with the model-output grouping in fitting_format().
    raw_burden = read.csv(o$respicompass$hospital_burden_agegroups, fileEncoding = "UTF-8-BOM") %>%
      normalise_country_name()   # "Czech Republic" -> "Czechia"

    data_hosp_burden = raw_burden %>%
      filter(country == opts$country_name) %>%
      remap_age_groups(o$respicompass_age_map) %>%
      transmute(date      = as.Date(NA),
                age_group,
                value     = total_rsv_hospitalisations,
                metric    = "hospital_admissions",
                data_freq = "total") %>%
      setDT()

    # Throw warning if no burden data for this country
    if (nrow(data_hosp_burden) == 0){
      warning("No RespiCompass hospital burden found for country ", opts$country)
    }

    # Combine admissions (weekly) and burden (seasonal totals per age group)
    combine_df = bind_rows(data_hosp_admissions, data_hosp_burden)
    fit$data = combine_df

    # Check no data is negative
    if (any(fit$data$value < 0, na.rm = TRUE)){
      stop("Negative data values identified")
    }

    return(fit)
  }
  
  if (opts$data_source$epi == "USER-DEFINED") {
    
    combined_df = readRDS(paste0(o$pth$data_epi, "/", o$analysis_name, ".RDS"))
    
    fit$data = combined_df
    
    # Check no data is negative
    if (any(fit$data$value < 0))
      stop("Negative data values identified")
    
    return(fit)
  }
}

# -------------------------------------------------------- -
# Filter to dates of interest ----
# -------------------------------------------------------- -
get_data_dates = function(opts) {

  # Calibrating to epi data: go back to start of burn in
  days = list(fit = opts$data_days, burn = opts$data_burn_in)
  
  # Dates of data we'll fit to
  date_end   = format_date(opts$data_end)
  date_start = date_end - days$fit - days$burn + 1
  
  # Sequence of dates for fitting and plotting
  all_dates  = seq(date_start,  date_end,  by = "day")
  
  # Datatable of fitting dates and day indices
  dates_df = data.table(day  = 1 : length(all_dates), 
                        date = all_dates) %>%
    filter(day > days$burn)
  
  return(dates_df)
}

# ---------------------------------------------------------------- -
# Change time period from day to week, month, quarter, or year ----
# ---------------------------------------------------------------- -
change_time_period = function(data, dates_df, data_period) {
  
  # A trivial process is looking at daily data
  if (data_period != "day") {
    
    # Changing daily data to different time period 
    data = data %>%
      mutate(date = ceiling_date(date, data_period)) %>%
      group_by(metric, date) %>%
      summarise(value = sum(value)) %>% 
      ungroup() %>%
      inner_join(dates_df, by = "date") %>%
      filter(!is.na(value)) %>%
      select(date, day, metric, value) %>%
      setDT()
  }
  
  return(data)
}
