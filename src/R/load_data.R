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
      filter(country == opts$country_name) %>%
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

    # --- weekly hospital admissions (only total age group) --- #
    # Load raw hospital admissions data from RespiCompass
    raw_data = read.csv(o$respicompass$hospital_admissions, fileEncoding = "UTF-8-BOM")
    
    # Map to connect country name with country ISO2
    country_map = read.csv("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/refs/heads/main/supporting-files/countries.csv")
    
    # Select only columns of interest and convert dates to R-interpretable
    data_hosp_admissions = raw_data %>%
      mutate(date  = as.Date(date_wk_floor) + 6,
             value = case_counts,
             country = opts$country_name) %>%
      select(date, country, value) %>%
      left_join(country_map, by=c("country")) %>%
      filter(iso2_code == opts$country) %>%
      select(-country, -iso2_code) %>%
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
      warning("No RespiCompass hospital data found for country ", opts$country)
    }
    
    
    # --- Total hospital burden (admissions) for each age group --- #
    # Load raw total hospital admissions data from RespiCompass
    raw_data = read.csv(o$respicompass$hospital_burden_agegroups, fileEncoding = "UTF-8-BOM")
    
    # Compute "Weekly totals -> 4-week totals" to be used to estimate values per age group below
    periods <- raw_data %>%               # your 2nd dataframe
      distinct(date_28days_floor) %>%
      mutate(
        period_start = as.Date(date_28days_floor) + 6,
        period_end   = as.Date(date_28days_floor) + 6 + weeks(4)
      )
    weekly_4wk <- data_hosp_admissions %>%
      crossing(periods) %>%
      filter(date >= period_start & date < period_end) %>%
      group_by(period_start) %>%
      summarise(
        total_4wk = sum(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(date = period_start)
    
    # Select only columns of interest
    data_hosp_burden = raw_data %>%
      mutate(date = as.Date(date_28days_floor) + 6,
             age_group = age_gp_modelling,
             value = NA,
             country = opts$country_name) %>%
      select(country, date, age_group, value, proportion) %>%
      left_join(country_map, by=c("country")) %>%
      filter(iso2_code == opts$country) %>%
      select(-country, -iso2_code) %>%
      mutate(date = format_date(date),
             metric = "hospital_admissions",
             data_freq = "4-weekly") %>%
      # From proportions to value
      left_join(weekly_4wk, by = "date") %>%
      mutate(
        value = total_4wk * proportion
      ) %>%
      filter(!is.na(value)) %>%
      # Change age group names
      mutate(age_group = case_when(age_group  == "< 3 months" ~ "0-3m",
                                   age_group  == "3-5 months" ~ "3-6m",
                                   age_group  == "6-11 months" ~ "6-12m",
                                   age_group  == "1-4 years" ~ "1-5y",
                                   age_group  == "5-64 years" ~ "5-65y",
                                   age_group == "65+ years" ~ "65+y",
                                   TRUE ~ age_group)) %>%
      select(date, age_group, value, metric, data_freq) %>%
      setDT()
    
    # Throw warning if no hospital data for this country
    if (nrow(data_hosp_burden) == 0){
      warning("No RespiCompass hospital data found for country ", opts$country)
    }
    
    
    
    # Filter dates of interest for fitting and plotting
    combine_df = bind_rows(data_hosp_admissions, data_hosp_burden)
    fit$data = combine_df
    
    # Check no data is negative
    if (any(fit$data$value < 0)){
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
