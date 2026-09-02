########################################################## #
# RESULTS EVALUATION
#
# Post-hoc evaluation of RSVinea model output for the
# RespiCompass 2025/2026 submission round.
#
# Loads scenario outputs produced by run_scenarios() and
# formats them for RespiCompass submission, then generates
# comparative plots against the static model and observed data.
#
# Run after launch.R has completed.
#
# STATIC MODEL COMPARISON
# The second half of this script compares results against the
# companion static RSV model. Its output file must be placed at:
#
#   output/3_results/IRL_respiCompass_2025_2026_results_staticModel.parquet
#
# The file can be obtained from the static model repository:
#   https://github.com/safinea-respiratory/safinea-static-RSV-model
# (output folder, file: IRL_respiCompass_2025_2026_results_staticModel.parquet)
#
# The comparison section is skipped automatically if the file
# is not present.
########################################################## #

# Clear global environment
rm(list = ls())

# Load all required packages and functions
source("R/dependencies.R")

# Load country list
country_list = read.csv("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/refs/heads/main/supporting-files/countries.csv")
pop_df = read.csv("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/refs/heads/main/Previous_Rounds/2025-2026_round_1/auxiliary-data/population/population_estimates.csv")
burden_df = read.csv("data/epidemiological/RSV_monthly_prop_age.csv")

# Age group mapping - in correct order
age_levels = data.frame( age_group = c("0-1m", "1-2m", "2-3m", "3-4m", "4-5m", "5-6m", 
                                       "6-7m", "7-8m", "8-9m", "9-10m", "10-11m", "11-12m", 
                                       "12-13m", "13-14m", "14-15m", "15-16m", "16-17m", "17-18m", 
                                       "18-19m", "19-20m", "20-21m", "21-22m", "22-23m", "23-24m", 
                                       "24-30m", "30-36m", "3-4y",  "4-5y", "5-18y", 
                                       "18-65y",
                                       "65+y"))
age_group_map = age_levels %>%
  mutate(age_group_respiCompass = c("0-3m", "0-3m", "0-3m",
                                    "3-6m", "3-6m", "3-6m",
                                    "6-12m", "6-12m", "6-12m", "6-12m", "6-12m", "6-12m",
                                    "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y", "1-5y",
                                    "5-65y", "5-65y",
                                    "65+y"))

# Define metrics to show
metric_levels = c("hospital_admissions", "cases")

# Load results for the configured country
country_code = "IE"

# Set options
o = set_options(do_step = c(1:3), analysis_name = country_code)

# Collate and interpret inputs so we know what to plot
#list[f, baseline] = fig_properties(o, list(...))
f = NULL
# Full list of scenarios as defined in yaml file 
f$scenarios = parse_yaml(o, "*read*") %>% names()

# Full scenario names as defined in yaml file 
f$scenario_names = parse_yaml(o, "*read*") %>% unname()

# ---- Extract model predictions ----
# Initiate the dataframe
results_list = list()

# Loop through scenarios to plot
for (scenario in f$scenarios) {
  
  #result = try_load(o$pth$scenarios, paste0(scenario, '_100k'))
  result = try_load(o$pth$scenarios, paste0(scenario, '_raw'))
  
  # Format model output and store in list to be concatenated
  results_list[[scenario]] = format_results(o, result)  
}

# Concatenate plotting dataframes for all scenarios
results_df = rbindlist(results_list) 
results_df$location = country_code

# ---- Prepare output ----

# Formatting the dates
duration = results_df %>% select(time) %>% unique() %>% nrow()
all_dates = seq(from = ymd(o$plot_start_date), by = "day", length.out = duration) 
dates_df  = data.table(date = all_dates, 
                       time  = 1 : length(all_dates))

# Get a dataframe with "incidence_prop_vacc" metric to be able to separate vacc and unvacc numbers for each metric
prop_vacc_df = results_df %>%
  filter(metric == "incidence_prop_vacc") %>% 
  rename(prop_outcomes_in_vacc = value) %>%
  select(-metric)

# Prepare output
results_df = results_df %>% inner_join(dates_df, by = "time") %>%
  # Add needed columns
  filter(metric %in% metric_levels) %>%
  left_join(prop_vacc_df, by = join_by(time, age_group, variant, param_id, scenario, location)) %>%
  left_join(age_group_map, by = join_by(age_group)) %>%
  mutate(date_week = floor_date(as.Date(date), "week", week_start = 1),
         scenario = factor(scenario, levels = f$scenarios),
         param_id = str_replace(param_id, "^((?:[^_]+_){2}).*$", "\\1")
  ) %>%
  select(-time) %>%
  # Compute outcome proportion in vacc vs unvacc
  mutate(value_vacc = prop_outcomes_in_vacc * value) %>%
  # Aggregate to weeks and lower number of age groups
  group_by(location, date_week, age_group_respiCompass, metric, variant, param_id, scenario) %>%
  summarise(value = sum(value),
            value_vacc = sum(value_vacc),
            value_unvacc = value - value_vacc,
            , .groups = "drop") %>%
  rename(date = date_week,
         age_group = age_group_respiCompass)




results_all_df_output = results_df

#write_parquet(results_all_df_output, "output/3_results/IRL_respiCompass_2025_2026_results.parquet")
#results_all_df_output = read_parquet("output/3_results/IRL_respiCompass_2025_2026_results.parquet")


### --- Prepare file for submission

respiCompass_df_submission_pre = results_all_df_output %>%
  mutate(target = ifelse(metric=="hospital_admissions", "rsv_hospitalisations",
                         ifelse(metric=="cases","rsv_infections", metric))) %>%
  # Rename age groups
  mutate(age_group = case_when(age_group  == "0-3m" ~ "0-2mo",
                               age_group  == "3-6m" ~ "3-5mo",
                               age_group  == "6-12m" ~ "6-11mo",
                               age_group  == "1-5y" ~ "1-4y",
                               age_group  == "5-65y" ~ "5-64y",
                               TRUE ~ age_group)) %>%
  #Rename baseline scenario
  mutate(scenario = as.character(scenario)) %>%
  #Rename other relevant variables
  mutate(target_end_date = date + 6,
         horizon = as.integer( (target_end_date - as.Date("2025-09-07")) / 7 ),
         round_id = '2025_2026_1_RSV',
         output_type='sample',
         #
         a = as.integer(str_extract(param_id, "(?<=s)\\d{4}(?=_)")),
         b = as.integer(str_extract(param_id, "(?<=_)(\\d{4})(?=_)")),
         a = as.integer(a),
         b = as.integer(b),
         output_type_id = as.character( (a - 1) * 30 + b )
  ) %>%
  filter(horizon<39) %>% # requirement by the RespiCompass
  rename(scenario_id = scenario) %>%
  select(-variant) %>%
  # Go from wide (value_vacc/value_unvacc/value) to long
  pivot_longer(
    cols = c(value_vacc, value_unvacc, value),
    names_to = "imm_raw",
    values_to = "value"
  ) %>%
  # Map column name -> imm status label
  mutate(
    imm_status = case_when(
      imm_raw == "value_vacc"   ~ "immYes",
      imm_raw == "value_unvacc" ~ "immNo",
      imm_raw == "value"        ~ "immTotal"
    ),
    # Combine age group + imm status
    pop_group = paste0(age_group, "_", imm_status)
  ) %>%
  # Select only necessary columns
  dplyr::select(round_id, scenario_id, location, target, 
                pop_group, horizon, target_end_date, output_type, output_type_id, value, imm_status) %>%
  # Add total age group
  { 
    df_long <- .
    
    totals <- df_long %>%
      group_by(
        round_id, scenario_id, location, target,
        horizon, target_end_date, output_type, output_type_id,
        imm_status
      ) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      mutate(pop_group = paste0("total_", imm_status))
    
    bind_rows(df_long, totals)
  } %>%
  select(-imm_status)

respiCompass_df_submission = respiCompass_df_submission_pre

#write_parquet(respiCompass_df_submission, "output/3_results/respiCompass_2025_2026_results_dynamicModel.parquet", compression = "gzip")
#respiCompass_df_submission = read_parquet("output/3_results/respiCompass_2025_2026_results_dynamicModel.parquet")



#### --- Plot infections over time

# Pop data
pop_df_total = pop_df %>% 
  group_by(country) %>%
  summarise(population = sum(population)) %>% 
  left_join(country_list) %>%
  rename(location = iso2_code)

# Time plot 
respiCompass_df_submission %>%
  filter(pop_group == "total_immTotal",
         target == "rsv_infections",
         scenario_id == "baseline") %>%
  left_join(pop_df_total) %>%
  mutate(value = value / (population/1e5)) %>%
  group_by(location, scenario_id, target, pop_group, horizon) %>%   # drop x if you only want 1 point per country
  summarise(
    median = median(value, na.rm = TRUE),
    lower  = quantile(value, 0.025, na.rm = TRUE),
    upper  = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Plot
  ggplot(aes(x = horizon, y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(size = 0.8) +
  facet_wrap(~location, scales = "free_x") + 
  theme_minimal() +
  labs(
    y = "Value per 1e5",
    x = "horizon"
  )



baseline_results = respiCompass_df_submission %>%
  filter(scenario_id == "baseline") %>%
  rename(value_baseline = value) %>%
  select(-scenario_id)

plot_df = left_join(respiCompass_df_submission, baseline_results) %>%
  filter(scenario_id != "baseline") %>%
  mutate(rel_change = (value-value_baseline)  / value_baseline)




# Plot burden by age group
respiCompass_df_submission %>%
  filter(target=="rsv_hospitalisations",
         str_ends(pop_group, "immTotal"),
         target_end_date >= as.Date("2025-09-01"),
         target_end_date <= as.Date("2026-03-01")
  ) %>%
  # Sum over all dates
  group_by(round_id, scenario_id, location, target, pop_group, output_type, output_type_id) %>%
  # Sum over all dates
  summarise(value = sum(value)
            , .groups = "drop") %>%
  # Plot
  ggplot(aes(x = scenario_id, y = value)) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  #geom_hline(yintercept = 0) +
  facet_wrap(~pop_group, scales = "free") + 
  theme_bw() +
  labs(x="Population group", y="", title="Hospitalisation burden by scenario")

# Plot burden by age group per week
respiCompass_df_submission %>%
  filter(target=="rsv_hospitalisations",
         str_ends(pop_group, "immTotal"),
         target_end_date >= as.Date("2025-09-01")) %>%
  # Sum over all dates
  group_by(round_id, target_end_date, horizon, location, target, pop_group, output_type, scenario_id) %>%
  # Sum over all dates
  summarise(
    median = median(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop") %>%
  # Plot
  ggplot(aes(x = horizon, y = median, color = scenario_id, fill = scenario_id)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  #geom_hline(yintercept = 0) +
  facet_wrap(~pop_group, scales = "free") + 
  theme_bw() +
  labs(x="Population group", y="", title="Hospitalisation burden by scenario")

# Plot burden by age group per week ONLY FOR THOSE BORN IN THE SEASON
cohort_start <- as.Date("2024-09-01")
cohort_end   <- as.Date("2025-02-28")
respiCompass_df_submission %>%
  # filter(horizon == -20,
  #        scenario_id == "baseline",
  #        target=="rsv_hospitalisations",
  #        output_type_id == 1,
  #        str_ends(pop_group, "immTotal")) %>%
  mutate(
    age_band = str_extract(pop_group, "^[^_]+"),
    age_start = as.numeric(str_extract(age_band, "^\\d+")),
    age_end   = as.numeric(str_extract(age_band, "(?<=-)\\d+")),
    age_unit  = str_extract(age_band, "mo|y")
  ) %>%
  mutate(
    birth_latest = case_when(
      age_unit == "mo" ~ target_end_date %m-% months(age_start),
      age_unit == "y"  ~ target_end_date %m-% years(age_start)
    ),
    birth_earliest = case_when(
      age_unit == "mo" ~ target_end_date %m-% months(age_end+1),
      age_unit == "y"  ~ target_end_date %m-% years(age_end+1)
    )
  ) %>%
  filter(
    birth_latest >= cohort_start,
    birth_earliest <= cohort_end
  ) %>%
  filter(target=="rsv_hospitalisations",
         str_ends(pop_group, "immTotal"),
         target_end_date >= as.Date("2024-09-01"),
         target_end_date <= as.Date("2025-03-01"),
  ) %>%
  #
  mutate(
    overlap_start = pmax(birth_earliest, cohort_start),
    overlap_end   = pmin(birth_latest, cohort_end),
    
    overlap_days = as.numeric(overlap_end - overlap_start + 1),
    window_days  = as.numeric(birth_latest - birth_earliest + 1),
    
    prop_in_cohort = case_when(
      overlap_end >= overlap_start ~ overlap_days / window_days,
      TRUE ~ 0
    ),
    
    value_adj = value * prop_in_cohort
  ) %>%
  # Sum over all relevant age groups
  group_by(round_id, target_end_date, horizon, location, target, output_type, scenario_id, output_type_id) %>%
  summarise(value_adj = sum(value_adj, na.rm= TRUE),
            .groups = "drop") %>%
  # Get median and quantiles
  group_by(round_id, target_end_date, horizon, location, target, output_type, scenario_id) %>%
  summarise(
    median = median(value_adj, na.rm = TRUE),
    lower = quantile(value_adj, 0.025, na.rm = TRUE),
    upper = quantile(value_adj, 0.975, na.rm = TRUE),
    .groups = "drop") %>%
  # Plot
  ggplot(aes(x = horizon, y = median, color = scenario_id, fill = scenario_id)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  #geom_hline(yintercept = 0) +
  #facet_wrap(~pop_group, scales = "free") + 
  theme_bw() +
  labs(x="Horizon", y="", title="Hospitalisation burden by scenario")





########################################################## #
# STATIC vs DYNAMIC MODEL COMPARISON
#
# Compares RSVinea (dynamic ODE model) against the companion
# static model. To run this section, place the static model
# output file at:
#
#   output/3_results/IRL_respiCompass_2025_2026_results_staticModel.parquet
#
# Obtain it from:
#   https://github.com/safinea-respiratory/safinea-static-RSV-model
#   (output/IRL_respiCompass_2025_2026_results_staticModel.parquet)
#
# This section is skipped automatically if the file is not present.
########################################################## #

# NB: output/ is namespaced by git branch (see set_dirs in R/directories.R), so
# this must be resolved via o$pth$results rather than hardcoded as "output/...".
static_model_path <- file.path(o$pth$results,
                               "IRL_respiCompass_2025_2026_results_staticModel.parquet")

if (!file.exists(static_model_path)) {
  message("Skipping static vs dynamic model comparison.")
  message("Place the static model output file at: ", static_model_path)
  message("Obtain it from: https://github.com/safinea-respiratory/safinea-static-RSV-model")
} else {
  
  static_df_load <- read_parquet(static_model_path) %>%
    mutate(Model = "Static", target_end_date = as.Date(target_end_date))
  
  ### LOAD DATA
  # Load weekly data
  hospital_admissions_df = read.csv("data/epidemiological/RSV_weekly_counts.csv") %>%
    mutate(target_end_date = as.Date(date_wk_floor) + 6,
           weekly_rsv_hospitalisations = case_counts) %>%
    select(target_end_date, season_name, weekly_rsv_hospitalisations) %>%
    setDT()
  # Load raw total hospital admissions data from RespiCompass
  raw_data = read.csv("data/epidemiological/RSV_monthly_prop_age.csv", fileEncoding = "UTF-8-BOM")
  # Compute "Weekly totals -> 4-week totals" to be used to estimate values per age group below
  periods <- raw_data %>%               # your 2nd dataframe
    distinct(date_28days_floor) %>%
    mutate(
      period_start = as.Date(date_28days_floor) + 6,
      period_end   = as.Date(date_28days_floor) + 6 + weeks(3)
    )
  weekly_4wk <- hospital_admissions_df %>%
    crossing(periods) %>%
    filter(target_end_date >= period_start & target_end_date <= period_end) %>%
    group_by(period_start) %>%
    summarise(
      total_4wk = sum(weekly_rsv_hospitalisations, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(date = period_start)
  
  # Select only columns of interest
  hospital_burden_df = raw_data %>%
    mutate(date = as.Date(date_28days_floor) + 6,
           age_group = age_gp_modelling,
           value = NA) %>%
    select(date, age_group, value, proportion) %>%
    # From proportions to value
    left_join(weekly_4wk, by = "date") %>%
    mutate(
      value = total_4wk * proportion,
      burden_start_date = date,
      burden_end_date = date + weeks(3),
      total_rsv_hospitalisations = value
    ) %>%
    filter(!is.na(value)) %>%
    # Change age group names
    mutate(age_group = case_when(age_group  == "< 3 months" ~ "0-2mo",
                                 age_group  == "3-5 months" ~ "3-5mo",
                                 age_group  == "6-11 months" ~ "6-11mo",
                                 age_group  == "1-4 years" ~ "1-4y",
                                 age_group  == "5-64 years" ~ "5-64y",
                                 age_group == "65+ years" ~ "65+y",
                                 TRUE ~ age_group),
           pop_group = paste0(age_group,"_immTotal"),
           season_start_year = if_else(month(burden_start_date) >= 8,
                                       year(burden_start_date),
                                       year(burden_start_date) - 1),
           season = paste0(season_start_year, "/", season_start_year + 1)) %>%
    select(burden_start_date, burden_end_date, age_group, pop_group, total_rsv_hospitalisations, season, season_start_year) %>%
    setDT()
  
  
  # static_df_baseline = static_df_load %>%
  #   filter(scenario_id == "E") %>%
  #   rename(value_baseline = value) %>%
  #   select(-scenario_id)
  # 
  # static_df = static_df_load %>% left_join(static_df_baseline) %>%
  #   filter(scenario_id != "E")
  
  
  # Merge both models together
  df_both_models = bind_rows(static_df_load  %>% mutate(location = "IE") %>% filter(scenario_id != "test_vacc"), 
                             respiCompass_df_submission %>% mutate(Model = "Dynamic"))
  
  age_labels <- c(
    "0-2mo_immTotal"   = "<3 months",
    "3-5mo_immTotal"  = "3-6 months",
    "6-11mo_immTotal"  = "6-12 months",
    "1-4y_immTotal"  = "1-5 years",
    "5-64y_immTotal" = "15-64 years",
    "65+y_immTotal"   = "65+ years",
    "total_immTotal"   = "all"
  )
  
  age_order <- c(
    "0-2mo_immTotal",
    "3-5mo_immTotal",
    "6-11mo_immTotal",
    "1-4y_immTotal",
    "5-64y_immTotal",
    "65+y_immTotal",
    "total_immTotal"
  )
  
  
  # Calibration of static model
  static_df_load %>% 
    filter(
      scenario_id == "baseline",
      pop_group == "total_immTotal"
    ) %>%
    group_by(target_end_date) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(hospital_admissions_df) %>%
    ggplot(aes(x = target_end_date)) +
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper,
        fill = "Model variation"
      ),
      alpha = 0.3
    ) +
    geom_line(
      aes(
        y = median_value,
        color = "Model median"
      ),
      linewidth = 1
    ) +
    geom_point(
      aes(
        y = weekly_rsv_hospitalisations,
        color = "Observed data"
      ),
      size = 2
    ) +
    scale_fill_manual(
      values = c("Model variation" = "steelblue")
    ) +
    scale_color_manual(
      values = c(
        "Model median" = "steelblue",
        "Observed data" = "darkred"
      )
    ) +
    labs(
      x = "Date",
      y = "Weekly hospitalisations",
      fill = "",
      color = ""
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  # Calibration of dynamic model
  df_both_models %>% 
    filter(
      Model == "Dynamic",
      scenario_id == "baseline",
      pop_group == "total_immTotal",
      target == "rsv_hospitalisations"
    ) %>%
    group_by(target_end_date) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(hospital_admissions_df) %>%
    ggplot(aes(x = target_end_date)) +
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper,
        fill = "Model variation"
      ),
      alpha = 0.3
    ) +
    geom_line(
      aes(
        y = median_value,
        color = "Model median"
      ),
      linewidth = 1
    ) +
    geom_point(
      aes(
        y = weekly_rsv_hospitalisations,
        color = "Observed data"
      ),
      size = 2
    ) +
    scale_fill_manual(
      values = c("Model variation" = "steelblue")
    ) +
    scale_color_manual(
      values = c(
        "Model median" = "steelblue",
        "Observed data" = "darkred"
      )
    ) +
    labs(
      x = "Date",
      y = "Weekly hospitalisations",
      fill = "",
      color = ""
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ################################################
  ############ 4-weekly age specific burden #######
  ################################################
  target_seasons <- c("2023/2024", "2024/2025", "2025/2026")
  
  # ── 4-weekly burden ────────────────────────────────────────────────────────────
  
  model_burden_4wk <- df_both_models %>%
    filter(
      scenario_id == "baseline",
      target == "rsv_hospitalisations",
      pop_group %in% names(age_labels),
      pop_group != "total_immTotal"
    ) %>%
    left_join(
      hospital_burden_df %>%
        select(burden_start_date, burden_end_date, pop_group) %>%
        distinct(),
      by = join_by(pop_group, target_end_date >= burden_start_date, target_end_date <= burden_end_date)
    ) %>%
    group_by(Model, pop_group, burden_start_date, burden_end_date, output_type_id) %>%
    summarise(burden = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Model, pop_group, burden_start_date, burden_end_date) %>%
    summarise(
      median_value = median(burden, na.rm = TRUE),
      lower        = quantile(burden, 0.025, na.rm = TRUE),
      upper        = quantile(burden, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mid_date  = burden_start_date + (burden_end_date - burden_start_date) / 2,
      age_label = factor(age_labels[pop_group], levels = age_labels[age_order])
    )
  
  obs_burden_4wk <- hospital_burden_df %>%
    filter(pop_group %in% names(age_labels), 
           pop_group != "total_immTotal",
           burden_start_date >= ymd("2023-08-01")) %>%
    mutate(
      mid_date  = burden_start_date + (burden_end_date - burden_start_date) / 2,
      age_label = factor(age_labels[pop_group], levels = age_labels[age_order])
    )
  
  plot_4wk <- function(model_name) {
    ggplot() +
      geom_point(
        data = obs_burden_4wk,
        aes(x = mid_date, y = total_rsv_hospitalisations, color = "Observed"),
        size = 2
        #position = position_nudge(x = 12.5)
      ) +
      geom_pointrange(
        data = model_burden_4wk %>% filter(Model == model_name, burden_start_date >= ymd("2023-08-01")),
        aes(x = mid_date, y = median_value, ymin = lower, ymax = upper, color = "Model"),
        size = 0.3,
        shape = 21,
        fill = NA
      ) +
      facet_wrap(~ age_label, scales = "free_y", ncol = 3) +
      scale_color_manual(
        values = c("Model" = "black", "Observed" = "darkred")
      ) +
      labs(
        title = paste(model_name, "model — 4-weekly burden"),
        x     = "Time",
        y     = "4-weekly burden",
        color = ""
      ) +
      coord_cartesian(ylim = c(0, NA)) +
      theme_bw() +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey85"),
        legend.position  = "bottom"
      )
  }
  
  plot_4wk("Static")
  plot_4wk("Dynamic")
  
  
  ################################################
  ############ Seasonal age specific burden #######
  ################################################
  
  # Build season lookup
  date_to_season <- hospital_burden_df %>%
    select(burden_start_date, burden_end_date, pop_group, season) %>%
    distinct()
  
  # Convert to data.table
  model_dt <- as.data.table(
    df_both_models %>%
      filter(
        scenario_id == "baseline",
        target == "rsv_hospitalisations",
        pop_group %in% names(age_labels),
        pop_group != "total_immTotal"
      )
  )
  
  burden_dt <- as.data.table(date_to_season)
  
  # Non-equi join to assign season
  model_with_season <- burden_dt[
    model_dt,
    on = .(
      pop_group == pop_group,
      burden_start_date <= target_end_date,
      burden_end_date >= target_end_date
    )
  ]
  
  # Aggregate to seasonal burden per sample, then summarise
  model_seasonal <- model_with_season %>%
    filter(!is.na(season)) %>%
    group_by(Model, pop_group, season, output_type_id) %>%
    summarise(seasonal_burden = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Model, pop_group, season) %>%
    summarise(
      median_value = median(seasonal_burden, na.rm = TRUE),
      lower = quantile(seasonal_burden, 0.05, na.rm = TRUE),
      upper = quantile(seasonal_burden, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Observed seasonal burden
  obs_seasonal <- hospital_burden_df %>%
    filter(pop_group %in% names(age_labels), pop_group != "total_immTotal") %>%
    group_by(pop_group, season) %>%
    summarise(seasonal_burden = sum(total_rsv_hospitalisations, na.rm = TRUE), .groups = "drop") %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Seasons to display
  target_seasons <- c("2023/2024", "2024/2025", "2025/2026")
  
  # Plot
  ggplot() +
    geom_pointrange(
      data = model_seasonal %>% filter(season %in% target_seasons),
      aes(x = season, y = median_value, ymin = lower, ymax = upper, color = Model),
      position = position_dodge(width = 0.4),
      size = 0.4
    ) +
    geom_point(
      data = obs_seasonal %>% filter(season %in% target_seasons),
      aes(x = season, y = seasonal_burden, color = "Observed"),
      size = 3,
      position = position_nudge(x = 0.)
    ) +
    facet_wrap(~ age_label, scales = "free_y", ncol = 3) +
    scale_color_manual(
      values = c("Static" = "lightblue", "Dynamic" = "black", "Observed" = "darkred")
    ) +
    labs(
      title = "Seasonal burden by model",
      x = "Season",
      y = "Seasonal burden",
      color = ""
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey85"),
      legend.position = "bottom"
    )
  
  
  
  #################################
  #################################
  
  
  
  
  
  
  
  ################################################
  ############ Seasonal age specific burden - NO VACC #######
  ################################################
  
  # Build season lookup
  date_to_season <- hospital_burden_df %>%
    select(burden_start_date, burden_end_date, pop_group, season) %>%
    distinct()
  
  # Convert to data.table
  model_dt <- as.data.table(
    df_both_models %>%
      filter(
        scenario_id == "no_vacc",
        target == "rsv_hospitalisations",
        pop_group %in% names(age_labels),
        pop_group != "total_immTotal"
      )
  )
  
  burden_dt <- as.data.table(date_to_season)
  
  # Non-equi join to assign season
  model_with_season <- burden_dt[
    model_dt,
    on = .(
      pop_group == pop_group,
      burden_start_date <= target_end_date,
      burden_end_date >= target_end_date
    )
  ]
  
  # Aggregate to seasonal burden per sample, then summarise
  model_seasonal <- model_with_season %>%
    filter(!is.na(season)) %>%
    group_by(Model, pop_group, season, output_type_id) %>%
    summarise(seasonal_burden = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Model, pop_group, season) %>%
    summarise(
      median_value = median(seasonal_burden, na.rm = TRUE),
      lower = quantile(seasonal_burden, 0.05, na.rm = TRUE),
      upper = quantile(seasonal_burden, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Observed seasonal burden
  obs_seasonal <- hospital_burden_df %>%
    filter(pop_group %in% names(age_labels), pop_group != "total_immTotal") %>%
    group_by(pop_group, season) %>%
    summarise(seasonal_burden = sum(total_rsv_hospitalisations, na.rm = TRUE), .groups = "drop") %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Seasons to display
  target_seasons <- c("2023/2024", "2024/2025", "2025/2026")
  
  # Plot
  ggplot() +
    geom_pointrange(
      data = model_seasonal %>% filter(season %in% target_seasons),
      aes(x = season, y = median_value, ymin = lower, ymax = upper, color = Model),
      position = position_dodge(width = 0.4),
      size = 0.4
    ) +
    # geom_point(
    #   data = obs_seasonal %>% filter(season %in% target_seasons),
    #   aes(x = season, y = seasonal_burden, color = "Observed"),
    #   size = 3,
    #   position = position_nudge(x = 0.)
    # ) +
    facet_wrap(~ age_label, scales = "free_y", ncol = 3) +
    scale_color_manual(
      values = c("Static" = "lightblue", "Dynamic" = "black", "Observed" = "darkred")
    ) +
    labs(
      title = "Seasonal burden by model",
      x = "Season",
      y = "Seasonal burden",
      color = ""
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey85"),
      legend.position = "bottom"
    )
  
  
  
  #################################
  #################################
  
  
  
  
  
  ################################################
  ############ Seasonal age specific burden - HIGH VACC #######
  ################################################
  
  # Build season lookup
  date_to_season <- hospital_burden_df %>%
    select(burden_start_date, burden_end_date, pop_group, season) %>%
    distinct()
  
  # Convert to data.table
  model_dt <- as.data.table(
    df_both_models %>%
      filter(
        scenario_id == "high_vacc",
        target == "rsv_hospitalisations",
        pop_group %in% names(age_labels),
        pop_group != "total_immTotal"
      )
  )
  
  burden_dt <- as.data.table(date_to_season)
  
  # Non-equi join to assign season
  model_with_season <- burden_dt[
    model_dt,
    on = .(
      pop_group == pop_group,
      burden_start_date <= target_end_date,
      burden_end_date >= target_end_date
    )
  ]
  
  # Aggregate to seasonal burden per sample, then summarise
  model_seasonal <- model_with_season %>%
    filter(!is.na(season)) %>%
    group_by(Model, pop_group, season, output_type_id) %>%
    summarise(seasonal_burden = sum(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Model, pop_group, season) %>%
    summarise(
      median_value = median(seasonal_burden, na.rm = TRUE),
      lower = quantile(seasonal_burden, 0.05, na.rm = TRUE),
      upper = quantile(seasonal_burden, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Observed seasonal burden
  obs_seasonal <- hospital_burden_df %>%
    filter(pop_group %in% names(age_labels), pop_group != "total_immTotal") %>%
    group_by(pop_group, season) %>%
    summarise(seasonal_burden = sum(total_rsv_hospitalisations, na.rm = TRUE), .groups = "drop") %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Seasons to display
  target_seasons <- c("2023/2024", "2024/2025", "2025/2026")
  
  # Plot
  ggplot() +
    geom_pointrange(
      data = model_seasonal %>% filter(season %in% target_seasons),
      aes(x = season, y = median_value, ymin = lower, ymax = upper, color = Model),
      position = position_dodge(width = 0.4),
      size = 0.4
    ) +
    # geom_point(
    #   data = obs_seasonal %>% filter(season %in% target_seasons),
    #   aes(x = season, y = seasonal_burden, color = "Observed"),
    #   size = 3,
    #   position = position_nudge(x = 0.)
    # ) +
    facet_wrap(~ age_label, scales = "free_y", ncol = 3) +
    scale_color_manual(
      values = c("Static" = "lightblue", "Dynamic" = "black", "Observed" = "darkred")
    ) +
    labs(
      title = "Seasonal burden by model",
      x = "Season",
      y = "Seasonal burden",
      color = ""
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey85"),
      legend.position = "bottom"
    )
  
  
  
  #################################
  #################################
  
  
  
  
  
  
  
  
  ################################################
  ############ Relative change in burden by season - NO VACC #######
  ################################################
  
  # Build season lookup
  date_to_season <- hospital_burden_df %>%
    select(burden_start_date, burden_end_date, pop_group, season) %>%
    distinct()
  
  # Convert to data.table - include all scenarios
  model_dt <- as.data.table(
    df_both_models %>%
      filter(
        scenario_id %in% c("baseline", "no_vacc", "high_vacc"),
        target == "rsv_hospitalisations",
        pop_group %in% names(age_labels),
        pop_group != "total_immTotal"
      )
  )
  
  burden_dt <- as.data.table(date_to_season)
  
  # Non-equi join to assign season
  model_with_season <- burden_dt[
    model_dt,
    on = .(
      pop_group         == pop_group,
      burden_start_date <= target_end_date,
      burden_end_date   >= target_end_date
    )
  ]
  
  # Aggregate to seasonal burden per sample per scenario
  model_seasonal_samples <- model_with_season %>%
    filter(!is.na(season)) %>%
    group_by(Model, pop_group, season, scenario_id, output_type_id) %>%
    summarise(seasonal_burden = sum(value, na.rm = TRUE), .groups = "drop")
  
  # Baseline samples
  baseline_samples <- model_seasonal_samples %>%
    filter(scenario_id == "baseline") %>%
    select(Model, pop_group, season, output_type_id, baseline_burden = seasonal_burden)
  
  # Compute relative change vs baseline per sample
  relative_change <- model_seasonal_samples %>%
    filter(scenario_id != "baseline") %>%
    left_join(baseline_samples, by = c("Model", "pop_group", "season", "output_type_id")) %>%
    mutate(rel_change = (seasonal_burden - baseline_burden) / baseline_burden * 100) %>%
    group_by(Model, pop_group, season, scenario_id) %>%
    summarise(
      median_value = median(rel_change, na.rm = TRUE),
      lower        = quantile(rel_change, 0.05, na.rm = TRUE),
      upper        = quantile(rel_change, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Seasons to display
  target_seasons <- c("2023/2024", "2024/2025", "2025/2026")
  
  # Plot function
  plot_relative <- function(scenario_name) {
    ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_pointrange(
        data = relative_change %>% filter(season %in% target_seasons, 
                                          scenario_id == scenario_name,
                                          pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal","6-11mo_immTotal")),
        aes(
          x     = season,
          y     = median_value,
          ymin  = lower,
          ymax  = upper,
          color = Model
        ),
        position = position_dodge(width = 0.4),
        size = 0.4
      ) +
      facet_wrap(~ age_label, scales = "free_y", ncol = 3) +
      scale_color_manual(
        values = c("Static" = "lightblue", "Dynamic" = "black")
      ) +
      labs(
        title = paste("Relative change in seasonal burden vs baseline —", scenario_name),
        x     = "Season",
        y     = "Relative change (%)",
        color = ""
      ) +
      theme_bw() +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey85"),
        legend.position  = "bottom"
      )
  }
  
  plot_relative("no_vacc")
  plot_relative("high_vacc")
  
  
  
  
  # Plot function
  plot_RR <- function(scenario_name) {
    ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_pointrange(
        data = relative_change %>% filter(season %in% target_seasons, 
                                          scenario_id == scenario_name,
                                          pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal","6-11mo_immTotal")),
        aes(
          x     = season,
          y     = median_value/100 + 1,
          ymin  = lower/100 + 1,
          ymax  = upper/100 + 1,
          color = Model
        ),
        position = position_dodge(width = 0.4),
        size = 0.4
      ) +
      facet_wrap(~ age_label, scales = "free_y", ncol = 3) +
      scale_color_manual(
        values = c("Static" = "lightblue", "Dynamic" = "black")
      ) +
      labs(
        title = paste("Relative Risk in seasonal burden vs baseline —", scenario_name),
        x     = "Season",
        y     = "Relative Risk",
        color = ""
      ) +
      theme_bw() +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey85"),
        legend.position  = "bottom"
      )
  }
  
  plot_RR("no_vacc")
  plot_RR("high_vacc")
  
  
  
  ########################################
  ### RR relative to 2023/2024 season ####
  ########################################
  
  # Build season lookup
  date_to_season <- hospital_burden_df %>%
    select(burden_start_date, burden_end_date, pop_group, season) %>%
    distinct()
  
  # Convert to data.table - include all scenarios
  model_dt <- as.data.table(
    df_both_models %>%
      filter(
        #scenario_id %in% c("baseline", "no_vacc", "high_vacc"),
        target == "rsv_hospitalisations",
        pop_group %in% names(age_labels),
        pop_group != "total_immTotal"
      )
  )
  
  burden_dt <- as.data.table(date_to_season)
  
  # Non-equi join to assign season
  model_with_season <- burden_dt[
    model_dt,
    on = .(
      pop_group         == pop_group,
      burden_start_date <= target_end_date,
      burden_end_date   >= target_end_date
    )
  ]
  
  # Aggregate to seasonal burden per sample per scenario
  model_seasonal_samples <- model_with_season %>%
    filter(!is.na(season)) %>%
    group_by(Model, pop_group, season, scenario_id, output_type_id) %>%
    summarise(seasonal_burden = sum(value, na.rm = TRUE), .groups = "drop")
  
  # Reference: 2023/2024 burden per sample per scenario
  reference_samples <- model_seasonal_samples %>%
    filter(season == "2023/2024") %>%
    select(Model, pop_group, scenario_id, output_type_id, reference_burden = seasonal_burden)
  
  # Compute relative risk vs 2023/2024 of the same scenario per sample
  relative_risk <- model_seasonal_samples %>%
    filter(season != "2023/2024") %>%
    left_join(reference_samples, by = c("Model", "pop_group", "scenario_id", "output_type_id")) %>%
    mutate(rr = seasonal_burden / reference_burden) %>%
    group_by(Model, pop_group, season, scenario_id) %>%
    summarise(
      median_value = median(rr, na.rm = TRUE),
      mean_value = mean(rr, na.rm = TRUE),
      lower        = quantile(rr, 0.025, na.rm = TRUE),
      upper        = quantile(rr, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(age_label = factor(age_labels[pop_group], levels = age_labels[age_order]))
  
  # Seasons to display (excluding reference season)
  target_seasons <- c("2024/2025", "2025/2026")
  
  # Plot function
  plot_rr <- function(scenario_name) {
    ggplot() +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      geom_pointrange(
        data = relative_risk %>% filter(season %in% target_seasons, scenario_id == scenario_name),
        aes(
          x     = season,
          y     = median_value,
          ymin  = lower,
          ymax  = upper,
          color = Model
        ),
        position = position_dodge(width = 0.4),
        size = 0.4
      ) +
      facet_wrap(~ age_label, scales = "free_y", ncol = 3) +
      scale_color_manual(
        values = c("Static" = "lightblue", "Dynamic" = "black")
      ) +
      labs(
        title = paste("Relative risk vs 2023/2024 —", scenario_name),
        x     = "Season",
        y     = "Relative risk",
        color = ""
      ) +
      theme_bw() +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey85"),
        legend.position  = "bottom"
      )
  }
  
  plot_rr("baseline")
  plot_rr("no_vacc")
  plot_rr("high_vacc")
  
  
  
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  #################################
  
  
  
  
  
  
  
  
  
  
  ### OTHER PLOTS ###
  
  df_season <- df_both_models %>%
    filter(str_ends(pop_group, "immTotal"),
           target == "rsv_hospitalisations",
           target_end_date <= ymd("2026-01-25")) %>%
    mutate(
      season_start_year = if_else(month(target_end_date) >= 8,
                                  year(target_end_date),
                                  year(target_end_date) - 1),
      season = paste0(season_start_year, "/", season_start_year + 1)
    ) %>%
    group_by(
      season_start_year, season,
      Model, scenario_id, pop_group, location, target, output_type_id
    ) %>%
    summarise(
      total_hosps = sum(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(season >= 2023)
  
  
  
  
  
  ## Total burden - BASELINE
  df_burden_baseline = df_season %>%
    filter(scenario_id == "baseline") %>%
    left_join(hospital_burden_df %>% 
                group_by(pop_group, season, season_start_year) %>% 
                summarise(total_data_hosps = sum(total_rsv_hospitalisations)) %>% 
                ungroup()) %>%
    filter(pop_group %in% age_order) %>%
    mutate(pop_group = factor(pop_group, levels = age_order))
  
  # Plot
  df_burden_baseline %>% 
    filter(pop_group != "total_immTotal") %>%
    ggplot(
      aes(x = season, y = total_hosps, fill = Model)
    ) +
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_point(
      aes(y = total_data_hosps, color = "Observed data"),
      size = 2
    ) +
    facet_wrap( ~ pop_group, scales = "free_y", ncol = 3, labeller = labeller(pop_group = age_labels)) +
    labs(
      x = "Season",
      y = "Total hospitalisations",
      fill = "Model",
      color = ""
    ) +
    scale_color_manual(
      values = c("Observed data" = "darkred")
    ) +
    theme_bw() +
    ylim(c(0, NA)) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  
  
  
  
  ## Total burden - NO VACC
  df_season %>%
    filter(scenario_id == "no_vacc",
           pop_group != "total_immTotal") %>%
    ggplot(
      aes(x = season, y = total_hosps, fill = Model)
    ) +
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    facet_wrap( ~ pop_group, scales = "free_y", ncol = 3, labeller = labeller(pop_group = age_labels)) +
    labs(
      x = "Season",
      y = "Total hospitalisations",
      fill = "Model",
      color = ""
    ) +
    theme_bw() +
    ylim(c(0, NA)) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  ## Total burden - HIGH VACC
  df_season %>%
    filter(scenario_id == "high_vacc",
           pop_group != "total_immTotal") %>%
    ggplot(
      aes(x = season, y = total_hosps, fill = Model)
    ) +
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    facet_wrap( ~ pop_group, scales = "free_y", ncol = 3, labeller = labeller(pop_group = age_labels)) +
    labs(
      x = "Season",
      y = "Total hospitalisations",
      fill = "Model",
      color = ""
    ) +
    theme_bw() +
    ylim(c(0, NA)) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  
  
  # To compare burden and change against an external estimate
  df_season %>% 
    filter(pop_group == "0-2mo_immTotal") %>%
    #filter(scenario_id == "no_vacc") %>%
    ggplot(
      aes(x = season, y = total_hosps, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap(scenario_id ~ pop_group, scales = "free_y", ncol = 7) +
    labs(
      x = "Season",
      y = "Total hospitalisations",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  # Plot relative decrease
  df_rel <- df_season %>%
    group_by(
      season_start_year, season,
      Model, pop_group, location, target,
      output_type_id
    ) %>%
    mutate(
      baseline_hosps = total_hosps[scenario_id == "baseline"][1]
    ) %>%
    ungroup() %>%
    filter(scenario_id != "baseline") %>%
    mutate(
      rel_decrease = total_hosps / baseline_hosps - 1
    )
  
  df_rel %>%
    ggplot(
      aes(x = season, y = rel_decrease, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap(scenario_id ~ pop_group, scales = "free_y", ncol = 7) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Season",
      y = "Relative decrease vs baseline",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Relative change - NO VACC
  df_rel %>%
    filter(scenario_id == "no_vacc") %>%
    filter(pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal","6-11mo_immTotal")) %>%
    ggplot(
      aes(x = season, y = rel_decrease, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap( ~ pop_group, scales = "free_y", ncol = 3, labeller = labeller(pop_group = age_labels)) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Season",
      y = "Relative decrease vs baseline",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Relative change - HIGH VACC
  df_rel %>%
    filter(scenario_id == "high_vacc") %>%
    filter(pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal","6-11mo_immTotal")) %>%
    ggplot(
      aes(x = season, y = rel_decrease, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap( ~ pop_group, scales = "free_y", ncol = 3, labeller = labeller(pop_group = age_labels)) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Season",
      y = "Relative decrease vs baseline",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  # Relative Risk - NO VACC
  df_rel %>%
    filter(scenario_id == "no_vacc") %>%
    filter(pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal","6-11mo_immTotal")) %>%
    ggplot(
      aes(x = season, y = rel_decrease+1, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap( ~ pop_group, scales = "free_y", ncol = 3, labeller = labeller(pop_group = age_labels)) +
    scale_y_continuous() +
    labs(
      x = "Season",
      y = "Relative Risk [RR]",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Relative Risk - HIGH VACC
  df_rel %>%
    filter(scenario_id == "high_vacc") %>%
    filter(pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal","6-11mo_immTotal")) %>%
    ggplot(
      aes(x = season, y = rel_decrease+1, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap( ~ pop_group, scales = "free_y", ncol = 3, labeller = labeller(pop_group = age_labels)) +
    scale_y_continuous() +
    labs(
      x = "Season",
      y = "Relative Risk [RR]",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  
  #
  df_rel %>%
    filter(pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal","6-11mo_immTotal")) %>%
    ggplot(
      aes(x = season, y = rel_decrease + 1, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap(scenario_id ~ pop_group, scales = "free_y", ncol = 3) +
    #scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Season",
      y = "Relative Risk [RR]",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  # Plot RR with comparative being 2023/2024 season (of the same scenario)
  df_rel_2324 <- df_season %>%
    group_by(
      Model, scenario_id, pop_group, location, target, output_type_id
    ) %>%
    mutate(
      hosps_2324 = total_hosps[season == "2023/2024"][1],
      rr_vs_2324 = total_hosps / hosps_2324
    ) %>%
    ungroup() %>%
    filter(!is.na(hosps_2324))
  
  df_rel_2324 %>%
    filter(pop_group %in% c("0-2mo_immTotal", "3-5mo_immTotal", "6-11mo_immTotal"),
           scenario_id == "baseline") %>%
    ggplot(
      aes(x = season, y = rr_vs_2324, fill = Model)
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    facet_wrap(scenario_id ~ pop_group, scales = "free_y", ncol = 3) +
    labs(
      x = "Season",
      y = "Relative risk vs 2023/2024 season",
      fill = "Model"
    ) +
    ylim(c(0,NA)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
} # end static model comparison section
