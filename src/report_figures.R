########################################################## #
# REPORT FIGURES
#
# Reproduces the figures comparing the dynamic (ODE) model in
# this repository against the companion static RSV model and
# the observed Irish hospitalisation data.
#
# Produces, in output/<branch>/3_results/<analysis>/report_figures/ :
#
#   7.1  Static and dynamic model fit to weekly hospitalisations
#   7.2  Static model, age-specific 28-day burden vs observed
#   7.3  Dynamic model, age-specific 28-day burden vs observed
#   7.4  Seasonal burden by model and age group, vs observed
#   7.5  Seasonal burden under the no vaccination scenario
#   7.6  Relative change under no vaccination, infants
#   7.7  Relative risk under no vaccination, infants
#   7.8  Seasonal burden under the high uptake scenario
#   7.9  Relative change under high uptake, infants
#   7.10 Relative risk under high uptake, infants
#   7.11 Seasonal burden under the catch-up scenarios
#   7.12 Relative change under the catch-up scenarios, infants
#   7.13 Relative risk under the catch-up scenarios, infants
#
# Usage:
#  Via RStudio: Open RSVinea.Rproj and source this file (Ctrl+Shift+S)
#
# PREREQUISITES
#  1. run_scenarios() must have completed for this analysis, so
#     that 2_scenarios/<analysis>/scenarios/*_raw.rds exist.
#  2. The static model output must be placed at the path in
#     STATIC_FILE below (see o$pth$output for the directory).
#     The static comparison is skipped if the file is absent.
#
# NOTE ON UNCERTAINTY
#  Points are the MEDIAN; bars are the 5th-95th percentile,
#  taken across the static model's samples and across the
#  dynamic model's parameter sets respectively.
#
#  Relative change and relative risk are computed PAIRED - each
#  scenario draw against its own baseline draw - then summarised.
#  Season 1 therefore collapses to exactly 0% / RR = 1, as no
#  vaccination programme is configured in that season.
########################################################## #

# Clear global environment
rm(list = ls())

# Load all required packages and functions
source("R/dependencies.R")

# Tidy up
if (interactive()) clf()  # Close figures
if (interactive()) clc()  # Clear console

# Set options
o = set_options(do_step = 1, analysis_name = "IE")

message("Producing report figures (", o$analysis_name, ")")

# ---- Settings ---------------------------------------------------------- ----

# Static model output (RespiCompass submission format). Placed by the user.
STATIC_FILE = file.path(o$pth$output, "2025_2026_1_RSV_IRL_staticModel.parquet")

# Scenario used as the 'high vaccination uptake' comparison
HIGH_UPTAKE = "newborn_high"

# Seasons to report on (season s runs Aug of year s to Jul of s+1)
SEASONS = 2023:2025

# Catch-up scenarios (figures 7.11 to 7.13). The campaign fires once, on
# vaccination_catch_up_date, so only that season differs between them -
# every earlier season is identical to baseline and carries no information.
CATCHUP_SCN    = c("catch_up_40", "catch_up_60", "catch_up_80")
CATCHUP_SEASON = 2025

# Output directory for these figures
FIGDIR = file.path(o$pth$results, "report_figures")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Plotted values ----------------------------------------------------- ----
# The values drawn in figures 7.4 to 7.13 are collected here, so those figures
# can be reproduced as a table. One CSV per figure, plus a combined
# long-format file. The time-series figures (7.1 to 7.3) are not tabulated.
# `median` is the plotted point; `lo` / `hi` are the whisker ends (the 5th and
# 95th percentile). Observed data has no uncertainty, so lo / hi are NA.
PLOTVALS = list()

collect = function(figure, quantity, df) {
  cols = c("model","scenario","age_band","season","date","median","lo","hi")
  for (nm in cols) if (!nm %in% names(df)) df[[nm]] = NA
  PLOTVALS[[length(PLOTVALS) + 1]] <<-
    df %>% mutate(figure = figure, quantity = quantity) %>%
    select(all_of(c("figure","quantity", cols))) %>%
    mutate(date = as.Date(date))
}

# ---- Age bands and appearance ------------------------------------------ ----

# Model reporting bands, and the labels used in the report
BANDS = c("0-3m", "3-6m", "6-12m", "1-5y", "5-65y", "65+y")
LAB   = c("0-3m" = "<3 months",  "3-6m"  = "3-5 months", "6-12m" = "6-11 months",
          "1-5y" = "1-4 years",  "5-65y" = "5-64 years", "65+y"  = "65+ years")

# Infant bands only - used for the relative change / relative risk figures
INFANT = c("0-3m", "3-6m", "6-12m")

# Static model age labels -> the model's reporting bands
STATIC_BAND = c("0-2mo" = "0-3m", "3-5mo" = "3-6m", "6-11mo" = "6-12m",
                "1-4y"  = "1-5y", "5-64y" = "5-65y", "65+y"   = "65+y")

COL_DYN = "black"; COL_STA = "#9ECAE1"; COL_OBS = "#8B1A1A"

seas = function(d) ifelse(month(d) >= 8, year(d), year(d) - 1)
slab = function(s) paste0(s, "/", s + 1)

thm = theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x      = element_text(angle = 45, hjust = 1))

# Median and 5th-95th percentile over whatever is not in `keyc`
qsum = function(df, keyc)
  df %>% group_by(across(all_of(keyc))) %>%
  summarise(med = median(value),
            lo  = quantile(value, o$quantiles[1]),
            hi  = quantile(value, o$quantiles[2]), .groups = "drop")

# ---- Load observed data ------------------------------------------------- ----

message(" - Loading observed data")

fit = setup_calibration(o)   # See calibration.R
fit = load_data(o, fit)      # See load_data.R

p          = parse_yaml(o, "baseline")$parsed
start_date = p$calibration_options$data_start
band_map   = unlist(p$age_group_map)

# Weekly total admissions, and the 28-day age-specific burden
obs_wk = fit$data %>%
  filter(metric == "hospital_admissions", age_group == "total") %>%
  select(date, value) %>% arrange(date)

obs_4w = fit$data %>%
  filter(metric == "hospital_admissions", data_freq == "4-weekly") %>%
  select(date, band = age_group, value)

obs_seas = obs_4w %>%
  mutate(season = seas(date)) %>% filter(season %in% SEASONS) %>%
  group_by(band, season) %>% summarise(value = sum(value), .groups = "drop")

# ---- Load static model -------------------------------------------------- ----

if (!file.exists(STATIC_FILE))
  stop("Static model output not found at:\n  ", STATIC_FILE,
       "\nPlace the file there, or update STATIC_FILE at the top of this script.")

message(" - Loading static model output")

# `_immTotal` combines immunised and non-immunised; the static model also
# reports each separately, and an administered_doses target the ODE model
# does not currently produce.
st_raw = nanoparquet::read_parquet(STATIC_FILE) %>%
  filter(target == "rsv_hospitalisations", grepl("_immTotal$", pop_group)) %>%
  mutate(age  = str_remove(pop_group, "_immTotal$"),
         band = ifelse(age == "total", "total", unname(STATIC_BAND[age])),
         sid  = output_type_id) %>%
  filter(!is.na(band)) %>%
  select(scenario = scenario_id, date = target_end_date, band, sid, value)

# ---- Load dynamic model ------------------------------------------------- ----

message(" - Loading dynamic model scenario output")

scen_files = list.files(o$pth$scenarios, pattern = "_raw\\.rds$", full.names = TRUE)
if (length(scen_files) == 0)
  stop("No scenario output found in:\n  ", o$pth$scenarios,
       "\nRun run_scenarios(o) first.")

dyn_daily = bind_rows(lapply(scen_files, function(f)
  readRDS(f) %>%
    filter(metric == "hospital_admissions") %>%
    mutate(band = unname(band_map[as.character(age_group)]),
           date = start_date + time - 1) %>%
    filter(!is.na(band)) %>%
    group_by(scenario, date, band, param_id) %>%
    summarise(v = sum(value, na.rm = TRUE), .groups = "drop")))

# Daily model output -> the observed week-end grid. param_id is
# s<sample>_<uncert>_<scenario>; stripping the scenario gives the underlying
# parameter set, so each scenario can be paired against its own baseline.
wk_dates = obs_wk$date
dyn_wk = dyn_daily %>%
  group_by(scenario, band, param_id) %>%
  group_modify(~{
    dn = as.numeric(.x$date); dv = .x$v
    tibble(date  = wk_dates,
           value = vapply(as.numeric(wk_dates),
                          function(x) sum(dv[dn > x - 7 & dn <= x]), numeric(1)))
  }) %>% ungroup() %>%
  mutate(sid = str_remove(param_id, paste0("_", scenario, "$"))) %>%
  select(-param_id)

# Weekly -> the observed 28-day periods. A burden value dated d covers
# [d, d + 28), which is how load_data() builds the observed burden.
to_4wk = function(df, keyc) {
  per = sort(unique(obs_4w$date))
  df %>% group_by(across(all_of(keyc))) %>%
    group_modify(~{
      dn = as.numeric(.x$date); dv = .x$value
      tibble(date  = per,
             value = vapply(as.numeric(per),
                            function(x) sum(dv[dn >= x & dn < x + 28]), numeric(1)))
    }) %>% ungroup()
}

st_4w  = to_4wk(st_raw %>% filter(scenario == "baseline", band != "total"), c("band","sid"))
dyn_4w = to_4wk(dyn_wk %>% filter(scenario == "baseline"),                  c("band","sid"))

# ---- Figure 7.1: model fit to weekly hospitalisations ------------------- ----

message(" - Figure 7.1")

f71 = bind_rows(
  st_raw %>% filter(scenario == "baseline", band == "total") %>%
    qsum("date") %>% mutate(model = "Static model"),
  dyn_wk %>% filter(scenario == "baseline") %>%
    group_by(date, sid) %>% summarise(value = sum(value), .groups = "drop") %>%
    qsum("date") %>% mutate(model = "Dynamic model")) %>%
  mutate(model = factor(model, levels = c("Static model", "Dynamic model")))

g = ggplot(f71, aes(date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "Model variation"), alpha = 0.9) +
  geom_line(aes(y = med, colour = "Model median"), linewidth = 0.6) +
  geom_point(data = obs_wk, aes(date, value, colour = "Observed data"), size = 1.3) +
  facet_wrap(~model, nrow = 1) +
  scale_fill_manual(NULL, values = c("Model variation" = "#C6DBEF")) +
  scale_colour_manual(NULL, values = c("Model median"  = "#2171B5",
                                       "Observed data" = COL_OBS)) +
  labs(x = "Date", y = "Weekly hospitalisations") + thm
ggsave(file.path(FIGDIR, "fig_7_1_model_fit.png"), g,
       width = 26, height = 11, units = "cm", dpi = o$save_resolution)

# ---- Figures 7.2 and 7.3: age-specific 28-day burden -------------------- ----

# open_circle: the static model reproduces the data by design, so its points
# are drawn unfilled to make the overlap with the observed points visible.
age_panel = function(mod, open_circle, file) {
  d  = qsum(mod, c("band","date")) %>% filter(band %in% BANDS) %>%
       mutate(band = factor(LAB[band], levels = LAB[BANDS]))
  ob = obs_4w %>% filter(band %in% BANDS) %>%
       mutate(band = factor(LAB[band], levels = LAB[BANDS]))
  g = ggplot(d, aes(date, med)) +
    { if (open_circle)
        geom_point(shape = 21, fill = NA, colour = COL_DYN, size = 1.9)
      else
        geom_pointrange(aes(ymin = lo, ymax = hi), colour = COL_DYN,
                        size = 0.25, fatten = 2) } +
    geom_point(data = ob, aes(date, value), colour = COL_OBS, size = 1.6) +
    facet_wrap(~band, scales = "free_y", nrow = 2) +
    labs(x = "Time", y = "4-weekly burden") + thm
  ggsave(file, g, width = 26, height = 15, units = "cm", dpi = o$save_resolution)
}

message(" - Figures 7.2 and 7.3")
age_panel(st_4w,  TRUE,  file.path(FIGDIR, "fig_7_2_static_age_fit.png"))
age_panel(dyn_4w, FALSE, file.path(FIGDIR, "fig_7_3_dynamic_age_fit.png"))

# ---- Figures 7.4, 7.5 and 7.8: seasonal burden -------------------------- ----

seasonal = function(scn) {
  bind_rows(
    st_raw %>% filter(scenario == scn, band != "total") %>%
      mutate(season = seas(date)) %>% filter(season %in% SEASONS) %>%
      group_by(band, season, sid) %>% summarise(value = sum(value), .groups = "drop") %>%
      qsum(c("band","season")) %>% mutate(model = "Static"),
    dyn_wk %>% filter(scenario == scn) %>%
      mutate(season = seas(date)) %>% filter(season %in% SEASONS) %>%
      group_by(band, season, sid) %>% summarise(value = sum(value), .groups = "drop") %>%
      qsum(c("band","season")) %>% mutate(model = "Dynamic")) %>%
    filter(band %in% BANDS) %>%
    mutate(band = factor(LAB[band], levels = LAB[BANDS]), season_lab = slab(season))
}

seas_plot = function(df, obs_layer, file, figure, scn) {
  g = ggplot(df, aes(season_lab, med, colour = model)) +
    geom_pointrange(aes(ymin = lo, ymax = hi),
                    position = position_dodge(width = 0.5), size = 0.35, fatten = 2.2)
  if (!is.null(obs_layer))
    g = g + geom_point(data = obs_layer, aes(season_lab, value, colour = "Observed"),
                       size = 2, inherit.aes = FALSE)
  g = g + facet_wrap(~band, scales = "free_y", nrow = 2) +
    scale_colour_manual(NULL, values = c("Dynamic"  = COL_DYN, "Static" = COL_STA,
                                         "Observed" = COL_OBS)) +
    expand_limits(y = 0) + labs(x = "Season", y = "Seasonal burden") + thm
  ggsave(file, g, width = 26, height = 15, units = "cm", dpi = o$save_resolution)

  collect(figure, "Seasonal burden",
          df %>% transmute(model, scenario = scn, age_band = as.character(band),
                           season = season_lab, median = med, lo, hi))
  if (!is.null(obs_layer))
    collect(figure, "Seasonal burden",
            obs_layer %>% transmute(model = "Observed", scenario = scn,
                                    age_band = as.character(band),
                                    season = season_lab, median = value))
}

message(" - Figures 7.4, 7.5 and 7.8")
obs_layer = obs_seas %>% filter(band %in% BANDS) %>%
  mutate(band = factor(LAB[band], levels = LAB[BANDS]), season_lab = slab(season))

seas_plot(seasonal("baseline"),  obs_layer,
          file.path(FIGDIR, "fig_7_4_seasonal_baseline.png"),    "7.4", "baseline")
seas_plot(seasonal("no_vacc"),   NULL,
          file.path(FIGDIR, "fig_7_5_seasonal_no_vacc.png"),     "7.5", "no_vacc")
seas_plot(seasonal(HIGH_UPTAKE), NULL,
          file.path(FIGDIR, "fig_7_8_seasonal_high_uptake.png"), "7.8", HIGH_UPTAKE)

# ---- Figures 7.6, 7.7, 7.9 and 7.10: relative change and relative risk --- ----

# Paired within draw: each scenario draw against its own baseline draw. This
# is much tighter than differencing two summarised means, because the shared
# parameter uncertainty cancels.
rel = function(scn) {
  paired = function(df, model_name)
    df %>% filter(scenario %in% c("baseline", scn), band %in% INFANT) %>%
      mutate(season = seas(date)) %>% filter(season %in% SEASONS) %>%
      group_by(scenario, band, season, sid) %>%
      summarise(v = sum(value), .groups = "drop") %>%
      pivot_wider(names_from = scenario, values_from = v) %>%
      mutate(rr = .data[[scn]] / baseline, model = model_name)

  bind_rows(paired(st_raw, "Static"), paired(dyn_wk, "Dynamic")) %>%
    group_by(model, band, season) %>%
    summarise(rr_med = median(rr),
              rr_lo  = quantile(rr, o$quantiles[1]),
              rr_hi  = quantile(rr, o$quantiles[2]),
              pc_med = median(100 * (rr - 1)),
              pc_lo  = quantile(100 * (rr - 1), o$quantiles[1]),
              pc_hi  = quantile(100 * (rr - 1), o$quantiles[2]), .groups = "drop") %>%
    mutate(band = factor(LAB[band], levels = LAB[INFANT]), season_lab = slab(season))
}

rel_plot = function(df, ycol, ylab, file, figure, scn) {
  g = ggplot(df, aes(season_lab, .data[[paste0(ycol, "_med")]], colour = model)) +
    # Null value: no change is 0% on the relative-change scale, but RR = 1
    geom_hline(yintercept = if (ycol == "rr") 1 else 0,
               linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_pointrange(aes(ymin = .data[[paste0(ycol, "_lo")]],
                        ymax = .data[[paste0(ycol, "_hi")]]),
                    position = position_dodge(width = 0.5), size = 0.35, fatten = 2.2) +
    facet_wrap(~band, scales = "free_y", nrow = 1) +
    scale_colour_manual(NULL, values = c("Dynamic" = COL_DYN, "Static" = COL_STA)) +
    labs(x = "Season", y = ylab) + thm
  if (ycol == "rr") g = g + expand_limits(y = 0)
  ggsave(file, g, width = 24, height = 11, units = "cm", dpi = o$save_resolution)

  collect(figure, ylab,
          df %>% transmute(model, scenario = paste0(scn, " vs baseline"),
                           age_band = as.character(band), season = season_lab,
                           median = .data[[paste0(ycol, "_med")]],
                           lo     = .data[[paste0(ycol, "_lo")]],
                           hi     = .data[[paste0(ycol, "_hi")]]))
}

message(" - Figures 7.6, 7.7, 7.9 and 7.10")
r_nv = rel("no_vacc")
r_hi = rel(HIGH_UPTAKE)

rel_plot(r_nv, "pc", "Relative change (%)",
         file.path(FIGDIR, "fig_7_6_relchange_no_vacc.png"),      "7.6",  "no_vacc")
rel_plot(r_nv, "rr", "Relative Risk",
         file.path(FIGDIR, "fig_7_7_relrisk_no_vacc.png"),        "7.7",  "no_vacc")
rel_plot(r_hi, "pc", "Relative change (%)",
         file.path(FIGDIR, "fig_7_9_relchange_high_uptake.png"),  "7.9",  HIGH_UPTAKE)
rel_plot(r_hi, "rr", "Relative Risk",
         file.path(FIGDIR, "fig_7_10_relrisk_high_uptake.png"),   "7.10", HIGH_UPTAKE)

# ---- Figures 7.11 to 7.13: the catch-up scenarios ----------------------- ----

message(" - Figures 7.11 to 7.13")

# Order the scenarios by their catch-up coverage, and put the coverage in the
# tick label. Read from the yaml so the axis cannot drift from the config.
# Axis labels are the scenario `name` from the yaml (e.g. "Catch up scenario
# 80%"), so the figure reads the same as the configuration. Coverage is read
# alongside it purely to order the axis left to right.
cu_meta = lapply(setNames(nm = c("baseline", CATCHUP_SCN)), function(x) {
  y  = parse_yaml(o, x)$parsed
  nm = y$.name
  list(cov  = y$vaccination_catch_up_coverage,
       name = if (is.null(nm)) x else nm)   # fall back to the id if unnamed
})
cu_cov = vapply(cu_meta, function(z) z$cov, numeric(1))
cu_ord = setdiff(names(sort(cu_cov)), "baseline")
cu_lab = setNames(vapply(cu_meta[cu_ord], function(z) z$name, character(1)), cu_ord)

cu_dat = bind_rows(
  st_raw %>% filter(scenario %in% c("baseline", CATCHUP_SCN), band != "total") %>%
    mutate(season = seas(date)) %>% filter(season == CATCHUP_SEASON) %>%
    group_by(model = "Static", scenario, band, sid) %>%
    summarise(v = sum(value), .groups = "drop"),
  dyn_wk %>% filter(scenario %in% c("baseline", CATCHUP_SCN)) %>%
    mutate(season = seas(date)) %>% filter(season == CATCHUP_SEASON) %>%
    group_by(model = "Dynamic", scenario, band, sid) %>%
    summarise(v = sum(value), .groups = "drop"))

cu_season_lab = slab(CATCHUP_SEASON)

# --- 7.11 seasonal burden, all age bands (baseline is shown in 7.4) ---
cu_burden = cu_dat %>% filter(scenario != "baseline", band %in% BANDS) %>%
  group_by(model, band, scenario) %>%
  summarise(med = median(v), lo = quantile(v, o$quantiles[1]),
            hi = quantile(v, o$quantiles[2]), .groups = "drop") %>%
  mutate(band = factor(LAB[band], levels = LAB[BANDS]),
         scenario = factor(scenario, levels = cu_ord))

g = ggplot(cu_burden, aes(scenario, med, colour = model)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.5), size = 0.35, fatten = 2.2) +
  facet_wrap(~band, scales = "free_y", nrow = 2) +
  scale_x_discrete(labels = cu_lab) +
  scale_colour_manual(NULL, values = c("Dynamic" = COL_DYN, "Static" = COL_STA)) +
  expand_limits(y = 0) +
  labs(x = "Scenario", y = "Seasonal burden") + thm
ggsave(file.path(FIGDIR, "fig_7_11_catchup_seasonal_burden.png"), g,
       width = 26, height = 16, units = "cm", dpi = o$save_resolution)

collect("7.11", "Seasonal burden",
        cu_burden %>% transmute(model, scenario = as.character(scenario),
                                age_band = as.character(band), season = cu_season_lab,
                                median = med, lo, hi))

# --- 7.12 / 7.13 relative change and relative risk, infants ---
# Paired: each scenario draw against its own baseline draw.
cu_rel = cu_dat %>% filter(band %in% INFANT) %>%
  pivot_wider(names_from = scenario, values_from = v) %>%
  pivot_longer(all_of(cu_ord), names_to = "scenario", values_to = "v") %>%
  mutate(rr = v / baseline) %>%
  group_by(model, band, scenario) %>%
  summarise(rr_med = median(rr), rr_lo = quantile(rr, o$quantiles[1]),
            rr_hi = quantile(rr, o$quantiles[2]),
            pc_med = median(100 * (rr - 1)),
            pc_lo  = quantile(100 * (rr - 1), o$quantiles[1]),
            pc_hi  = quantile(100 * (rr - 1), o$quantiles[2]), .groups = "drop") %>%
  mutate(band = factor(LAB[band], levels = LAB[INFANT]),
         scenario = factor(scenario, levels = cu_ord))

cu_plot = function(ycol, ylab, href, file, figure) {
  g = ggplot(cu_rel, aes(scenario, .data[[paste0(ycol, "_med")]], colour = model)) +
    geom_hline(yintercept = href, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_pointrange(aes(ymin = .data[[paste0(ycol, "_lo")]],
                        ymax = .data[[paste0(ycol, "_hi")]]),
                    position = position_dodge(width = 0.5), size = 0.35, fatten = 2.2) +
    facet_wrap(~band, scales = "free_y", nrow = 1) +
    scale_x_discrete(labels = cu_lab) +
    scale_colour_manual(NULL, values = c("Dynamic" = COL_DYN, "Static" = COL_STA)) +
    labs(x = "Scenario", y = ylab) + thm
  ggsave(file.path(FIGDIR, file), g, width = 24, height = 12,
         units = "cm", dpi = o$save_resolution)

  collect(figure, ylab,
          cu_rel %>% transmute(model,
                               scenario = paste0(as.character(scenario), " vs baseline"),
                               age_band = as.character(band), season = cu_season_lab,
                               median = .data[[paste0(ycol, "_med")]],
                               lo     = .data[[paste0(ycol, "_lo")]],
                               hi     = .data[[paste0(ycol, "_hi")]]))
}

cu_plot("pc", "Relative change (%)", 0,
        "fig_7_12_catchup_relchange.png", "7.12")
cu_plot("rr", "Relative Risk",       1,
        "fig_7_13_catchup_relrisk.png",  "7.13")

# ---- Summary tables ----------------------------------------------------- ----

show_rel = function(df, title) {
  cat("\n=== ", title, " (relative change %, median [5-95%]) ===\n", sep = "")
  print(as.data.frame(
    df %>% select(model, band, season_lab, pc_med, pc_lo, pc_hi) %>%
      mutate(across(starts_with("pc"), ~round(.x, 1))) %>%
      arrange(band, season_lab, model)), row.names = FALSE)
}
show_rel(r_nv, "No vaccination vs baseline")
show_rel(r_hi, paste0(HIGH_UPTAKE, " vs baseline"))

cat("\n=== Catch-up scenarios vs baseline, ", cu_season_lab,
    " (relative change %, median [5-95%]) ===\n", sep = "")
print(as.data.frame(
  cu_rel %>% select(model, band, scenario, pc_med, pc_lo, pc_hi) %>%
    mutate(across(starts_with("pc"), ~round(.x, 1))) %>%
    arrange(band, scenario, model)), row.names = FALSE)

# ---- Write the plotted values ------------------------------------------- ----

vals = bind_rows(PLOTVALS) %>%
  mutate(across(c(median, lo, hi), ~round(.x, 4))) %>%
  select(-date) %>%                     # seasonal figures only; no time axis
  arrange(figure, age_band, season, model)

# Combined long-format table, and one file per figure
write.csv(vals, file.path(FIGDIR, "figure_values_all.csv"), row.names = FALSE, na = "")

# Figures 7.4 to 7.10 only - the weekly / 28-day time series in 7.1 to 7.3
# are not tabulated.
fig_file = c("7.4" = "fig_7_4_seasonal_baseline",
             "7.5" = "fig_7_5_seasonal_no_vacc",
             "7.6" = "fig_7_6_relchange_no_vacc",
             "7.7" = "fig_7_7_relrisk_no_vacc",
             "7.8" = "fig_7_8_seasonal_high_uptake",
             "7.9" = "fig_7_9_relchange_high_uptake",
             "7.10" = "fig_7_10_relrisk_high_uptake",
             "7.11" = "fig_7_11_catchup_seasonal_burden",
             "7.12" = "fig_7_12_catchup_relchange",
             "7.13" = "fig_7_13_catchup_relrisk")
for (fg in names(fig_file)) {
  d = vals %>% filter(figure == fg) %>%
    select(where(~!all(is.na(.x))))   # drop columns that do not apply to this figure
  write.csv(d, file.path(FIGDIR, paste0(fig_file[[fg]], ".csv")), row.names = FALSE, na = "")
}

cat("\n=== rows written per figure ===\n")
print(as.data.frame(vals %>% count(figure, quantity, name = "rows")), row.names = FALSE)

# ---- Report-ready wide tables ------------------------------------------- ----
# One row per age band and season, one column per model, each cell formatted as
# "median (lo-hi)". This is the shape to paste into a document; the long file
# above is the one to compute from.

# Decimals differ by quantity: counts to the unit, percentages to 0.1, RR to 0.01
digits_for = function(q) ifelse(q == "Seasonal burden", 0,
                        ifelse(q == "Relative Risk", 2, 1))

# formatC() takes a scalar `digits`, so format element-wise: the number of
# decimals varies by row because the quantities do.
fmt_num = function(x, d) mapply(function(xi, di)
  if (is.na(xi)) NA_character_ else formatC(xi, format = "f", digits = di),
  x, d, USE.NAMES = FALSE)

fmt_ci = function(med, lo, hi, d) {
  m = fmt_num(med, d)
  ifelse(is.na(lo) | is.na(hi), m,
         paste0(m, " (", fmt_num(lo, d), "-", fmt_num(hi, d), ")"))
}

wide = vals %>%
  mutate(cell = fmt_ci(median, lo, hi, digits_for(quantity))) %>%
  select(figure, quantity, scenario, age_band, season, model, cell) %>%
  pivot_wider(names_from = model, values_from = cell) %>%
  # keep the figure order, and the age bands in reporting order
  mutate(fig_num = as.numeric(figure),
         age_band = factor(age_band, levels = LAB)) %>%
  arrange(fig_num, age_band, season) %>%
  select(-fig_num) %>%
  select(figure, quantity, scenario, age_band, season,
         any_of(c("Dynamic", "Static", "Observed")))

write.csv(wide, file.path(FIGDIR, "figure_tables_wide.csv"), row.names = FALSE, na = "")
cat("\n=== report-ready wide table (first rows) ===\n")
print(as.data.frame(head(wide, 4)), row.names = FALSE)

message("* Figures and value tables written to ", FIGDIR)
