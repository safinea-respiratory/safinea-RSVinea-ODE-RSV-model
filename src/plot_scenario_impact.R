########################################################## #
# SCENARIO IMPACT ACROSS COUNTRIES
#
# Impact of the adult RSV vaccination scenarios on hospitalisation burden,
# compared against each country's own baseline, for every country that has
# scenario output. Produces eight figures:
#
#   1  pct_averted_scenario_ages   % averted, scenario's own eligible ages
#   2  pct_averted_union_ages      % averted, union of eligible ages per country
#   3  abs_averted                 absolute averted hospitalisations (count)
#   4a doses_per100k_total         doses per 100k total population
#   4b doses_per100k_eligible      doses per 100k eligible population
#   5  averted_per100k_total       averted hosps per 100k total population
#   6  averted_per100k_eligible    averted hosps per 100k eligible (scenario)
#   7  averted_per100k_union       averted hosps per 100k eligible (union)
#
# All: x = country, median + 90% interval ACROSS PAIRED SAMPLES, facet =
# scenario, colour = season (except the dose plots, which have no season
# dimension because there is a single campaign date).
#
# Countries with no scenario output, or with baseline only, are SKIPPED with a
# message rather than erroring, so this can be run mid-sweep.
#
# Reads output/2_scenarios/<ISO>/scenarios/<scenario>_raw.rds, written by
# run_scenarios() (do_step 2). Standalone; does not depend on
# results_evaluation.R.
########################################################## #

rm(list = ls())
source("R/dependencies.R")

# ---- Settings ----------------------------------------------------------------
OUT_DIR  <- "output/3_results"
PREFIX   <- "scenario_impact"
BASELINE <- "baseline"
METRIC   <- "hospital_admissions"   # the outcome averted
DOSES    <- "n_doses"
QUANT    <- c(0.05, 0.95)

# Seasons are bounded on 1 August (the model's own convention). A 730-day run
# from 2026-09-01 ends 2028-08-31, spilling 31 days into 2028/2029 - a sliver
# outside the RSV season whose burden is ~0, where a relative change is
# meaningless. Drop anything shorter, and report what was dropped.
MIN_SEASON_DAYS <- 180

scen_root <- "output/2_scenarios"

# ---- Helpers -----------------------------------------------------------------
season_label <- function(d) {
  y <- ifelse(month(d) >= 8, year(d), year(d) - 1)
  paste0(y, "/", y + 1)
}

scen_files <- function(iso) {
  d <- file.path(scen_root, iso, "scenarios")
  if (!dir.exists(d)) return(character(0))
  f <- list.files(d, pattern = "_raw[.]rds$")
  setNames(file.path(d, f), sub("_raw[.]rds$", "", f))
}

# Population by RespiCompass reporting band -> the model's fine age-group labels
# used in the scenario output. Needed for the per-100k denominators.
pop_by_band <- function(o, iso) {
  p <- read.csv(o$pop_url, fileEncoding = "UTF-8-BOM") %>% normalise_iso2()
  p <- p[p$country == iso, ]
  if (!nrow(p)) return(NULL)
  setNames(p$population, remap_lookup(o)[p$age_group])
}
# respicompass_age_map gives data-label -> reporting band; the eligible age
# groups in the yaml are FINE model labels (65-69y ...), which coincide with the
# reporting bands for adults, so map through and match on those.
remap_lookup <- function(o) {
  m <- parse_yaml(o, "baseline")$parsed$respicompass_age_map
  setNames(unlist(m), names(m))
}

# ---- Collect every quantity for one country ----------------------------------
collect_one <- function(iso) {

  files <- scen_files(iso)
  if (!length(files))                { message("  - ", iso, ": no scenario output, skipped"); return(NULL) }
  if (!BASELINE %in% names(files))   { message("  - ", iso, ": no baseline, skipped");        return(NULL) }
  others <- setdiff(names(files), BASELINE)
  if (!length(others))               { message("  - ", iso, ": baseline only, skipped");      return(NULL) }

  o     <- set_options(do_step = 2, analysis_name = iso, quiet = TRUE)
  fit   <- load_data(o, setup_calibration(o))
  start <- min(fit$dates_model$date)

  # Eligible age groups PER SCENARIO. Scenarios can override any yaml field, so
  # a scenario that widened eligibility would otherwise be summed over the
  # baseline's groups and its effect understated with no error.
  elig <- lapply(names(files), function(s)
    unlist(parse_yaml(o, s)$parsed$adult_vaccination_agegroups))
  names(elig) <- names(files)
  # Union across the scenarios this country actually has, so plots 2 and 7
  # compare every scenario on the SAME population.
  elig_union <- sort(unique(unlist(elig[others])))
  if (!length(elig_union)) { message("  - ", iso, ": no eligible age groups, skipped"); return(NULL) }

  # Population denominators
  popv <- pop_by_band(o, iso)
  if (is.null(popv)) { message("  - ", iso, ": no population data, skipped"); return(NULL) }
  pop_total <- sum(popv, na.rm = TRUE)

  # adult_vaccination_agegroups are FINE model labels ("65-69y"), but popv is
  # keyed by REPORTING bands ("65-70y"). Map through age_group_map before
  # looking population up. Without this only "80+y" matches (the one label that
  # happens to be identical in both schemes) and the eligible denominator
  # collapses to the 80+ population alone, inflating every per-eligible figure
  # by ~3.3x - which is exactly what the doses-per-100k-eligible plot showed.
  agm <- unlist(parse_yaml(o, "baseline")$parsed$age_group_map)
  pop_of <- function(ages) {
    bands <- unique(agm[ages])
    if (anyNA(bands))
      warning(iso, ": age group(s) ", paste(ages[is.na(agm[ages])], collapse = ", "),
              " not found in age_group_map - excluded from the population denominator")
    bands <- bands[!is.na(bands)]
    if (!length(bands)) return(NA_real_)
    sum(popv[bands], na.rm = TRUE)
  }

  # Seasons the run actually covers
  days_ref <- data.table(season = season_label(start + 0:730))[, .N, by = season]
  keep_seasons <- days_ref[N >= MIN_SEASON_DAYS, season]
  dropped <- days_ref[N < MIN_SEASON_DAYS]
  if (nrow(dropped))
    message("    (", iso, ": dropping partial season(s) ",
            paste(sprintf("%s = %dd", dropped$season, dropped$N), collapse = ", "), ")")

  # ---- burden per (sample, season) for a given scenario and age set ----------
  # hospital_admissions is an incidence FLOW, so summing over time within a
  # season is the seasonal total. Grouping by `sim` keeps samples separate -
  # they are uncertainty and must never be summed together.
  burden_of <- function(scen, ages) {
    r <- as.data.table(readRDS(files[[scen]]))
    r <- r[metric == METRIC & age_group %in% ages]
    if (!nrow(r)) return(NULL)
    r[, sim    := sub(paste0("_", scen, "$"), "", param_id)]
    r[, season := season_label(start + time - 1L)]
    r <- r[season %in% keep_seasons]
    r[, .(burden = sum(value, na.rm = TRUE)), by = .(sim, season)]
  }

  # ---- doses: CUMULATIVE counter, so take the FINAL value, never the sum -----
  doses_of <- function(scen, ages) {
    r <- as.data.table(readRDS(files[[scen]]))
    r <- r[metric == DOSES & age_group %in% ages]
    if (!nrow(r)) return(NULL)
    r[, sim := sub(paste0("_", scen, "$"), "", param_id)]
    r[time == max(time), .(doses = sum(value, na.rm = TRUE)), by = sim]
  }

  out <- list()
  for (scen in others) {

    ages_s <- elig[[scen]]
    if (!length(ages_s)) next

    # --- averted hospitalisations, on two age bases -------------------------
    for (basis in c("scenario", "union")) {
      ages <- if (basis == "scenario") ages_s else elig_union
      b_sc <- burden_of(scen,     ages)
      b_bl <- burden_of(BASELINE, ages)
      if (is.null(b_sc) || is.null(b_bl)) next
      m <- merge(b_sc, b_bl, by = c("sim", "season"), suffixes = c("", "_bl"))
      if (!nrow(m)) next
      # Paired per sample: each scenario draw against its OWN baseline draw, so
      # shared parameter uncertainty cancels instead of inflating the interval.
      m[, `:=`(pct     = 100 * (burden - burden_bl) / burden_bl,
               averted = burden_bl - burden)]
      m[, `:=`(av_per100k_tot  = 1e5 * averted / pop_total,
               av_per100k_elig = 1e5 * averted / pop_of(ages))]
      out[[length(out) + 1]] <- m[, .(
        iso = iso, scen = scen, season, basis, sim,
        pct, averted, av_per100k_tot, av_per100k_elig)]
    }

    # --- doses (no season dimension: a single campaign date) ----------------
    d <- doses_of(scen, ages_s)
    if (!is.null(d) && nrow(d)) {
      # Sanity check on the eligible-population denominator. With a single
      # campaign, doses per eligible person cannot exceed the coverage, and
      # should sit just below it (only susceptibles are vaccinated, so those in
      # E/I/H/R on the campaign date are missed). A value above 1 means the
      # denominator is wrong - which is how the fine-vs-reporting age label
      # mismatch was caught. Report rather than silently plotting nonsense.
      cov  <- parse_yaml(o, scen)$parsed$adult_vacc_coverage
      per  <- median(d$doses) / pop_of(ages_s)
      if (is.finite(per) && is.finite(cov) && cov > 0) {
        if (per > 1)
          warning(iso, "/", scen, ": ", round(per, 2), " doses per eligible person",
                  " (> 1) - eligible population denominator looks wrong")
        else if (per > cov * 1.02)
          warning(iso, "/", scen, ": doses/eligible = ", round(per, 3),
                  " exceeds coverage ", cov, " - check the denominator")
      }
      out[[length(out) + 1]] <- data.table(
        iso = iso, scen = scen, season = NA_character_, basis = "doses", sim = d$sim,
        pct = NA_real_, averted = NA_real_,
        av_per100k_tot  = 1e5 * d$doses / pop_total,
        av_per100k_elig = 1e5 * d$doses / pop_of(ages_s))
    }
  }

  if (!length(out)) { message("  - ", iso, ": nothing usable, skipped"); return(NULL) }
  res <- rbindlist(out, fill = TRUE)
  message("  + ", iso, ": ", paste(others, collapse = ", "),
          " | eligible ", paste(elig_union, collapse = "/"),
          " | n = ", length(unique(res$sim)))
  res
}

# ---- Run over all countries --------------------------------------------------
isos <- sort(list.dirs(scen_root, recursive = FALSE, full.names = FALSE))
message("* Scanning ", length(isos), " countries")
raw <- rbindlist(lapply(isos, function(i)
  tryCatch(collect_one(i),
           error = function(e) { message("  ! ", i, ": ", conditionMessage(e), " - skipped"); NULL })),
  fill = TRUE)

if (is.null(raw) || !nrow(raw)) stop("No country produced usable scenario output.")
message("\n* ", length(unique(raw$iso)), " countries, ",
        length(unique(raw$scen)), " scenarios")

# ---- Summarise across samples ------------------------------------------------
summ <- function(d, col) {
  d[!is.na(get(col)), .(median = median(get(col), na.rm = TRUE),
                        lo     = quantile(get(col), QUANT[1], na.rm = TRUE),
                        hi     = quantile(get(col), QUANT[2], na.rm = TRUE),
                        n_sim  = .N),
    by = .(iso, scen, season)]
}

# ---- Plot builder ------------------------------------------------------------
make_plot <- function(d, title, ylab, by_season = TRUE, log_y = FALSE) {
  d <- copy(d)
  d[, iso := factor(iso, levels = sort(unique(iso)))]
  g <- ggplot(d, aes(x = iso, y = median))
  # A zero reference line is meaningless on a log axis (log10(0) = -Inf) and
  # would drop the panel, so only draw it on linear axes.
  if (!log_y)
    g <- g + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")
  if (by_season) {
    d[, season := factor(season, levels = sort(unique(season)))]
    g <- g + geom_pointrange(aes(ymin = lo, ymax = hi, colour = season, group = season),
                             position = position_dodge(width = 0.6), size = 0.35) +
      scale_colour_viridis_d(name = "Season", option = "D", end = 0.8)
  } else {
    g <- g + geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.35, colour = "grey20")
  }
  if (log_y) {
    # Guard: log10 needs strictly positive values. Anything <= 0 (a scenario
    # that INCREASED burden) cannot be shown on a log axis - drop it loudly
    # rather than letting ggplot silently discard the row.
    bad <- d[!is.finite(median) | median <= 0 | !is.finite(lo) | lo <= 0]
    if (nrow(bad)) {
      warning("log scale: dropping ", nrow(bad), " non-positive point(s): ",
              paste(unique(paste(bad$iso, bad$scen, bad$season)), collapse = "; "))
      d <- d[is.finite(median) & median > 0 & is.finite(lo) & lo > 0]
      g <- g %+% d
    }
    g <- g + scale_y_log10(labels = scales::label_comma()) +
      annotation_logticks(sides = "l", outside = FALSE, alpha = 0.4)
  }

  g + facet_wrap(~ scen, ncol = 1, scales = "free_y") +
    labs(title = title,
         subtitle = paste0("Median and ", 100 * (QUANT[2] - QUANT[1]),
                           "% interval across paired samples",
                           if (log_y) "  -  log10 y axis" else ""),
         x = NULL, y = ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey85"),
          legend.position = "bottom",
          panel.grid.major.x = element_blank())
}

save_plot <- function(g, name, n_facet, tab) {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(OUT_DIR, paste0(PREFIX, "_", name))
  ggsave(paste0(f, ".png"), g, width = 12, height = 3 + 2.4 * n_facet,
         dpi = 200, limitsize = FALSE)
  fwrite(tab[order(scen, season, iso)], paste0(f, ".csv"))
  message("  saved ", basename(f), ".png / .csv")
}

sc <- raw[basis == "scenario"]; un <- raw[basis == "union"]; dz <- raw[basis == "doses"]
nf <- length(unique(raw$scen))

message("\n* Writing figures")
specs <- list(
  list(d = sc, col = "pct",             n = "1_pct_averted_scenario_ages",
       t = "Relative change in seasonal RSV hospitalisations vs baseline (scenario's eligible ages)",
       y = "Relative change (%)  -  negative = averted", s = TRUE),
  list(d = un, col = "pct",             n = "2_pct_averted_union_ages",
       t = "Relative change in seasonal RSV hospitalisations vs baseline (union of eligible ages)",
       y = "Relative change (%)  -  negative = averted", s = TRUE),
  # Log10: absolute counts span ~2 orders of magnitude across countries
  # (CY ~150 to DE ~13,000), so a linear axis is dominated by the largest
  # countries and the smaller ones are unreadable.
  list(d = sc, col = "averted",         n = "3_abs_averted",
       t = "Hospitalisations averted vs baseline (scenario's eligible ages)",
       y = "Hospitalisations averted (count)", s = TRUE, log = TRUE),
  list(d = dz, col = "av_per100k_tot",  n = "4a_doses_per100k_total",
       t = "Vaccine doses administered, per 100k total population",
       y = "Doses per 100k total population", s = FALSE),
  list(d = dz, col = "av_per100k_elig", n = "4b_doses_per100k_eligible",
       t = "Vaccine doses administered, per 100k eligible population",
       y = "Doses per 100k eligible population", s = FALSE),
  list(d = sc, col = "av_per100k_tot",  n = "5_averted_per100k_total",
       t = "Hospitalisations averted per 100k total population",
       y = "Averted per 100k total population", s = TRUE),
  list(d = sc, col = "av_per100k_elig", n = "6_averted_per100k_eligible",
       t = "Hospitalisations averted per 100k eligible population (scenario's eligible ages)",
       y = "Averted per 100k eligible population", s = TRUE),
  list(d = un, col = "av_per100k_elig", n = "7_averted_per100k_union",
       t = "Hospitalisations averted per 100k eligible population (union of eligible ages)",
       y = "Averted per 100k eligible population", s = TRUE)
)

for (sp in specs) {
  if (is.null(sp$d) || !nrow(sp$d)) { message("  ! ", sp$n, ": no data, skipped"); next }
  tab <- summ(sp$d, sp$col)
  if (!nrow(tab)) { message("  ! ", sp$n, ": no data, skipped"); next }
  save_plot(make_plot(tab, sp$t, sp$y, by_season = sp$s,
                      log_y = isTRUE(sp$log)), sp$n, nf, tab)
}

message("\n* Done - ", OUT_DIR)
