########################################################## #
# EVALUATE
#
# Single entry point for post-processing the adult RSV vaccination
# scenarios. Replaces results_evaluation.R and plot_scenario_impact.R.
#
# Does four things, in order:
#
#   1. LOAD    every country's scenario output (this model), plus the
#              companion static model's parquet.
#   2. PLOT    scenario impact, one figure per (season x model) - see
#              "Figure layout" below.
#   3. SHAPE   both models into the project's submission format (see
#              "Submission format" below - this is NOT the RespiCompass
#              wiki schema).
#   4. SAVE    one parquet per model.
#
# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
# The scenario grid has TWO axes - eligibility age threshold (60/65/70/75/80+)
# and uptake (25/50/75%) - on top of country, season and model. That is five
# categorical dimensions and only x, y, colour and facet to spend, so something
# has to move to separate FILES.
#
#   x         country
#   y         the figure's quantity
#   facet     uptake (3 panels: 25 / 50 / 75%)
#   colour    eligibility age threshold (sequential - within each country the
#             five thresholds should fan out in order, so a crossing stands out)
#   file      season x model
#
# Uptake as the facet also makes the uptake x VE reference band cleaner: each
# panel has a single uptake and therefore a single band, rather than three
# translucent ones overlapping.
#
# Splitting model across files (rather than shape, as the old script did) keeps
# each panel to five points per country. The cost is that the model comparison
# is no longer side-by-side; figure 8 (the scenario surface) carries it instead,
# with fill = the difference between the two models. Figure 8 is a tile plot
# (x = age, y = uptake) and is unaffected by this choice.
#
# ---------------------------------------------------------------------------
# Submission format
# ---------------------------------------------------------------------------
# Based on RespiCompass round-1 2026/2027 (round1_2627_rsv.md) with these
# PROJECT-SPECIFIC departures - do not "correct" them back to the wiki:
#   - columns: scenario_id, location, target, pop_group, target_end_date,
#              output_type_id, value.  NO round_id / horizon / output_type.
#   - immYes / immNo are dropped; only the _immTotal groups are kept.
#   - age groups: <1year, 1-17, 18-59, 60-64, 65-69, 70-74, 75-79, 80+, total
#     (so the hub's 0-2mo/3-5mo/6-11mo collapse to <1year, and 1-4 + 5-17
#      collapse to 1-17).
#   - rsv_infections is not produced; targets are rsv_hospitalisations and
#     administered_doses only.
#   - administered_doses carries the ACTUAL campaign date as target_end_date,
#     read from adult_vaccination_dates in the config - not the hub's
#     hard-coded scenario-start date.
#
# Both models are put through the SAME shaping function, so they are identical
# in structure by construction.
#
# ---------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------
# * Output is namespaced by git branch - see set_dirs() in R/directories.R.
# * Each <scenario>_raw.rds is read EXACTLY ONCE and reduced immediately to
#   (a) daily admissions by reporting band and (b) final dose counts. Those
#   files are ~128 MB each and the old scripts re-read them per figure basis.
#   Dropping immYes/immNo also means incidence_prop_vacc / hosp_prop_vacc are
#   no longer needed at all.
# * Peak memory is dominated by the accumulated submission frame
#   (weeks x bands x samples x scenarios x countries). Expect ~1 GB at 100
#   samples and 28 countries; see SHARD_SUBMISSION if that is too much.
# * Countries with no scenario output are skipped with a message, so this can
#   be run mid-sweep.
########################################################## #

rm(list = ls())
source("R/dependencies.R")

if (interactive()) clf()
if (interactive()) clc()

# ---- Settings ---------------------------------------------------------- ----

# Bootstrap one options list purely to resolve the branch-namespaced output
# root. The per-country loop builds its own `o` for each ISO. This only picks
# paths - it does not restrict which countries are processed - but it must name
# a country that has an input/<ISO>.yaml.
BOOTSTRAP_ISO <- "AT"
o_paths  <- set_options(do_step = 2, analysis_name = BOOTSTRAP_ISO, quiet = TRUE)
OUT_ROOT <- o_paths$pth$output

SCEN_ROOT <- file.path(OUT_ROOT, "2_scenarios")
OUT_DIR   <- file.path(OUT_ROOT, "3_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BASELINE <- "baseline"          # internal id of the comparator scenario
METRIC   <- "hospital_admissions"
DOSES    <- "n_doses"
QUANT    <- c(0.05, 0.95)

# Seasons run 1 August to 31 July. A 730-day run from 2026-09-01 spills a few
# weeks into 2028/2029; a relative change over a sliver of a season is
# meaningless, so drop anything shorter and say so.
MIN_SEASON_DAYS <- 180

DYNAMIC_LABEL <- "Dynamic (this model)"
STATIC_LABEL  <- "Static"

# Companion static model, in ITS OWN format (round_id / horizon / output_type,
# hub age bands, imm splits). Reshaped by this script exactly like our own
# output. Missing file is reported and skipped, not an error.
STATIC_FILE <- file.path(OUT_ROOT, "2026_2027_1_RSV_staticModel.parquet")

# Set TRUE to write one submission parquet per country instead of holding the
# whole frame in memory and combining at the end.
SHARD_SUBMISSION <- FALSE

# Point / interval styling. In geom_pointrange `size` is the interval LINE
# width and `fatten` multiplies the point on top of it, so raise POINT_FATTEN
# to enlarge markers without thickening the bars.
POINT_SIZE   <- 0.45
POINT_FATTEN <- 5

# ---- Submission age bands ---------------------------------------------- ----
# Model reporting band (age_group_map values in default.yaml) -> submission
# label. The three infant bands collapse to <1year and 1-5y + 5-18y collapse
# to 1-17, per the project format.
BAND_TO_SUB <- c("0-3m"   = "<1year",
                 "3-6m"   = "<1year",
                 "6-12m"  = "<1year",
                 "1-5y"   = "1-17",
                 "5-18y"  = "1-17",
                 "18-60y" = "18-59",
                 "60-65y" = "60-64",
                 "65-70y" = "65-69",
                 "70-75y" = "70-74",
                 "75-80y" = "75-79",
                 "80+y"   = "80+")

# Static model's own band labels -> the same submission labels.
STATIC_TO_SUB <- c("0-2mo"  = "<1year",
                   "3-5mo"  = "<1year",
                   "6-11mo" = "<1year",
                   "1-4"    = "1-17",
                   "5-17"   = "1-17",
                   "18-59"  = "18-59",
                   "60-64"  = "60-64",
                   "65-69"  = "65-69",
                   "70-74"  = "70-74",
                   "75-79"  = "75-79",
                   "80+"    = "80+")

SUB_LEVELS <- c("<1year", "1-17", "18-59", "60-64", "65-69",
                "70-74", "75-79", "80+", "total")

# Reporting window for rsv_hospitalisations: the round asks for the two RSV
# seasons from 2026-09-01 to 2028-05-28, i.e. the 91 weeks ending 2026-09-06
# .. 2028-05-28. Dropping the `horizon` COLUMN from the project format does not
# drop the reporting PERIOD, so the window still applies. The model runs on to
# 2028-08-30; that tail is simulated (and plotted) but not submitted.
# NB the first week is a 6-day sum: the hub's week 36 starts 2026-08-31 but the
# simulation starts 2026-09-01.
SUB_WEEK_FIRST <- as.Date("2026-09-06")
SUB_WEEK_LAST  <- as.Date("2028-05-28")

# ---- Helpers ------------------------------------------------------------ ----

season_label <- function(d) {
  y <- ifelse(month(d) >= 8, year(d), year(d) - 1)
  paste0(y, "/", y + 1)
}

scen_files <- function(iso) {
  d <- file.path(SCEN_ROOT, iso, "scenarios")
  if (!dir.exists(d)) return(character(0))
  f <- list.files(d, pattern = "_raw[.]rds$")
  setNames(file.path(d, f), sub("_raw[.]rds$", "", f))
}

# Scenario id -> the two plot axes. Ids follow the round's convention
# <letter>.<n>-<coverage>, where n = 1..5 is the eligibility threshold
# (60/65/70/75/80) and the suffix is the uptake percentage.
#
# This is the ONE parser used for both models: the static file arrives as a
# parquet with no yaml beside it, and keying both models off the shared
# scenario id is what guarantees they are scored on the same axes. For our own
# output the result is cross-checked against the config (see check_axes).
AGE_BY_INDEX <- c("1" = 60, "2" = 65, "3" = 70, "4" = 75, "5" = 80)

scenario_axes <- function(id) {
  m <- str_match(id, "^([A-Z])\\.([1-5])-(\\d+)$")
  data.table(scen      = id,
             family    = m[, 2],
             elig_age  = unname(AGE_BY_INDEX[m[, 3]]),
             uptake    = as.numeric(m[, 4]) / 100)
}

# Cross-check the id-derived axes against what the config actually configures,
# so a mislabelled scenario is caught rather than silently plotted in the wrong
# panel. Returns invisibly; warns on disagreement.
check_axes <- function(o, scen, ax) {
  p <- parse_yaml(o, scen)$parsed
  cov_yaml <- p$adult_vacc_coverage
  ages     <- unlist(p$adult_vaccination_agegroups)
  age_yaml <- suppressWarnings(min(as.numeric(str_extract(ages, "^\\d+")), na.rm = TRUE))
  if (is.finite(cov_yaml) && !is.na(ax$uptake) && abs(cov_yaml - ax$uptake) > 1e-9)
    warning(scen, ": id implies uptake ", ax$uptake, " but config sets ", cov_yaml)
  if (is.finite(age_yaml) && !is.na(ax$elig_age) && age_yaml != ax$elig_age)
    warning(scen, ": id implies ", ax$elig_age, "+ but config's youngest eligible band is ", age_yaml)
  invisible(NULL)
}

# Population by submission band, for the per-100k denominators.
#
# Built ONCE for every country in the hub file, not per country we happened to
# simulate. Population is reference data: keying it to the dynamic run meant a
# country present in the static file but not (yet) in our own output was
# silently dropped from every figure with a misleading "no population data".
# Reading the CSV once instead of once per country is also a good deal faster.
pop_all_countries <- function(o) {
  p <- read.csv(o$pop_url, fileEncoding = "UTF-8-BOM") %>% normalise_iso2()
  p$sub <- STATIC_TO_SUB[p$age_group]        # hub data uses the static labels
  p <- p[!is.na(p$sub), ]
  if (!nrow(p)) return(list())
  split(p, p$country) |>
    lapply(function(x) tapply(x$population, x$sub, sum))
}

# ---- 1. LOAD - this model ----------------------------------------------- ----
# Reads each scenario file once and returns the two reduced tables everything
# else is built from:
#   hosp : iso, scen, sim, date, sub_band, value   (daily admissions)
#   dose : iso, scen, sim, sub_band, doses         (final cumulative value)
load_dynamic_country <- function(iso) {

  files <- scen_files(iso)
  if (!length(files))              { message("  - ", iso, ": no scenario output, skipped"); return(NULL) }
  if (!BASELINE %in% names(files)) { message("  - ", iso, ": no baseline, skipped");        return(NULL) }
  others <- setdiff(names(files), BASELINE)
  if (!length(others))             { message("  - ", iso, ": baseline only, skipped");      return(NULL) }

  o     <- set_options(do_step = 2, analysis_name = iso, quiet = TRUE)
  fit   <- load_data(o, setup_calibration(o))
  start <- min(fit$dates_model$date)

  # Fine model age group -> reporting band -> submission label, from the config
  # rather than re-hardcoded here.
  agm      <- unlist(parse_yaml(o, BASELINE)$parsed$age_group_map)
  fine2sub <- BAND_TO_SUB[agm]
  names(fine2sub) <- names(agm)
  if (anyNA(fine2sub))
    stop(iso, ": no submission band for model age group(s): ",
         paste(names(fine2sub)[is.na(fine2sub)], collapse = ", "))

  # Campaign date, for the administered_doses target_end_date. Read from the
  # config so it follows adult_vaccination_dates instead of being hard-coded.
  camp <- ymd(unlist(parse_yaml(o, BASELINE)$parsed$adult_vaccination_dates))
  camp <- min(camp, na.rm = TRUE)

  hosp <- dose <- vector("list", length(files))

  for (i in seq_along(files)) {
    scen <- names(files)[i]
    r <- as.data.table(readRDS(files[[i]]))
    r <- r[metric %in% c(METRIC, DOSES) & !is.na(age_group)]
    if (!nrow(r)) next

    # sim = the trajectory, i.e. param_id with the trailing "_<scenario>"
    # removed. param_id is "s<param_set>_<fitting_set>_<scenario>".
    r[, sim      := sub(paste0("_", scen, "$"), "", param_id)]
    r[, sub_band := fine2sub[age_group]]
    r <- r[!is.na(sub_band)]

    h <- r[metric == METRIC,
           .(value = sum(value, na.rm = TRUE)), by = .(sim, time, sub_band)]
    if (nrow(h)) {
      h[, `:=`(iso = iso, scen = scen, date = start + time - 1L)]
      h[, time := NULL]
      hosp[[i]] <- h
    }

    # Doses are a CUMULATIVE counter: take the final value, never the sum.
    d <- r[metric == DOSES]
    if (nrow(d)) {
      d <- d[time == max(time), .(doses = sum(value, na.rm = TRUE)), by = .(sim, sub_band)]
      d[, `:=`(iso = iso, scen = scen)]
      dose[[i]] <- d
    }
    rm(r); invisible(gc(FALSE))
  }

  hosp <- rbindlist(hosp, fill = TRUE)
  dose <- rbindlist(dose, fill = TRUE)
  if (!nrow(hosp)) { message("  - ", iso, ": no admissions output, skipped"); return(NULL) }

  # Axes for every non-baseline scenario, validated against the config.
  ax <- rbindlist(lapply(others, scenario_axes))
  for (s in others) {
    a <- ax[scen == s]
    if (nrow(a) && !is.na(a$elig_age)) check_axes(o, s, a) else
      warning(iso, ": cannot parse scenario id '", s, "' into (age, uptake) axes")
  }

  message("  + ", iso, ": ", length(others), " scenarios, ",
          uniqueN(hosp$sim), " trajectories")

  list(hosp = hosp, dose = dose, axes = ax, start = start, camp = camp,
       scenarios = others)
}

# ---- 2. LOAD - static model --------------------------------------------- ----
# Returns the same two reduced tables, so both models flow through identical
# downstream code. The static file is WEEKLY, not daily, so `date` here is the
# week-ending date; season assignment and submission weeks both work off it.
load_static <- function() {

  if (!file.exists(STATIC_FILE)) {
    message("* Static model file not found (", STATIC_FILE, ") - single-model output")
    return(NULL)
  }
  read_pq <- if (requireNamespace("nanoparquet", quietly = TRUE)) nanoparquet::read_parquet else
             if (requireNamespace("arrow", quietly = TRUE))       arrow::read_parquet else NULL
  if (is.null(read_pq)) {
    message("* Neither nanoparquet nor arrow installed - cannot read the static parquet")
    return(NULL)
  }

  d <- as.data.table(read_pq(STATIC_FILE))
  d <- normalise_iso2(d, col = "location")        # static uses EL for Greece
  d[, target_end_date := as.Date(target_end_date)]

  # Only the _immTotal rows are the whole band; immYes/immNo split it and would
  # double-count if summed alongside. 'total_immTotal' is the file's own total
  # and must not be folded into the age bands.
  hosp <- d[target == "rsv_hospitalisations" &
            grepl("_immTotal$", pop_group) &
            !grepl("^total_", pop_group)]
  hosp[, band     := sub("_immTotal$", "", pop_group)]
  hosp[, sub_band := STATIC_TO_SUB[band]]
  if (anyNA(hosp$sub_band))
    warning("static: unmapped age band(s) dropped: ",
            paste(unique(hosp$band[is.na(hosp$sub_band)]), collapse = ", "))
  hosp <- hosp[!is.na(sub_band),
               .(value = sum(value, na.rm = TRUE)),
               by = .(iso = location, scen = scenario_id, sim = output_type_id,
                      date = target_end_date, sub_band)]

  # Doses: the file reports these BOTH per eligible age band AND as a national
  # 'undefined' total, and the bands sum exactly to that total - so summing every
  # administered_doses row double-counts. Take the total when it is present and
  # fall back to the bands when it is not (older static files carried only the
  # 'undefined' row, which is how this slipped through on the earlier data).
  dz <- d[target == "administered_doses"]
  if (nrow(dz) && "undefined" %in% dz$pop_group) {
    dz <- dz[pop_group == "undefined"]
  } else if (nrow(dz)) {
    message("  (static: no 'undefined' dose total, summing the per-band rows)")
  }
  # NB a constant cannot go in `by` (data.table requires every `by` element to
  # be the same length as the subset), so sub_band is attached afterwards.
  dose <- if (nrow(dz)) {
    tmp <- dz[, .(doses = sum(value, na.rm = TRUE)),
              by = .(iso = location, scen = scenario_id, sim = output_type_id)]
    tmp[, sub_band := "undefined"][]
  } else {
    data.table(iso = character(), scen = character(), sim = character(),
               doses = numeric(), sub_band = character())
  }

  # The static file states its own dose date; keep it rather than imposing ours.
  # NB it emits an administered_doses row for EVERY week, zero except in the
  # campaign week, so the campaign date is the first week with a non-zero value -
  # not min(target_end_date), which is just the start of the reporting period.
  nz   <- dz[value > 0]
  camp <- if (nrow(nz)) min(nz$target_end_date, na.rm = TRUE)
          else if (nrow(dz)) min(dz$target_end_date, na.rm = TRUE)
          else as.Date(NA)

  message("* Static model: ", uniqueN(hosp$iso), " countries, ",
          uniqueN(hosp$scen), " scenarios, ", uniqueN(hosp$sim), " trajectories")

  list(hosp = hosp, dose = dose, camp = camp)
}

# ---- 3. Scenario impact -------------------------------------------------- ----
# Per (iso, scenario, season, sample): burden vs the SAME sample's baseline, on
# two age bases. Pairing per sample means shared parameter uncertainty cancels
# instead of inflating the interval.
impact_table <- function(hosp, dose, pop_by_iso, axes, model_label) {

  if (is.null(hosp) || !nrow(hosp)) return(NULL)
  hosp <- copy(hosp)
  hosp[, season := season_label(date)]

  # Seasons this run actually covers, per country.
  # Our output is daily, the static file weekly, so the minimum-length test is
  # expressed in the units each one actually carries.
  min_pts <- if (model_label == STATIC_LABEL) MIN_SEASON_DAYS / 7 else MIN_SEASON_DAYS
  span <- hosp[, .(days = uniqueN(date)), by = .(iso, season)]
  drop <- span[days < min_pts, .(iso, season)]
  if (nrow(drop))
    message("    (", model_label, ": dropping partial season(s) ",
            paste(unique(sprintf("%s %s", drop$iso, drop$season)), collapse = ", "), ")")
  if (nrow(drop)) hosp <- hosp[!drop, on = .(iso, season)]
  if (!nrow(hosp)) return(NULL)

  out <- list()

  for (this_iso in unique(hosp$iso)) {

    h   <- hosp[iso %in% this_iso]
    pop <- pop_by_iso[[this_iso]]
    if (is.null(pop)) { message("  ! ", this_iso, ": no population data, skipped"); next }
    pop_total <- sum(pop, na.rm = TRUE)

    scens <- setdiff(unique(h$scen), BASELINE)
    if (!length(scens)) next
    ax_i <- axes[scen %in% scens][!is.na(elig_age)]

    # Eligible submission bands per scenario, and their union across the
    # scenarios this country actually has, so every scenario in figures 2 and 7
    # is scored on the SAME population.
    bands_for <- function(min_age) {
      adult <- c("60-64" = 60, "65-69" = 65, "70-74" = 70, "75-79" = 75, "80+" = 80)
      names(adult)[adult >= min_age]
    }
    elig <- setNames(lapply(ax_i$elig_age, bands_for), ax_i$scen)
    elig_union <- sort(unique(unlist(elig)))
    if (!length(elig_union)) next

    # NB argument names must NOT match column names: inside `i` a data.table
    # column masks a same-named variable, and the `..` prefix only works in `j`.
    burden <- function(want_scen, want_bands)
      h[scen %in% want_scen & sub_band %in% want_bands,
        .(burden = sum(value, na.rm = TRUE)), by = .(sim, season)]

    for (s in scens) {
      for (basis in c("scenario", "union")) {
        bands <- if (basis == "scenario") elig[[s]] else elig_union
        if (!length(bands)) next
        b_sc <- burden(s, bands); b_bl <- burden(BASELINE, bands)
        if (!nrow(b_sc) || !nrow(b_bl)) next
        m <- merge(b_sc, b_bl, by = c("sim", "season"), suffixes = c("", "_bl"))
        if (!nrow(m)) next
        pop_elig <- sum(pop[bands], na.rm = TRUE)
        m[, `:=`(pct     = 100 * (burden - burden_bl) / burden_bl,
                 averted = burden_bl - burden)]
        out[[length(out) + 1]] <- m[, .(
          iso = this_iso, scen = s, season, basis, sim,
          pct, averted,
          av_per100k_tot  = 1e5 * averted / pop_total,
          av_per100k_elig = 1e5 * averted / pop_elig)]
      }

      # Doses: no season dimension (a single campaign).
      d <- dose[iso == this_iso & scen %in% s]
      if (nrow(d)) {
        bands <- elig[[s]]
        # Static doses arrive as a single national total (see load_static);
        # ours are split by band. Use whatever is there.
        dd <- if (all(d$sub_band == "undefined"))
                d[, .(doses = sum(doses, na.rm = TRUE)), by = sim] else
                d[sub_band %in% bands, .(doses = sum(doses, na.rm = TRUE)), by = sim]
        if (nrow(dd)) {
          pop_elig <- sum(pop[bands], na.rm = TRUE)
          out[[length(out) + 1]] <- data.table(
            iso = this_iso, scen = s, season = NA_character_, basis = "doses",
            sim = dd$sim, pct = NA_real_, averted = NA_real_,
            av_per100k_tot  = 1e5 * dd$doses / pop_total,
            av_per100k_elig = 1e5 * dd$doses / pop_elig)
        }
      }
    }
  }

  if (!length(out)) return(NULL)
  res <- rbindlist(out, fill = TRUE)
  res[, model := model_label]
  merge(res, axes, by = "scen", all.x = TRUE)
}

# ---- 4. Plotting --------------------------------------------------------- ----

summ <- function(d, col) {
  d[!is.na(get(col)),
    .(median = median(get(col), na.rm = TRUE),
      lo     = quantile(get(col), QUANT[1], na.rm = TRUE),
      hi     = quantile(get(col), QUANT[2], na.rm = TRUE),
      n_sim  = .N),
    by = .(iso, scen, season, model, elig_age, uptake)]
}

# Expected-reduction reference. For a directly-protected population the
# seasonal reduction is roughly UPTAKE x VE, so each uptake level gets a band
# spanning VE at the season's start and end, plus a midline.
#
# VE is the mean across all 500 replicate waning curves. Season s (1-indexed
# from the campaign) spans months (s-1)*12 to s*12. Drawn negative to match the
# plot's negative-is-averted axis.
#
# These are a reference, not a fit: points on the midline mean the realised
# reduction matches uptake x VE exactly, and departures show where transmission
# dynamics add to or subtract from it (indirect protection, susceptible
# depletion, epidemic timing relative to the campaign).
ve_bands <- function(uptakes, season, camp, ve_col = "VE_inf") {
  wf <- file.path("data", "respicompass_cache", "waning_curves.csv")
  if (!file.exists(wf)) { message("  ! waning curves not found, no VE bands"); return(NULL) }
  w  <- as.data.table(read.csv(wf))

  # Months since the campaign at the start and end of this season. Measured
  # from the ACTUAL campaign date rather than assuming it coincides with the
  # season boundary: with a 1 November campaign, season 1 covers months 0-9
  # since vaccination, not 0-12, so the season-index shortcut would quote a VE
  # that is three months too waned.
  y  <- as.integer(substr(season, 1, 4))
  s0 <- max(as.Date(paste0(y, "-08-01")), camp)
  s1 <- as.Date(paste0(y + 1L, "-07-31"))
  if (s1 <= camp) return(NULL)
  m0 <- max(0L, as.integer(floor(as.numeric(s0 - camp) / 30.44)))
  m1 <- as.integer(floor(as.numeric(s1 - camp) / 30.44))
  if (!all(c(m0, m1) %in% w$month)) return(NULL)
  v0 <- mean(w[month == m0][[ve_col]]); v1 <- mean(w[month == m1][[ve_col]])
  rbindlist(lapply(uptakes, function(cv) data.table(
    uptake = cv,
    ymin = -100 * cv * v0, ymax = -100 * cv * v1,
    ymid = -100 * cv * (v0 + v1) / 2)))
}

uptake_lab <- function(x) paste0(round(100 * x), "%")

# Layout B: x = country, facet = eligibility age, colour = uptake.
make_layout_b <- function(tab, title, subtitle, ylab, log_y = FALSE, bands = NULL) {

  d <- copy(tab)
  d[, iso := factor(iso, levels = sort(unique(iso)))]
  d[, uptake_f := factor(uptake_lab(uptake), levels = uptake_lab(sort(unique(uptake))))]
  d[, elig_f := factor(paste0(elig_age, "+"),
                       levels = paste0(sort(unique(elig_age)), "+"))]

  g <- ggplot(d, aes(x = iso, y = median))

  # Bands first, so point ranges draw over them. Uptake is the FACET, so each
  # panel has exactly one uptake and therefore exactly one band - it needs no
  # colour of its own, and a neutral grey keeps the colour scale free for the
  # eligibility threshold. (When uptake was the colour aesthetic this had to be
  # three overlapping translucent bands per panel.) `uptake_f` is the facet
  # variable, so ggplot drops each band into its own panel automatically.
  if (!is.null(bands) && nrow(bands)) {
    b <- copy(bands)
    b[, uptake_f := factor(uptake_lab(uptake), levels = levels(d$uptake_f))]
    g <- g +
      geom_rect(data = b, inherit.aes = FALSE,
                aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
                fill = "grey30", alpha = 0.12) +
      geom_hline(data = b, aes(yintercept = ymid),
                 linetype = "dashed", linewidth = 0.4, colour = "grey25")
  }

  if (!log_y) g <- g + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")

  if (log_y) {
    bad <- d[!is.finite(median) | median <= 0 | !is.finite(lo) | lo <= 0]
    if (nrow(bad)) {
      warning("log scale: dropping ", nrow(bad), " non-positive point(s): ",
              paste(unique(paste(bad$iso, bad$scen)), collapse = "; "))
      d <- d[is.finite(median) & median > 0 & is.finite(lo) & lo > 0]
      g <- g %+% d
    }
  }

  g <- g +
    geom_pointrange(aes(ymin = lo, ymax = hi, colour = elig_f, group = elig_f),
                    position = position_dodge(width = 0.8),
                    size = POINT_SIZE, fatten = POINT_FATTEN) +
    scale_colour_viridis_d(name = "Eligibility age", option = "C", end = 0.85) +
    # Colour = eligibility age, facet = uptake. Both are ordered, so the
    # sequential palette still reads correctly: within a country the five
    # thresholds should fan out in order, and a crossing is worth a look.
    # ncol = 1 so panels stack: 28 countries on x makes each panel wide.
    # Fixed y (not free) - comparing 25% against 75% is the point of the facet,
    # and a per-panel axis would silently rescale that away.
    facet_wrap(~uptake_f, ncol = 1) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey85"),
          legend.position = "bottom",
          panel.grid.major.x = element_blank())

  if (log_y) g <- g + scale_y_log10(labels = scales::label_comma())
  g
}

# Layout C: the scenario grid IS a 2-D surface, so tile it. This is the only
# view that shows the INTERACTION between the two axes - whether extra uptake
# buys more at 60+ or at 80+ - and it scales to every country at once.
# Uncertainty is not shown here; layout B carries that.
make_layout_c <- function(tab, title, subtitle, fill_lab, diverging = FALSE) {

  d <- copy(tab)
  d[, uptake_f := factor(uptake_lab(uptake), levels = uptake_lab(sort(unique(uptake))))]
  d[, elig_f := factor(paste0(elig_age, "+"),
                       levels = paste0(sort(unique(elig_age)), "+"))]

  g <- ggplot(d, aes(x = elig_f, y = uptake_f, fill = median)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_wrap(~iso) +
    labs(title = title, subtitle = subtitle,
         x = "Eligibility age threshold", y = "Uptake", fill = fill_lab) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey85"),
          legend.position = "bottom",
          panel.grid = element_blank())

  if (diverging)
    g + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0)
  else
    g + scale_fill_viridis_c(option = "C", direction = -1, end = 0.92)
}

save_fig <- function(g, name, tab, n_facet = 5, wide = FALSE) {
  f <- file.path(OUT_DIR, name)
  ggsave(paste0(f, ".png"), g,
         width  = if (wide) 14 else 12,
         height = if (wide) 11 else 2.6 + 2.2 * n_facet,
         dpi = 200, limitsize = FALSE)
  fwrite(tab[order(scen, iso)], paste0(f, ".csv"))
  message("  saved ", basename(f), ".png / .csv")
}

# ---- 5. Submission shaping ----------------------------------------------- ----
# One function, applied to BOTH models, so they are identical by construction.
#
#   hosp : iso, scen, sim, date, sub_band, value
#          `date` is DAILY for our model and WEEKLY (week-ending) for static;
#          both are aggregated onto the same Monday-start weeks below.
#   dose : iso, scen, sim, sub_band, doses
#   camp : campaign date, used as target_end_date for administered_doses
to_submission <- function(hosp, dose, camp) {

  if (is.null(hosp) || !nrow(hosp)) return(NULL)

  h <- copy(hosp)
  # Week-ending Sunday, matching the round's target_end_date convention.
  h[, target_end_date := floor_date(date, "week", week_start = 1) + 6L]
  wk <- h[, .(value = sum(value, na.rm = TRUE)),
          by = .(location = iso, scenario_id = scen, output_type_id = sim,
                 target_end_date, sub_band)]

  # Age totals. Summed from the bands, which are mutually exclusive by
  # construction (the static file's own 'total_*' rows were excluded on read).
  tot <- wk[, .(value = sum(value, na.rm = TRUE), sub_band = "total"),
            by = .(location, scenario_id, output_type_id, target_end_date)]

  hosp_rows <- rbind(wk, tot, fill = TRUE)[
    target_end_date >= SUB_WEEK_FIRST & target_end_date <= SUB_WEEK_LAST
    , .(scenario_id, location, target = "rsv_hospitalisations",
        pop_group = paste0(sub_band, "_immTotal"),
        target_end_date, output_type_id, value)]

  dose_rows <- NULL
  if (!is.null(dose) && nrow(dose)) {
    dz <- dose[, .(value = sum(doses, na.rm = TRUE)),
               by = .(location = iso, scenario_id = scen, output_type_id = sim)]
    dose_rows <- dz[, .(scenario_id, location, target = "administered_doses",
                        pop_group = "undefined",
                        target_end_date = camp, output_type_id, value)]
  }

  out <- rbind(hosp_rows, dose_rows, fill = TRUE)
  setorder(out, scenario_id, location, target, pop_group, target_end_date, output_type_id)
  out[]
}

# Report anything that would make the file wrong, rather than writing quietly.
check_submission <- function(d, label) {
  if (is.null(d) || !nrow(d)) { message("  ! ", label, ": nothing to write"); return(invisible()) }
  bad_pg <- setdiff(sub("_immTotal$", "", unique(d$pop_group)), c(SUB_LEVELS, "undefined"))
  if (length(bad_pg))
    warning(label, ": unexpected pop_group(s): ", paste(bad_pg, collapse = ", "))
  n_traj <- uniqueN(d$output_type_id)
  if (n_traj < 100 || n_traj > 300)
    warning(label, ": ", n_traj, " trajectories (RespiCompass expects 100-300)")
  if (!"administered_doses" %in% d$target)
    warning(label, ": administered_doses target is missing")
  h <- d[target == "rsv_hospitalisations"]
  if (nrow(h)) {
    nw <- uniqueN(h$target_end_date)
    if (nw != 91)
      warning(label, ": ", nw, " reporting weeks, expected 91 (",
              format(min(h$target_end_date)), " .. ", format(max(h$target_end_date)), ")")
  }
  if (anyNA(d$value)) warning(label, ": ", sum(is.na(d$value)), " NA values")
  message("  ", label, ": ", format(nrow(d), big.mark = ","), " rows, ",
          uniqueN(d$location), " countries, ", uniqueN(d$scenario_id), " scenarios, ",
          n_traj, " trajectories")
}

# ========================================================================= ====
# RUN
# ========================================================================= ====

# ---- Load ---------------------------------------------------------------- ----
isos <- sort(list.dirs(SCEN_ROOT, full.names = FALSE, recursive = FALSE))
message("* Scanning ", length(isos), " countries in ", SCEN_ROOT)

dyn <- lapply(isos, load_dynamic_country)
names(dyn) <- isos
dyn <- dyn[!vapply(dyn, is.null, logical(1))]
if (!length(dyn)) stop("No usable scenario output found under ", SCEN_ROOT)

dyn_hosp <- rbindlist(lapply(dyn, `[[`, "hosp"), fill = TRUE)
dyn_dose <- rbindlist(lapply(dyn, `[[`, "dose"), fill = TRUE)
dyn_axes <- unique(rbindlist(lapply(dyn, `[[`, "axes"), fill = TRUE))
dyn_camp <- min(do.call(c, lapply(dyn, `[[`, "camp")), na.rm = TRUE)

# Shared across both models, so they are scored on identical denominators.
POP <- pop_all_countries(o_paths)

message("* This model: ", length(dyn), " countries, ",
        uniqueN(dyn_hosp$scen), " scenarios, ", uniqueN(dyn_hosp$sim), " trajectories")
message("* Campaign date (from config): ", format(dyn_camp))

sta <- load_static()

# ---- Impact tables ------------------------------------------------------- ----
message("\n* Building impact tables")
imp <- list(impact_table(dyn_hosp, dyn_dose, POP, dyn_axes, DYNAMIC_LABEL))

if (!is.null(sta)) {
  sta_axes <- unique(rbindlist(lapply(setdiff(unique(sta$hosp$scen), BASELINE),
                                      scenario_axes), fill = TRUE))
  imp[[2]] <- impact_table(sta$hosp, sta$dose, POP, sta_axes, STATIC_LABEL)
  bad <- sta_axes[is.na(elig_age) | is.na(uptake)]$scen
  if (length(bad))
    message("  ! Static scenario id(s) not in <letter>.<1-5>-<coverage> form, so they",
            " cannot be placed on the (age, uptake) axes and are EXCLUDED from every",
            " figure: ", paste(bad, collapse = ", "),
            "
    (they are still reshaped into the submission file)")
}
imp <- rbindlist(imp[!vapply(imp, is.null, logical(1))], fill = TRUE)
imp <- imp[!is.na(elig_age) & !is.na(uptake)]

# ---- Figures ------------------------------------------------------------- ----
message("\n* Writing figures")

seasons <- sort(unique(imp$season[!is.na(imp$season)]))
models  <- unique(imp$model)

# (season-bearing figures) x (season) x (model)
specs <- list(
  list(basis = "scenario", col = "pct", n = "1_pct_averted_scenario_ages", bands = TRUE,
       t = "Relative change in seasonal RSV hospitalisations vs baseline (scenario's eligible ages)",
       y = "Relative change (%)  -  negative = averted"),
  list(basis = "union", col = "pct", n = "2_pct_averted_union_ages", bands = TRUE,
       t = "Relative change in seasonal RSV hospitalisations vs baseline (union of eligible ages)",
       y = "Relative change (%)  -  negative = averted"),
  # Log10: absolute counts span ~2 orders of magnitude across countries, so a
  # linear axis is dominated by the largest and the smallest are unreadable.
  list(basis = "scenario", col = "averted", n = "3_abs_averted", log = TRUE,
       t = "Hospitalisations averted vs baseline (scenario's eligible ages)",
       y = "Hospitalisations averted (count)"),
  list(basis = "scenario", col = "av_per100k_tot", n = "5_averted_per100k_total",
       t = "Hospitalisations averted per 100k total population",
       y = "Averted per 100k total population"),
  list(basis = "scenario", col = "av_per100k_elig", n = "6_averted_per100k_eligible",
       t = "Hospitalisations averted per 100k eligible population (scenario's eligible ages)",
       y = "Averted per 100k eligible population"),
  list(basis = "union", col = "av_per100k_elig", n = "7_averted_per100k_union",
       t = "Hospitalisations averted per 100k eligible population (union of eligible ages)",
       y = "Averted per 100k eligible population")
)

for (sp in specs) {
  for (mdl in models) {
    for (i in seq_along(seasons)) {
      sn  <- seasons[i]
      d   <- imp[basis == sp$basis & model == mdl & season == sn]
      if (!nrow(d)) next
      tab <- summ(d, sp$col)
      if (!nrow(tab)) next
      bnd <- if (isTRUE(sp$bands)) ve_bands(sort(unique(tab$uptake)), sn, dyn_camp) else NULL
      nm  <- paste0("scenario_impact_", sp$n, "_",
                    gsub("/", "-", sn), "_",
                    if (mdl == STATIC_LABEL) "static" else "dynamic")
      save_fig(make_layout_b(tab, sp$t,
                             paste0(mdl, "  |  season ", sn,
                                    "  |  median and ", 100 * (QUANT[2] - QUANT[1]),
                                    "% interval across paired samples",
                                    if (isTRUE(sp$bands))
                                      "\nShaded band = expected reduction (uptake x VE), dashed line = its midpoint"
                                    else ""),
                             sp$y, log_y = isTRUE(sp$log), bands = bnd),
               nm, tab, n_facet = uniqueN(tab$uptake))
    }
  }
}

# Dose figures: no season dimension (a single campaign).
dose_specs <- list(
  list(col = "av_per100k_tot",  n = "4a_doses_per100k_total",
       t = "Vaccine doses administered, per 100k total population",
       y = "Doses per 100k total population"),
  # Kept as a SANITY CHECK, not a result. Doses are counted as uptake x the
  # whole eligible population (see ageing_event in R/model.R), so this should
  # sit flat at 100000 x uptake for every country, give or take demographic
  # drift between the census denominator and the campaign date. A country off
  # its line means the eligible-population denominator is wrong - which is how
  # the fine-vs-reporting age label mismatch was caught previously.
  list(col = "av_per100k_elig", n = "4b_doses_per100k_eligible",
       t = "Vaccine doses administered, per 100k eligible population  (SANITY CHECK - expect flat)",
       y = "Doses per 100k eligible population")
)

for (sp in dose_specs) {
  for (mdl in models) {
    d <- imp[basis == "doses" & model == mdl]
    if (!nrow(d)) next
    tab <- summ(d, sp$col)
    if (!nrow(tab)) next
    nm <- paste0("scenario_impact_", sp$n, "_",
                 if (mdl == STATIC_LABEL) "static" else "dynamic")
    save_fig(make_layout_b(tab, sp$t,
                           paste0(mdl, "  |  median and ", 100 * (QUANT[2] - QUANT[1]),
                                  "% interval across paired samples"),
                           sp$y),
             nm, tab, n_facet = uniqueN(tab$uptake))
  }
}

# ---- Figure 8: scenario surface (layout C) ------------------------------- ----
for (mdl in models) {
  for (sn in seasons) {
    d <- imp[basis == "scenario" & model == mdl & season == sn]
    if (!nrow(d)) next
    tab <- summ(d, "pct")
    if (!nrow(tab)) next
    nm <- paste0("scenario_impact_8_surface_pct_averted_", gsub("/", "-", sn), "_",
                 if (mdl == STATIC_LABEL) "static" else "dynamic")
    save_fig(make_layout_c(tab,
                           "Relative change in seasonal RSV hospitalisations vs baseline",
                           paste0(mdl, "  |  season ", sn,
                                  "  |  median across paired samples (no interval shown - see figure 1)"),
                           "Relative change (%)"),
             nm, tab, wide = TRUE)
  }
}

# Model difference, if both are present: the one view that puts the two models
# side by side now that they are in separate files elsewhere.
if (length(models) > 1) {
  for (sn in seasons) {
    d <- imp[basis == "scenario" & season == sn]
    tab <- summ(d, "pct")
    w <- dcast(tab, iso + scen + elig_age + uptake ~ model, value.var = "median")
    if (!all(c(DYNAMIC_LABEL, STATIC_LABEL) %in% names(w))) next
    w[, median := get(DYNAMIC_LABEL) - get(STATIC_LABEL)]
    w <- w[!is.na(median)]
    if (!nrow(w)) next
    w[, `:=`(season = sn, model = "difference", scen = scen)]
    nm <- paste0("scenario_impact_8_surface_model_difference_", gsub("/", "-", sn))
    save_fig(make_layout_c(w,
                           "Model difference in relative change (dynamic - static)",
                           paste0("Season ", sn,
                                  "  |  negative = this model averts more than the static model"),
                           "Difference (pp)", diverging = TRUE),
             nm, w, wide = TRUE)
  }
}

# ---- Submission files ---------------------------------------------------- ----
message("\n* Building submission files")

sub_dyn <- to_submission(dyn_hosp, dyn_dose, dyn_camp)
check_submission(sub_dyn, DYNAMIC_LABEL)

sub_sta <- if (!is.null(sta)) to_submission(sta$hosp, sta$dose, sta$camp) else NULL
if (!is.null(sub_sta)) check_submission(sub_sta, STATIC_LABEL)

write_sub <- function(d, file) {
  if (is.null(d) || !nrow(d)) return(invisible())
  p <- file.path(OUT_DIR, file)
  write_parquet(d, p, compression = "gzip")
  message("  wrote ", basename(p), " (", format(nrow(d), big.mark = ","), " rows)")
}

write_sub(sub_dyn, "2026_2027_1_RSV_dynamicModel.parquet")
write_sub(sub_sta, "2026_2027_1_RSV_staticModel_reshaped.parquet")

message("\n* Done - ", OUT_DIR)
