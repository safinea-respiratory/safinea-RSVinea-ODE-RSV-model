########################################################## #
# MODEL
#
# Simulate transmission model for a given parameter set
# (defined by p). Output is a long datatable, denoted m, and
# several other model properties of interest.
#
########################################################## #

# -------------------------------------------------------- -
# Main function: The guts of the model ----
# -------------------------------------------------------- -
model = function(o, scenario, fit = NULL, uncert = NULL, do_plot = TRUE, verbose = "date") {
  
  # Initiate model timer
  tic("model")
  
  # ---- Generate model parameters ----
  if (verbose != "none") message(" - Parsing input")
  data = fit$data
  dates_df = fit$dates_model
  fit$data = NULL
  fit$dates_model = NULL
  
  # Load and parse user-defined inputs for this scenario
  yaml = parse_yaml(o, scenario, fit = fit, uncert = uncert)
  
  # Shorthand for model parameters
  p = yaml$parsed
  
  #---- Contact matrix ----
  # Import synthetic contact matrices (Prem et al. 2021)
  contact_matrices = load_contacts(o)
  
  # Select contact matrix for country specified in .yaml
  country = p$contact_matrix_countries
  p$contact_matrix = contact_matrices[[country]]

  p$age_groups  = unlist(p$age_groups)
  p$age_breaks  = unlist(p$age_breaks_months) / 12
  
  
  # Change age bands
  p$contact_matrix <- reband_contact_matrix(mat = p$contact_matrix,
                                            original_breaks = seq(0, 80, 5),
                                            target_breaks = p$age_breaks)
  row.names(p$contact_matrix) = p$age_groups
  colnames(p$contact_matrix) = p$age_groups
  p$n_age = length(p$age_groups)
  
  # Age group mapping to match data and risk vectors (defined in default.yaml)
  p$age_group_map = data.frame(
    smaller_group = p$age_groups,
    larger_group  = unlist(p$age_group_map[p$age_groups])
  )
  
  # Get full country name from the RespiCompass country list.
  # The YAML stores the ISO-2 code under calibration_options$country;
  # use that to look up the name consistently with how load_data.R does it.
  country_iso2 = p$calibration_options$country
  country_name = o$countries_df %>%
    filter(iso2_code == country_iso2) %>%
    pull(country)

  # Load the population data and adjust the group size
  population_github_df = read.csv(o$pop_url, fileEncoding = "UTF-8-BOM") %>%
    filter(country == country_name) %>%
    remap_age_groups(o$respicompass_age_map) 
  # Adjust the population group sizes assuming uniform distribution (taking into account the age band width)
  p$population <- redistribute_population(
    coarse_df = population_github_df,
    fine_age_groups = p$age_groups,
    fine_age_breaks = p$age_breaks,
    country_name = country_name
  )
  
  p$births = o$births %>%
    filter(country == country_name) 
  
  # The burden data loaded in load_data.R is already filtered to the configured
  # country, so we only need to select the 4-weekly frequency here.
  p$burden = data %>%
    filter(data_freq == "4-weekly")
  
  # Account for relative susceptibility, hospitalisation and mortality rates by age
  p = age_relativity(p)

  # ---- Model set up ---
  if (verbose != "none") message(" - Running model")
  
  # List of compartments.
  # Naming convention: <state><tier>[_stage], where state is one of
  #   V  = vaccinated, not yet infected (waning vaccine immunity)
  #   S  = susceptible
  #   E  = latent (exposed, pre-infectious)
  #   Ev = latent, came from vaccinated stream (for VE-against-severity tracking)
  #   I  = infectious
  #   Iv = infectious, came from vaccinated stream
  #   H  = hospitalised
  #   R  = recovered (temporary immunity, wanes to next susceptible tier)
  #   D  = dead (cumulative, absorbing)
  # Tier = number of PRIOR infections: 0 = naive, 1 = one prior, 2 = two, 3 = three+.
  #
  # Infant vaccination (V0): single compartment; waning encoded via age-group
  #   position in infant_vaccine_rel_protection (age ≈ time since vaccination).
  #
  # Adult vaccination (V1_j, V2_j, V3_j): staged waning chain of length W.
  #   Each month the ageing event advances individuals one step (V_k_j → V_k_{j+1}).
  #   At stage W individuals return to S_k (fully waned). Adult protection at
  #   stage j is given by adult_vaccine_rel_protection[j].
  #
  # E0v / I0v track infants infected while in V0 (for infant VE-against-severity).
  # E1v–E3v / I1v–I3v track adults infected while in any V_k stage (same purpose).
  p$compartments = c(
    "V0",
    paste0("V1_", seq_len(p$W)),   # adult waning chain, tier 1
    paste0("V2_", seq_len(p$W)),   # adult waning chain, tier 2
    paste0("V3_", seq_len(p$W)),   # adult waning chain, tier 3
    "S0", "S1", "S2", "S3",
    "E0", "E0v", "E1", "E1v", "E2", "E2v", "E3", "E3v",
    "I0", "I0v", "I1", "I1v", "I2", "I2v", "I3", "I3v",
    "H0", "H1", "H2", "H3",
    "R0", "R1", "R2", "R3",
    "D0", "D1", "D2", "D3"
  )
  
  # Seed the initial infection
  states = initiate_epidemic(p, verbose)
  
  # make sure all event days are integers
  p[grepl('_day',names(p))] <- lapply(p[grepl('_day',names(p))],round)
  
  # Define times to have an ageing event (first day of each month within the
  # simulation window). Events are computed here from the data start date and
  # n_days, so they follow the simulation period automatically.
  p$start_date = min(dates_df$date)
  start_date <- p$start_date
  end_date   <- p$start_date + p$n_days
  from <- if_else(day(start_date) == 1, # if the first day is start of the month
                  floor_date(start_date, "month"),
                  ceiling_date(start_date, "month"))
  to <- floor_date(end_date, "month")
  event_dates <- seq(from, to, by = "1 month")
  event_times <- as.numeric(event_dates - p$start_date)

  # ---- Sanity check: season_effect coverage ----
  # The ODE applies a per-season scalar on beta, indexed by the season number
  # returned by get_season_number(). The vector season_effect has length 3
  # (season1 = 1, season2_effect, season3_effect). If the simulation spans
  # more seasons than that, later seasons fall back to scalar = 1 in the ODE,
  # which is likely not what the user intended — warn so they can extend the
  # configuration (e.g. add season4_effect) or shorten n_days.
  n_seasons_simulated = get_season_number(end_date, start_date)
  n_seasons_configured = length(c(1, p$season2_effect, p$season3_effect))
  if (n_seasons_simulated > n_seasons_configured)
    warning("Simulation spans ", n_seasons_simulated, " seasons but only ",
            n_seasons_configured, " season_effect values are configured ",
            "(season1 = 1, season2_effect, season3_effect). Seasons ",
            n_seasons_configured + 1, "-", n_seasons_simulated,
            " will use scalar = 1. Either shorten n_days or extend ",
            "season_effect in the yaml.")
  
  # Solve ODE model
  out = deSolve::ode(y = states,
                     times = seq(1, p$n_days, by = 1),
                     func = rsv_model,
                     p = p,
                     events = list(func = ageing_event, time = event_times),
                     method = "vode",
                     atol = 1e-4,    # absolute tolerance
                     rtol = 1e-4)    # relative tolerance
  
  # Convert output to data frame
  out_df = as_tibble(out) 
  
  # Coerce all columns to numeric (not deSolve format)
  out_df[] = lapply(out_df, as.numeric)
  
  # Format output
  m = format_output(out_df, p)

  # Model runtime (in seconds)
  time_clock = toc(quiet = TRUE)
  time_taken = round(time_clock$toc - time_clock$tic) %>%
    seconds_to_period()
  
  # Option to constrain the model output
  if(!is.null(p$metrics_to_save)){
    
    # check if requested metrics are valid
    if(!all(p$metrics_to_save %in% unique(m$metric))){
      message("  > Requested output metrics do not align with model metrics: \n - ", paste(unique(m$metric),collapse='\n - '))
    }
    
    # filter metrics
    m = m %>%
      filter(!is.na(value),
             metric %in% p$metrics_to_save)
  }
  
  # Display model run time, if appropriate
  if (verbose != "none")
    message("  > Model runtime: ", time_taken)
  
  # ---- Close up ----
  # Combine inputs and outputs into one list
  results = list(yaml    = yaml$raw, 
                 input   = yaml$parsed, 
                 output  = m)
  
  return(results)
}

# ---------------------------------------------------------
# Extract state variables for each age group
# ---------------------------------------------------------
extract_state_variables = function(state, states, age_groups) {
  
  sapply(states, function(comp) {
    state[paste0(comp, "_", age_groups)]
  }, simplify = FALSE)
  
}

# ---------------------------------------------------------
# ODE model
# ---------------------------------------------------------
rsv_model = function(t, y, p){
  
  # Assert non-negativity (which may be caused by numerical errors in solver)
  states = pmax(y, 0)

  extract_states = p$compartments

  states  <- as.data.frame(matrix(states,nrow=p$n_age,ncol=length(extract_states),byrow=FALSE))
  names(states) <- extract_states

  # Pre-extract adult V-stage matrices (n_age × W) before entering 'with'.
  # These are accessible inside 'with' via the parent environment.
  # Element [a, j] = number of individuals in age group a at waning stage j.
  V1_mat <- as.matrix(states[, paste0("V1_", seq_len(p$W)), drop=FALSE])
  V2_mat <- as.matrix(states[, paste0("V2_", seq_len(p$W)), drop=FALSE])
  V3_mat <- as.matrix(states[, paste0("V3_", seq_len(p$W)), drop=FALSE])

  with(states, {
    
    # Seasonality: cosine wave with peak at `peak_day` (day of year),
    # raised to `seasonality_exponent` while preserving sign to allow
    # asymmetric peaks/troughs. Amplitude=0 disables seasonality.
    cval <- cos(2 * pi * (t - p$peak_day) / 365)
    seasonality_factor <- 1 + p$amplitude * sign(cval) * abs(cval)^p$seasonality_exponent

    # Which RSV season are we in (1 = season containing start_date)?
    season_nr = get_season_number(p$start_date + t - 1, p$start_date)

    # Per-season scalar on transmission rate.
    season_effect = c(1, p$season2_effect, p$season3_effect)
    season_scalar = if (season_nr >= 1 && season_nr <= length(season_effect))
                      season_effect[season_nr] else 1

    # ---- Force of infection ----
    # Effective infectious population: vaccinated streams (I_kv) contribute equally
    # to unvaccinated streams — vaccination does not reduce onward transmission.
    infectious_A = p$first_infection_infectiousness  * (I0 + I0v) +
                   p$second_infection_infectiousness * (I1 + I1v) +
                   p$third_infection_infectiousness  * (I2 + I2v) +
                   p$third_infection_infectiousness  * (I3 + I3v)

    FoI_tmp  = p$beta_A * season_scalar * seasonality_factor * p$contact_scalar *
               (infectious_A / p$population$population)
    lambda_A = rowSums(matrix(rep(FoI_tmp, p$n_age), ncol = p$n_age, byrow = TRUE) * p$contact_matrix)

    # ---- Vaccine immunity ----
    # Infant: residual susceptibility indexed by age group (age ≈ time since vaccination).
    # 1 = no protection, 0 = full protection.
    infant_vaccine_immunity <- 1 - p$infant_vacc_IE * unlist(p$infant_vaccine_rel_protection)

    # Adult: residual susceptibility indexed by waning stage j = 1..W.
    adult_vacc_immunity <- 1 - p$adult_vacc_IE * unlist(p$adult_vaccine_rel_protection)

    # VE against severity — convert overall (trial-reported) to conditional on infection.
    # VE_sev_cond = 1 - (1 - VE_hosp_overall) / (1 - VE_acq)
    infant_vacc_IE_hosp_cond <- 1 - (1 - p$infant_vacc_IE_hosp) / (1 - p$infant_vacc_IE)
    adult_vacc_IE_hosp_cond  <- 1 - (1 - p$adult_vacc_IE_hosp)  / (1 - p$adult_vacc_IE)

    # ---- Adult V-stage infection flows ----
    # For each tier k and waning stage j, infection flow from V_k_j[a] to E_kv[a]:
    #   flow[a,j] = susceptibility * prior_prot_k * adult_vacc_immunity[j] * lambda_A[a] * V_k_mat[a,j]
    # Computed as n_age × W matrices: sweep scales column j by adult_vacc_immunity[j],
    # then row-multiplication by the age-specific infection rate handles lambda_A.
    V1_infection_flow <- (p$susceptibility * p$prior_infection_protection  * lambda_A) *
                           sweep(V1_mat, 2, adult_vacc_immunity, "*")
    V2_infection_flow <- (p$susceptibility * p$prior_2infection_protection * lambda_A) *
                           sweep(V2_mat, 2, adult_vacc_immunity, "*")
    V3_infection_flow <- (p$susceptibility * p$prior_3infection_protection * lambda_A) *
                           sweep(V3_mat, 2, adult_vacc_immunity, "*")

    # Assign per-stage derivatives so the generic assembly loop below can find them
    for (.j in seq_len(p$W)) {
      assign(paste0("dV1_", .j), -V1_infection_flow[, .j])
      assign(paste0("dV2_", .j), -V2_infection_flow[, .j])
      assign(paste0("dV3_", .j), -V3_infection_flow[, .j])
    }

    # Total daily inflow into vaccinated-exposed compartments (sum over all stages)
    E1v_inflow <- rowSums(V1_infection_flow)
    E2v_inflow <- rowSums(V2_infection_flow)
    E3v_inflow <- rowSums(V3_infection_flow)

    # ---- Per-age hospitalisation risk for the current RSV season ----
    # Step function on the season-start dates (Aug 1 of each season_year).
    # Each row of p$p_hosp_A_by_season is one season's age-resolved p_hosp_A
    # vector, derived from that season's pooled hospital-burden split. The
    # biological motivation is strain turnover between seasons (severity
    # by age band can differ when the circulating virus changes); within a
    # season p_hosp_A is held constant. See age_relativity() for details.
    #
    # findInterval returns:
    #   0                                -> simulation date is before the first
    #                                        observed season; use the pooled-
    #                                        across-all-seasons fallback.
    #   k in [1, nrow(by_season)]        -> row k (most recent season-start <=
    #                                        current date). For dates after the
    #                                        last observed season this holds
    #                                        the most recent row (no extrapola-
    #                                        tion beyond observed data).
    current_date_num = as.numeric(p$start_date) + floor(t) - 1
    .s_idx = findInterval(current_date_num, p$p_hosp_A_season_starts_num)
    p_hosp_A_t = if (.s_idx >= 1) p$p_hosp_A_by_season[.s_idx, ] else p$p_hosp_A_fallback

    # ---- Incidence ----
    incidence_vacc_infant <- (p$susceptibility * infant_vaccine_immunity) * lambda_A * V0
    incidence_vacc_adult  <- E1v_inflow + E2v_inflow + E3v_inflow
    incidence_A <- p$susceptibility * lambda_A * (S0 + S1 + S2 + S3) +
                   incidence_vacc_infant + incidence_vacc_adult

    cases_A = p$p_confirm_A * incidence_A

    # Hospital admissions: vaccinated streams use reduced p_hosp (VE against severity)
    hospital_admissions_A =
      p$p_confirm_hosp_A * (
        p_hosp_A_t * (1/p$theta) * (I0 + I1 + I2 + I3) +                              # unvaccinated
        p_hosp_A_t * (1 - infant_vacc_IE_hosp_cond) * (1/p$theta) * I0v +             # infant vaccinated
        p_hosp_A_t * (1 - adult_vacc_IE_hosp_cond)  * (1/p$theta) * (I1v + I2v + I3v) # adult vaccinated
      )

    deaths = p$p_confirm_death * p$p_death * 1/p$mu * (H0 + H1 + H2 + H3)

    # Share of new infections in vaccinated individuals (all streams combined)
    incidence_prop_vacc = ifelse(incidence_A > 0,
                                 (incidence_vacc_infant + incidence_vacc_adult) / incidence_A,
                                 0)

    #---- Ordinary differential equations ----

    # Naive tier — infant vaccination stream (V0 → E0v, S0 → E0)
    dV0  <- -(p$susceptibility * infant_vaccine_immunity) * lambda_A * V0
    dS0  <- -p$susceptibility * lambda_A * S0
    dE0  <-  p$susceptibility * lambda_A * S0 - (1/p$gamma_A) * E0
    dE0v <-  incidence_vacc_infant             - (1/p$gamma_A) * E0v
    dI0  <-  (1/p$gamma_A) * E0  - (1/p$theta) * I0
    dI0v <-  (1/p$gamma_A) * E0v - (1/p$theta) * I0v
    # p_hosp_A_t already includes age effect; see age_relativity()
    dH0  <-  p_hosp_A_t * (1/p$theta) * I0 +
             p_hosp_A_t * (1 - infant_vacc_IE_hosp_cond) * (1/p$theta) * I0v -
             (1 - p$p_death) * (1/p$delta) * H0 - p$p_death * (1/p$mu) * H0
    dR0  <-  (1 - p_hosp_A_t) * (1/p$theta) * I0 +
             (1 - p_hosp_A_t * (1 - infant_vacc_IE_hosp_cond)) * (1/p$theta) * I0v +
             (1 - p$p_death) * (1/p$delta) * H0 - (1/p$omega_1) * R0
    dD0  <-  p$p_death * (1/p$mu) * H0

    # Tier 1 — adult vaccination stream (V1_j → E1v, S1 → E1)
    dS1  <- -p$susceptibility * p$prior_infection_protection * lambda_A * S1 + (1/p$omega_1) * R0
    dE1  <-  p$susceptibility * p$prior_infection_protection * lambda_A * S1 - (1/p$gamma_A) * E1
    dE1v <-  E1v_inflow - (1/p$gamma_A) * E1v
    dI1  <-  (1/p$gamma_A) * E1  - (1/p$theta) * I1
    dI1v <-  (1/p$gamma_A) * E1v - (1/p$theta) * I1v
    dH1  <-  p_hosp_A_t * (1/p$theta) * I1 +
             p_hosp_A_t * (1 - adult_vacc_IE_hosp_cond) * (1/p$theta) * I1v -
             (1 - p$p_death) * (1/p$delta) * H1 - p$p_death * (1/p$mu) * H1
    dR1  <-  (1 - p_hosp_A_t) * (1/p$theta) * I1 +
             (1 - p_hosp_A_t * (1 - adult_vacc_IE_hosp_cond)) * (1/p$theta) * I1v +
             (1 - p$p_death) * (1/p$delta) * H1 - (1/p$omega_2) * R1
    dD1  <-  p$p_death * (1/p$mu) * H1

    # Tier 2
    dS2  <- -p$susceptibility * p$prior_2infection_protection * lambda_A * S2 + (1/p$omega_2) * R1
    dE2  <-  p$susceptibility * p$prior_2infection_protection * lambda_A * S2 - (1/p$gamma_A) * E2
    dE2v <-  E2v_inflow - (1/p$gamma_A) * E2v
    dI2  <-  (1/p$gamma_A) * E2  - (1/p$theta) * I2
    dI2v <-  (1/p$gamma_A) * E2v - (1/p$theta) * I2v
    dH2  <-  p_hosp_A_t * (1/p$theta) * I2 +
             p_hosp_A_t * (1 - adult_vacc_IE_hosp_cond) * (1/p$theta) * I2v -
             (1 - p$p_death) * (1/p$delta) * H2 - p$p_death * (1/p$mu) * H2
    dR2  <-  (1 - p_hosp_A_t) * (1/p$theta) * I2 +
             (1 - p_hosp_A_t * (1 - adult_vacc_IE_hosp_cond)) * (1/p$theta) * I2v +
             (1 - p$p_death) * (1/p$delta) * H2 - (1/p$omega_3) * R2
    dD2  <-  p$p_death * (1/p$mu) * H2

    # Tier 3+
    dS3  <- -p$susceptibility * p$prior_3infection_protection * lambda_A * S3 + (1/p$omega_3) * R2 + (1/p$omega_4) * R3
    dE3  <-  p$susceptibility * p$prior_3infection_protection * lambda_A * S3 - (1/p$gamma_A) * E3
    dE3v <-  E3v_inflow - (1/p$gamma_A) * E3v
    dI3  <-  (1/p$gamma_A) * E3  - (1/p$theta) * I3
    dI3v <-  (1/p$gamma_A) * E3v - (1/p$theta) * I3v
    dH3  <-  p_hosp_A_t * (1/p$theta) * I3 +
             p_hosp_A_t * (1 - adult_vacc_IE_hosp_cond) * (1/p$theta) * I3v -
             (1 - p$p_death) * (1/p$delta) * H3 - p$p_death * (1/p$mu) * H3
    dR3  <-  (1 - p_hosp_A_t) * (1/p$theta) * I3 +
             (1 - p_hosp_A_t * (1 - adult_vacc_IE_hosp_cond)) * (1/p$theta) * I3v +
             (1 - p$p_death) * (1/p$delta) * H3 - (1/p$omega_4) * R3
    dD3  <-  p$p_death * (1/p$mu) * H3

    # Combine derivatives into a named vector.
    # The loop finds dV0, dV1_1..dV1_W, dV2_1..dV2_W, dV3_1..dV3_W,
    # dS0..dD3, dE0v..dI3v via get() in the current environment.
    derivatives <- matrix(0, p$n_age, length(extract_states))
    for (i_state in seq_along(extract_states)) {
      derivatives[, i_state] <- get(paste0("d", extract_states[i_state]))
    }

    derivatives <- c(derivatives)
    names(derivatives) <- paste0("d", rep(extract_states, each = p$n_age), "_", p$age_groups)

    names(incidence_A)          = paste0("incidence_A",          "_", p$age_groups)
    names(cases_A)              = paste0("cases_A",              "_", p$age_groups)
    names(hospital_admissions_A)= paste0("hospital_admissions_A","_", p$age_groups)
    names(deaths)               = paste0("deaths",               "_", p$age_groups)
    names(incidence_prop_vacc)  = paste0("incidence_prop_vacc",  "_", p$age_groups)
    names(seasonality_factor)   = "seasonality_factor"

    return(list(derivatives,
                incidence_A,
                cases_A,
                hospital_admissions_A,
                incidence_prop_vacc,
                deaths,
                seasonality_factor))

  }) # end with
  
  
}

# ---------------------------------------------------------
# Monthly ageing event
# ---------------------------------------------------------
# Fired by deSolve on the first day of each month. Three things happen:
#
# 1. DEMOGRAPHIC AGEING: for every compartment we move a fraction
#    `1/width_months` of each age bin into the next-older bin, so an
#    n-month-wide bin empties on a roughly n-month timescale. The oldest
#    bin (e.g. "80+y") has infinite width and never empties.
#
# 2. ADULT V-STAGE WANING: the adult vaccination waning chain is advanced
#    one step. Individuals in V_k_j move to V_k_{j+1} (j = 1..W-1);
#    those in V_k_W (fully waned) return to S_k and become eligible for
#    re-vaccination. This is done after demographic ageing so age-group
#    movement and stage advancement are independent.
#
# 3. VACCINATION:
#    a. Infant routine: newborns (0-1m) are split between S0 and V0
#       according to infant_vacc_coverage, which is non-zero only during
#       the configured infant_vaccination_start/end windows.
#    b. Infant catch-up: on infant_vaccination_catch_up_date, a fraction
#       infant_vaccination_catch_up_coverage of S0 in the listed age groups
#       is moved into V0.
#    c. Adult campaign: on each adult_vaccination_dates, a fraction
#       adult_vacc_coverage of S1/S2/S3 in adult_vaccination_agegroups
#       is moved into V1_1/V2_1/V3_1 (entering the waning chain at stage 1).
#
# Note: D* compartments are not age-shifted, so cumulative deaths retain
# the age group at time of death. Mortality is currently disabled in the
# configs (p_death = 0); revisit if age-at-observation tracking is needed.
ageing_event <- function(t, y, parms) {

  # Extract names and split into prefix and age group
  comp_names <- names(y)
  parts      <- str_match(comp_names, "^(.*)_(.+)$")
  prefixes   <- parts[, 2]
  age_labels <- parts[, 3]

  widths <- sapply(age_labels, bin_width_months)
  new_y  <- y  # copy to modify

  # Current event date (origin is parms$start_date)
  current_date <- parms$start_date + round(t)
  if (day(current_date) != 1) return(y)

  # ---- 1. Demographic ageing (all compartments) ----
  # Infant vaccination coverage for this month (non-zero only in season window)
  if (any(current_date >= ymd(parms$infant_vaccination_start) &
          current_date <= ymd(parms$infant_vaccination_end))) {
    infant_vacc_coverage <- parms$infant_vacc_coverage
  } else {
    infant_vacc_coverage <- 0
  }

  # Monthly births
  births_val <- with(parms$births, {
    idx <- match(current_date, date)
    if (!is.na(idx)) births[idx] else 0
  })

  for (pref in unique(prefixes)) {

    idx <- which(prefixes == pref)

    # Shift individuals from younger to older age bins
    for (k in seq_along(idx)) {
      i <- idx[k]
      if (k == length(idx)) next  # last bin never empties (infinite width)
      w    <- widths[i]
      frac <- if (is.finite(w)) 1/w else 0
      move        <- y[i] * frac
      new_y[i]         <- new_y[i]         - move
      new_y[idx[k+1]]  <- new_y[idx[k+1]] + move
    }

    # Replenish / zero the 0-1m bin
    first_label <- age_labels[idx[1]]
    if (first_label == "0-1m") {
      first_bin <- idx[1]
      if (pref == "S0") {
        new_y[first_bin] <- births_val * (1 - infant_vacc_coverage)
      } else if (pref == "V0") {
        new_y[first_bin] <- births_val * infant_vacc_coverage
      } else {
        new_y[first_bin] <- 0
      }
    }
  }

  # ---- 2. Adult V-stage waning advancement ----
  # Advance each tier's waning chain by one month.
  # Process stages from last to first to avoid overwriting values mid-loop.
  for (k in 1:3) {
    idx_Sk <- which(prefixes == paste0("S", k))

    for (j in parms$W:1) {
      idx_j <- which(prefixes == paste0("V", k, "_", j))
      if (j == parms$W) {
        # Last stage: waned individuals return to susceptible pool
        new_y[idx_Sk] <- new_y[idx_Sk] + new_y[idx_j]
      } else {
        # All other stages: advance to next stage
        idx_j1        <- which(prefixes == paste0("V", k, "_", j + 1))
        new_y[idx_j1] <- new_y[idx_j1] + new_y[idx_j]
      }
      new_y[idx_j] <- 0
    }
  }

  # ---- 3a. Infant catch-up campaign ----
  if (current_date == ymd(parms$infant_vaccination_catch_up_date)) {
    new_y2  <- new_y
    ind_S0  <- which(age_labels %in% parms$infant_vaccination_catch_up_agegroup & prefixes == "S0")
    ind_V0  <- which(age_labels %in% parms$infant_vaccination_catch_up_agegroup & prefixes == "V0")
    new_y[ind_S0] <- new_y2[ind_S0] * (1 - parms$infant_vaccination_catch_up_coverage)
    new_y[ind_V0] <- new_y[ind_V0]  + new_y2[ind_S0] * parms$infant_vaccination_catch_up_coverage
  }

  # ---- 3b. Adult vaccination campaign ----
  # On each campaign date, move adult_vacc_coverage fraction of S1/S2/S3
  # in the eligible age groups into V1_1/V2_1/V3_1 (waning chain stage 1).
  if (any(current_date == ymd(parms$adult_vaccination_dates))) {
    for (k in 1:3) {
      idx_Sk  <- which(prefixes == paste0("S", k) &
                       age_labels %in% parms$adult_vaccination_agegroups)
      idx_Vk1 <- which(prefixes == paste0("V", k, "_1") &
                       age_labels %in% parms$adult_vaccination_agegroups)
      to_vacc          <- new_y[idx_Sk] * parms$adult_vacc_coverage
      new_y[idx_Sk]    <- new_y[idx_Sk]  - to_vacc
      new_y[idx_Vk1]   <- new_y[idx_Vk1] + to_vacc
    }
  }

  return(new_y)
}

# Helper: compute bin width (months)
bin_width_months <- function(age_label) {
  # age_label is only the suffix like "7-8m", not the full name!
  if (str_detect(age_label, "\\+y$")) return(Inf)
  if (str_detect(age_label, "m$")) {
    # Capture only the two numbers around the dash
    nums <- as.numeric(unlist(str_match(age_label, "^(\\d+)-(\\d+)m$")[,2:3]))
    return(nums[2] - nums[1])
  }
  if (str_detect(age_label, "y$")) {
    nums <- as.numeric(unlist(str_match(age_label, "^(\\d+)-(\\d+)y$")[,2:3]))
    return((nums[2] - nums[1]) * 12)
  }
  stop(paste("Cannot parse age group:", age_label))
}


# ---------------------------------------------------------
# Seed initial state at t = 0 from yaml-defined per-age-group conditions
# ---------------------------------------------------------
# For each age group `g` (with population N_g), the yaml block
# `initial_conditions` specifies a tier (0, 1 or 3) and the proportions to
# place in S, E and I of that tier. The compartment populations are then:
#
#   S{tier}_g = prop_S * N_g
#   E{tier}_g = prop_E * N_g
#   I{tier}_g = init_inf * prop_I_unit * N_g
#   R{tier}_g = N_g - (S + E + I)        # remainder goes to recovered
#   all other compartments = 0
#
# Notes:
#   - Tier 2 (S2/E2/I2/R2) is intentionally not seeded; it fills dynamically
#     via flow R1 -> S2 during simulation.
#   - V0 is not seeded; newborns enter V0 via the monthly ageing event when
#     infant vaccination is active.
#   - V1_j/V2_j/V3_j (adult waning chain) all start at zero; they are filled
#     by the adult campaign logic in ageing_event on adult_vaccination_dates.
#   - E0v/I0v and E1v-E3v/I1v-I3v (vaccinated exposed/infectious streams)
#     start at zero and fill via ODE infection flows once V compartments are
#     populated.
# ---------------------------------------------------------
initiate_epidemic = function(p, verbose){

  # Total model population (across all age groups)
  n = p$population$population %>% sum()

  # Throw an error if trying to create trivial number of people
  if (is.null(n) || n == 0)
    stop("Attemping to create zero people")

  # Display (or not) how many people we are creating
  if (verbose != "none")
    message("  > Initiating population of ", thou_sep(n))

  # All compartments start at zero...
  initial_conditions = setNames(
    lapply(p$compartments, function(comp) rep(0, p$n_age)),
    p$compartments
  )

  # ...then populate per the yaml-defined initial conditions
  ic = rbindlist(p$initial_conditions)
  ic = ic[match(p$age_groups, age_group)]      # align ordering with p$age_groups

  # Sanity check: every age group must have an initial condition row
  if (any(is.na(ic$tier)))
    stop("Missing initial_conditions entries for age group(s): ",
         paste(p$age_groups[is.na(ic$tier)], collapse = ", "))

  # Apply the seeding tier-by-tier
  pop = p$population$population
  for (i in seq_len(p$n_age)) {
    tier  = ic$tier[i]
    S_amt = ic$prop_S[i]                  * pop[i]
    E_amt = ic$prop_E[i]                  * pop[i]
    I_amt = p$init_inf * ic$prop_I_unit[i] * pop[i]
    R_amt = pop[i] - S_amt - E_amt - I_amt

    initial_conditions[[paste0("S", tier)]][i] = S_amt
    initial_conditions[[paste0("E", tier)]][i] = E_amt
    initial_conditions[[paste0("I", tier)]][i] = I_amt
    initial_conditions[[paste0("R", tier)]][i] = R_amt
  }

  if (any(unlist(initial_conditions) < 0)){
    stop("Negative initial values. Stopping the simulation!")
  }

  # Combine initial state variables into a named vector
  states = unlist(lapply(p$compartments, function(comp) {
    setNames(initial_conditions[[comp]], paste0(comp, "_", p$age_groups))
  }))

  return(states)
}

# ---------------------------------------------------------
# Prepare final model output
# ---------------------------------------------------------
format_output = function(out_df, p) {
  
  #---- Susceptibles ----
  # Select susceptible columns
  S_cols = c(outer(c("S0_", "S1_", "S2_", "S3_"), p$age_groups, paste0)) %>% as.vector()
  
  # Summarise susceptibles by age group and time
  S_df = out_df %>% pivot_longer(cols = all_of(S_cols),
                                 names_to = "compartment",
                                 values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val, na.rm =TRUE),
              .groups = "drop") %>%
    mutate(metric = "susceptibles",
           variant = NA_character_)  # Susceptible compartment not disaggregated by variant
  
  #---- New infections (incidence) ----
  # Select incidence columns
  new_A_cols = c(outer("incidence_A_", p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise incidence by variant, age group and time
  new_A_df = out_df %>% pivot_longer(cols = all_of(new_A_cols),
                                     names_to = "compartment",
                                     values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "new_infections",
           variant = "A") 
  
  #---- Proportion of incidence in vaccinated individuals ----
  # The per-age fraction of new infections occurring in vaccinated individuals.
  # Used in results_evaluation to split burden into vaccinated / unvaccinated
  # streams for the RespiCompass submission format.
  incidence_prop_vacc_cols = c(outer("incidence_prop_vacc_", p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise incidence in vaccinated by variant, age group and time
  incidence_prop_vacc_df = out_df %>% pivot_longer(cols = all_of(incidence_prop_vacc_cols),
                                       names_to = "compartment",
                                       values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "incidence_prop_vacc",
           variant = "A") 
  
  #---- Confirmed cases (ILI+) ----
  # Select incidence columns
  cases_A_cols = c(outer("cases_A_", p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise incidence by variant, age group and time
  cases_A_df = out_df %>% pivot_longer(cols = all_of(cases_A_cols),
                                       names_to = "compartment",
                                       values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "cases",
           variant = "A") 
  
  
  #---- Latent (pre-infectious) ----
  # Select latent (pre-infectious) columns — includes vaccinated exposed streams (E_kv)
  E_A_cols = c(outer(c("E0_", "E0v_", "E1_", "E1v_", "E2_", "E2v_", "E3_", "E3v_"), p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise latent (pre-infectious) by age group and time
  E_A_df = out_df %>%  pivot_longer(cols = all_of(E_A_cols),
                                    names_to = "compartment",
                                    values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "latent",
           variant = "A") 
  
  #---- Infectious ----
  # Select infectious columns — includes vaccinated infectious streams (I_kv)
  I_A_cols = c(outer(c("I0_", "I0v_", "I1_", "I1v_", "I2_", "I2v_", "I3_", "I3v_"), p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise infectious by age group and time
  I_A_df = out_df %>%  pivot_longer(cols = all_of(I_A_cols),
                                    names_to = "compartment",
                                    values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "infectious",
           variant = "A") 
  
  #---- Hospital occupancy----
  # Select hospitalised columns
  H_cols = c(outer(c("H0_", "H1_", "H2_", "H3_"), p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise hospitalised by age group and time
  H_df = out_df %>% pivot_longer(cols = all_of(H_cols),
                                 names_to = "compartment",
                                 values_to = "val") %>%
    
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "hospital_occupancy",
           variant = NA_character_) # Hospitalised compartment not disaggregated by variant
  
  #---- Hospital admissions ----
  # Select hospital admissions columns
  admit_A_cols = c(outer("hospital_admissions_A_", p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise hospital admissions by variant, age group and time
  admit_A_df = out_df %>%  pivot_longer(cols = all_of(admit_A_cols),
                                        names_to = "compartment",
                                        values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "hospital_admissions",
           variant = "A") 
  
  #---- Recovered ----
  # Select recovered columns
  R_cols = c(outer(c("R0_", "R1_", "R2_", "R3_"), p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise receovered by age group and time
  R_df = out_df %>% pivot_longer(cols = all_of(R_cols),
                                 names_to = "compartment",
                                 values_to = "val") %>%
    
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "recovered",
           variant = NA_character_) # Recovered compartment not disaggregated by variant
  
  #---- Deceased ----
  # Select deceased columns (cumulative deaths)
  D_cols = c(outer(c("D0_", "D1_", "D2_", "D3_"), p$age_groups, paste0)) %>%
    as.vector()
  
  D_df = out_df %>% pivot_longer(cols = all_of(D_cols),
                                 names_to = "compartment",
                                 values_to = "val") %>%
    
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "deceased",
           variant = NA_character_) # Deceased compartment not disaggregated by variant
  
  #---- Vaccinated ----
  # Select vaccinated columns: infant stream (V0) + all adult waning stages (V1_j, V2_j, V3_j)
  V_cols = c(
    outer("V0_", p$age_groups, paste0),
    outer(paste0("V1_", seq_len(p$W), "_"), p$age_groups, paste0),
    outer(paste0("V2_", seq_len(p$W), "_"), p$age_groups, paste0),
    outer(paste0("V3_", seq_len(p$W), "_"), p$age_groups, paste0)
  ) %>% as.vector()
  
  V_df = out_df %>% pivot_longer(cols = all_of(V_cols),
                                 names_to = "compartment",
                                 values_to = "val") %>%
    
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "vaccinated",
           variant = NA_character_) # Deceased compartment not disaggregated by variant
  
  
  #---- Deaths ----
  # Select death columns
  death_cols = c(outer("deaths_", p$age_groups, paste0)) %>%
    as.vector()
  
  # Summarise hospital admissions by variant, age group and time
  deaths_df = out_df %>%  pivot_longer(cols = all_of(death_cols),
                                       names_to = "compartment",
                                       values_to = "val") %>%
    mutate(age_group = sub(".*_", "", compartment),
           age_group = factor(age_group, levels = p$age_groups)) %>%  # Extract age group number from compartment name
    group_by(time, age_group) %>%
    summarise(value = sum(val),
              .groups = "drop") %>%
    mutate(metric = "deaths")
  
  
  #---- Seasonality ----
  # NB: This is an input, not an output, but is helpful for visualisation
  seasonality = out_df %>% select(time, seasonality_factor) %>%
    rename(value = seasonality_factor) %>%
    mutate(metric = "seasonality",
           age_group = NA_character_,
           variant = NA_character_) %>%
    select(time, age_group, value, metric, variant)
  
  # Compile output
  m = bind_rows(S_df,
                new_A_df,
                incidence_prop_vacc_df,
                cases_A_df,
                E_A_df,
                I_A_df,
                admit_A_df,
                H_df,
                R_df,
                D_df,
                V_df,
                deaths_df,
                seasonality)
  
  #---- Total living population ----
  # Define the subset of metrics to be summed
  pop_metrics = c("susceptibles", "latent", "infectious", "hospital_occupancy", "recovered", "vaccinated")  
  
  # Compute 'total' as the sum of the selected metrics within each age group
  pop_total = m %>%
    filter(metric %in% pop_metrics) %>%  # Keep only the relevant metrics
    group_by(time, age_group) %>%  # Sum within each age group
    summarise(metric = "total", value = sum(value, na.rm = TRUE), .groups = "drop")
  
  overall_total = pop_total %>% group_by(time) %>% summarise(value = sum(value))
  
  m = bind_rows(m, pop_total) %>%
    setDT()
  
  return(m)
}

# ---------------------------------------------------------
# Load contact matrices
# ---------------------------------------------------------
load_contacts = function(o) {
  load(o$contact_matrices)
  return(contact_all)
}

# ---------------------------------------------------------
# Account for relative susceptibility by age
# ---------------------------------------------------------
age_relativity = function(p){
  
  # Initialise susceptibility by age group
  p$susceptibility = p$age_group_map %>% 
    mutate(susceptibility = 1) %>% 
    select(-smaller_group) %>%
    group_by(larger_group) %>%
    slice(1) %>% # Keep just 'larger groups'
    ungroup() %>%         
    # Read in relative susceptibility by age group 
    mutate(susceptibility = case_when(larger_group %in% c("0-1m", "1-2m", "2-3m") ~ p$rel_sus_a,
                                      larger_group %in% c("3-4m", "4-5m", "5-6m") ~ p$rel_sus_b,
                                      larger_group %in% c("65+y") ~ p$rel_sus_c,
                                      TRUE ~ susceptibility)) %>%
    # Map susceptibility to smaller age groups
    right_join(p$age_group_map, by = "larger_group") %>%
    rename(age_group = smaller_group) %>%
    mutate(age_group = factor(age_group, levels = p$age_groups)) %>%
    arrange(age_group) %>%
    select(susceptibility) %>%
    unlist()
  
  names(p$susceptibility) = NULL

  # ---- Season-specific age relativity of p_hosp_A ----
  #
  # p_hosp_A is the per-infection probability of being hospitalised, by age
  # group. It varies by season as different RSV strains can circulate and 
  # the relative severity across age bands can shift. 
  # We therefore compute ONE p_hosp_A vector PER SEASON, using that season's 
  # full burden pooled together.
  #
  # Season convention (matches get_season_number()): an RSV season runs
  # from Aug 1 of year Y to Jul 31 of year Y+1, labelled here by Y (the
  # "season_year"). 

  # Aggregate burden by (season_year, age_group) and pivot wide.
  # season_year = calendar year if month >= 8, else calendar year - 1.
  burden_wide = p$burden %>%
    mutate(season_year = if_else(month(date) >= 8, year(date), year(date) - 1)) %>%
    group_by(season_year, age_group) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = age_group, values_from = value, values_fill = 0) %>%
    arrange(season_year)

  # Make sure the three reporting age groups we need exist as columns.
  # If a band is entirely absent from p$burden (data source dropped it, or
  # the band mapping is wrong) pivot_wider will simply not create the
  # column and the arithmetic below would silently break on NULL. Fill
  # missing bands with zero so the downstream fallback path triggers.
  for (band in c("0-3m", "3-6m", "6-12m")) {
    if (!band %in% names(burden_wide)) {
      warning("age_relativity(): reporting band '", band, "' missing from ",
              "p$burden — filling with zeros and falling back to ratio = 1 ",
              "(no relative adjustment).")
      burden_wide[[band]] = 0
    }
  }

  # Per-season ratios
  ratio1_s = burden_wide$`0-3m` / burden_wide$`3-6m`
  ratio2_s = 2 * burden_wide$`0-3m` / burden_wide$`6-12m`

  # Pooled (across all seasons) ratios — used as the fallback for any
  # season whose own ratio is degenerate, and as the simulation-time
  # fallback for dates before the first observed season.
  total_0_3m  = sum(burden_wide$`0-3m`,  na.rm = TRUE)
  total_3_6m  = sum(burden_wide$`3-6m`,  na.rm = TRUE)
  total_6_12m = sum(burden_wide$`6-12m`, na.rm = TRUE)
  ratio1_fallback = total_0_3m / total_3_6m
  ratio2_fallback = 2 * total_0_3m / total_6_12m

  # A ratio is "bad" if dividing by it later would produce a non-finite
  # value. That covers NaN (0/0 from a missing season), Inf (x/0 from a
  # missing denominator band), and 0 (which is finite but would become
  # Inf inside compute_p_hosp_A_row when we do `... / ratio`).
  ratio_is_bad = function(r) !is.finite(r) | r == 0

  # If the pooled fallback is itself bad (entire denominator band empty
  # across every season, or numerator '0-3m' empty too), fall back to
  # ratio = 1 — the 3-6m / 6-12m bands then use the same formula as 0-3m
  # (no relative adjustment). Warn so the user knows the data is degenerate.
  if (ratio_is_bad(ratio1_fallback)) {
    warning("age_relativity(): pooled '0-3m'/'3-6m' burden ratio is ",
            "degenerate (", ratio1_fallback, ") across all seasons in ",
            "p$burden — falling back to ratio1 = 1.")
    ratio1_fallback = 1
  }
  if (ratio_is_bad(ratio2_fallback)) {
    warning("age_relativity(): pooled '0-3m'/'6-12m' burden ratio is ",
            "degenerate (", ratio2_fallback, ") across all seasons in ",
            "p$burden — falling back to ratio2 = 1.")
    ratio2_fallback = 1
  }

  # Replace per-season bad ratios (NaN, Inf, or 0) with the pooled fallback
  ratio1_s[ratio_is_bad(ratio1_s)] = ratio1_fallback
  ratio2_s[ratio_is_bad(ratio2_s)] = ratio2_fallback
  
  # NOTE: The two lines below collapse per-season ratios to the pooled fallback
  # (x/x = 1 for any finite non-zero x, so ratio1_s = ratio1_fallback everywhere).
  # This is intentional: per-season variation is retained structurally but the
  # current default uses the pooled estimate across all seasons for robustness.
  # Remove these two lines to activate per-season age-relative hospitalisation.
  ratio1_s = ratio1_fallback * ratio1_s/ratio1_s
  ratio2_s = ratio2_fallback * ratio2_s/ratio2_s

  # Helper: build a length-n_age p_hosp_A vector from a (ratio1, ratio2)
  # pair. Draws fresh rnorm() values per call — see TODO below.
  compute_p_hosp_A_row = function(ratio1, ratio2) {
    p$age_group_map %>%
      mutate(p_hosp_A = p$p_hosp_A) %>% # Baseline = oldest age-group hosp risk
      select(-smaller_group) %>%
      group_by(larger_group) %>%
      slice(1) %>%   # Keep just 'larger groups'
      ungroup() %>%
      # TODO: the rnorm() draws below inject fresh randomness into p_hosp_A
      # on every call, so calibration is not reproducible for the same fit.
      # Either make this deterministic (use mean = 0.25) or route through
      # the yaml `uncertainty:` block so the draw is logged and seeded.
      mutate(p_hosp_A = case_when(
        larger_group %in% c("0-3m")  ~ p_hosp_A * p$rel_hosp_a_A,
        larger_group %in% c("3-6m")  ~ p_hosp_A * p$rel_hosp_a_A * (1 + rnorm(1, 0.25, 0.03)) / ratio1,
        larger_group %in% c("6-12m") ~ p_hosp_A * p$rel_hosp_a_A * (1 + rnorm(1, 0.25, 0.03)) / ratio2,
        larger_group %in% c("1-5y")  ~ p_hosp_A * p$rel_hosp_b_A,
        larger_group %in% c("65+y")  ~ p_hosp_A * p$rel_hosp_c_A,
        TRUE ~ p_hosp_A)) %>%
      # Map hospitalisation risk back to fine age groups
      right_join(p$age_group_map, by = "larger_group") %>%
      rename(age_group = smaller_group) %>%
      mutate(age_group = factor(age_group, levels = p$age_groups)) %>%
      arrange(age_group) %>%
      pull(p_hosp_A)
  }

  # Per-season matrix: rows = seasons (sorted by season_year ascending),
  # cols = age groups. Rownames are the season_year (e.g. "2024" for the
  # 2024-25 season) for traceability in any [0, 1] check error messages.
  n_s = nrow(burden_wide)
  p$p_hosp_A_by_season = matrix(0, nrow = n_s, ncol = p$n_age,
                                dimnames = list(burden_wide$season_year, p$age_groups))
  for (s in seq_len(n_s)) {
    p$p_hosp_A_by_season[s, ] = compute_p_hosp_A_row(ratio1_s[s], ratio2_s[s])
  }
  
  # Numeric season-start dates (Aug 1 of season_year, in days-since-epoch)
  # parallel to the rows of the matrix, for fast findInterval lookups in
  # the ODE.
  p$p_hosp_A_season_starts_num = as.numeric(
    as.Date(paste0(burden_wide$season_year, "-08-01")))

  # Pooled-across-all-seasons fallback, used when the simulated date is
  # before the first observed season (e.g. burn-in).
  p$p_hosp_A_fallback = compute_p_hosp_A_row(ratio1_fallback, ratio2_fallback)

  # Hard check: every p_hosp_A entry must be a valid probability in [0, 1].
  # With the fallback ladder above the ratios going into compute_p_hosp_A_row
  # are always finite and non-zero, so the matrix can't contain Inf/NaN.
  # The remaining failure modes are configuration-driven: e.g. the yaml
  # combination p_hosp_A * rel_hosp_*_A exceeds 1, or rnorm() flips sign.
  # Stop with a pointer to the offending (season, age_group) cell so the
  # user can fix the yaml rather than silently feed a non-probability into
  # the ODE.
  bad_cells = which(p$p_hosp_A_by_season < 0 | p$p_hosp_A_by_season > 1,
                    arr.ind = TRUE)
  if (nrow(bad_cells) > 0) {
    bad_seasons = rownames(p$p_hosp_A_by_season)[bad_cells[, "row"]]
    bad_ages    = p$age_groups[bad_cells[, "col"]]
    bad_vals    = p$p_hosp_A_by_season[bad_cells]
    stop("age_relativity(): p_hosp_A outside [0, 1] in ", nrow(bad_cells),
         " (season, age_group) cell(s). First offender: season ",
         bad_seasons[1], ", age_group ", bad_ages[1],
         ", p_hosp_A = ", signif(bad_vals[1], 4),
         ". Check yaml: p_hosp_A * rel_hosp_*_A must produce probabilities ",
         "in [0, 1] (e.g. for 0-3m, p_hosp_A * rel_hosp_a_A <= 1).")
  }
  if (any(p$p_hosp_A_fallback < 0 | p$p_hosp_A_fallback > 1)) {
    bad = which(p$p_hosp_A_fallback < 0 | p$p_hosp_A_fallback > 1)
    stop("age_relativity(): pooled-fallback p_hosp_A outside [0, 1] for ",
         "age_group(s): ", paste(p$age_groups[bad], collapse = ", "),
         ". Check baseline p_hosp_A and rel_hosp_*_A in the yaml.")
  }
  
  # Keep p$p_hosp_A as a scalar-vector default (the pooled fallback) for
  # any code path that reads it directly. The ODE no longer reads this; it
  # looks up p$p_hosp_A_by_season / p$p_hosp_A_fallback per time step.
  p$p_hosp_A = p$p_hosp_A_fallback
  names(p$p_hosp_A) = NULL


  # Format age-dependent mortality risk, given hospitalisation
  # Initialise mortality risk by age group
  p$p_death = p$age_group_map %>% 
    mutate(p_death = p$p_death) %>% # Baseline is mortality risk of oldest age group
    select(-smaller_group) %>%
    group_by(larger_group) %>%
    slice(1) %>% # Keep just 'larger groups'
    ungroup() %>%         
    # Read in relative mortality risk  by age group
    mutate(p_death = case_when(larger_group %in% c("0-3m") ~ p_death * p$rel_death_a,
                               larger_group %in% c("3-6m") ~ p_death * p$rel_death_b,
                               larger_group %in% c("65+y") ~ p_death * p$rel_death_c,
                               TRUE ~ p_death)) %>%
    # Map mortality risk to smaller age groups
    right_join(p$age_group_map, by = "larger_group") %>%
    rename(age_group = smaller_group) %>%
    mutate(age_group = factor(age_group, levels = p$age_groups)) %>%
    arrange(age_group) %>%
    select(p_death) %>%
    unlist()
  
  names(p$p_death) = NULL
  
  return(p)
}


# ------------------------------------------------------------
# get_season_number()
#
# Determines which epidemiological season a date belongs to.
# Seasons are defined as running from August 1 to July 31.
# The season containing `start_date` is defined as Season 1.
#
# Arguments:
#   date        - A Date (or vector of Dates) to classify
#   start_date  - A Date defining Season 1
#
# Returns:
#   Numeric season number (can be vectorised)
# ------------------------------------------------------------
get_season_number <- function(date, start_date) {
  
  date <- as.Date(date)
  start_date <- as.Date(start_date)
  
  # Season year starts on August 1
  season_year <- function(d) {
    ifelse(month(d) >= 8, year(d), year(d) - 1)
  }
  
  base_season_year <- season_year(start_date)
  current_season_year <- season_year(date)
  
  return(current_season_year - base_season_year + 1)
}
