########################################################## #
# CALIBRATION
#
# Fit parameters to a set of epidemiological targets, over time 
# 
# Uses adaptive sampling by comparing likelihood of model
# output given the observed target data
#
# Function output is a set of best-fitting parameter sets and 
# associated model output. 
########################################################## #

# -------------------------------------------------------- -
# Parent function for model calibration process  ----
# -------------------------------------------------------- -
run_calibration = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Calibrating model")
  
  # ---- Load input and data ----
  
  # Initiate fit list and perform a few checks on input .yaml
  fit = setup_calibration(o)
  
  # Load data - see load_data.R
  fit = load_data(o, fit)
  
  # specify reference period, to limit the time horizon
  o$limit_n_days <- as.numeric(diff(range(fit$data$date, na.rm = TRUE))+1)
  
  # specify reference metrics, to limit the model output
  o$fit_metrics <- get_fit_metrics(parse_yaml(o, scenario = "baseline")$parsed)
  
  # Iterate through adaptive sampling rounds (including initial 'r0' step)
  for (r_val in 0 : fit$input$adaptive_sampling$rounds) {
    message(" - Adaptive sampling round ", r_val)
    
    # (Re)sample parameter sets for this sampling round
    param_ids = sample_parameters(o, fit, r_val)
    
    # Update round index string
    r_idx = paste0("r", r_val)
    
    # Check if previous results exists
    sets_df = try_load(o$pth$fitting, paste0(r_idx, "_samples"), throw_error = FALSE)
    
    # Simulate parameters, if not existing or explicitly requested
    if (is.null(sets_df) || !o$plot_only){
      
      # Simulate these parameter samples
      simulate_parameters(o, fit, r_idx, param_ids)
      
      # Calculate likelihood
      likelihood(o, fit, r_idx)
      
      # Compile parameter sets
      compile_samples(o, r_val)
      
    } else {
      message("  > Reusing previous fit")
    }
    
    message("  > Plotting performance")
    
    # Fitting plots inc. by age group (see plotting.R)
    plot_best_samples(o, fit, fig_name = paste0("Model_fit_", fit$opts$country), round_idx = r_idx)
    
  }
  
  # Save final result
  save_calibration(o, fit)
  
}

# ---------------------------------------------------------------- -
# Initiate fit list and perform a few checks on input yaml  ----
# ---------------------------------------------------------------- -
setup_calibration = function(o) {
  
  # Load parsed parameters from yaml file
  p = parse_yaml(o, scenario = "baseline")$parsed
  
  # Shorthand for calibration options
  opts = p$calibration_options
  
  # ---- Sanity checks on input yaml ----
  
  # Metrics to be fitted over some time frame
  fit_metrics = get_fit_metrics(p)
  fit_days    = get_fit_days(p)
  
  # Throw error if n_days from yaml file is not at least max of fit_days
  if (p$n_days < fit_days)
    stop("Model calibration requires n_days >= ", fit_days)
  
  # ---- Construct list for fitting details ----
  
  # Initiate fit list (with easy access to calibration options)
  fit = list(input = p, opts = opts)
  
  # Extract parameters to calibrate
  fit_df = list2dt(fit$input$calibration_parameters) %>%
    mutate(lower = as.numeric(lower), 
           upper = as.numeric(upper))
  
  # Shorthand for fitted parameters
  fit$params = fit_df$param
  
  # Bounds of all parameters to be fitted, as matrix
  fit$bounds = as.matrix(fit_df[, .(lower, upper)])
  
  # Population scalar
  if(!opts$country == "none"){
    # Look up full country name from the RespiCompass country list
    opts$country_name = o$countries_df %>%
      filter(iso2_code == opts$country) %>%
      pull(country)
    pop_data = read.csv(o$pop_url, fileEncoding = "UTF-8-BOM") %>%
      normalise_iso2() %>%                 # Eurostat 'EL' -> ISO-2 'GR' (Greece)
      filter(country == opts$country) %>%
      remap_age_groups(o$respicompass_age_map)
    
  } 
  
  return(fit)
}

# -------------------------------------------------------- -
# Save final file with best parameter sets ----
# -------------------------------------------------------- -
save_calibration = function(o, fit) {
  
  # ---- Best simulated parameter set ----
  # Load all simulated samples
  all_samples = try_load(o$pth$fitting, "rx_samples")
  
  # Sort parameter set according to the likelihood
  samples_sorted = all_samples %>%
    arrange(desc(likelihood))  %>% # sort by likelihood
    unique()
  
  fit$best_simulated = samples_sorted %>%
    select(all_of(fit$params))
  
  fit$likelihoods = samples_sorted %>%
    select(-all_of(fit$params))
  
  # Store best of all as list for consistency 
  fit$best_one = as.list(fit$best_simulated[1,])
  
  # Save this final result
  saveRDS(fit, file = paste0(o$pth$fitting, "fit_result.rds"))
}

# -------------------------------------------------------- -
# Load final fit file and select 'best' parameter set ----
# -------------------------------------------------------- -
load_calibration = function(o, index, ...) {
  
  # Load model fitting result - throw an error if it doesn't exist
  err_msg = " > Cannot find a fitting file for this analysis - have you run step 1?"
  fit_result = try_load(o$pth$fitting, "fit_result", msg = err_msg, throw_error = FALSE, ...)
  
  if(is.null(fit_result)){
    message(err_msg)
    message(" > Use default values")
  } else{
    fit_result <- fit_result$best_simulated[index,] %>% 
      as.list()
  }
  
  return(fit_result)
}

# -------------------------------------------------------- -
# Generate a set of parameter sets to simulate ----
# -------------------------------------------------------- -
sample_parameters = function(o, fit, r_val) {
  
  opts = fit$input$adaptive_sampling
  
  # First iteration is simple: sample initial points
  if (r_val == 0){
    
    message("  > Sampling initial parameters")
    log_bounds = log10(fit$bounds)
    if (any( is.na(log_bounds) | is.infinite(log_bounds) | is.nan(log_bounds) )){
      stop("Logarithm of bounds is NA, NaN, or Inf. Stopping the simulations!")
    }
    samples = tgp::lhs(opts$init_samples, log_bounds) %>%
      10**. %>% # exponentiate back from log space to lin space
      as_named_dt(fit$params)
  }
  # On subsequent iterations, weight, resample and perturb
  if (r_val > 0) {
    message(paste0("  > Sampling parameters for round ", r_val))
    
    # Set bounds
    bounds_wide = bind_cols(param = fit$params, fit$bounds) %>%
      pivot_wider(names_from = param, values_from = c(lower, upper),
                  names_glue = "{param}_{.value}") %>%
      setDT()
    
    # Load samples from previous round
    round_path    = paste0("r", r_val-1, "_samples")
    round_samples = try_load(o$pth$fitting, round_path)
    
    # ---- Weight by likelihood ----
    #
    # Parameter sets are resampled (with replacement) in proportion to how well
    # they fitted, so the next round concentrates near the better ones. Two
    # things can go wrong:
    #   * weights nearly EQUAL  -> no selection; the round is just random search
    #   * one weight DOMINATES  -> all draws are copies of a single set, the
    #                              population collapses to one point and the
    #                              remaining rounds cannot explore
    # The effective sample size, ESS = 1 / sum(w^2), measures which regime we are
    # in: ESS = n when all weights are equal, ESS = 1 when one dominates. (It is
    # the reciprocal of the probability that two independent draws pick the same
    # set.) We temper the weights, w ∝ exp(alpha * L), choosing alpha in (0, 1]
    # so ESS meets a target: alpha = 1 uses the likelihood as-is, smaller alpha
    # softens the differences to preserve diversity. alpha is never raised above
    # 1, so we never claim more information than the likelihood actually holds.
    #
    # NB subtracting max(L) before exponentiating is essential: log-likelihoods
    # are large and negative, and exp() of them underflows to 0 for every sample,
    # which would make all weights 0 and the resampling weights NaN. Subtracting
    # the max divides every weight by the same constant, so ratios are unchanged.
    L <- round_samples$likelihood
    if (all(!is.finite(L)))
      stop("All likelihoods are non-finite in round ", r_val - 1,
           " - cannot resample. Check the model output and the likelihood.")
    # Demote failed runs to well below the worst surviving set (kept, but
    # effectively never resampled) rather than to an arbitrary fixed constant.
    L[!is.finite(L)] <- min(L[is.finite(L)], na.rm = TRUE) - 10

    ess_at     <- function(a) { w <- exp(a * (L - max(L))); w <- w / sum(w); 1 / sum(w^2) }
    ess_target <- 0.5 * length(L)   # conventional SMC choice: keep half the sets effective

    alpha <- if (ess_at(1) >= ess_target) 1 else
             uniroot(function(a) ess_at(a) - ess_target, c(1e-8, 1))$root

    w <- exp(alpha * (L - max(L)))

    # Report the diagnostic: ESS at alpha = 1 says how much the data actually
    # discriminated between these parameter sets. ESS ~ n means this round
    # learned nothing (e.g. an observation model too dispersed to be informative).
    message(sprintf("  > ESS(untempered) = %.1f / %d | alpha = %.3g | ESS(used) = %.1f",
                    ess_at(1), length(L), alpha, ess_at(alpha)))

    # Resample the required number of sets (with replacement, so better-fitting
    # sets are drawn repeatedly and the population concentrates on them).
    resampled = round_samples %>%
      mutate(weight_norm = w / sum(w)) %>%
      slice_sample(n = opts$init_samples,
                   replace = TRUE,
                   weight_by = weight_norm) %>%
      select(-c(param_id, round, likelihood, weight_norm))

    # ---- Adaptive perturbation kernel ----
    #
    # Resampling produces DUPLICATES, so the jitter is what turns them back into
    # distinct candidates: it is the only step that generates new points to test,
    # and its width is the search radius.
    #
    # A fixed percentage step (previously sd = 0.02) is the wrong radius, because
    # a fixed % of a parameter's VALUE says nothing about how uncertain we are
    # about that parameter. Measured on real output, 2% was ~30% of the plausible
    # range for rel_hosp_a_A but only ~0.4% of it for k - leaving k effectively
    # frozen while rel_hosp_a_A regularly overshot its bounds.
    #
    # Instead take the step size from the SPREAD of the resampled particles for
    # each parameter: that spread is the current estimate of its uncertainty, so
    # it is the right scale to explore at. This also anneals automatically (wide
    # while particles are scattered, narrow once they concentrate) and correctly
    # gives an already-converged parameter a small step.
    #
    # Spread is measured on the LOG scale (the kernel is multiplicative and every
    # fitted parameter is positive), with a robust IQR-based estimator so a couple
    # of stray particles cannot inflate the step for the whole population.
    kernel_scale = if (!is.null(o$kernel_scale)) o$kernel_scale else 0.5
    kernel_floor = if (!is.null(o$kernel_floor)) o$kernel_floor else 0.01

    kernel_sd = vapply(fit$params, function(prm) {
      v = resampled[[prm]]
      v = v[is.finite(v) & v > 0]
      if (length(v) < 2) return(kernel_floor)
      s = stats::IQR(log(v)) / 1.349          # robust; equals sd for a normal
      if (!is.finite(s) || s <= 0) s = stats::sd(log(v))
      if (!is.finite(s)) s = 0
      max(kernel_scale * s, kernel_floor)
    }, numeric(1))

    samples = resampled %>%

      # Perturb: multiplicative log-normal step, sized per parameter
      mutate(across(.cols = all_of(fit$params),
                    .fns  = ~ .x * exp(rnorm(n(), mean = 0, sd = kernel_sd[[cur_column()]])))) %>%

      # Keep values inside the prior by REFLECTING rather than clamping. Clamping
      # parked every overshoot exactly on the bound, creating a spike there that
      # is easily misread as "the prior is too narrow" - it accounted for 3.6% of
      # all sampled values, and 8.3% of peak_day. See reflect_into() in auxiliary.R.
      mutate(across(.cols = all_of(fit$params),
                    .fns  = ~ reflect_into(.x,
                                           bounds_wide[[paste0(cur_column(), "_lower")]],
                                           bounds_wide[[paste0(cur_column(), "_upper")]]))) %>%
      as_named_dt(fit$params)
  }
  
  # Number of samples to simulate 
  n_samples = nrow(samples)
  
  # Prepare for simulation
  params_df = samples %>%
    as_named_dt(fit$params) %>%
    mutate(param_idx = 1 : n_samples, 
           round     = r_val) %>%
    mutate(param_id = get_param_id(round, param_idx)) %>%
    select(param_id, round, all_of(fit$params)) %>%
    setDT()
  
  # Extract parameter IDs
  param_ids = params_df$param_id
  
  # Save simulation df for reference 
  saveRDS(params_df, file = paste0(o$pth$fitting, "r", r_val, "_params.rds"))
  
  return(param_ids)
}

# -------------------------------------------------------- -
# Simulate newly-generated parameter sets ----
# -------------------------------------------------------- -
simulate_parameters = function(o, fit, r_idx, param_ids) {
  
  # Full path to all output files to be produced
  output_files = paste0(o$pth$fit_samples, param_ids, ".rds")
  
  # If we don't want to overwrite, remove reference to files that already exist
  if (o$overwrite_samples == FALSE){
    param_ids = param_ids[!file.exists(output_files)]
  }
  
  # Number of simulations to run 
  n_simulations = length(param_ids)
  
  # ---- Simulate model ----
  
  message("  > Simulating ", thou_sep(n_simulations), " samples")
  
  # Skip this process if nothing to run
  if (n_simulations > 0) {
    
    # ---- Run fitting samples -----
    # Round of adaptive sampling
    r_val = as.numeric(str_remove(r_idx, "r"))
    
    # Load parameter set samples
    sim_df = try_load(o$pth$fitting, paste0(r_idx, "_params"))
    
    # Create a cluster 
    cl = makeCluster(o$parallel)
    
    # Source dependencies to each worker node
    clusterEvalQ(cl, source("R/dependencies.R"))
    
    # Export necessary info to each worker node
    clusterExport(cl, c("o")) 
    
    #for(task_id in 1:nrow(sim_df)){ # Keep this for debugging without parallelisation
    process_task = function(task_id) {
      print(Sys.time())
      # Select parameter set associated with this task ID
      param_df = sim_df[task_id, ]
      param_id = param_df$param_id

      # Parameter values in list format (for input into model)
      param_list = param_df %>%
        select(-round, -param_id) %>%
        unique() %>%
        as.list()

      # Append flag that we want to perform fit
      fit_list = list.append(param_list, .perform_fit = TRUE)
      fit_list$data = fit$data
      fit_list$dates_model = fit$dates_model

      # Fix uncertainty parameters at the central (median) value of their
      # distributions during calibration, so the likelihood reflects the
      # fitted parameters rather than uncertainty-sampling noise. Full
      # parameter uncertainty is propagated later in the scenario runs.
      uncert_list = sample_average(o)  # See uncertainty.R

      message(" - Running model")

      result = model(o,
                     scenario = "baseline",
                     fit      = fit_list,
                     uncert   = uncert_list,
                     verbose  = "none")

      # include round, and filter to fit metrics if specified
      output_df = result$output %>%
        { if (!is.null(o$fit_metrics)) filter(., !is.na(value), metric %in% c(o$fit_metrics, "total")) else . } %>%
        mutate(param_id = !!param_id,
               round    = r_val,
               .before  = 1)

      # Save sample output as an RDS file
      saveRDS(output_df, paste0(o$pth$fit_samples, param_id, ".rds"))

      return(output_df)
    }

    res_list = pblapply(X = 1:nrow(sim_df), FUN = process_task, cl=cl)
    stopCluster(cl)
  }
  
  # ---- Load and summarise output ----
  message("  > Concatenating simulation outcomes")
  
  # # Load all files into one long datatable
  output_matrix = foreach(i_file = output_files,
                          .combine = rbind) %do% { readRDS(i_file) } 
  output_df = data.table(output_matrix)                    
  
  # Save aggregated output to file
  saveRDS(output_df, file = paste0(o$pth$fitting, r_idx, "_output.rds"))
}

# -------------------------------------------------------- -
# Calculate likelihood ----
# -------------------------------------------------------- -
likelihood = function(o, fit, r_idx, do_plot = FALSE) {
  message("  > Calculating likelihood")
  
  # Load simulated parameter sets and associated model output
  param_df = try_load(o$pth$fitting, paste0(r_idx, "_params"))
  output_df = try_load(o$pth$fitting, paste0(r_idx, "_output"))
  
  # Filter model output to fitting targets and summarise 
  model_df = output_df %>% 
    fitting_format(fit, r_idx)
  
  # Compute daily, monthly, and total metric
  model_df = aggregate_model_output(model_df, fit$data, date_col = "date")
  
  # Append weights to data
  data_df = fit$data %>%
    format_weights(fit$input) %>%
    setDT()
  
  # ---- Observation-model over-dispersion ----
  # `k` is the negative-binomial SIZE: Var = mu + mu^2/k, so LARGER k means LESS
  # over-dispersion (k -> Inf is Poisson). It behaves like every other model
  # parameter: if `k` is listed under calibration_parameters it is FITTED, and
  # param_df carries one value per sample; otherwise the fixed value from the
  # yaml is used for every sample.
  if ("k" %in% names(param_df)) {
    disp_df = param_df %>% select(param_id, .k = k)
  } else {
    k_fixed = fit$input$k
    if (is.null(k_fixed))
      stop("Observation-model dispersion 'k' not found in the parsed yaml. Add ",
           "`k:` to config/default.yaml, or list it under calibration_parameters.")
    disp_df = param_df %>% select(param_id) %>% mutate(.k = k_fixed)
  }

  likelihood_df = model_df %>%
    inner_join(data_df, by = c("age_group", "metric", "date", "data_freq"), relationship = "many-to-many") %>%
    left_join(disp_df, by = "param_id") %>%
    # Normalise both target and value to ensure comparable scales
    mutate(value = pmax(value, 0)) %>%
    mutate(
      this_likelihood = weight + dnbinom(round(target), size = .k, mu = value, log = TRUE)
    ) %>%
    # Note: log-likelihood is averaged (mean) within each metric/data_freq,
    # then summed across metrics. The mean step weights each metric equally
    # regardless of its number of data points; switch the mean to sum if you
    # want observation-count weighting instead.
    group_by(param_id, round, metric, data_freq) %>%
    summarise(likelihood_metric = mean(this_likelihood)) %>%
    ungroup() %>%
    group_by(param_id, round) %>%
    summarise(likelihood = sum(likelihood_metric)) %>%
    ungroup()
  
  # Join parameter sets to likelihood
  samples_df = param_df %>%
    inner_join(likelihood_df, by = c("round", "param_id")) %>%
    filter(!is.na(likelihood))
  
  # Save simulation df for reference 
  saveRDS(samples_df, file = paste0(o$pth$fitting, r_idx, "_samples.rds"))
  
}

# -------------------------------------------------------- -
# Compile parameter sets ----
# -------------------------------------------------------- -
compile_samples = function(o, r_val) {
  
  # Initiate list of samples
  samples_list = list()
  
  # Iterate through all already simulated rounds
  for (i in 0 : r_val) {
    
    # Load samples from each round
    round_path    = paste0("r", i, "_samples")
    round_samples = try_load(o$pth$fitting, round_path)
    
    # Store in list as we iterate
    samples_list[[i + 1]] = round_samples
  }
  
  # Compile all samples into a datatable
  all_samples = rbindlist(samples_list)
  
  # Save all these samples with a slightly different naming convention
  saveRDS(all_samples, file = paste0(o$pth$fitting, "rx_samples.rds"))
  
}

# -------------------------------------------------------- -
# Format model output for fitting ----
# -------------------------------------------------------- -
fitting_format = function(model_output, model_input, r_idx) {
  # Formatting the dates
  all_dates = seq(model_input$opts$data_start, model_input$opts$data_end, by = "day")
  dates_df  = data.table(date = all_dates, 
                         day  = 1 : length(all_dates)) %>%
    filter(day > model_input$opts$data_burn_in)
  
  # Map age groups
  output_group = model_output %>% select(age_group) %>% unique() %>% filter(!is.na(age_group))
  
  # Map fine age groups to reporting bands using the canonical mapping carried
  # in the parsed config (default.yaml age_group_map), rather than re-hardcoding
  band_map = unlist(model_input$input$age_group_map)
  age_group_map = output_group %>%
    mutate(data_group = unname(band_map[as.character(age_group)]))
  
  # Remove burn in phase in model output and map to data dates
  model_df = model_output %>%
    # Filter dates to match target data
    filter(time > model_input$opts$data_burn_in) %>%
    mutate(day = as.integer(time)) %>%
    inner_join(dates_df, by = "day") %>%
    select(-c(time, day)) %>%
    group_by(param_id, metric, variant) %>%
    # Group age groups to match data
    left_join(age_group_map, by = "age_group") %>%
    group_by(param_id, round, metric, date, data_group) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    rename(age_group = data_group) %>%
    ungroup()
  
  
  # Remove 'seasonality' metrics before scaling
  model_df = model_df %>%
    filter(metric != "seasonality")
  
  # Separate 'total' age group to avoid affecting scaling
  total_group = model_df %>%
    group_by(param_id, round, metric, date) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(age_group = "total") %>%
    ungroup()
  
  # Combine age-group and total rows
  model_df = bind_rows(model_df, total_group) %>%
    setDT()
  
  # Shorthand for calibration weights list
  w = model_input$input$calibration_weights
  
  # Weights for each metric
  weight_df = as.data.table(w$metric) %>% 
    pivot_longer(cols = everything(), 
                 names_to  = "metric", 
                 values_to = "weight")
  
  # Summarise by variant, where appropriate
  model_df = model_df %>%
    left_join(weight_df, by = "metric") %>%
    filter(weight > 0) %>%
    group_by(param_id, round, age_group, metric, date) %>%
    summarise(value = sum(value)) %>%
    ungroup()
  
  # Write file for plotting purposes
  saveRDS(model_df, file = paste0(o$pth$fitting, r_idx, "plot_output.rds"))
  
  return(model_df)
}

# -------------------------------------------------------- -
# Extract and select calibration weightings ----
# -------------------------------------------------------- -
format_weights = function(data, model_input) {
  
  # Shorthand for calibration weights list
  w = model_input$calibration_weights
  
  # Weights for each metric
  weight = as.data.table(w$metric) %>% 
    pivot_longer(cols = everything(), 
                 names_to  = "metric", 
                 values_to = "weight")
  
  # All possible times
  weight_df = data %>%
    filter(!is.na(value)) %>%
    select(date, age_group, metric, target = value, data_freq) %>%
    left_join(weight, by = "metric") %>%
    setDT()
  
  return(weight_df)
}

# ---------------------------------------------------------------- -
# Summarise where a round's parameter sets sit within their priors ----
# ---------------------------------------------------------------- -
# Diagnostic for "are the prior bounds sensible?". Reports the DISTRIBUTION of
# the sampled sets, not the single best-fitting one: the best set is one draw
# from a noisy search and routinely sits in the tail of its own distribution, so
# judging a prior by whether the argmax landed near a bound is unreliable.
#
# Returns, per parameter, the position of the sample median/quartiles within the
# prior range (0 = lower bound, 1 = upper) plus the share of samples in the
# outer 5% at each end. A prior is only worth widening if a real share of the
# MASS piles up against one end - not because one sample did.
posterior_summary = function(o, cc, r_idx = NULL) {

  pth = paste0("output/1_calibration/", cc, "/")
  fit = try_load(pth, "fit_result", throw_error = FALSE)
  if (is.null(fit)) return(NULL)

  # Default to the last round ACTUALLY RUN for this fit. NB do not just take the
  # highest-numbered r*_samples file on disk: a previous run with more rounds
  # leaves stale files behind, and those would silently be reported instead.
  if (is.null(r_idx)) {
    if (is.null(fit$likelihoods) || !"round" %in% names(fit$likelihoods)) return(NULL)
    r_idx = paste0("r", max(fit$likelihoods$round, na.rm = TRUE))
  }
  s = try_load(pth, paste0(r_idx, "_samples"), throw_error = FALSE)
  if (is.null(s)) return(NULL)

  bd = data.frame(param = fit$params, fit$bounds)
  do.call(rbind, lapply(fit$params, function(prm) {
    v  = s[[prm]]
    lo = bd$lower[bd$param == prm]; hi = bd$upper[bd$param == prm]

    # Position is measured on the LOG scale, because sample_parameters draws the
    # initial Latin hypercube in log10 space. On a wide prior such as k = [1, 100]
    # a linear position is misleading: the log-uniform median (k = 10) sits at a
    # linear position of only 0.09, which looks like "piled up against the floor"
    # when it is simply the middle of the prior.
    use_log = is.finite(lo) && lo > 0 && is.finite(hi) && hi > lo
    pos = if (use_log) (log(v) - log(lo)) / (log(hi) - log(lo)) else (v - lo) / (hi - lo)

    data.frame(iso = cc, round = r_idx, param = prm, scale = ifelse(use_log, "log", "lin"),
               q25 = round(quantile(pos, .25, names = FALSE), 3),
               med = round(median(pos), 3),
               q75 = round(quantile(pos, .75, names = FALSE), 3),
               frac_low  = round(mean(pos < 0.05), 3),   # share of MASS at each end
               frac_high = round(mean(pos > 0.95), 3),
               at_bound  = sum(v == lo | v == hi))       # exact hits => clamping
  }))
}

# -------------------------------------------------------- -
# Create consistent IDs for parameter sets ----
# -------------------------------------------------------- -
get_param_id = function(round, param) {
  
  # Pad values for file name consistency
  param_id = paste0("r", sprintf("%02i", round),
                    "p", sprintf("%04i", param))
  
  return(param_id)
}

# -------------------------------------------------------- -
# Get days we need to simulate for fitting purposes ----
# -------------------------------------------------------- -
get_fit_days = function(p) {
  
  fit_days = p$calibration_options$data_days + 
    p$calibration_options$data_burn_in
  
  return(fit_days)
}

# -------------------------------------------------------- -
# Which metrics are to be reported for fitting purposes ----
# -------------------------------------------------------- -
get_fit_metrics = function(p) {
  
  # All non-trivial epi metrics required
  fit_metrics = p$calibration_weights$metrics
  fit_metrics = names(fit_metrics[fit_metrics > 0])
  
  return(fit_metrics) 
}


# ------------------------------------------------------------ -
# Alter key parameters in yaml file for fitting simulation ----
# ------------------------------------------------------------ -
apply_fit = function(y, fit_list) {
  
  # Only needed if fitting OR using fitted parameters
  if (!is.null(fit_list)) {
    
    # Check if we are fitting here - indicated by .perform_fit
    if (isTRUE(fit_list$.perform_fit)) {
      
      # In this case, we only need to run the bare essentials...
      
      # Bare essential metrics
      fit_metrics = get_fit_metrics(y)
      
      # Function for removing any disaggregation of metrics 
      reduce_fn = function(x) {if (!is.null(x$by)) {x$by = "none"}; return(x)}
      
      # Now safely remove flag, it's done it's job
      fit_list$.perform_fit = NULL
    }
    
    # ---- Alter parameters ----
    
    # Loop through parameters to alter
    for (param in names(fit_list)) {
      
      # Check whether parameter we want to overwrite already exists
      eval_str("param_valid = is.numeric(y$", param, ")")
      
      # Check parameter is valid - throw error if not
      if (!param_valid)
        stop("Attempting to alter invalid parameter '", param, "' during model fitting")
      
      # Re-define paramter value (using eval to enable change of listed items)
      eval_str("y$", param, " = fit_list[[param]]")
    }
  }
  
  return(y)
}



# --------------------------------------------------------------------------------- -
# Aggregate model output for different time freq: daily, weekly, monthly, total ----
# --------------------------------------------------------------------------------- -

aggregate_model_output <- function(data_model, data_reported, date_col = "date") {
  # store name of the date column
  date_name <- date_col
  
  
  # get all grouping columns (everything except date and value)
  group_cols <- setdiff(names(data_model), c(date_name, "value"))
  
  # DAILY
  daily <- data_model %>%
    mutate(data_freq = "daily")
  
  # WEEKLY
  weekly <- data_model %>%
    mutate(period = ceiling_date(data_model[[date_name]], "week", , change_on_boundary = FALSE)) %>%
    group_by(across(all_of(group_cols)), period) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    rename(!!date_name := period) %>%
    mutate(data_freq = "weekly")
  
  # 4-WEEKLY
  # Get starting date from data
  if (length(data_reported) == 1 && is.na(data_reported)){
    four_week = NULL
    
  } else{
    dates_4weekly = intersect(data_model$date %>% unique(),
                              data_reported %>% filter(data_freq == "4-weekly") %>%
                                pull(date)) %>% as.Date()

    if (length(dates_4weekly) == 0) {

      # No 4-weekly targets in the reported data, so there is no grid to anchor
      # to and nothing to aggregate. NB the RespiCompass 2026/2027 round supplies
      # the age-stratified burden as SEASONAL TOTALS (data_freq = "total"), not
      # as a 4-weekly series, so this is the normal path for that round.
      four_week = NULL

    } else {

      start_date_4weekly = min(dates_4weekly)
      anchor_monday <- floor_date(start_date_4weekly, "week", week_start = 1)
      four_week <- data_model %>%
        filter(date >= anchor_monday) %>%
        mutate(
          wk_monday = floor_date(.data[[date_name]], "week", week_start = 1),
          period = anchor_monday + 6 + weeks(4) * (as.integer(wk_monday - anchor_monday) %/% 28)
        ) %>%
        group_by(across(all_of(group_cols)), period) %>%
        summarise(value = sum(value), .groups = "drop") %>%
        rename(!!date_name := period) %>%
        mutate(
          data_freq = "4-weekly"
        )
    }
  }
  
  # MONTHLY
  monthly <- data_model %>%
    mutate(period = floor_date(data_model[[date_name]], "month")) %>%
    group_by(across(all_of(group_cols)), period) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    rename(!!date_name := period) %>%
    mutate(data_freq = "monthly")
  
  # TOTAL (overall)
  total <- data_model %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(!!date_name := as.Date(NA),
           data_freq = "total") %>%
    select(all_of(names(daily)))
  
  
  return( bind_rows(daily, weekly, four_week, monthly, total) )
  
}
