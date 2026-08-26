########################################################## #
# SCENARIOS
#
# Simulates all scenarios defined in yaml file for this 
# analysis, always including baseline.
#
########################################################## #

# -------------------------------------------------------- -
# Parent function for running scenarios ----
# -------------------------------------------------------- -
run_scenarios = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Running all scenarios")
  
  # Create and read scenarios to simulate (see parse_input.R)
  all_scenarios = names(parse_yaml(o, "*create*"))
  
  # ---- Set up full set of simulations ----
  
  # Sample values from all parameter uncertainty distributions
  uncert_df = sample_uncertainty(o)  # See uncertainty.R
  
  # Generate set of sim IDs (depends on uncertainty options)
  sim_df = create_sim_id(o, uncert_df, all_scenarios) 
  
  # Save to file for later use
  saveRDS(sim_df, file = paste0(o$pth$simulations, "all_simulations.rds"))
  
  # ---- Simulate with local parallelisation ----
  
  # Create a cluster 
  cl = makeCluster(o$parallel)
  
  # Source dependencies to each worker node
  clusterEvalQ(cl, source("R/dependencies.R"))
  
  # Export necessary info to each worker node
  clusterExport(cl, c("o")) 
  
  #for(task_id in 1:nrow(sim_df)){ # Keep to simplify debugging without cluster
  process_task = function(task_id) {

    # Select analysis to run based on job ID
    this_sim = sim_df[task_id, ]
    
    # If yes, load uncertainty datatable
    uncert_df = try_load(o$pth$uncertainty, "uncertainty")
    
    # Select values relevant for this simulation
    this_uncert = uncert_df[param_set == this_sim$param_set, ]
    
    # Deal with uncertaties that come from scenarios (aka differ between scenarios)
    scen <- this_sim$scenario   # e.g. "A"
    pat  <- paste0("^scenarios\\$", scen, "\\$")
    this_uncert <- this_uncert[
      # keep rows that don't start with "scenarios$"
      !grepl("^scenarios\\$", param) | grepl(pat, param)
    ][
      # for scenario rows, strip the prefix
      grepl("^scenarios\\$", param),
      param := sub("^scenarios\\$[^$]+\\$", "", param)
    ]
    
    
    # Convert to list
    uncert_list = this_uncert$value %>% 
      setNames(this_uncert$param) %>%
      as.list()
    
    # Loaded fitted model parameters relevant for this simulation (if available)
    fit_list = load_calibration(o, index = this_sim$fitting_set)
    
    # Read in target data
    fit_list$data = get_target_data(o)
    # Load data - see load_data.R
    fit1 = load_data(o, setup_calibration(o))
    fit_list$dates_model = fit1$dates_model
    
    message(" - Running model")

    # Adult vaccine waning replicate for this simulation, mapped 1:1 from the
    # FITTING sample (sample 1 -> rep 1, sample 2 -> rep 2, ...) so the 500
    # RespiCompass waning curves propagate into the projection intervals.
    # Keyed on fitting_set rather than param_set because n_best_samples is the
    # index that is routinely > 1 (n_parameter_sets is often left at 1, which
    # would pin every simulation to rep 1 and erase the waning uncertainty).
    # Deterministic by construction: re-running a sim_id reproduces its curve.
    # get_waning_curve() wraps with modulo if there are more samples than reps.
    o$waning_rep = this_sim$fitting_set

    # Simulate model with this parameter set
    result = model(o, this_sim$scenario,
                   fit     = fit_list,
                   uncert  = uncert_list,
                   verbose = "none")
    
    # Extract metrics of interest over time points of interest
    output_df = result$output %>%
      filter(!is.na(value)) %>%
      mutate(param_id = this_sim$sim_id)
    
    # Save sample output as an RDS file
    saveRDS(output_df, paste0(o$pth$simulations, this_sim$sim_id, ".rds"))
  }
  
  res_list = pblapply(X = 1:nrow(sim_df), FUN = process_task, cl=cl)
  stopCluster(cl)
  
  message(" - Quantifying parametric uncertainty")
  
  # ---- Load and summarise output ---- 
  
  # Load full set of simulations to summarise
  sim_df = try_load(o$pth$simulations, "all_simulations")
  
  n_scenarios = length(unique(sim_df$scenario))
  
  for(task_id in 1:n_scenarios){
    
    # Scenarios to summarise
    scenario_name = unique(sim_df$scenario)[task_id]
    
    message("   * Loading scenario ", scenario_name)
    
    # All simulations to summarise over (calibration and uncertainty)
    sim_ids = sim_df %>% 
      filter(scenario == scenario_name) %>%
      pull(sim_id)
    
    # Iterate over simulations and rbind model output
    output_matrix = foreach(sim_id = sim_ids,
                            .combine = rbind) %do%{
                              sim_result = try_load(o$pth$simulations, sim_id)
                            }
    
    
    # Throw an error if no results found
    if (nrow(output_matrix) == 0)
      stop("No results available for scenario '", scenario_name, "'")
    
    # Initiate analysis list - we'll add summarised model output to this
    result = list(analysis_name = o$analysis_name,
                  scenario_name = scenario_name)
    
    # Bind model output into single datatable
    raw_df = output_matrix[!is.na(value), ] %>%
      mutate(scenario = scenario_name)
    
    # Save this raw (ie unsummarised) form of model output
    saveRDS(raw_df, file = paste0(o$pth$scenarios, scenario_name, "_raw.rds"))
    
    # Summarise model output
    result = process_results(o, result, raw_df)
    
    # Store summarised result
    saveRDS(result, file = paste0(o$pth$scenarios, scenario_name, "_summarised.rds"))
  }
  
}

# -------------------------------------------------------------------- -
# Generate set of simulation IDs (depends on uncertainty options) ----
# -------------------------------------------------------------------- -
create_sim_id = function(o, uncert_df, all_scenarios) {
  
  # Check we have a positive number of samples to generate
  if (o$n_parameter_sets < 1)
    stop("The value of 'n_parameter_sets' must be a positive integer")
  
  # Unique parameter sets
  if (!is.null(uncert_df)) 
    param_sets = unique(uncert_df$param_set)
  
  # Trivial if not simulating parameter uncertainty
  if (is.null(uncert_df)) 
    param_sets = 1
  
  # IDs of unique uncertainty sets (padded by zeros)
  param_id = str_pad(param_sets, 4, pad = "0")
  
  # Number of fitting samples
  fitting_sets = seq(1:o$n_best_samples)
  
  # IDs of unique fitting sets (padded by zeros)
  fitting_id = str_pad(fitting_sets, 4, pad = "0")
  
  all_sims = expand_grid(param_id = param_id, 
                         fitting_id = fitting_id,
                         scenario = all_scenarios)
  
  # Use all of this to create unique simulation identifiers
  sim_df = all_sims %>% 
    mutate(sim_id = unite(all_sims, "x", sep = "_")$x, 
           sim_id = paste0("s", sim_id), 
           param_set = as.numeric(param_id),
           fitting_set =  as.numeric(fitting_id)) %>%
    select(sim_id, scenario, fitting_set, param_set) %>% 
    arrange(scenario, fitting_set, param_set) %>% 
    setDT()
  
  return(sim_df)
}

