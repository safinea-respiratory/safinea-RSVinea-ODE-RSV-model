########################################################## #
# RESULTS
#
# Plot standard set of epi. results.
#
########################################################## #

# -------------------------------------------------------- -
# Parent function for plotting standard results set ----
# -------------------------------------------------------- -
run_results = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(3, o$do_step)) return()
  
  message("* Producing results and figures")

  all_scenarios = names(parse_yaml(o, "*read*"))

  # Check plotting flag 
  if (o$plot_scenarios == TRUE && length(all_scenarios) > 0) {
     
    message(" - Plotting scenarios")
     
    # Alternative scenarios - all metrics
    plot_scenarios(o,  fig_name = paste0("scenario_", o$analysis_name), scenarios = all_scenarios) 
  }
  
  # Generate parameter table
  message(" - Generate parameter table")
  write_parameter_table(o)
}

# ------------------------------------------------------------ -
# Post process raw model output and append to 'result' list ----
# ------------------------------------------------------------ -
process_results = function(o, result, raw_output) {

  # Map fine model age groups to reporting bands using the canonical
  # mapping in default.yaml (see age_band_lookup in auxiliary.R)
  band_map = age_band_lookup(o)
  age_group_map = raw_output %>%
    select(age_group) %>%
    unique() %>%
    filter(!is.na(age_group)) %>%
    mutate(data_group = unname(band_map[as.character(age_group)]))
  
  group_output = raw_output %>% left_join(age_group_map, by = "age_group") %>%
    group_by(scenario, param_id, metric, time, data_group) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    rename(age_group = data_group) %>%
    
    filter(!metric == "seasonality") %>%
    group_by(scenario, param_id, time, age_group) %>%
    select(scenario, param_id, metric, time, age_group, value) %>%
    ungroup() 
  
  pop_output = raw_output %>% left_join(age_group_map, by = "age_group") %>%
    group_by(scenario, param_id, metric, time, data_group) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    rename(age_group = data_group) %>%
    
    filter(!metric == "seasonality") %>%
    group_by(scenario, param_id, time, metric) %>%
    summarise(value = sum(value)) %>%
    mutate(age_group = "total") %>%
    select(scenario, param_id, metric, time, age_group, value) %>%
    ungroup() 
  
  output_df = bind_rows(group_output, pop_output)
  
  # Create key summary statistics across parameter uncertainty
  result$output =  output_df  %>% 
    group_by(scenario, metric, age_group, time) %>%
    summarise(mean   = mean(value),
              median = quantile(value, 0.5, na.rm = TRUE),
              lower  = quantile(value, o$quantiles[1], na.rm = TRUE),
              upper  = quantile(value, o$quantiles[2], na.rm = TRUE),
              .groups = "drop") 
  
  return(result)
}

# -------------------------------------------------------- -
# Aggregate model parameters and generate summary table ----
# -------------------------------------------------------- -
write_parameter_table = function(o) {
  
  # Load data
  fit <- as.data.frame(load_calibration(o))
 
  # Load and parse user-defined inputs for this scenario
  yaml = parse_yaml(o, scenario = "baseline", fit = fit)
  
  # Shorthand for model parameters
  param_data = yaml$parsed
  
  # Adjust for age-specific parameters
  param_data = adjust_age_specific_param(param_data)
    
  # Initiate output table
  param_tbl_out <- data.frame(matrix(NA,ncol=4,nrow=length(param_data)))
  names(param_tbl_out) <- c("parameter",'n_best_samples','best_sample','prior')
  param_tbl_out$parameter = names(param_data)
  
  # Include posteriors
  for(i in 1:length(param_data)){
    values <- param_data[[i]]
    if(is.numeric(values) || length(unlist(values))==1){
      param_tbl_out[i,2] = get_param_summary(values,o$n_best_samples)
      param_tbl_out[i,3] = get_param_summary(values,1)
    } 
  }
  
  # Obtain priors
  calib_param <- yaml$parsed$calibration_parameters
  param_prior <- param_data
  for(i_param in 1:length(calib_param)){
    calib_subset <- calib_param[[i_param]]
    param_prior[calib_subset[[1]]] <- list(unlist(calib_subset[-1]))
  }
  
  # Adjust for age-specific parameters
  param_prior = adjust_age_specific_param(param_prior)
  
  # Include priors in Table
  # note: select param with lower and upper value 
  param_tbl_out$prior[] <- ''
  for(i_param in 1:length(param_prior)){
    calib_subset <- param_prior[[i_param]]
    if(length(calib_subset) > 1 && !is.null(names(calib_subset)) && all(names(calib_subset) == c('lower','upper'))){
      sel_row <- param_tbl_out$parameter == names(param_prior)[i_param]
      param_tbl_out$prior[sel_row] <- paste0("[",paste(calib_subset,collapse = ';'),"]")
    }
  }
  
  # Include number of samples at first NA row
  row_empty = which(is.na(param_tbl_out$best_sample))[1]
  param_tbl_out[row_empty,] = c("samples",o$n_best_samples,1,'') 
  
  # remove other NA's
  param_tbl_out = param_tbl_out[!is.na(param_tbl_out$best_sample),]
    
  fig_name = paste0("scenario_", str_to_title(o$analysis_name))
  
  # Store as csv 
  write.table(param_tbl_out, file = paste0(o$pth$results, "parameters_tbl_",o$analysis_name,".csv"), sep = ',', row.names = FALSE)
}

# -------------------------------------------------------- -
# Combine general and age-specific parameters ----
# -------------------------------------------------------- -
adjust_age_specific_param = function(param_data){
  param_data$qi_a <- param_data$beta_A * param_data$rel_sus_a
  param_data$qi_b <- param_data$beta_A * param_data$rel_sus_b
  param_data$qi_c <- param_data$beta_A * param_data$rel_sus_c
  param_data$qi_d <- param_data$beta_A 
  
  # These are now absolute probabilities set directly by their own parameter,
  # so there is no longer a p_hosp_A multiplier to apply (see age_relativity()
  # in model.R). p_hosp_d remains the 5-18y / 18-60y baseline.
  param_data$p_hosp_a <- param_data$p_hosp_a_A
  param_data$p_hosp_b <- param_data$p_hosp_b_A
  param_data$p_hosp_c <- param_data$p_hosp_c_A
  param_data$p_hosp_d <- param_data$p_hosp_A
  
  param_data$p_death_a <- param_data$p_death * param_data$rel_death_a
  param_data$p_death_b <- param_data$p_death * param_data$rel_death_b
  param_data$p_death_c <- param_data$p_death * param_data$rel_death_c
  param_data$p_death_d <- param_data$p_death
  
  return(param_data)
}

# -------------------------------------------------------- -
# Get parameter summary in terms of the mean and 95%CrI ----
# -------------------------------------------------------- -
get_param_summary = function(values, n_best_samples){
  
  # if only one value, return
  if(length(values) == 1 || n_best_samples == 1){
    return(values[1])
  }
  
  # select
  values = values[1:min(length(values),n_best_samples)]
  
  # format
  param_txt = paste0(
         format(mean(values),digits = 3, scientific = FALSE), ' [',
         format(quantile(values,0.025),digits = 2, scientific = FALSE),';',
         format(quantile(values,0.975),digits = 2, scientific = FALSE), ']')
     
  return(param_txt)
}
