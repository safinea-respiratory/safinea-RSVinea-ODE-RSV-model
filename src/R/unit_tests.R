########################################################## #
# UNIT TESTS
#
# A series of self-contained functions to test out different
# aspects of model functionality.
#
########################################################## #

# ------------------------------------------------------------ -
# Run a single simulation and produce a few standard outputs ----
# ------------------------------------------------------------ -
run_model_test = function(o, scenario = "baseline", rerun = TRUE, benchmark = FALSE) {

  # Only continue if specified by do_step
  if (!is.element(0, o$do_step)) return()
  
  # Load parsed parameters from yaml file
  p = parse_yaml(o, scenario = scenario)$parsed
  opts = c(type = p$calibration_type, p$calibration_options)
  
  message("* Testing single model simulation")

  # ---- Load the epidemiological data ----
  # model() is not self-contained: age_relativity() derives the age-relativity
  # ratios from the 4-weekly hospital burden that arrives as fit$data (see
  # `p$burden` in model.R). Passing fit = NULL therefore fails with
  # "no applicable method for 'filter' applied to an object of class NULL"
  # before the ODE is ever reached. Load the data the same way run_calibration()
  # does.
  fit = setup_calibration(o)   # See calibration.R
  fit = load_data(o, fit)      # See load_data.R

  # ---- Choose the parameter set to test at ----
  # Run at the CENTRE OF THE CALIBRATION PRIORS rather than at the bare yaml
  # scalars. The scalar defaults are not a runnable parameterisation on their
  # own: p_hosp_A (0.0458) and rel_hosp_a_A (50) are inherited from an older
  # parameterisation and multiply to 2.29, so age_relativity() correctly rejects
  # them as a probability. Every country fits these parameters, so the scalars
  # are never used in a real run. Testing at the prior midpoint therefore
  # exercises a parameter set the model is actually meant to run at, and keeps
  # this test honest without silently retuning model defaults.
  prior_mid = rowMeans(fit$bounds)
  names(prior_mid) = fit$params

  test_fit = c(as.list(prior_mid),
               list(data        = fit$data,
                    dates_model = fit$dates_model))

  # ---- Simulate model ----
 
  # File to save after a simulation
  test_file = paste0(o$pth$testing, "model_test.rds")
  
  # Check whether we want to simulate (or simply plot a previously created file)
  if (rerun == TRUE || !file.exists(test_file)) {

   set.seed(123)
    
   # If uncertainty defined, take the average over the distribution(s)
   uncert_list = sample_average(o)  # See uncertainty.R
    
    # Run model for the defined scenario (see model.R)
    result = model(o, 
                   scenario = scenario,
                   fit     = test_fit,
                   uncert  = uncert_list,
                   do_plot = TRUE,
                   verbose = "date")
    saveRDS(result, test_file)
  }
  
  # Load simulated file - we may have skipped simulating
  result = readRDS(test_file)
  
  # Fitting plots inc. by age group (see plotting.R)
  plot_simple_output(o, result, fig_name = paste0("Unit_test_", o$analysis_name))
  
  # Option to benchmark current output with previous
  if(benchmark == TRUE){
    
    # Check if there exists a benchmark file
    benchmark_test_file = gsub('.rds','_benchmark.rds',test_file)
    if(!file.exists(benchmark_test_file)){
      message("* No benchmark file to validate model simulation output")
      message("* Call `rebase_benchmark()` to store reference values")
    } else{
      
      # load benchmark results and compare with current results
      benchmark_result <- readRDS(benchmark_test_file)
      boolean_equal <- all.equal(result, benchmark_result)
      
      
      # Print comparison to terminal
      if(is_bare_logical(boolean_equal) && boolean_equal){
        message("* Unit-test model results did not change")
      } else{
        message("* Unit-test model results changed!")
        message(paste(boolean_equal,collapse = '\n'))
      }
    }
  } 
}
  
# -------------------------------------------------------- -
# Rename 'model_test.rds' as benchmark to enable validation
# -------------------------------------------------------- -
rebase_benchmark = function(){
  test_file = paste0(o$pth$testing, "model_test.rds")
  saveRDS(readRDS(test_file),file=gsub('.rds','_benchmark.rds',test_file))
  message(" * Rebased benchmark")
}

