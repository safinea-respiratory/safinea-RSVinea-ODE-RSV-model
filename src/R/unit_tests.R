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
                   fit     = NULL,
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

