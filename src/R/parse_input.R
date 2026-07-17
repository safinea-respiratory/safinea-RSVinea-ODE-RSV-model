########################################################## #
# PARSE INPUT 
#
# Load yaml input file of interest and parse all inputs such
# that they can be interpreted by the transmission model.
#
########################################################## #

# ------------------------------------------------------------ -
# Read yaml input file and perform some additional parsing ----
# ------------------------------------------------------------ -
parse_yaml = function(o, scenario, fit = NULL, uncert = NULL) {
  
  # NOTE: The initial parsing work is done by read_yaml from yaml package
  
  # Check that the user-defined file actually exists 
  if (!file.exists(o$pth$params_user))
    stop("\n You are attempting to run analysis: ", o$analysis_name, 
         "\n  but no yaml file was found with this name", 
         "\n  (missing file: ", o$pth$params_user, ")")
  
  read_yaml_fn = "read_yaml" # YAML is read in R

  # Load default model parameters
  y = get(read_yaml_fn)(o$pth$params_default)
  
  # Load user-defined model parameters
  y_user = get(read_yaml_fn)(o$pth$params_user)
  
  # Overwrite any parameter defaults for which we have user-defined values 
  # NOTE: A few sanity checks on the user-defined inputs are performed here
  list[y, u] = overwrite_defaults(y, y_user)

  # Read in uncertainty
  list[y, u] = parse_uncertainty(o, y, u = u)

  # If we only need uncertainty details, return out now
  if (!is.null(uncert) && uncert[[1]] == "*read*")
    return(u)
  
  #  ---- Apply fitting/fitted parameters if provided ----

  # Apply fitted parameters (or simulate parameter sets when calibrating)
  y = apply_fit(y, fit)  # See calibration.R
  
  # Apply parameter uncertainty values

  y = apply_uncertainty(y, uncert) # See uncertainty.R
  
  # ---- Apply scenario of choice ----
  
  # Apply values for specified scenario
  # NOTE: Some general checks on the content of all scenarios are performed here
  y = parse_scenarios(o, y, scenario)
  
  # If we only need to read scenario names, return out here with named vector
  if (scenario %in% c("*read*", "*create*"))
    return(y$scenario_names)

  # ---- Network ----
  # Convert comma-separated string to vector of strings
  y = str2vec(y, "contact_matrix_countries")
  
  # Throw an error if country codes are in not in ISO-3 format
  if (any(sapply(y$contact_matrix_countries, nchar) != 3))
   stop("Use ISO Alpha-3 country codes for 'contact_matrix_countries' parameters")

  # ---- Calibration options ----
  
  # Easy access calibration options
  opts = y$calibration_options

  # Convert data end to date
  opts$data_end = format_date(opts$data_end)
  
  # Calculate data start date
  opts$data_start = opts$data_end - opts$data_days - opts$data_burn_in + 1
  
  # Overwrite calibration options list and append calibration period
  y$calibration_options = c(opts, data_period = y$calibration_time_period)
  
  # Remove redundant items
  y[c("calibration_time_period")] = NULL

  # ---- Final formatting ----
  
  # Order elements alphabetically to make it easy to search
  y = y[order(names(y))]
  
  # Provide both the raw yaml file and the parsed parameters used in the model
  yaml = list(raw = as.yaml(y), parsed = y)
  
  return(yaml)
}

# ------------------------------------------------------------------------ -
# Overwrite any parameter defaults for which we have user-defined values ---- 
# ------------------------------------------------------------------------ -
overwrite_defaults = function(y, y_overwrite) {

  # Scenario block field redundant once files have been loaded
  y_overwrite$scenario_block = NULL
  
  # Sanity checks on user-defined inputs
  do_checks(y1 = y, y2 = y_overwrite)
  
  # The key exception for overwriting is scenarios, which we'll concatenate
  # except if 'baseline' is redefined or renamed
  if(any(unlist(y_overwrite$scenarios) == "baseline")){ # if y_overwite also contains "baseline"
    y$scenarios = y_overwrite$scenarios
  } else {
    y$scenarios = c(y$scenarios, y_overwrite$scenarios)
  }
  
  # User defined scenarios can now be removed from this list
  y_overwrite$scenarios = NULL
  
  # ---- Name all unnamed lists with unique IDs ----
  
  # Start by naming scenarios with their unique IDs
  scenario_ids = unlist(lapply(y$scenarios, function(x) x$id))
  y$scenarios  = setNames(y$scenarios, scenario_ids)
  
  # Then for the lists within each scenario (as scenarios can have nested unnamed lists)
  for (scenario in names(y$scenarios))
    y$scenarios[[scenario]] = name_lists(y$scenarios[[scenario]])
  
  # Then do the same for all other unnamed lists
  y = name_lists(y, y_overwrite = y_overwrite)
  
  # Also for what we still need to overwrite with
  y_overwrite = name_lists(y_overwrite)
  
  # ---- Calibration and uncertainty parameters ----
  
  # If user has defined a calibration block, use this directly
  if (!is.null(y_overwrite$calibration_parameters))
    y$calibration_parameters = y_overwrite$calibration_parameters
  
  # Any user defined calibration can now be removed from this list
  y_overwrite$calibration_parameters = NULL
  
  # Remove uncertainty parameters (we apply these after sampling from distributions)
  list[y_overwrite, u] = parse_uncertainty(o, y_overwrite)
  list[y, u]           = parse_uncertainty(o, y, u = u)
  
  # ---- Functional forms ----
  
  # All individual parameter items (see auxiliary.R)
  overwrite_items = names(unlist_format(y_overwrite))

  # Which of the these are functions - these are special cases
  fn_call_idx = grepl("\\$fn$", overwrite_items)
  fn_items = str_remove(overwrite_items[fn_call_idx], "\\$fn")
  
  # Skip this if no function items
  if (length(fn_items) > 0) {
    
    # Apply function using eval 
    #
    # NOTE: No checks needed as key-pairs can differ depending on function defined
    for (fn in fn_items)
      eval_str("y$", fn, " = y_overwrite$", fn)
    
    # All items that pertain to a function call
    fn_item_idx = paste0("^") %>%
      paste0(fn_items, collapse = "|") %>%
      str_replace_all("\\$", "\\\\$") %>%
      grepl(overwrite_items)
    
    # Drop these function items from recursive overwriting
    overwrite_items = overwrite_items[!fn_item_idx]
  }

  # ---- Whole-sequence (array) items ----
  #
  # unlist() flattens a multi-element yaml sequence into INDEXED names
  # ("item1", "item2", ...). Those names cannot be resolved with `$` on either
  # side, so if left in the loop below both `y$item1` and `y_overwrite$item1`
  # evaluate to NULL, the class check passes trivially ("NULL" == "NULL"), and
  # the class-preserving step builds get("as.NULL") and errors out.
  #
  # Overwrite such items whole instead, then drop their flattened names from the
  # loop. This is what allows a country yaml to override vector parameters, e.g.
  # adult_vaccination_dates, infant_vaccination_start, adult_vaccination_agegroups
  # or a per-age vector such as background_mortality_rate.
  is_seq = function(v)
    (is.atomic(v) || is.list(v)) && is.null(names(v)) && length(v) > 1

  seq_items = names(y_overwrite)[vapply(y_overwrite, is_seq, logical(1))]

  if (length(seq_items) > 0) {

    for (seq_item in seq_items) {

      # Must exist in the defaults (every parameter is initialised there)
      if (is.null(y[[seq_item]]))
        stop(" ! Unrecognised item in input yaml file: ", seq_item)

      y[[seq_item]] = y_overwrite[[seq_item]]
    }

    # Drop the flattened (indexed) names of these items from the loop below
    overwrite_items = setdiff(overwrite_items,
                              names(unlist_format(y_overwrite[seq_items])))
  }

  # ---- Overwrite defaults ----
  
  # Loop through whats left: values to overwrite with
  for (overwrite_item in overwrite_items) {
    
    # The value to overwrite with, and the original (to check for consistency)
    overwrite_value = eval_str("y_overwrite$", overwrite_item)
    default_value   = eval_str("y$", overwrite_item)
    
    # Normally we will be replacing some non-trivial value
    if (length(default_value) > 0) {
      
      # Check class consistency
      do_checks(y1 = setNames(default_value, overwrite_item), 
                y2 = setNames(overwrite_value, overwrite_item))
      
    } else {  # Otherwise value has been trivialised within name_lists
      
      # However class should still be consistent
      if (class(overwrite_value) != class(default_value))
        stop("Inconsistent data class in input yaml file: \n", 
             paste0(" ! ", overwrite_item, ": ", class(default_value), 
                    " -> ", class(overwrite_value), "\n"))
    }
    
    # If a string, wrap in quotes so eval knows it's not a variable
    if (is.character(overwrite_value))
      overwrite_value = paste0("'", overwrite_value, "'")
    
    # Ensure class is preserved by using as.xxx
    as_class = paste0("get('as.", class(default_value), "')(", overwrite_value, ")")

    # Apply the assignment using eval for n-nested lists
    eval_str("y$", overwrite_item, " = ", as_class)
  }
  
  return(list(y, u)) 
}

# -------------------------------------------------------- -
# Sanity checks on user-defined inputs ----
# -------------------------------------------------------- -
do_checks = function(y1 = NULL, y2 = NULL) {
  
  # Names of paramters defined by user
  user_params = names(y2)
  
  # Check if any parameters have been repeated
  duplicate_items = user_params[duplicated(user_params)]
  
  # Throw an error if this is the case
  if (length(duplicate_items) > 0)
    stop(" ! Duplicated items in input yaml file: ", paste0(duplicate_items, collapse = ", "))
  
  # Check if any unknown items been defined
  unknown_items = setdiff(user_params, names(y1))
  
  # Throw an error if this is the case
  if (length(unknown_items) > 0)
    stop(" ! Unrecognised items in input yaml file: ", paste0(unknown_items, collapse = ", "))
  
  # Data class of all user-defined parameters
  class_user = unlist(lapply(y2, class))
  
  # Default data class of these same parameters
  class_def = unlist(lapply(y1, class))[user_params]
  
  # Parameters that will cause problems as they do not have the same data class
  errors = names(class_def)[!compare.list(class_def, class_user)]
  
  # This isn't a problem if these are 'function lists'
  for (err_param in errors) {
    
    # Check if either class is a list
    is_list = c(class_def[err_param], class_user[err_param]) == "list"
    
    # If yes, continue
    if (any(is_list)) {
      
      # Extract the elements from the list of interest
      list_items = get(paste0("y", which(is_list)))[[err_param]]
      
      # If this is a 'function list', remove this parameter from vector of errors
      if (names(list_items)[[1]] %in% c("fn", "uncertainty"))
        errors = setdiff(errors, err_param)
    }
  }

  # Any inconsistencies remain?
  if (length(errors) > 0) {
    
    # Throw an error reporting which parameters, what it should be, and what it is
    stop("Inconsistent data class in input yaml file: \n", 
         paste0(" ! ", errors, ": ", class_def[errors], " -> ", class_user[errors], "\n"))
  }
}

# -------------------------------------------------------- -
# Apply names to un-named lists with unique IDs ----
# -------------------------------------------------------- -
name_lists = function(y, y_overwrite = NULL) {
  
  # Items that are themselves lists
  y_lists = y[unlist(lapply(y, is.list))]
  
  # Of those lists, the ones which have NULL names (ie unnamed lists)
  y_null = y_lists[unlist(lapply(lapply(y_lists, names), is.null))]
  
  # Of those unnamed lists, the ones which have unique IDs
  y_id = y_null[unlist(lapply(y_null, function(x) names(x[[1]])[1] == "id"))]
  
  # Apply the IDs as the name of the previously unnamed list
  y_named = lapply(y_id, function(x) setNames(x, list2dt(x, fill = TRUE)$id))
  
  # Which of this now named lists will be overwritten
  for (param in intersect(names(y_overwrite), names(y_named))) {
    
    # IDs of this parameters before and after the overwrite
    param_ids_default = unlist(lapply(y[[param]], function(x) x$id))
    param_ids = unlist(lapply(y_overwrite[[param]], function(x) x$id))
    
    # We'll want to reset all values if not identical
    if (!identical(param_ids, param_ids_default)) {
      
      # Store the basic structure of each item
      param_format = y_named[[param]][1]
      
      # Loop through elements are extract class
      for (param_el in names(param_format[[1]])) {
        param_class = class(param_format[[1]][[param_el]])
        
        # Reset value by trivialising, but retraining class
        eval_str("param_format[[1]]$", param_el, " = ", param_class, "(0)")
      }
      
      # Set up the necessary number of items and apply the basic structure
      y_named[[param]] = NULL
      y_named[[param]][param_ids] = param_format
    }
  }
  
  # Store these newly named lists
  y[names(y_named)] = y_named
  
  return(y)
}

# -------------------------------------------------------- -
# Construct function call as a string to be evaluated ----
# -------------------------------------------------------- -
parse_fn = function(fn_args, along = NULL, evaluate = TRUE) {
  
  # Extract function name from input list
  fn = fn_args$fn
  
  # Remove this function field to leave only function arguments
  fn_args$fn = NULL
  
  # We may want to append the values to evaluate the function 'along'
  if (!is.null(along))
    fn_args = append(along, fn_args)
  
  # Collapse all key-value arguments into a comma-seperated string
  args_string = paste0(names(fn_args), " = ", fn_args, collapse = ", ")
  
  # Concatenate function name with arguments
  fn_call = paste0(fn, "(", args_string, ")")
  
  # Either return that string evaluated...
  if (evaluate == TRUE) {
    fn_eval = eval_str(fn_call)
    
    return(fn_eval)
    
    # ... or simply the function string
  } else {
    return(fn_call)
  }
}

# -------------------------------------------------------------------------- -
# Parse user-defined scenarios - in a separate function for readability ----
# -------------------------------------------------------------------------- -
parse_scenarios = function(o, y, scenario) {

  # Append scenario ID to main y list, simply for reference
  y$.id = scenario
  
  # Trivial process for the baseline scenario - just store 'name'
  if (scenario == "baseline") {
    y$.name = y$scenarios$baseline$name
    
  } else {  # Otherwise, work to do...
    # Extract short and long names of all scenarios we've defined
    scenarios_id   = names(y$scenarios)
    scenarios_name = unlist(lapply(y$scenarios, function(x) x$name))
    
    # ---- Sanity checks on scenario definitions ----
    # Any scenarios for which we do not have a unique 'name' item
    unset_names = scenarios_id[!scenarios_id %in% names(scenarios_name)]
    duplicate_names = scenarios_id[duplicated(scenarios_name)]
    
    # Throw an error if any names are missing
    if (length(unset_names) > 0)
      stop("! All scenarios must have a 'name' item: ", paste0(unset_names, collapse = ", "))

    # Throw an error if any names are duplicated
    if (length(duplicate_names) > 0)
      stop("! All scenarios must have a unique 'name': ", paste0(duplicate_names, collapse = ", "))
    
    # If we only want to read scenario names, return out here
    if (scenario %in% c("*read*", "*create*"))
      return(list(scenario_names = scenarios_name))
    
    # Check the selected scenario actually exists
    if (!scenario %in% scenarios_id) {

      # Report that such a scenario isn't defined
      stop("Scenario '", scenario, "' not recognised")
    }
    
    # ---- Apply scenario values ----
    
    # Trivial process for the baseline scenario
    if (scenario != "baseline") {
     
      # Index of the alternative scenario we want to apply
      scenario_idx = which(scenarios_id == scenario)
      
      # Remove id and name items to leave only parameters to overwrite with
      y_scenario = list.remove(y$scenarios[[scenario_idx]], c("id", "name"))

      # Use overwrite function to check and apply item values
      list[y, ] = overwrite_defaults(y, y_scenario)
    }
    
    # Keep reference to this scenario 'name'
    y$.name = scenarios_name[[scenario]]
  }
  
  # Can now safely remove scenarios field
  y$scenarios = NULL
  
  return(y)
}
