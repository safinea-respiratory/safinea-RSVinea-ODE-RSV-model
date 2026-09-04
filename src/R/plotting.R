########################################################## #
# PLOTTING
#
# All plotting functions in one place.
#
########################################################## #

# -------------------------------------------------------- -
# Plot simple output ----
# -------------------------------------------------------- -
plot_simple_output = function(o, result, fig_name) {
  
  message("* Plotting simple output")
  
  output_df = result$output %>% setDT() 
  
  metric_levels = c("cases", "hospital_admissions", "deaths")
  
  # Plot 1: Model output, by age group, plotted against data
  plot1_df = output_df %>% 
    filter(metric %in% metric_levels) %>%
    
    # Summarise over age group and variant
    group_by(time, metric) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%
    
    # Factorise metric to set the correct order
    mutate(metric = factor(metric, metric_levels)) 
  
  # Define colour palette
  o$palette_simple = viridis(length(metric_levels))
  
  g1 = ggplot(plot1_df, aes(x = time, y = value, colour = metric)) + 
    
    # Line plot
    geom_line(show.legend = FALSE, alpha = 0.6) +
    
    # Facet the plot
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +  
    
    # Minimal theme
    theme_minimal() +
    
    # Define colour for data source
    scale_colour_manual(values = o$palette_simple) +  
    
    # Adjust legend title and size
    labs(colour = NULL) +  # Change legend title
    theme(legend.key.size = unit(0.1, "cm"),   # Reduce legend key size
          legend.text = element_text(size = 4),  # Reduce legend text size
          legend.title = element_text(size = 4)) # Adjust legend title size
  
  #---- Improve aesthetics ----
  # Prettify plot 1
  g1 = g1 + theme_classic() + 
    theme(strip.text   = element_text(size = 10), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 8), 
          axis.text    = element_text(size = 8), 
          axis.line    = element_blank(), 
          panel.border = element_rect(linewidth = 1, colour = "black", fill = NA), 
          strip.background = element_blank())
  
  if (!is.null(fig_name))
    fig_save(o, g1, paste0(fig_name, "_test"), height = o$save_height/4)
  
  return()
}

# ------------------------------------------------------------------------ -
# Plot metrics over time for multiple metrics, groups, and/or scenarios ----
# ------------------------------------------------------------------------ -
plot_scenarios = function(o, fig_name, ...) {
  
  # Collate and interpret inputs so we know what to plot
  #list[f, baseline] = fig_properties(o, list(...))
  f = NULL
  # Full list of scenarios as defined in yaml file 
  f$scenarios = parse_yaml(o, "*read*") %>% names()
  
  # Full scenario names as defined in yaml file 
  f$scenario_names = parse_yaml(o, "*read*") %>% unname()
  
  # ---- Extract model predictions ----
  # Initiate plotting dataframe
  plot_list = list()
  
  # Loop through scenarios to plot
  for (scenario in f$scenarios) {
    
    #result = try_load(o$pth$scenarios, paste0(scenario, '_100k'))
    result = try_load(o$pth$scenarios, paste0(scenario, '_raw'))
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, result)  
  }
  
  # Concatenate plotting dataframes for all scenarios
  plot_df = rbindlist(plot_list) 
  
  # ---- Apply dates ----
  
  # Formatting the dates
  duration = plot_df %>% select(time) %>% unique() %>% nrow()
  all_dates = seq(from = ymd(o$plot_start_date), by = "day", length.out = duration) 
  dates_df  = data.table(date = all_dates, 
                         time  = 1 : length(all_dates))
  
  # Define correct order of age_group and metric
  age_levels = c("0-1m", "1-2m", "2-3m", "3-4m", "4-5m", "5-6m", 
                 "6-7m", "7-8m", "8-9m", "9-10m", "10-11m", "11-12m", 
                 "12-13m", "13-14m", "14-15m", "15-16m", "16-17m", "17-18m", 
                 "18-19m", "19-20m", "20-21m", "21-22m", "22-23m", "23-24m", 
                 "24-30m", "30-36m", "3-4y",  "4-5y", "5-18y", 
                 "18-65y",
                 "65+y")
  metric_levels = c("hospital_admissions")
  
  
  # Add dates
  plot_df = plot_df %>% inner_join(dates_df, by = "time") %>%
    mutate(date = as.Date(date),
           scenario = factor(scenario, levels = f$scenarios)) %>% # Sort scenarios in legend
    select(-time)
  
  # Compute daily, monthly, and total metric
  plot_df = aggregate_model_output(plot_df, NA, date_col = "date")
  
  # Compute "total" age group
  plot_df_total = plot_df %>% 
    group_by(metric, variant, param_id, scenario, date, data_freq) %>%
    summarise(value = sum(value),
              age_group = "total") %>%
    ungroup()
  
  # Merge "total" age group with others
  plot_df = bind_rows(plot_df, plot_df_total)
  
  # Compute ribbon
  plot_df = plot_df %>% 
    group_by(age_group, metric, variant, scenario, date, data_freq) %>%
    summarise(mean   = mean(value),
              #median = quantile(value, 0.5, na.rm = TRUE),
              lower  = quantile(value, o$quantiles[1], na.rm = TRUE),
              upper  = quantile(value, o$quantiles[2], na.rm = TRUE),
              .groups = "drop")
  
  age_group_df = plot_df %>% filter(age_group %in% age_levels) %>%
    filter(metric %in% metric_levels,
           data_freq == "weekly") %>%
    # Factorise age_group and metric to set the correct order
    mutate(age_group = factor(age_group, age_levels),  
           metric = factor(metric, metric_levels)) %>%
    # Create facet_label
    mutate(facet_label = paste(age_group, metric, sep = " ")) %>%
    arrange(age_group, metric) %>%
    mutate(facet_label = factor(facet_label, levels = unique(facet_label))) %>%
    setDT()
  
  total_df = plot_df %>%
    filter(metric %in% metric_levels,
           data_freq == "weekly",
           age_group == "total") %>%
    # Create facet_label
    mutate(facet_label = paste(age_group, metric, sep = " ")) %>%
    mutate(facet_label = factor(facet_label, levels = unique(facet_label))) %>%
    # Factorise metric to set the correct order
    mutate(metric = factor(metric, metric_levels)) %>%
    setDT()
  
  scenario_labels = setNames(f$scenario_names, f$scenarios)
  
  # Check if any names are still missing - throw an error if so
  missing_names = f$scenarios[is.na(f$scenario_names)]
  if (length(missing_names) > 0)
    stop("Scenario names not recognised: \n", paste(missing_names, collapse = "\n"))
  
  
  
  # Plot 1 (Scenarios by age group)
  plot1_df = age_group_df %>%
    mutate(scenario = factor(scenario, levels = f$scenarios))
  
  # Plot 2 (Scenarios summarised over all age groups)
  plot2_df = total_df %>%
    arrange(metric)
  
  # Plot 4 (Baseline projection)
  plot4_df = age_group_df %>% 
    filter(scenario == "baseline")
  
  # Plot 5 (also age_group 'all')
  age_levels_full <- c(age_levels,'total')
  plot5_df = plot_df %>%
    filter(metric %in% metric_levels) %>%
    # Factorise age_group and metric to set the correct order
    mutate(age_group = factor(age_group, age_levels_full),
           metric = factor(metric, metric_levels)) %>%
    # Create facet_label
    mutate(facet_label = paste(age_group, metric, sep = " ")) %>%
    arrange(age_group, metric) %>%
    mutate(facet_label = factor(facet_label, levels = unique(facet_label))) %>%
    setDT()
  
  # ---- Load fitting data ----
  data_df  = get_target_data(o)
  
  # Currently only the all-age-group aggregate (g2) is generated by default.
  # The additional plots below are available but not saved on each run —
  # uncomment as needed:
  #   g1  = scenarios by fine age group
  #   g1b = g1 zoomed to o$plot_zoom_date onwards
  #   g3  = deaths (if deaths metric is enabled)
  #   g4  = baseline-only projection

  # Plot 2: Scenarios summarised over all age groups
  g2 = get_scenario_ggplot(plot2_df, data_df %>% filter(data_freq == "weekly"), scenario_labels)

  # Save to output
  if (!is.null(fig_name)) {
    fig_save(o, g2, paste0(fig_name, "_aggregated"), height = o$save_height / 4)
  }
  
}

# -------------------------------------------------------- -
# Calibration plots ----
# -------------------------------------------------------- -
plot_best_samples = function(o, fit, fig_name, round_idx) {
  # Percentage of best parameter sets to plot
  p_sets = c(1)
  
  # Whether to only plot density for supported limits
  density_trim = FALSE
  
  # ---- Load samples, model output, and data ----
  # Load all samples, their likelihood, and associated model output
  sets_df = try_load(o$pth$fitting, paste0(round_idx, "_samples"))
  model_df  = try_load(o$pth$fitting, paste0(round_idx, "plot_output"))
  
  # # Read in target data
  data_df <- get_target_data(o)
  
  # Compute daily, monthly, and total metric
  model_df = aggregate_model_output(model_df, data_df, date_col = "date")
  
  # Extract parameter bounds
  param_df = data.table(param = fit$params, fit$bounds)
  
  # Define correct order of metric
  metric_levels = c("hospital_admissions")
  
  # ---- Create plot for each set of samples ----
  
  # Likelihood limits (used for colouring)
  likelihood_limits = c(0, max(sets_df$likelihood))
  
  # Initiate plotting list
  plot_list = list()
  
  # Iterate through number of sample sets to plot
  for (i in seq_along(p_sets)) {
    
    # Number of parameter sets to plot on this iteration
    n_sets =  max(10, round(nrow(sets_df) * p_sets[i] / 100))
    
    # IDs of parameter sets to plot on this iteration
    plot_id = sets_df %>%
      slice_max(likelihood, n = n_sets) %>%
      pull(param_id)
    
    # Plot 1: Model output, total age group, plotted against data
    plot1_df = model_df %>%
      ungroup() %>%
      mutate(date = as.Date(date)) %>%
      filter(param_id %in% plot_id) %>%
      left_join(sets_df, by = c("param_id", "round")) %>%
      select(param_id, likelihood, age_group, metric, date, value, data_freq) %>%
      # Filter total age group and correct metric levels
      filter(age_group == "total") %>%
      filter(data_freq == "weekly") %>%
      filter(metric %in% metric_levels) %>%      
      # Factorise age_group and metric to set the correct order
      mutate(metric = factor(metric, levels = metric_levels)) %>%
      # Create facet_label
      mutate(facet_label = paste(age_group, metric, sep = " ")) %>%
      arrange(age_group, metric) %>%
      mutate(facet_label = factor(facet_label, levels = unique(facet_label))) %>%
      mutate(mean = value,
             scenario = param_id)
    
    g1 = get_scenario_ggplot(plotX_df = plot1_df,
                             dataX_df = data_df %>%
                               filter(age_group == "total",
                                      data_freq == "weekly"))
    
    # Plot total burden per age group
    plot1A_df = model_df %>%
      ungroup() %>%
      filter(data_freq == "4-weekly",
             param_id %in% plot_id) %>%
      left_join(sets_df, by = c("param_id", "round")) %>%
      select(param_id, likelihood, age_group, metric, value, data_freq, date) %>%
      # Filter total age group and correct metric levels
      filter(metric %in% metric_levels) %>%
      # Factorise age_group and metric to set the correct order
      mutate(metric = factor(metric, levels = metric_levels)) %>%
      # Create facet_label
      mutate(facet_label = paste(age_group, metric, date, sep = " ")) %>%
      arrange(age_group, metric) %>%
      mutate(facet_label = factor(facet_label, levels = unique(facet_label))) %>%
      mutate(mean = value,
             scenario = param_id)
    
    # --- (optional) keep only the metric/frequency you want ---
    model_use <- plot1A_df 
    
    data_use <- data_df %>%
      filter(data_freq == "4-weekly")
    
    # Order age groups (adjust if your set differs)
    age_order <- c("0-3m","3-6m","6-12m","1-5y","5-65y","65+y")
    
    # 1) Summarise model to one value per param_id x age_group (mean here)
    model_sum <- model_use %>%
      group_by(param_id, age_group, date) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      mutate(source = "Model")
    
    # 2) Replicate observed data for each param_id so it shows in every facet
    observed_rep <- data_use %>%
      select(age_group, value, date) %>%
      mutate(source = "Observed") %>%
      crossing(model_use %>% distinct(param_id))
    
    # 3) Combine for plotting
    plot_df <- bind_rows(model_sum, observed_rep) %>%
      filter(age_group %in% age_order) %>%
      mutate(age_group = factor(age_group, levels = age_order))
    
    # 4) Plot: side-by-side (tight) bars, faceted by param_id
    ggplot(plot_df, aes(x = age_group, y = value, fill = source)) +
      geom_col(position = position_dodge(width = 0.6), width = 0.5) +
      facet_wrap(~ param_id) +                 # add scales="free_y" if needed
      labs(x = "Age group", y = "Number", fill = "Source") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 5) a second plot
    # Compute ratio per param_id & age_group
    ratio_df <- plot_df %>%
      filter(source %in% c("Model", "Observed")) %>%
      group_by(age_group, date) %>%
      summarise(
        mean_model = mean(value[source == "Model"]),
        observed   = value[source == "Observed"][1],
        .groups = "drop"
      ) %>%
      mutate(ratio = observed / mean_model )
    
    summ <- plot_df %>%
      filter(source == "Model") %>%
      group_by(age_group, date) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        sd_value   = sd(value, na.rm = TRUE),
        quant_50 = quantile(value,probs = 0.5),
        quant_5 = quantile(value,probs = 0.05),
        quant_95 = quantile(value,probs = 0.95),
        .groups = "drop"
      )
    
    age_labels <- c(
      "0-3m"   = "<3 months",
      "3-6m"  = "3-6 months",
      "6-12m"  = "6-12 months",
      "1-5y"  = "1-5 years",
      "5-65y" = "15-64 years",
      "65+y"   = "65+ years"
    )
    
    # Plot
    g1a = 
      ggplot(summ, aes(x = date, y = quant_50)) +
      geom_point() +
      geom_errorbar(aes(ymin = quant_5,
                        ymax = quant_95),
                    width = 0.2) +
      facet_wrap(~ age_group, scales = "free", labeller = labeller(age_group = age_labels)) +
      # Points for Observed
      geom_point(
        data = plot_df %>% filter(source == "Observed") %>% group_by(age_group, date) %>% slice(1),
        aes(x = date, y = value),
        color = "darkred", size = 2,
        position = position_dodge(width = 0.8)
      ) +
      # # Ratio as text
      # geom_text(
      #   data = ratio_df,
      #   aes(y = observed, label = sprintf("%.2f", ratio)),
      #   hjust = -0.5, vjust = 0, color = "black", size = 3
      # ) +
      #ylim(c(0,NA)) +
      coord_cartesian(ylim = c(0, NA)) + 
      labs(x = "Time", y = "4-weekly burden", fill = "Source") +
      #theme_minimal() +
      #scale_y_log10() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    # Plot seasonal burden
    summ_season = plot_df %>%
      filter(source == "Model") %>%
      mutate(season_year = if_else(month(date) >= 8, year(date), year(date) - 1)) %>%
      # Sum across dates in a given season
      group_by(age_group, season_year, param_id) %>%
      summarise(
        season_value = sum(value),
        .groups = "drop"
      ) %>%
      # Mean and sd over samples
      group_by(age_group, season_year) %>%
      summarise(
        mean_value = mean(season_value, na.rm = TRUE),
        sd_value   = sd(season_value, na.rm = TRUE),
        quant_50 = quantile(season_value,probs = 0.5),
        quant_5 = quantile(season_value,probs = 0.05),
        quant_95 = quantile(season_value,probs = 0.95),
        .groups = "drop"
      ) %>%
      mutate(season = paste0(season_year,'/',season_year+1))
      
    
    g1b = ggplot(summ_season, aes(x = season, y = quant_50)) +
      geom_point() +
      geom_errorbar(aes(ymin = quant_5,
                        ymax = quant_95),
                    width = 0.2) +
      facet_wrap(~ age_group, scales = "free", labeller = labeller(age_group = age_labels)) +
      # Points for Observed
      geom_point(
        data = plot_df %>% 
          filter(source == "Observed") %>% 
          mutate(season_year = if_else(month(date) >= 8, year(date), year(date) - 1),
                 season = paste0(season_year,'/',season_year+1)) %>%
          group_by(age_group, param_id, season) %>% 
          summarise(
            value = sum(value, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          group_by(age_group, season) %>% 
          slice(1),
        aes(x = season, y = value),
        color = "darkred", size = 3,
        position = position_dodge(width = 0.8)
      ) +
      # # Ratio as text
      # geom_text(
      #   data = ratio_df,
      #   aes(y = observed, label = sprintf("%.2f", ratio)),
      #   hjust = -0.5, vjust = 0, color = "black", size = 3
      # ) +
      ylim(c(0,NA)) +
      labs(x = "Season", y = "Seasonal burden", fill = "Source") +
      #theme_minimal() +
      #scale_y_log10() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ###############################################################
    ###############################################################
    
    
    # Plot 2: Model output, total age group, plotted against data with geom_ribbon
    plot2_df = model_df %>% ungroup() %>%
      mutate(date = as.Date(date)) %>%
      
      # Select best samples
      filter(param_id %in% plot_id) %>%
      left_join(sets_df, by = c("param_id", "round")) %>%
      select(param_id, likelihood, age_group, metric, date, value, data_freq) %>%
      filter(data_freq == "weekly") %>%
      group_by(age_group, metric, date) %>%
      
      # Find quantiles
      summarise(mean   = mean(value),
                median = quantile(value, 0.5, na.rm = TRUE),
                lower  = quantile(value, o$quantiles[1], na.rm = TRUE),
                upper  = quantile(value, o$quantiles[2], na.rm = TRUE),
                .groups = "drop") %>%
      
      filter(age_group == "total") %>%
      filter(metric %in% metric_levels) %>%
      
      # Create facet_label
      mutate(facet_label = paste(age_group, metric, sep = " ")) %>%
      arrange(age_group, metric) %>%
      mutate(facet_label = factor(facet_label, levels = unique(facet_label)))  %>%
      
      # Add scenario label
      mutate(scenario = "calibration")
    
    g2 = get_scenario_ggplot(plotX_df = plot2_df, 
                             dataX_df = data_df %>% 
                               filter(age_group == "total",
                                      data_freq == "weekly"),
                             scenario_labels = c(calibration="calibration"))
    
    
    # Plot 3: Posterior distributions
    plot3_df = sets_df %>%
      filter(param_id %in% plot_id) %>%
      select(likelihood, all_of(fit$params)) %>%
      unique() %>%
      pivot_longer(cols = -likelihood, 
                   names_to = "param") %>%
      setDT()
    
    g3 = ggplot(plot3_df) + 
      geom_segment(data    = param_df, 
                   mapping = aes(x = lower, 
                                 y = 0,     
                                 xend = upper, 
                                 yend = 0), 
                   alpha = 0) +
      geom_density(mapping = aes(x = value, 
                                 y = after_stat(scaled), 
                                 weight = likelihood, 
                                 fill   = param, 
                                 colour = param), 
                   alpha = 0.5, trim = density_trim, show.legend = FALSE) + 
      facet_wrap(~param, nrow = 5, scales = "free_x") + 
      scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    
    
    
    
    #---- Improve aesthetics ----
    # Prettify plot 3
    g3 = g3 + theme_classic() + 
      theme(strip.text   = element_text(size = 12), 
            axis.title   = element_blank(), 
            axis.text.x  = element_text(size = 8), 
            axis.text.y  = element_text(size = 8), 
            axis.line    = element_blank(), 
            panel.border = element_rect(linewidth = 1, colour = "black", fill = NA), 
            strip.background = element_blank())
    
  }
  
  #if (fit$input$adaptive_sampling$rounds == as.integer(str_split(round_idx,"r")[[1]][[2]])) browser()
  
  if (!is.null(fig_name))
    fig_save(o, g1, paste0(fig_name, "_fit"), round_idx)
  
  if (!is.null(fig_name))
    fig_save(o, g2, paste0(fig_name,"_ribbon"), round_idx)
  
  if (!is.null(fig_name))
    fig_save(o, g3, paste0(fig_name,"_posteriors"), round_idx)
  
  if (!is.null(fig_name))
    fig_save(o, g1a, paste0(fig_name,"_age_burden"), round_idx)
  
  if (!is.null(fig_name))
    fig_save(o, g1b, paste0(fig_name,"_seasonal_age_burden"), round_idx)
  
  return()
}

# --------------------------------------------------------------------- -
# Generate a ggplot object based on the provided plot_df and data_df ----
# --------------------------------------------------------------------- -
get_scenario_ggplot = function(plotX_df, dataX_df, scenario_labels = NULL, x_min = NULL){
  
  # Option to adjust x-limits
  if(!is.null(x_min)){
    plotX_df = plotX_df %>% filter(date > ymd(x_min))
    if(!is.null(dataX_df)) dataX_df = dataX_df %>% filter(date > ymd(x_min))
  }
  
  # check if scenario is part of plotX_df, and resolve if not
  if(!"scenario" %in% names(plotX_df)){
    plotX_df= plotX_df %>% mutate(scenario = NA)
  }
  
  # define whether a ribbon needs to be included
  bool_add_ribbon <- ("lower" %in% names(plotX_df) && "upper" %in% names(plotX_df))
  
  # Define number of age groups (rows)
  num_age <- length(unique(plotX_df$age_group))
  num_scenario = length(unique(plotX_df$scenario))
  
  # improve facet notation
  plotX_df$facet_label <- improve_facet_notation(plotX_df$facet_label)
  dataX_df$facet_label <- improve_facet_notation(dataX_df$facet_label)
  
  # Check common age groups
  if(!any(dataX_df$age_group %in% plotX_df$age_group)){
    dataX_df <- NULL
  }
  
  # init color param
  scale_colour_values <- c()
  scale_colour_labels <- c()
  
  # Start figure
  gX = ggplot(plotX_df, aes(x = date, y = mean, colour = scenario, fill = scenario)) 
  
  # Option to add reference data
  if(!is.null(dataX_df)){
    # Points with colour based on data source (observed / synthetic)
    gX = gX + geom_point(data = dataX_df, 
                         aes(x = date, y = value, colour = factor(source)),  # explicitly redefine x
                         size = 1,
                         inherit.aes = FALSE)  # avoid inheriting global aesthetics
    # Define colour for data source
    ###scale_colour_values = c(scale_colour_values,observed="gray60", synthetic="darkred")
    ###scale_colour_labels = c(scale_colour_labels, observed="observed", synthetic="synthetic")
  }
  
  # Add scenario projections
  if(bool_add_ribbon){
    gX = gX +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA, linewidth = 0.0) +
      geom_line(linewidth = 0.75)
  } else{
    gX = gX +
      ###geom_line(linewidth = 0.50, color = "black")
      geom_line(linewidth = 0.50)
  } 
  
  # Facet and layout the plot
  gX = gX +
    facet_wrap(~ facet_label,
               scales = "free_y", nrow = num_age) +  # Force factor usage here
    
    # Minimal theme
    theme_minimal() +
    
    # Adjust x and y axes
    scale_x_date(date_labels = "%d/%m/%y") + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), 
                       name = "per 100,000 population") +
    
    # Adjust legend title and size
    theme(legend.text = element_text(size = 4),  # Reduce legend text size
          legend.title = element_text(size = 4))  # Adjust legend title size
  
  # Note: scenario colours currently use ggplot2 defaults. For explicit control,
  # use scale_color_manual() / scale_fill_manual() with turbo() for multiple
  # scenarios or "black" for a single scenario.
  
  
  # resolve issue when data points are there
  if(!is.null(dataX_df)){
    gX = gX + guides(fill = "none") # Hide the colour legend
  }
  # Option to remove legend when scenario_labels is NULL
  if(length(unique(plotX_df$scenario)) == 1 || is.null(scenario_labels)){
    gX = gX + labs(colour = NULL) 
  }
  # Option to remove legend overall
  if(length(unique(plotX_df$scenario)) == 1){
    gX = gX + guides(colour = "none", fill = "none")
  }
  
  # Prettify plot
  gX = gX + theme_classic() + 
    theme(strip.text   = element_text(size = 8), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 8), 
          axis.text    = element_text(size = 7), 
          axis.line    = element_blank(), 
          panel.grid.major = element_line(color = "gray80", linewidth = 0.2),
          panel.grid.minor = element_line(color = "gray80", linewidth = 0.2),
          panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA), 
          strip.background = element_blank()) 
  
  return(gX)
}  


# --------------------------------------------------------- -
# Perform some textual enhancements to the face labels ----
# --------------------------------------------------------- -
improve_facet_notation = function(facet_label){
  
  facet_label <- gsub(' ','y ',facet_label)
  facet_label <- gsub('_%',' (%)',facet_label)
  facet_label <- gsub('n_','number of ',facet_label)
  facet_label <- gsub('_',' ',facet_label)
  
  return(factor(facet_label, levels = unique(facet_label)))
}

# ---------------------------------------------------------------- -
# Load target data while accounting for calibration weights ----
# ---------------------------------------------------------------- -
get_target_data = function(o){
  
  # First try to load pre-generated data
  data_df <- try_load(o$pth$fitting, "plot_data", throw_error = FALSE)
  if(!is.null(data_df)){
    return(data_df)
  }
  
  # Else, generate and store
  
  # Initiate fit list and perform a few checks on input .yaml
  fit = setup_calibration(o)
  
  # Load data - see load_data.R
  fit = load_data(o, fit)
  
  # Shorthand for calibration weights list
  w = fit$input$calibration_weights
  
  # Weights for each metric
  weight_df = as.data.table(w$metric) %>% 
    pivot_longer(cols = everything(), 
                 names_to  = "metric", 
                 values_to = "weight")
  
  # Labelling as observed / synthetic
  source_df = as.data.table(w$metric) %>% 
    pivot_longer(cols = everything(), 
                 names_to  = "metric", 
                 values_to = "weight") %>%
    filter(!weight == 0) %>%
    rename(source = weight) %>%
    mutate(source = case_when(metric == "cases" ~ fit$opts$cases,
                              metric == "hospital_admissions" ~ fit$opts$hospital_admissions,
                              metric == "deaths" ~ fit$opts$deaths,
                              TRUE ~ "missing"))
  # Read in target data
  data_df = fit$data %>%
    filter(!is.na(value)) %>%
    
    # # Filter out metrics not included in fitting
    left_join(weight_df, by = "metric") %>%
    filter(weight > 0) %>%
    group_by(age_group, metric, date, data_freq) %>%
    summarise(value = sum(value)) %>%

    # Label observed / synthetic
    left_join(source_df, by = "metric") %>%
    
    # Match labelling of model output
    mutate(facet_label = paste(age_group, metric, sep = " ")) %>%
    setDT()
  
  # Write file for later plotting purposes
  saveRDS(data_df, file = paste0(o$pth$fitting, "plot_data.rds"))
  
  return(data_df)
  
}

# -------------------------------------------------------- -
# Save a ggplot figure to file with default settings ----
# -------------------------------------------------------- -
fig_save = function(o, g, ..., path = "results", width = o$save_width, height = o$save_height) {
  
  # Collapse inputs into vector of strings
  fig_name_parts = unlist(list(...))
  
  # Construct file name to concatenate with file path
  save_name = paste(fig_name_parts, collapse = " - ")
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in o$figure_format) {
    save_pth  = paste0(o$pth[[path]], save_name, ".", fig_format)
    
    # Save figure (size specified in options.R)
    ggsave(save_pth, 
           plot   = g, 
           device = fig_format, 
           dpi    = o$save_resolution, 
           width  = width, 
           height = height, 
           units  = o$save_units)
  }
}

# -------------------------------------------------------- -
# Save a plotting dataframe to file ----
# -------------------------------------------------------- -
save_dataframe = function(o, plot_df, fig_name) {
  
  # Give a generic name if otherwise unnamed figure
  if (is.null(fig_name))
    fig_name = "Untitled"
  
  # A bit of formatting to remove illegal chars
  save_name = fig_name %>%
    paste(collapse = "_") %>%
    str_replace_all(" ", "_") %>%
    paste0(".rds")
  
  # Save plotting dataframe in figure folder
  saveRDS(plot_df, file = paste0(o$pth$results, save_name))
}

#------------------------------------------------------ -
# Format model outcomes ready for plotting ----
# ----------------------------------------------------- -
format_results = function(o, results) {
  
  # Reduce model output down to what we're interested in
  output_df = results %>%
    filter(is.na(time) | time >= o$plot_from, 
           is.na(time) | time <= o$plot_to)
  
  # ---- Format output ----
  # Final formatting touches
  output_df = output_df %>%
    arrange(scenario, metric, age_group, time) 
  
  return(output_df)
}



