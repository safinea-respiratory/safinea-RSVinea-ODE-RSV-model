########################################################## #
# AUXILIARY FUNCTIONS
#
# A series of helpful R functions.
#
########################################################## #

# -------------------------------------------------------- -
# Remap age group labels using a named mapping ----
# -------------------------------------------------------- -
remap_age_groups = function(df, mapping, col = "age_group") {
  df %>%
    mutate(!!col := coalesce(unlist(mapping)[.data[[col]]], .data[[col]]))
}

# -------------------------------------------------------- -
# Canonical fine -> reporting age-group lookup ----
# Single source of truth: age_group_map in default.yaml, so the
# fine-to-reporting mapping is never hard-coded per script.
# Returns a named character vector (names = fine groups, values = bands).
# -------------------------------------------------------- -
age_band_lookup = function(o) {
  unlist(yaml::read_yaml(o$pth$params_default)$age_group_map)
}

# -------------------------------------------------------- -
# Set as datatable and rename columns in one line ----
# -------------------------------------------------------- -
as_named_dt = function(x, new_names) {
  
  # Convert to datatable
  dt = as.data.table(x)
  
  # Check new names are correct length
  old_names = names(dt)
  if (length(old_names) != length(new_names))
    stop("Inconsistent number of column names provided")
  
  # Set new column names
  named_dt = setnames(dt, old_names, new_names)
  
  return(named_dt)
}

# -------------------------------------------------------- -
# Clear the console ----
# -------------------------------------------------------- -
clc = function() cat("\014")

# -------------------------------------------------------- -
# Clear all figures ----
# -------------------------------------------------------- -
clf = function() graphics.off()

# -------------------------------------------------------- -
# Create colour scheme ----
# -------------------------------------------------------- -
colour_scheme = function(map, pal = NULL, n = 1, ...) {
  
  # Has colour palette been defined
  if (is.null(pal)) {
    
    # That's ok as long as it's defined within the map argument
    if (!grepl("::", map))
      stop("Palette not defined - Use 'pal = my_pal' or 'map = my_map::my_pal'")
    
    # Separate out the map and the palette
    pal = str_remove(map, ".*\\::")
    map = str_remove(map, "\\::.*")
  }
  
  # Initiate colours variable
  colours = NULL
  
  # Built in colour schemes
  if (map == "base")
    colours = get(pal)(n, ...)	
  
  # A load of colour maps from the pals package
  #
  # See: https://www.rdocumentation.org/packages/pals/versions/1.6
  if (map == "pals")
    colours = get(pal)(n, ...)
  
  # Stylish HCL-based colour maps
  #
  # See: https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
  if (grepl("_hcl$", map))
    colours = get(map)(palette = pal, n = n, ...)
  
  # Colour Brewer colour schemes
  if (map == "brewer")
    colours = brewer_pal(palette = first_cap(pal), ...)(n)
  
  # Viridis colour schemes
  if (map == "viridis")
    colours = viridis_pal(option = pal, ...)(n)
  
  # Throw an error if colours not yet defined
  if (is.null(colours))
    stop("Colour map '", map, "' not recognised (supported: base, pals, hcl, brewer, viridis)")
  
  return(colours)
}


# -------------------------------------------------------------------- -
# Evaluate a string (in calling function environment) using eval ----
# -------------------------------------------------------------------- -
eval_str = function(...)
  eval(parse(text = paste0(...)), envir = parent.frame(n = 1))


# -------------------------------------------------------- -
# Format heterogeneous styles of dates ----
# -------------------------------------------------------- -
format_date = function(dates, convert = "ymd") {
  styles = c("dmy", "dmY", "ymd", "Ymd")
  dates = parse_date_time(dates, styles)
  dates = get(convert)(dates)
  return(dates)
}

# -------------------------------------------------------- -
# Convert list to datatable ----
# -------------------------------------------------------- -
list2dt = function(x, ...) {
  dt = rbindlist(lapply(x, as.data.table), ...)
  return(dt)
}

# -------------------------------------------------------- -
# Suppress output from a function call ----
# -------------------------------------------------------- -
quiet = function(x) { 
  sink_con = file("sink.txt")
  sink(sink_con, type = "output")
  sink(sink_con, type = "message")
  on.exit(sink(type   = "output"))
  on.exit(sink(type   = "message"), add = TRUE)
  on.exit(file.remove("sink.txt"),  add = TRUE)
  invisible(force(x)) 
}

# ---------------------------------------------------------------------------- -
# Wrapper for consistent behaviour of base::sample when length(x) is one ----
# ---------------------------------------------------------------------------- -
sample_vec = function(x, ...) x[sample.int(length(x), ...)]

# -------------------------------------------------------- -
# Initiate progress bar with normal-use options ----
# -------------------------------------------------------- -
start_progress_bar = function(n_tasks) {
  pb = txtProgressBar(min = 0, max = n_tasks,
                      initial = 0, width = 100, style = 3)
  return(pb)
}

# ---------------------------------------------------------- -
# Convert comma-separated string to vector of elements ----
# ---------------------------------------------------------- -
str2vec = function(x, v) {
  x[[v]] = x[[v]] %>% 
    str_split(",", simplify = TRUE) %>% 
    str_remove_all(" ")
  return(x)
}

# ---------------------------------------------------------- -
# Bi-directional setdiff - elements not in both x and y ----
# ---------------------------------------------------------- -
symdiff = function(x, y) setdiff(union(x, y), intersect(x, y))

# -------------------------------------------------------- -
# Format a number with thousand mark separators ----
# -------------------------------------------------------- -
thou_sep = function(val) {
  format_val = format(val, scientific = FALSE,
                      trim = TRUE, 
                      drop0trailing = TRUE, 
                      big.mark = ",")
  return(format_val)
}

# -------------------------------------------------------- -
# Load a file if it exists, throw an error if not ----
# -------------------------------------------------------- -
try_load = function(pth, file, msg = NULL, type = "rds", throw_error = TRUE, sep = FALSE) {
  
  # Initiate trivial output
  file_contents = NULL
  
  # Set default error message
  if (is.null(msg))
    msg = "Cannot load file"
  
  # Switch case for loading function
  loading_fnc = switch(
    tolower(type), 
    
    # Support both RDS and CSV
    "rds" = "readRDS", 
    "csv" = "read.csv",
    
    # Throw an error if anything else requested
    stop("File type '", type, "' not supported")
  )
  
  # Concatenate path and file name
  file_name = paste0(pth, ifelse(sep, file_sep(), ""), file, ".", type)
  
  # If file doesn't exist, throw an error if desired
  if (!file.exists(file_name) && throw_error == TRUE)
    stop(msg, " [missing: ", file_name, "]")
  
  # If file exists, try to load it
  if (file.exists(file_name)) {
    
    # Get the loading function and attempt to load file
    file_contents = tryCatch(
      get(loading_fnc)(file_name),
      
      # Catch the error - we may not want to throw it
      error = function(e) {
        
        # Throw descriptive error if desired
        if (throw_error == TRUE) 
          stop(msg, " [unreadable: ", file_name, "]")
      }
    )
  }
  
  return(file_contents)
}

# -------------------------------------------------------- -
# Un-list and return names with separator of choice ----
# -------------------------------------------------------- -
unlist_format = function(x, sep = "$", ...) {
  y = unlist(x, ...)
  names(y) = gsub("\\.", sep, names(y))
  return(y)
}


# ---------------------------------------------- -
# Change age brackets of the contact matrix ----
# ---------------------------------------------- -
#-- New bands for contact matrix -## 
reband_contact_matrix <- function(mat, original_breaks, target_breaks, pop = NULL) {
  
  # Change years to months
  original_breaks = 12*original_breaks
  target_breaks = 12*target_breaks
  
  # Expand to 1-month age bands
  n_orig <- nrow(mat)
  age_bins <- data.frame(
    group = 1:n_orig,
    lower = original_breaks[-length(original_breaks)],
    upper = original_breaks[-1]
  )
  
  # Build a full month-level version of the matrix
  expand_to_1m <- function(mat, age_bins) {
    n <- sum(age_bins$upper - age_bins$lower)
    full_mat <- matrix(0, nrow = n, ncol = n)
    
    row_idx <- 1
    for (i in 1:nrow(age_bins)) {
      r_len <- age_bins$upper[i] - age_bins$lower[i]
      col_idx <- 1
      for (j in 1:nrow(age_bins)) {
        c_len <- age_bins$upper[j] - age_bins$lower[j]
        full_mat[row_idx:(row_idx + r_len - 1), col_idx:(col_idx + c_len - 1)] <- mat[i, j]/r_len
        col_idx <- col_idx + c_len
      }
      row_idx <- row_idx + r_len
    }
    
    return(full_mat)
  }
  
  # Collapse from 1-month to target bands
  collapse_from_1m <- function(mat_1m, target_breaks) {
    n_new <- length(target_breaks) - 1
    collapsed <- matrix(0, n_new, n_new)
    
    month_to_group <- function(ages, breaks) findInterval(ages, breaks, rightmost.closed = TRUE)
    months <- 0:(nrow(mat_1m) - 1)
    row_map <- month_to_group(months, target_breaks)
    
    for (i in 1:n_new) {
      for (j in 1:n_new) {
        rows <- which(row_map == i)
        cols <- which(row_map == j)
        
        submat <- matrix(mat_1m[rows, cols], nrow = length(rows), ncol = length(cols))
        collapsed[i, j] <- mean(rowSums(submat)) 
        
      }
    }
    
    return(collapsed)
  }
  
  mat_1m <- expand_to_1m(mat, age_bins)
  
  new_mat <- collapse_from_1m(mat_1m, target_breaks)
  return(new_mat)
}


# ------------------------------------------------------------------ -
# Helper to convert age group label to numeric age range in years ----
# ------------------------------------------------------------------ -

# 
parse_age_range <- function(label) {
  if (label == "65+y") return(c(65, 100))
  
  parts <- str_split(label, "-")[[1]]
  # Determine the unit from the second part
  unit <- ifelse(str_detect(parts[2], "mo|m"), "m",
                 ifelse(str_detect(parts[2], "yr|y"), "y", NA))
  
  parse_part <- function(x, unit_hint) {
    if (str_detect(x, "mo|m")) as.numeric(str_remove(x, "mo|m")) / 12
    else if (str_detect(x, "yr|y")) as.numeric(str_remove(x, "yr|y"))
    else if (unit_hint == "m") as.numeric(x) / 12
    else if (unit_hint == "y") as.numeric(x)
    else as.numeric(x)  # fallback
  }
  
  c(parse_part(parts[1], unit), parse_part(parts[2], unit))
}


# ------------------------------------------------------------------ -
# Redistribute population across fine-grained age groups ----
# ------------------------------------------------------------------ -

redistribute_population <- function(coarse_df, fine_age_groups, fine_age_breaks, country_name) {
  # Filter for the country
  coarse_data <- coarse_df %>% filter(country == country_name)
  
  # Parse coarse group age ranges
  coarse_ranges <- t(sapply(coarse_data$age_group, parse_age_range))
  colnames(coarse_ranges) <- c("start", "end")
  coarse_data <- bind_cols(coarse_data, as.data.frame(coarse_ranges))
  
  # Build fine-grained age group dataframe
  fine_df <- data.frame(
    age_group = fine_age_groups,
    start = fine_age_breaks[-length(fine_age_breaks)],
    end = fine_age_breaks[-1]
  )
  
  # Redistribute
  redistributed <- lapply(1:nrow(fine_df), function(i) {
    f_start <- fine_df$start[i]
    f_end <- fine_df$end[i]
    
    overlaps <- coarse_data %>%
      filter(start < f_end & end > f_start)
    
    total_pop <- 0
    for (j in 1:nrow(overlaps)) {
      o_start <- overlaps$start[j]
      o_end <- overlaps$end[j]
      o_width <- o_end - o_start
      overlap_start <- max(f_start, o_start)
      overlap_end <- min(f_end, o_end)
      overlap_duration <- overlap_end - overlap_start
      
      if (overlap_duration > 0) {
        proportion <- overlap_duration / o_width
        contribution <- overlaps$population[j] * proportion
        total_pop <- total_pop + contribution
      }
    }
    
    return(data.frame(
      country = country_name,
      age_group = fine_df$age_group[i],
      population = (total_pop)
    ))
  })
  
  return( bind_rows(redistributed) )
}

