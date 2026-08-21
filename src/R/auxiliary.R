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
# Normalise country keys across RespiCompass data files ----
# -------------------------------------------------------- -
# The RespiCompass files are NOT internally consistent in how they key countries,
# and the model joins on the convention in supporting-files/countries.csv:
#
#   * GREECE   - population_estimates.csv and births_by_month.csv use Eurostat's
#                'EL'; countries.csv (and hence the model) uses ISO-2 'GR'.
#   * CZECHIA  - the target-data files (hospitaladmissions.csv,
#                hospitalburden_agegroups.csv) name it "Czech Republic";
#                countries.csv calls it "Czechia".
#
# Left unhandled, the filters silently match ZERO rows for those countries:
# Greece gets an empty population (model failure), Czechia an empty calibration
# target (nothing to fit to). Both helpers are idempotent, so it is safe to apply
# them to any RespiCompass table defensively.
normalise_iso2 = function(df, col = "country") {
  if (!is.null(df) && col %in% names(df)) {
    v = as.character(df[[col]])
    v[v == "EL"] = "GR"
    df[[col]] = v
  }
  return(df)
}

normalise_country_name = function(df, col = "country") {
  if (!is.null(df) && col %in% names(df)) {
    v = as.character(df[[col]])
    v[v == "Czech Republic"] = "Czechia"
    df[[col]] = v
  }
  return(df)
}

# -------------------------------------------------------- -
# Canonical fine -> reporting age-group lookup ----
# Single source of truth: age_group_map in default.yaml, so the
# fine-to-reporting mapping is never hard-coded per script.
# Returns a named character vector (names = fine groups, values = bands).
# -------------------------------------------------------- -
# Returns named character vector: fine group → reporting band
age_band_lookup = function(o) {
  unlist(yaml::read_yaml(o$pth$params_default)$age_group_map)
}

# Returns ordered character vector of all fine model age groups (from age_groups in YAML)
age_group_levels = function(o) {
  unlist(yaml::read_yaml(o$pth$params_default)$age_groups)
}

# Returns ordered character vector of unique RespiCompass reporting bands (from age_group_map in YAML)
respicompass_band_order = function(o) {
  unique(unlist(yaml::read_yaml(o$pth$params_default)$age_group_map))
}

# Returns data frame mapping fine age groups to RespiCompass reporting bands (for left_join)
age_group_map_df = function(o) {
  lkp <- age_band_lookup(o)
  data.frame(age_group            = names(lkp),
             age_group_respiCompass = unname(lkp),
             stringsAsFactors = FALSE)
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
  #
  # NB: the contact rate is divided by the COLUMN bin width (c_len), not the row
  # width. mat[i,j] is the mean number of contacts a person in ego-group i has
  # with alters in group j. Splitting a group into finer bands therefore behaves
  # differently by dimension:
  #   - COLUMNS (alters) must be SPLIT proportionally, since the alters are
  #     divided between the finer bands -> divide by c_len.
  #   - ROWS (egos) must be REPLICATED, since each person in a finer ego band has
  #     the same per-person contact rate -> no division.
  # Combined with collapse_from_1m() this gives
  #     collapsed[i,j] = mat[I,J] * len_target(j) / len_original(J),
  # which round-trips exactly (collapsed == mat) for ANY set of bin widths.
  #
  # This previously divided by r_len, which is equivalent ONLY when all original
  # bins share the same width (as in the plain 5-year Prem bands). It silently
  # broke once the last band was widened to 75-100y to cover the 80+ model group:
  # every 75+ ego row was scaled by 60/300 = 1/5, cutting their force of infection
  # fivefold and collapsing modelled burden in 75-79y and 80+y.
  expand_to_1m <- function(mat, age_bins) {
    n <- sum(age_bins$upper - age_bins$lower)
    full_mat <- matrix(0, nrow = n, ncol = n)

    row_idx <- 1
    for (i in 1:nrow(age_bins)) {
      r_len <- age_bins$upper[i] - age_bins$lower[i]
      col_idx <- 1
      for (j in 1:nrow(age_bins)) {
        c_len <- age_bins$upper[j] - age_bins$lower[j]
        full_mat[row_idx:(row_idx + r_len - 1), col_idx:(col_idx + c_len - 1)] <- mat[i, j]/c_len
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
  if (str_detect(label, "\\+y$")) {
    lower <- as.numeric(str_remove(label, "\\+y$"))
    return(c(lower, 100))
  }
  
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


# ------------------------------------------------------------------ -
# Derive age-specific background mortality RATE from RespiCompass counts ----
# ------------------------------------------------------------------ -
# RespiCompass provides annual all-cause DEATH COUNTS by (coarse) age band.
# We convert these to a per-capita annual RATE (deaths / population) and map the
# rate onto the model's fine age groups. The denominator comes from the model's
# already-redistributed fine population, summed back up to the mortality bands,
# so numerator and denominator are internally consistent.
#
# ASSUMPTION / LIMITATION: a single annual rate per band is applied uniformly
# across the year (see the flat monthly-rate note in ageing_event()); RespiCompass
# also ships mortality_month.csv (seasonal share) which is NOT used here.
#
# Returns a named numeric vector (one annual rate per fine age group), or NULL
# if the country is absent from the mortality data (caller then falls back to
# the yaml `background_mortality_rate` override).
compute_background_mortality <- function(mortality_df, population_fine,
                                         age_group_map, age_groups, iso2) {

  # Deaths for this country, keyed by RespiCompass mortality band. Guard against
  # multiple reference years by keeping only the most recent.
  mort <- mortality_df %>% filter(iso2_code == iso2)
  if (nrow(mort) == 0) return(NULL)
  if ("reference_year" %in% names(mort))
    mort <- mort %>% filter(reference_year == max(reference_year))
  deaths_band <- setNames(mort$total_deaths, mort$age_group)

  # RSVinea reporting bands -> RespiCompass mortality bands.
  # (All infant reporting bands collapse into the single '<1' mortality band.)
  # Keys are the RSVinea reporting bands (age_group_map values); values are the
  # raw RespiCompass mortality-file band labels.
  reporting_to_mortality <- c(
    "0-3m"   = "<1",    "3-6m"  = "<1",    "6-12m" = "<1",
    "1-5y"   = "1-4",   "5-18y" = "5-17",
    "18-60y" = "18-59", "60-65y" = "60-64", "65-70y" = "65-69",
    "70-75y" = "70-74", "75-80y" = "75-79", "80+y"  = "80+"
  )

  # Fine age group -> reporting band -> mortality band (length n_age, ordered)
  fine_to_reporting <- setNames(age_group_map$larger_group, age_group_map$smaller_group)
  band_of_age <- unname(reporting_to_mortality[fine_to_reporting[age_groups]])

  # Population per mortality band = sum of fine populations falling in the band
  pop_fine <- setNames(population_fine$population, population_fine$age_group)
  pop_band <- tapply(pop_fine[age_groups], band_of_age, sum)

  # Annual per-capita rate per band, then broadcast back to the fine age groups
  rate_band <- deaths_band[names(pop_band)] / pop_band
  rate_fine <- unname(rate_band[band_of_age])
  names(rate_fine) <- age_groups

  if (any(is.na(rate_fine)))
    warning("compute_background_mortality(): NA rate for age group(s): ",
            paste(age_groups[is.na(rate_fine)], collapse = ", "),
            " (country ", iso2, "). Check mortality-band coverage.")

  return(rate_fine)
}

