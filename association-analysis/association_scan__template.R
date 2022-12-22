library(tidyverse)

# Running from terminal
# "Usage: <script>.R <dataset> <endpoint type> <endpoint location>"

# Monitor time
start_time <- Sys.time()

# Install required packages -----------------------------------------------

message("Installing required packages...")

# Install newer versions of rlang
if (utils::packageVersion("rlang") <= "0.4.0") {
  install.packages("rlang", repos = "http://cran.r-project.org", quiet = TRUE)
}

# Install required packages
purrr::walk(
  .x = c("tidyverse", "argparse", "survival", "stats", "aws.s3"),
  .f = ~ if (!requireNamespace(.x, quietly = TRUE)) {
    install.packages(.x, quiet = TRUE)
  }
)

# Source functions 
source("utils.R")

# Set options for the future package --------------------------------------

future::plan("multiprocess")
options(future.plan = "multiprocess", mc.cores = 16L)
options(future.globals.maxSize = 960*1024^2)  # 730*1024^2
future::availableCores("mc.cores")

# Command line interface --------------------------------------------------
args <- commandArgs(TRUE)

if (length(args) != 4L) {
  stop(
    "Usage: <script>.R <dataset> <endpoint type> <endpoint location>",
    call. = FALSE
  )
}

dataset <- args[1L]
endpoint_type <- args[2L]
endpoint_location <- args[3L]

# Set up everything -------------------------------------------------------

# Path on s3 where the endpoints stored
dir_endpoints <- ...

# Path on s3 where the results should be saved
dir_results <- ...

# List the endpoint files in dir_endpoints
endpoint_files <- list_files_s3(bucket = ... , dir = dir_endpoints, file_format = ...)


# Read data ---------------------------------------------------------------

# Here one should read and merge the data for biomarkers, related spectrometer
# information and the variables needed for adjusting the analyses (age, sex, UKB
# assessment center)

# Read a data frame of the variables needed for adjusting the analyses
df_adj <- ...

# Read a data frame of the biomarker NMR data
df_nmr <- ...


# Preprocess the NMR data -------------------------------------------------

# Remove outliers 4*IQR from median
df_nmr <- df_nmr %>% 
  dplyr::group_by(.data$biomarker_id) %>% 
  dplyr::filter(
    abs(.data$concentration - median(.data$concentration, na.rm = TRUE)) < 4 * IQR(.data$concentration, na.rm = TRUE)
  ) %>% 
  dplyr::ungroup() 


# Adjust NMR data for the spectrometer
df_nmr <- df_nmr %>% 
  dplyr::group_by(.data$biomarker_id) %>% 
  dplyr::mutate(concentration = .data$concentration %>% log1p()) %>% 
  dplyr::ungroup() %>%
  correct_concentrations_for_spectrometer() %>% 
  dplyr::select(.data$ID, .data$biomarker_id, concentration = .data$residuals)


# Scale biomarker concentrations
df <- df %>%
  dplyr::group_by(.data$biomarker_id) %>%
  dplyr::mutate(
    concentration_scaled = .data$concentration %>% scale() 
  ) %>% 
  dplyr::ungroup() 


# Run associations --------------------------------------------------------

purrr::walk(
  .x = endpoint_files,
  .f = function(file){
    
    # Read the endpoint data
    # The endpoint file should have columns PREVALENT_X, INCIDENT_X or X_DEATH
    # and X_AGE denoting the prevalent, incident or mortality events and age at
    # the event, respectively.
    df_outcome <- 
      aws.s3::s3read_using(
        FUN = readRDS,
        object = file.path(dir_endpoints, file)
      ) %>% 
      dplyr::mutate(ID = as.character(.data$ID)) 
    
    # Extract name of the event and event age based on the endpoint type
    event <- 
      switch(
        endpoint_type,
        "PREVALENT" = df_outcome %>% dplyr::select(dplyr::starts_with("PREVAL_")) %>% names(),
        "INCIDENT" = df_outcome %>% dplyr::select(dplyr::starts_with("INCIDENT_")) %>% names(),
        "DEATH" = df_outcome %>% dplyr::select(dplyr::ends_with("_DEATH")) %>% names()
      )
    
    event_age <- 
      switch(
        endpoint_type,
        "PREVALENT" = NULL,
        "INCIDENT" = df_outcome %>% dplyr::select(dplyr::ends_with("_AGE")) %>% names(),
        "DEATH" = df_outcome %>% dplyr::select(dplyr::ends_with("_AGE")) %>% names()
      )
  
    
    # Define formula 
    # (logistic regression for prevalent disease endpoints, cox for incident disease and mortality endpoints)
    formula <- 
      if (endpoint_type == "PREVALENT") {
        formula(paste0(event, "~", "concentration_scaled + age + factor(sex) + factor(region)"))
      } else {
        formula(paste0("survival::Surv(time = BL_AGE, time2 = ", event_age,", event = ",  event, ") ~",
                       "concentration_scaled + factor(sex) + factor(region)"))
      }
    
    df <- df %>% 
      dplyr::inner_join(
        df_outcome %>% dplyr::select(.data$ID, !!event, !!event_age),
        by = "ID"
      ) %>% 
      tidyr::drop_na()
    
    message("Running associations with formula: ", c(formula))
    message("Number of samples: ",  df %>% pull(.data$ID) %>% unique() %>% length(), ".")
    
    if (endpoint_type != "PREVALENT") {
      # For incident and mortality endpoints, remove prevalent cases
      df <- df %>% 
        dplyr::rename("BL_AGE" = .data$age) %>%
        dplyr::filter(.data$BL_AGE < .data[[event_age]]) 
    }
    
    # Compute associations
    df_assoc <- df %>%
      # Factorize variables
      dplyr::mutate(
        region = factor(.data$region),
        sex = factor(.data$sex)
      ) %>% 
      tidyr::nest(data = -.data$biomarker_id) %>%
      dplyr::mutate(
        # purrr::map can be used here also
        model = furrr::future_map(
          .x = .data$data,
          # Compute associations if n_events >= 50
          .f = ~ if (sum(.x %>% dplyr::select(!!event), na.rm = T) >= 50) {
            {if (endpoint_type %in% c("INCIDENT", "DEATH")) {
              # Cox regression for incident and mortality events
              survival::coxph(
                formula = formula,
                data = .x
              ) 
            } else if (endpoint_type == "PREVALENT") {
              # Logistic regerssion for prevalent events
              stats::glm(
                formula = formula,
                data = .x,
                family = "binomial"
              ) 
            }} %>% 
              broom::tidy() %>%
              dplyr::filter(term == "concentration_scaled") %>%
              dplyr::select(
                estimate,
                se = std.error,
                pvalue = p.value
              ) %>% 
              dplyr::mutate(
                n = nrow(.x),
                n_events = sum(.x %>% dplyr::select(!!event))
              ) 
          } else {
            tibble::tibble(
              estimate = NA_real_,
              se = NA_real_,
              pvalue = NA_real_,
              n = nrow(.x),
              n_events = sum(.x %>% dplyr::select(!!event), na.rm = T)
            )
          },
          .progress = TRUE
        )
      ) %>%
      dplyr::select(-.data$data) %>%
      tidyr::unnest(.data$model) %>%
      tibble::add_column(endpoint = event) %>% 
      dplyr::select(
        .data$biomarker_id,
        .data$se,
        .data$pvalue,
        .data$n,
        .data$n_events,
        .data$endpoint
      )
    
    # Save results on s3 ds-ukb-results bucket
    message(paste0("Saving results in ", dir_results, "...\n"))
    
    aws.s3::s3write_using(
      df_assoc, 
      FUN = write_csv, 
      object = paste0(dir_results, "/", event, ".csv")
    )
  } 
)


# Print time --------------------------------------------------------------

message(paste0("Done! Time taken ", round(difftime(Sys.time(), start_time, units = "mins"), 2), " mins."))



