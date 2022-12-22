# Utility functions needed by the association scan


list_files_s3 <- function(bucket, dir, file_format) {
  
  files <- aws.s3::get_bucket_df(
    bucket = bucket,
    prefix = paste0(
      stringr::str_remove(
        string = dir, 
        pattern = paste0("s3://", bucket, "/")), 
      "/"),
    max = Inf # This is super important - otherwise will return only 1000 keys
  ) %>%
    dplyr::mutate(
      file_name = stringr::str_remove_all(
        string = .data$Key,
        pattern = paste0(
          stringr::str_remove(
            string = dir, 
            pattern = paste0("s3://", bucket, "/")), 
          "/")
      )
    ) %>%
    dplyr::filter(
      stringr::str_detect(
        string = .data$file_name,
        pattern = file_format)
    ) %>%
    dplyr::pull(.data$file_name)
  
  return(files)
}


correct_concentrations_spectrometer <- function(df) {
  
  # Fit linear regression to model the effect of the technical variables on biomarker concentrations 
  # and compute residuals to be used with further analysis
  
  df <- df %>%
    tidyr::drop_na() %>% 
    dplyr::mutate(
      spectrometer_serial_number = as.factor(.data$spectrometer_serial_number)
    ) %>% 
    tidyr::nest(data = -.data$biomarker_id) %>%
    dplyr::mutate(
      # purrr::map can be used here also
      residuals = furrr::future_map(
        .x = .data$data,
        .f =  ~ stats::lm(
          formula = formula("concentration ~ factor(spectrometer_serial_number)"),
          data = .x
        ) %>% 
          resid()
      )
    ) %>% 
    tidyr::unnest(c(.data$data, .data$residuals)) 
  
  return(df)
  
}
