# Utility function for reading and preprocessing summary statistics for plotting

read_summary_stats <- function(significance_thresh = 0.05/1000) {
  
  # Read list of clinically validated biomarkers ----------------------------
  
  validated_bmrs <- 
    readr::read_csv(
      "data/validated_biomarkers.csv",
      col_types = readr::cols()
    ) %>% 
    dplyr::pull(biomarker_id)
  
  # Include chapters A-N for the main analyses
  included_chapters <- LETTERS[1:14] %>% paste(collapse = "|")
  
  # Read summary statistics -------------------------------------------------
  
  df <- 
    readr::read_csv(
      "data/ukb_nightingale_biomarker_summary_statistics.csv",
      col_types = readr::cols()
    ) %>% 
    dplyr::filter(
      # Include incident events in the main analyses
      endpoint_type == "incident",
      biomarker_id %in% validated_bmrs,
      str_detect(.data$icd10, included_chapters)
    ) %>% 
    dplyr::select(
      age_group,
      biomarker_id, biomarker_name, group_name,
      estimate, se, pvalue,
      icd10, icd10_chapter, icd10_desc,
      n, n_events
    ) 
  

  # Add confidence intervals, hazard ratios and statistical significance 
  
  df <- df %>% 
    dplyr::mutate(
      hazard_ratio = exp(estimate),
      ci_upper = exp(estimate + 1.96 * se),
      ci_lower = exp(estimate - 1.96 * se),
      significance = ifelse(.data$pvalue < significance_thresh, TRUE, FALSE)
    )
  

  # Wrangle names to fit in figures -----------------------------------------

  df <- df %>% 
    dplyr::mutate(
      icd10_desc = icd10_desc %>% 
        str_replace_all(
          c(
            "classified elsewhere" = "cl. elsw.",
            "not elsewhere classified" = "NEC",
            "Other disorders of fluid, electrolyte and acid-base balance" = "Disorders of fluid, electrolyte and acid-base balance",
            "and other lipidemias" = "",
            "Osteoporosis without" = "Osteoporosis w/o"
            )
        ),
      icd10_chapter = .data$icd10_chapter %>% 
        dplyr::recode(
          "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism" =
            "Diseases of the blood and blood-forming organs"
        )
    ) %>% 
    dplyr::rename(disease = icd10_desc, chapter = icd10_chapter)
  
  return(df)
  
}