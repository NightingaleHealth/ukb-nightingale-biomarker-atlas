library(tidyverse)

# Setup -------------------------------------------------------------------

# Install pacman
if (suppressWarnings(!require("pacman"))){
  install.packages("pacman")
}

# Load and install required packages if not already installed
pacman::p_load(
  tidyverse, RColorBrewer, pheatmap, ggpubr,
  ggrepel, cowplot, egg, grDevices, devtools
)

# Source utility functions
purrr::walk(
  .x = list.files("utils", full.names = T),
  .f = ~source(.x)
)

dir_figures <- "figures"

# List incidence endpoints for plotting
selected_incidence_endpoints <- 
  tibble::tribble(
    ~icd10, ~disease_name, 
    "A41", "Sepsis",
    "C34", "Lung cancer",
    "F32", "Depression",
    "I21", "Myocardial infarction",
    "G47", "Sleep disorders",
    "M81", "Osteoporosis"
  ) %>% 
  dplyr::arrange(.data$icd10)

# List biomarkers shared by clinical chemistry and NMR
shared_biomarkers <- 
  c(
    "Total-C", 
    "Clinical LDL-C",
    "HDL-C",
    "Total triglycerides",
    "ApoB", 
    "ApoA1",
    "Albumin", 
    "Glucose"
  )

# Read summary statistics -------------------------------------------------

# Nightingale NMR biomarker results
df_assoc_nmr <- read_summary_stats() %>% 
  dplyr::filter(age_group == "Full population") %>% 
  dplyr::mutate(method = "Nightingale NMR")

# Clinical chemistry replication results
df_assoc_cc <- 
  readr::read_csv(
    "data/ukb_clinical_chemistry_summary_statistics.csv", 
    show_col_types = FALSE
  ) %>% 
  dplyr::mutate(method = "Clinical chemistry")

# Combine data ------------------------------------------------------------

df_plot <-  
  dplyr::bind_rows(
    df_assoc_nmr, df_assoc_cc
  ) %>% 
  dplyr::filter(
    icd10 %in% selected_incidence_endpoints$icd10,
    biomarker_name %in% shared_biomarkers
  ) %>% 
  dplyr::left_join(selected_incidence_endpoints, by = "icd10") %>% 
  dplyr::mutate(
    biomarker_name = factor(biomarker_name, levels = shared_biomarkers),
    disease_name = paste(.data$icd10, disease_name)
  ) %>% 
  dplyr::arrange(biomarker_name) %>% 
  dplyr::select(
    biomarker_name, disease_name, method, 
    estimate, se, pvalue, 
    n, n_events
  )

# Write source data -------------------------------------------------------

readr::write_tsv(
  df_plot %>% 
    dplyr::select(
      biomarker_name, disease = disease_name, method,
      estimate, se, pvalue
    ),
  "source-data/Figure5_source_data.tsv"
)

# Plotting ----------------------------------------------------------------

pdf(
  paste0(dir_figures, "/Figure5.pdf"),
  width = 18*0.393701,
  height = 10*0.393701
)

ggforestplot::forestplot(
  df = df_plot,
  name = biomarker_name,
  estimate = estimate,
  se = se, 
  colour = method,
  pvalue = pvalue,
  logodds = TRUE,
  psignif = 0.05/1000
) + 
  facet_wrap(~disease_name, ncol = 3, scales = "free_x") +
  labs(x = "Hazard ratio (95% CI)") +
  guides(colour = guide_legend(ncol = 1, reverse = T)) +
  theme(
    text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    axis.title.x = element_text(size = 7)
  ) 

dev.off()

