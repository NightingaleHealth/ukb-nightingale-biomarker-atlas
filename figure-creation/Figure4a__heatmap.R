library(tidyverse)

# Setup -------------------------------------------------------------------

# Install pacman
if (suppressWarnings(!require("pacman"))){
  install.packages("pacman")
}

# Load and install required packages if not already installed
pacman::p_load(
  tidyverse, RColorBrewer, pheatmap, ggpubr,
  ggrepel, cowplot, egg, grDevices
)

# Source utility functions
purrr::walk(
  .x = list.files("utils", full.names = T),
  .f = ~source(.x)
)

dir_figures <- "figures"

# Read summary statistics -------------------------------------------------

df_summary_stats <- read_summary_stats() %>% 
  dplyr::filter(age_group == "Full population")


# Extract data frame for plotting -----------------------------------------

# Get 3 most common diseases by each chapter
df_most_common <- df_summary_stats %>% 
  dplyr::filter(significance == TRUE) %>% 
  dplyr::group_by(chapter, icd10) %>% 
  dplyr::summarise(n_significant_associations = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(n_significant_associations)) %>% 
  dplyr::group_by(chapter) %>% 
  dplyr:: mutate(rank = row_number(-n_significant_associations)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(rank <= 3)

df_plot <- df_summary_stats %>% 
  dplyr::filter(icd10 %in% df_most_common$icd10) %>% 
  dplyr::mutate(
    disease = purrr::map_chr(
      .x = .data$icd10,
      .f = ~ukbtools::ukb_icd_code_meaning(as.character(.x))$meaning
    )
  )


# Convert to matrix format ------------------------------------------------

# Convert to matrix (for heatmap)
plot_mat <- df_plot %>% 
  dplyr::select(biomarker_name, disease, estimate) %>% 
  tidyr::spread(key = biomarker_name, value = estimate)   %>% 
  tibble::column_to_rownames("disease") %>% 
  as.matrix()

# Matrix for marking statistical significance
mat_signif_pval <- df_plot %>% 
  dplyr::mutate(significant_pval = dplyr::if_else(significance == TRUE, "*", "")) %>% 
  select(biomarker_name, disease, significant_pval) %>% 
  tidyr::spread(key = biomarker_name, value = significant_pval) %>% 
  tibble::column_to_rownames("disease") %>% 
  as.matrix() 

# Checks
identical(rownames(plot_mat), rownames(mat_signif_pval))
identical(colnames(plot_mat), colnames(mat_signif_pval))



# Plotting specs ----------------------------------------------------------

# Make a custom palette and breaks
custom_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100)

custom_breaks  <-
  c(
    seq(-0.7, 0, length.out = length(custom_palette)/2),
    seq(0.01, 0.7, length.out = length(custom_palette)/2)
  )

# Extract annotations for bioamarkers (groups) and diseases (ICD-10 chapters)
biomarker_annotations <- df_plot %>% 
  dplyr::select(biomarker_name, "Biomarker group" = group_name) %>% 
  dplyr::distinct() %>% 
  tibble::column_to_rownames("biomarker_name")

chapter_annotations <- df_plot %>% 
  dplyr::select(disease, "Chapter" = chapter) %>% 
  dplyr::distinct() %>% 
  tibble::column_to_rownames("disease")


# Specify colors for the annotations
annotation_colors <-
  list(
    `Biomarker group` = c(
      "Amino acids" = "#DEB4FCFF",
      "Aromatic amino acids" = "#FDA2D7FF",
      "Branched-chain amino acids" = "#FEBACDFF",
      "Apolipoproteins" = "#B2FFEFFF",
      "Cholesterol" = "#6D91E1FF",
      "Fatty acids" = "#ABC2F2FF",
      "Fatty acid ratios" = "#C9CFFEFF",
      "Fluid balance" = "#FDE988FF",
      "Glycolysis related metabolites" = "#8EFBAEFF",
      "Inflammation" = "#FC7E85FF",
      "Triglycerides" = "#B3DAF6FF"
    ),
    `Chapter` = default_color_palette() %>% 
      set_names(unique(df_plot %>% arrange(icd10) %>% pull(chapter)))
  )


# Save source data --------------------------------------------------------

readr::write_tsv(
  plot_mat %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("disease"),
  "source-data/Figure4a_source_data.tsv"
)

# Plotting ----------------------------------------------------------------

plot <- 
  pheatmap::pheatmap(
    mat = plot_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    gaps_row = 5,
    clustering_distance_cols = "correlation",
    clustering_distance_rows = "correlation",
    treeheight_row = 40,
    annotation_col = biomarker_annotations,
    annotation_row = chapter_annotations,
    annotation_colors = annotation_colors,
    color = custom_palette,
    breaks = custom_breaks,
    na_col = "#F3F3F3",
    border_color = NA, 
    cellwidth = 11,
    cellheight = 11,
    angle_col = 45,
    display_numbers = mat_signif_pval
  )

ggsave(
  filename =  file.path(dir_figures, "Figure4a.pdf"),
  plot = plot,
  height = 12,
  width =  22
)
