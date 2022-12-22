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

# Threshold for statistical significance
signif_thresh <- 0.05/1000


# Read summary statistics -------------------------------------------------

df_summary_stats <- read_summary_stats() %>% 
  dplyr::filter(age_group == "Full population")


# Extract data frame for plotting -----------------------------------------

df_plot <- df_summary_stats %>% 
  dplyr::filter(icd10 %in% selected_incidence_endpoints$icd10) %>% 
  dplyr::left_join(selected_incidence_endpoints, by = "icd10") %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    disease_name = paste(.data$icd10, disease_name),
    disease_name = factor(disease_name, levels = unique(.data$disease_name))
  ) %>% 
  dplyr::arrange(disease_name)


# Plotting ----------------------------------------------------------------

groups <- list(
  "lhs" = c("Amino acids", "Aromatic amino acids",  "Branched-chain amino acids",
            "Fluid balance","Glycolysis related metabolites", "Cholesterol"),
  "rhs" = c("Inflammation", "Apolipoproteins", "Fatty acid ratios", 
            "Fatty acids", "Triglycerides")
)


# Write source data -------------------------------------------------------

readr::write_tsv(
  df_plot %>% 
    dplyr::filter(group_name %in% c(groups$lhs, groups$rhs)) %>% 
    dplyr::select(biomarker_name, group_name, disease = disease_name, estimate, se, pvalue),
  "source-data/Figure3_source_data.tsv"
)


# Plot --------------------------------------------------------------------

plot_list <- 
  purrr::map2(
    .x = groups, .y = c(1, 2),
    .f = ~df_plot %>%  
      filter(group_name %in% .x) %>% 
      arrange(abs(estimate)) %>% 
      dplyr::mutate(group_name = factor(group_name, levels = .x)) %>% 
      ggforestplot::forestplot(
        df = .,
        name = biomarker_name,
        estimate = estimate,
        pvalue = pvalue,
        logodds = TRUE,
        psignif = signif_thresh,
        colour = disease_name
      ) + 
      ggforce::facet_col(
        facets = ~group_name,
        scales = "free",
        space = "free"
      ) +
      scale_color_manual(values = default_color_palette()) +
      labs(x = "Hazard ratio (95% CI), per 1-SD increment") +
      scale_x_continuous(
        trans = 'log10',
        limits = c(min(df_plot$ci_lower, na.rm = T), max(df_plot$ci_upper, na.rm = T)),
        breaks = c(0.6, 1, 1.4, 1.8, 2.2)
      ) +
      theme(
        text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        legend.position = ifelse(.y == 2, "right", "none"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 7)
      ) 
  )


# Extract legend
legend <- 
  cowplot::plot_grid(
    plotlist = list(NULL, cowplot::get_legend(plot_list[[2]])),
    ncol = 10
  ) 

# Remove legend from original plot
plot_list[[2]] <- plot_list[[2]] + theme(legend.position = "none")


# Combine left and right hand side
plot_main <- 
  cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 2,
    rel_widths = c(1, 1.1)
  )


# Export
pdf(paste0(dir_figures, "/Figure3.pdf"), width = 18*0.393701, height = 19*0.393701)
cowplot::plot_grid(
  plotlist = list(legend, plot_main),
  ncol = 1,
  rel_heights = c(1, 6)
)
dev.off()
