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

# Threshold for statistical significance
signif_thresh <- 0.05/1000

# List overlapping diseases in UKB and THL biobanks
overlapping_diseases <- 
  c(
    "All-cause mortality",
    "Major adverse cardiovascular event", 
    "Diabetes", 
    "COPD", 
    "Chronic kidney failure", 
    "Liver diseases"
  )

# Read associations -------------------------------------------------------

df_assoc <- readr::read_csv("data/ukb_thl_replication_summary_statistics.csv")

# Format ------------------------------------------------------------------

df_plot <- df_assoc %>% 
  dplyr::rename(disease = disease_name) %>% 
  dplyr::mutate(cohort_no = as.integer(factor(cohort))) %>% 
  dplyr::group_by(cohort) %>% 
  dplyr::mutate(color = c("#FB6467FF", "#4DBBD5FF", "#3C5488FF")[cohort_no]) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    hazard_ratio = exp(estimate),
    ci_upper = exp(estimate + 1.96 * se),
    ci_lower = exp(estimate - 1.96 * se),
    significance = ifelse(.data$pvalue > signif_thresh, TRUE, FALSE),
    cohort_significance = paste0(cohort, significance),
    fill = ifelse(significance, "white", .data$color),
    biomarker_name = factor(biomarker_name, levels = unique(masterdata::master_data$biomarkers$short_name))
  ) %>% 
  dplyr::arrange(biomarker_name) %>% 
  dplyr::arrange(group_name) %>% 
  dplyr::mutate(
    biomarker_name = factor(biomarker_name, levels = unique(.data$biomarker_name))
  )


# Define colors for potting
colors <- as.character(df_plot$color)
names(colors) <- as.character(df_plot$cohort)

fills <- as.character(df_plot$fill)
names(fills) <- as.character(df_plot$cohort_significance)


# Define biomarker annotation colors
biomarker_annotation_colors <- 
  c(
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
  )


# Plotting ----------------------------------------------------------------

plots <- purrr::map(
  # For supplement, plot remaining diseases
  .x = overlapping_diseases,
  .f = function(selected_disease) {
    
    df_subplot <-  df_plot %>% 
      dplyr::filter(.data$disease == !!selected_disease) 
    
    # Write source data
    readr::write_tsv(
      df_subplot %>% 
        dplyr::select(
          disease, cohort,
          biomarker_name, group_name,
          estimate, se, pvalue
        ),
      paste0("source-data/Figure6", letters[which(overlapping_diseases == selected_disease)], "_source_data.tsv")
    )
    
    plot <- df_subplot %>% 
      ggplot(
        aes(
          x = biomarker_name, 
          y = hazard_ratio, 
          color = cohort)
      ) +
      geom_pointrange(
        aes(
          ymin = ci_lower, 
          ymax = ci_upper, 
          fill = cohort_significance
        ),
        shape = 21,
        size = 0.4,
        fatten = 2,
        stroke = 0.5
      ) +
      geom_hline(aes(yintercept = 1)) +
      scale_y_continuous(
        trans = 'log10',
        breaks = seq(0.6, 1.8, 0.2)
      ) +
      scale_fill_manual(values = fills) +
      scale_colour_manual(values = colors) +
      guides(fill = "none") +
      labs(y = "Hazard ratio (95% CI)") +
      ggtitle(selected_disease) +
      theme_bw() +
      theme(
        text = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 7),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1, color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.2, colour = "black")
      )
  }
)


plot_legend <- df_plot %>% 
  ggplot(aes(x = biomarker_name, y = hazard_ratio, color = cohort)) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper)) +
  scale_colour_manual(values = colors) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 0.1))) +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm")
  )

legend <-  
  cowplot::plot_grid(
    plotlist = list(cowplot::get_legend(plot_legend)), 
    ncol = 2
  ) 

plots_full <- 
  cowplot::plot_grid(
    plotlist = plots,
    ncol = 2,
    labels = c("a", "b", "c", "d", "e", "f"),
    label_y = 1.02,
    label_size = 7
  )


# Export ------------------------------------------------------------------

pdf(file.path(dir_figures, "Figure6.pdf"), width = 18*0.393701, height = 18.5*0.393701)

cowplot::plot_grid(
  plotlist = list(legend, plots_full),
  ncol = 1,
  rel_heights = c(1.5, 10)
)

dev.off()



