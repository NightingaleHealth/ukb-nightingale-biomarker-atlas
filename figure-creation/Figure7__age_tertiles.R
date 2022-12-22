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

# Biomarker groups to plot for the main figure
selected_bmr_groups <- 
  c("Cholesterol", 
    "Apolipoproteins",
    "Triglycerides",
    "Fatty acid ratios",  
    "Inflammation")

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


# Read summary stats ------------------------------------------------------

df_summary_stats <- read_summary_stats() %>% 
  # Keep only age tertile results (i.e. exclude full population results)
  dplyr::filter(.data$age_group != "Full population")

# Extract data frame for plotting -----------------------------------------

df_plot <- df_summary_stats %>% 
  dplyr::filter(
    icd10 %in% selected_incidence_endpoints$icd10,
    group_name %in% selected_bmr_groups
  ) %>% 
  dplyr::left_join(selected_incidence_endpoints, by = "icd10") %>% 
  dplyr::mutate(
    disease_name = paste(.data$icd10, disease_name),
    disease_name = factor(disease_name, levels = unique(.data$disease_name)),
    age_group = factor(.data$age_group, levels = rev(unique(.data$age_group))),
    group_name = factor(.data$group_name, levels = selected_bmr_groups),
  ) %>% 
 dplyr::arrange(disease_name, group_name)


# Write source data -------------------------------------------------------

readr::write_tsv(
  df_plot %>%
    dplyr::select(
      age_group, disease = disease_name,
      biomarker_name, group_name, 
      estimate, se, pvalue
    ),
  paste0("source-data/Figure7_source_data.tsv")
)

# Plotting ----------------------------------------------------------------

plot_list <- 
  purrr::map2(
    .x = unique(selected_incidence_endpoints$icd10),
    .y = selected_incidence_endpoints$disease_name,
    .f = function(icd10, disease_name) {
 
      df_subplot <- df_plot %>% 
        dplyr::filter(icd10 == !!icd10) 
    
      plot <- 
        ggforestplot::forestplot(
          df_subplot,
          name = biomarker_name,
          estimate = estimate,
          pvalue = pvalue,
          logodds = TRUE,
          psignif = signif_thresh,
          colour = age_group
        ) + 
        ggforce::facet_col(
          facets = ~group_name, 
          scales = "free_y", 
          space = "free", 
          strip.position = "top"
        ) +
        masterdata::scale_color_ng_d(reverse = T) +
        scale_x_continuous(
          trans = 'log10',
          limits = c(min(df_plot$ci_lower, na.rm = T), max(df_plot$ci_upper, na.rm = T)),
          breaks = c(0.6,1, 1.4, 1.8, 2.2)
        ) +
        ggtitle(paste0(icd10, " ", disease_name)) +
        theme(
          text = element_text(size = 7),
          axis.text.y = {if (icd10 %in% selected_incidence_endpoints$icd10[c(1,4)]) {element_text(color = "black", size = 6.5)} else { element_blank()}},
          axis.title.x = element_blank(),
          legend.position = ifelse(icd10 == selected_incidence_endpoints$icd10[[1]], "right", "none"),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(size = 7),
          strip.text = element_text(size = 6.5),
          panel.spacing = unit(0.1,'cm')
        ) 
      
      return(plot)
    }
  )


# Extract legend
legend <- 
  cowplot::plot_grid(
    plotlist = list(cowplot::get_legend(plot_list[[1]])),
    ncol = 2
  ) 

# Remove legend from original plot
plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "none")


# Combine left and right hand side
plot_main <-
  gridExtra::grid.arrange(
    grobs = plot_list,
    widths = c(1.4,1,1),
    ncol = 3,
    bottom = grid::textGrob(label = "Hazard ratio (95% CI), per 1-SD increment", x = 0.58, gp = grid::gpar(fontsize = 7))
  )


# Export
pdf(paste0(dir_figures, "/Figure7.pdf"),  width = 18*0.393701, height = 21*0.393701)
cowplot::plot_grid(
  plotlist = list(legend, plot_main),
  ncol = 1,
  rel_heights = c(1, 15)
)
dev.off()

