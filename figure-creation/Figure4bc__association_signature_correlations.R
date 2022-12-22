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

# Selected ICD10 code pairs to plot
selected_pairs <-
  tibble::tribble(
    ~icd_code1, ~icd_code2,
    "G62", "K76",
    "I21", "I50",
  )

# Read summary statistics -------------------------------------------------

df_summary_stats <- read_summary_stats() %>% 
  dplyr::filter(age_group == "Full population") 


# Scatterplots of the association signature correlations ------------------

plot_list<-
  purrr::map2(
    .x = selected_pairs$icd_code1,
    .y = selected_pairs$icd_code2,
    .f = function(icd_code1, icd_code2) {

      df_selected_pair <- df_summary_stats %>%
        dplyr::filter(icd10 %in% c(icd_code1, icd_code2)) 
      
      # Write source data
      readr::write_tsv(
        df_selected_pair %>% 
          dplyr::select(biomarker_name, group_name, disease, chapter, estimate, se, pvalue),
        paste0("source-data/Figure4", letters[1 + which(selected_pairs$icd_code1 == icd_code1)], "_source_data.tsv")
      )
      
      # Format for plotting
      df_subplot <- df_selected_pair %>%
        dplyr::select(
          biomarker_name, icd10, hazard_ratio, ci_upper, ci_lower, significance
        )  %>%
        tidyr::pivot_wider(
          names_from = "icd10",
          values_from = c("hazard_ratio", "ci_upper", "ci_lower", "significance")
        ) %>% 
        dplyr::mutate(
          significance =
            dplyr::case_when(
              .data[[!!paste0("significance_", icd_code1)]] & .data[[!!paste0("significance_", icd_code2)]] ~ "Both *p* < 5e-5",
              .data[[!!paste0("significance_", icd_code1)]]  | .data[[!!paste0("significance_", icd_code2)]] ~ "One *p* < 5e-5",
              !.data[[!!paste0("significance_", icd_code1)]] & !.data[[!!paste0("significance_", icd_code2)]] ~ "Both *p* \u2265 5e-5",
              TRUE ~ NA_character_
            )
        )
      
      
      # Define colors
      colors <- c("black", "grey40", "gray80") %>%
        set_names(c("Both *p* < 5e-5", "One *p* < 5e-5"), "Both *p* \u2265 5e-5")
      

      # Extract limits
      min_limit <- min(
        c(df_subplot %>% pull(paste0("ci_lower_", icd_code1)),
          df_subplot %>% pull(paste0("ci_lower_", icd_code2))),
        na.rm = T)

      max_limit <- max(
        c(df_subplot %>% pull(paste0("ci_upper_", icd_code1)),
          df_subplot %>% pull(paste0("ci_upper_", icd_code2))),
        na.rm = T)
      
      ggplot(
        data = df_subplot,
        aes_string(
          x = paste0("hazard_ratio_", icd_code1), 
          y = paste0("hazard_ratio_", icd_code2), 
          color = "significance")
      ) +
        geom_hline(yintercept = 1, color = "#E64B35FF") +
        geom_vline(xintercept = 1, color = "#E64B35FF") +
        geom_abline(color = "darkgrey", linetype = "dashed") +
        geom_pointrange(
          aes_string(ymin = paste0("ci_lower_", icd_code2), ymax = paste0("ci_upper_", icd_code2)),
          size = 0.1
        ) + 
        geom_errorbarh(
          aes_string(xmin = paste0("ci_lower_", icd_code1), xmax = paste0("ci_upper_", icd_code1)),
          size = 0.1, height = 0
        ) + 
        geom_point(pch = 20, size = 3) +
        labs(
          x = paste0(ukbtools::ukb_icd_code_meaning(as.character(icd_code1))$meaning, ", hazard ratio (95% CI)"),
          y = paste0(ukbtools::ukb_icd_code_meaning(as.character(icd_code2))$meaning, ", hazard ratio (95% CI)"),
          color = "Significance"
        ) +
        scale_color_manual(values = colors) +
        ggrepel::geom_text_repel(aes(label = biomarker_name), show.legend = FALSE) + 
        scale_x_continuous(limits = c(min_limit, max_limit)) +
        scale_y_continuous(limits = c(min_limit, max_limit)) +
        ggpubr::stat_cor(
          aes(label = paste("R:", ..r..)), 
          label.y.npc = 0.99, 
          show.legend = FALSE,
          color = "black"
        ) +
        theme_bw() + 
        theme(
          legend.position = c(0.2, 0.78),
          legend.key.size = unit(0.5, "cm"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(size = 0.2, colour = "black"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          legend.text = ggtext::element_markdown(size = 10),
          legend.background = element_blank()
        )
      
    }
  )


# Save plot ---------------------------------------------------------------

grDevices::cairo_pdf(
  file.path(dir_figures, "Figure4bc.pdf"), width = 12, height = 5
)

cowplot::plot_grid(
  plotlist = plot_list %>%
    lapply(egg::set_panel_size, width = unit(3.5, "in"), height = unit(3.5, "in")),
  ncol = 2
)

dev.off()




