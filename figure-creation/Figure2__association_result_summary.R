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

# Directory to store figures 
dir_figures <- "figures"
dir.create(dir_figures)

# List biomarkers included in the main plot
selected_biomarkers <- 
  c(
    "Glycoprotein acetyls", "PUFA/MUFA",
    "Alanine", "Total BCAA"
  )

# Read summary statistics -------------------------------------------------

df_summary_stats <- read_summary_stats()

df_summary_stats <- df_summary_stats %>% 
  dplyr::filter(age_group == "Full population")


# Set up default colors ---------------------------------------------------

chapter_annotations <- df_summary_stats %>% 
  dplyr::arrange(disease) %>% 
  dplyr::pull(chapter) %>% 
  unique() 

chapter_annotation_colors <- default_color_palette()[1:length(chapter_annotations)] %>% 
  purrr::set_names(chapter_annotations)


# Panel a: Summary barplot to be displayed on top -------------------------

df_barplot <- df_summary_stats %>%
  dplyr::filter(significance == TRUE) %>% 
  dplyr::group_by(biomarker_name, group_name, chapter) %>%
  dplyr::summarise(signif_n_chapter = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(biomarker_name) %>% 
  dplyr::mutate(total_signif_n = sum(signif_n_chapter)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(total_signif_n)) %>% 
  dplyr::mutate(
    biomarker_name = factor(.data$biomarker_name, levels = unique(.data$biomarker_name)),
    chapter = factor(chapter, levels = chapter_annotations)
  )


# Write source data -------------------------------------------------------

readr::write_tsv(df_barplot, "source-data/Figure2a_source_data.tsv")

# Barplot -----------------------------------------------------------------

barplot <- 
  ggplot(
    data = df_barplot,
    aes(x = biomarker_name, y = signif_n_chapter, fill = chapter)
  ) +
  geom_col() +
  labs(
    x = "", y = "Number of significant\nassociations",
    fill = "ICD-10 chapter"
  ) +
  scale_fill_manual(values = chapter_annotation_colors) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme_bw(base_size = 7) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.2, "cm"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 6.5, face = "bold"),
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1, color = "black"),
    axis.title.y = element_text(size = 6.5),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.1, colour = "black"),
    plot.margin = margin(b=-0.5, l = 0.7, unit = "cm")
  )


# Panels b-e: examples ----------------------------------------------------

df_plot <- df_summary_stats %>% 
  dplyr::filter(biomarker_name %in% selected_biomarkers) %>% 
  dplyr::arrange(pvalue) %>% 
  dplyr::group_by(biomarker_name) %>% 
  dplyr::mutate(rank = row_number(pvalue)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    # Replace 0 p-values with the smallest non-zero normalized floating-point
    pvalue = ifelse(.data$pvalue == 0, .Machine$double.xmin, .data$pvalue),
    pvalue_label = formatC(.data$pvalue, format = "e", digits = 1)
  )

# Pad disease names such that they have the same length in all displayed plots
pad_width <- df_plot %>% 
  dplyr::filter(rank <= 20) %>% 
  pull(disease) %>% 
  nchar() %>% 
  max()

df_plot <- df_plot %>% 
  dplyr::mutate(
    disease = 
      stringr::str_pad(
        .data$disease, 
        width = pad_width, 
        side = "left"
      )
  )

# Create list of the plots
ls_plots <- 
  purrr::map(
    .x = selected_biomarkers,
    .f = function(biomarker) {
      
      
      # Filter data to be plotted
      df_biomarker_plot <- df_plot %>% 
        dplyr::filter(biomarker_name == !!biomarker) %>% 
        dplyr::arrange(icd10) %>% 
        dplyr::mutate(disease = factor(disease, levels = unique(.data$disease)))
 
      df_plot_assoc <- df_biomarker_plot %>% 
        dplyr::filter(rank <= 20) %>%
        dplyr::arrange(abs(estimate)) %>% 
        dplyr::mutate(disease = factor(disease, levels = unique(.data$disease)))
      
      # Write as a source data file
      readr::write_tsv(
        df_plot_assoc %>% 
          dplyr::arrange(desc(abs(estimate))) %>% 
          dplyr::select(biomarker_name, group_name, disease, chapter, estimate, se, pvalue),
        paste0("source-data/Figure2", letters[1 + which(selected_biomarkers == biomarker)], "_source_data.tsv")
      )
      
      custom_breaks <- 
        switch(
          biomarker,
          "Alanine" = seq(0.6, 2.2, 0.4),
          "PUFA/MUFA" = seq(0.6, 1.8, 0.2),
          "Glycoprotein acetyls" = seq(0.6, 1.8, 0.2),
          "Total BCAA" = seq(0.8, 2, 0.4)
        )
      
      plot_assoc <- df_plot_assoc %>% 
        dplyr::mutate(disease = fct_rev(disease)) %>% 
        dplyr::arrange(disease) %>% 
        ggforestplot::forestplot(
          name = disease, 
          logodds = TRUE,
          size = 0.1
        ) +
        geom_text(
          aes(label = pvalue_label),
          nudge_y = 0.2,
          parse = TRUE,
          size = 1.5,
          color = "grey40"
        ) +
        xlab("Hazard ratio (95% CI)") +
        scale_x_continuous(
          trans = 'log10',
          breaks = custom_breaks,
          expand = c(0.05, 0.05)
        ) +
        theme(
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 7),
          axis.text.y = element_text(size = 6.5)
        )
      
      plot_assoc_annotations <-  df_plot_assoc %>% 
        ggplot(aes(y = disease, x = 1, fill = chapter)) +
        geom_tile() +
        theme_void() +
        scale_y_discrete(position = "right") + 
        scale_fill_manual(values = chapter_annotation_colors) +
        theme(legend.position = 'none')
      
      plot <- 
        ggpubr::ggarrange(
          plot_assoc %>% egg::set_panel_size(width = unit(1, "in"), height = unit(2.1, "in")),
          plot_assoc_annotations %>% egg::set_panel_size(width = unit(2/20, "in"), height = unit(2.1, "in")), 
          nrow = 1, 
          widths = c(2+2, 3/20) ,
          align = 'h'
        )
      
      plot <- 
        gridExtra::grid.arrange(
          plot,
          top = grid::textGrob(biomarker, hjust = 0, x = 0.1, y = 0.5, gp=grid::gpar(fontsize = 7, fontface = "bold"))
        )
      
      return(plot)
    }
  )


# Combine all plots into a single figure ----------------------------------

plot_top <- 
  cowplot::plot_grid(
    plotlist =  list(barplot %>% egg::set_panel_size(width = unit(3.7, "in"), height = unit(0.9, "in"))),
    ncol = 1,
    labels = c("a"),
    label_size = 7,
    label_y = 0.99
  )

plot_bottom <- 
  cowplot::plot_grid(
    plotlist = ls_plots,
    ncol = 2,
    labels = c("b", "c", "d", "e"),
    vjust = 0,
    label_y = 0.99,
    label_size = 7
  )

grDevices::cairo_pdf(
  filename = paste0(dir_figures, "/Figure2.pdf"), 
  width = 18*0.393701, 
  height = 19*0.393701
) 

gridExtra::grid.arrange(
  grobs = list(plot_top, plot_bottom),
  ncol = 1,
  heights = c(1, 2.7)
)

dev.off()


