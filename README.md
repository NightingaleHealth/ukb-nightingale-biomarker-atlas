## Data analysis scripts for "Atlas of plasma NMR biomarkers for health and disease in 118,461 individuals from the UK Biobank"

### Description

This repository contains scripts related to paper "Atlas of plasma NMR biomarkers for health and disease in 118,461 individuals from the UK Biobank". We provide an example template for running the association analyses of NMR biomarkers against disease endpoints and scripts for the downstream analyses based on the summary statistics from the association analyses. 

A webtool for visualizing the results is available at https://nightingalehealth.com/atlas. Summary statistics from the
association analyses are also made publicly available for download via the webtool. 

### Requirements

The code has been tested with R version 4.1.1. 

The scripts will take care of installing the required packages, which include
the following:

- tidyverse
- rlang
- argparse
- future
- aws.s3 (if using s3-compatible data storage)
- survival
- stats
- RColorBrewer
- ggforestplot
- pheatmap
- ggpubr
- ggrepel
- cowplot
- egg

### Content
#### Association analyses

Example template for running the association scan of the individual NMR
biomarkers across disease endpoints is located under `association-analysis`.

In order to run the analyses, you need to implement the data reading and writing as
required by the available data infrastructure and location of the data. The script assumes you have available the
biomarker data, variables required for adjusting the analyses (age, sex, UK
biobank assessment center) and pre-defined endpoint data (as separate files for
each of the endpoints analysed).

The raw UK biobank data used for the analyses cannot be shared, but is available
for approved researchers through UK biobank. Summary statistics from the
association analyses in UK biobank are made publicly available for download in a
free webtool: https://nightingalehealth.com/atlas.


#### Figures

Scripts for making the figures presented in the manuscript, based on summary
statistics from the association analyses, are located under `figure-creation`.

For making the figures, you will need to have summary statistics from the association analysis available under `data` directory. The summary statistics from the association analyses are publicly available for download from
https://nightingalehealth.com/atlas, which we anticipate will help future
studies for replicating and extending the present results.

