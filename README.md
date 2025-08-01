# GEO DNAm Data

This repository proivdes a workflow to identify and access all available DNAm datasets provided by the Gene Expression Omnibus (GEO) https://www.ncbi.nlm.nih.gov/geo/.

# GEO Wrangling Workflow

The main workflow script is:

* [`get-geo-dnam.r`](https://github.com/yousefi138/geo-dnam/blob/main/get-geo-dnam.r)

This script is sourced from `run-all.r`. Details on the worklfow inputs and outputs are both described there and below:

## in: paths

## out: 

1. [paths$output]/blood.gses.csv:

  * lists datasets with DNAm derived from blood
  * column 'dnam.file' has the file path to methylation data for each dataset
    - these files have csv.gz format columns=samples, rows=cpg sites
    - these will have small numbers of missing values encoded as 'NA'
  * column 'counts.file' has the file path to cell count estimates for each dataset
    (these files have .csv.gz format columns=cell types, rows=samples)

2. [paths$output]/common-cpgs.txt: list of cpg sites common to 75% of the datasets


# Workflow output report

The main summary of the results is available in the below report file:

* [`report.html`](https://github.com/yousefi138/geo-dnam/blob/main/report.html)

