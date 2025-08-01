args <- commandArgs(trailingOnly=TRUE)

config.name <- "default"
if (length(args) > 0)
    config.name <- args[1]

paths <- config::get(config=config.name)
print(paths)

paths$data <- file.path(paths$project, "data")
paths$output <- file.path(paths$project, "results")
paths$cache <- file.path(paths$project, "results", "analysis-cache")
print(paths)

## dataset wrangling
## 
## in: paths
## out: 
## 1. [paths$output]/blood.gses.csv:
##   * lists datasets with DNAm derived from blood
##   * column 'dnam.file' has the file path to methylation data for each dataset
##     - these files have csv.gz format columns=samples, rows=cpg sites
##     - these will have small numbers of missing values encoded as 'NA'
##   * column 'counts.file' has the file path to cell count estimates for each dataset
##     (these files have .csv.gz format columns=cell types, rows=samples)
## 2. [paths$output]/common-cpgs.txt: list of cpg sites common to 75% of the datasets
#system("R CMD BATCH --vanilla example.r")
source("get-geo-dnam.r", echo=T, max.deparse.length = 500)


## run analysis looking at relationship between
## methylation predicted proteins and 
## tumor vs. normal tissue type. 
## render an html summary
packages <- c("rmarkdown", "knitr", "kableExtra", "tidyverse")
lapply(packages, require, character.only=T)

render("report.rmd", 
	output_format = "all")
