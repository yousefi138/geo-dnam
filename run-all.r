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


## updates
##
## Dec 12, 2025
## validation datasets
## * GSE107351 -- 116 of which 12 MLH1 epimutation carriers, 61 Lynch syndrome patients and 41 healthy controls
## * GSE183920 -- 601 bladder cancer cases of which 371 bad outcomes (recurrence or death)
## * GSE281691 -- 588 of which 272 metabolic liver cancer cases
## * GSE51032 -- 845 of which 235 later developed breast cancer, 166 colorectal cancer and 20 developed other primary cancers.
## omitted two blood datasets:
## * GSE97362, all (235) CHARGE and Kabuki syndrome patients
## * GSE188949, all (107) cord blood from preterm infants
##
## June 9, 2026
## some datasets had missing values for hundreds of DNAm scores
## this was due to the many small Waterfield et al. models
fnames = list.files(file.path(dir$output, "methylation-scores"),full.names=T)
counts = sapply(fnames, function(fname) {
  dat = fread(fname,data.table=F)
  rownames(dat) = dat[[1]]
  dat = dat[,-1]
  dat = as.matrix(dat)
  counts = apply(dat, 2, function(x) mean(!is.na(x)))
  dat = dat[!grepl("waterfield", rownames(dat)),,drop=F]
  counts.after = apply(dat,2, function(x) mean(!is.na(x)))
  cat(date(), fname, median(counts), median(counts.after), "\n")
  c(before=quantile(counts), after=quantile(counts.after))
})

colnames(counts) = sub(".csv.gz","",basename(colnames(counts)))

## before removing Waterfield
sum(counts["before.50%",] < 0.95)
sum(counts["before.25%",] < 0.95)
##[1] 47

## after removing Waterfield
sum(counts["after.50%",] < 0.95)
sum(counts["after.25%",] < 0.95)
#[1] 1
