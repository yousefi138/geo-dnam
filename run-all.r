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

## 
## in:
## out:
#system("R CMD BATCH --vanilla example.r")
source("get-geo-dnam.r", echo=T, max.deparse.length = 500)

## run analysis looking at relationship between
## methylation predicted proteins and 
## tumor vs. normal tissue type. 
## render an html summary
packages <- c("rmarkdown", "knitr", "kableExtra")
lapply(packages, require, character.only=T)

render("report.rmd", 
	output_format = "all")
