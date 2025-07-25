args <- commandArgs(trailingOnly=TRUE)

config.name <- "default"
if (length(args) > 0)
    config.name <- args[1]

paths <- config::get(config=config.name)
print(paths)

paths$data <- file.path(paths$project, "data")
paths$output <- file.path(paths$project, "results")
paths$cache <- file.path(paths$project, "results", "analysis-cache")
paths$scripts <- file.path(paths$project, "scripts")
print(paths)

## 
## in:
## out:
#source("example.r", echo=T, max.deparse.length = 500)
#system("R CMD BATCH --vanilla example.r")


## run analysis looking at relationship between
## methylation predicted proteins and 
## tumor vs. normal tissue type. 
## render an html summary
packages <- c("rmarkdown", "knitr")
lapply(packages, require, character.only=T)

render("analysis.rmd", 
	output_format = "all",
    output_dir = "docs")
