## ----globals -------------------------------------------------------------
dir <- paths

## ----load.data -------------------------------------------------------------
gses <- as.data.frame(
            data.table::fread(file.path(dir$output, "gses.csv"), 
                            na.strings = "NA"))

