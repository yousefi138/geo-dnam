## ----globals -------------------------------------------------------------
packages <- c(
  "geograbi", # https://github.com/yousefi138/geograbi
  "eval.save",
  "purrr", # for map2)
  "easyPubMed",
  "dplyr",
  "parallel")
lapply(packages, require, character.only=T)

options(mc.cores=20)

# set dirs  
dir <- paths
if(!dir.exists(dir$cache)) dir.create(dir$cache)
eval.save.dir(dir$cache)

## check or make all output dirs
c("series.files", "dnam") |>
  sapply(function(.x){
    if (!dir.exists(file.path(dir$output,.x)))
      dir.create(file.path(dir$output,.x), recursive = T)
    dir.exists(file.path(dir$output,.x))		
  })

## source py written functions
#source(file.path(dir$scripts, "src/get.characteristics.r"))
source(file.path(dir$scripts, "src/retrieve.papers.r"))

## ----gses -------------------------------------------------------------
eval.save({
  ## get gses and basic info on all methylation array records
  gses <- geograbi.retrieve.datasets(
    gpl = c(
      "GPL13534", # 450k
      "GPL21145", "GPL21145", # EPIC
      "GPL33022")) |> ## 5sec
    subset(organism =="Homo sapiens") |>
    subset(samples >100)
  
}, "gses", redo=F)
gses <- eval.ret("gses")

## ----series.files -------------------------------------------------------------
eval.save({
  ## download the series files, keep filename
  ## ## 2 hours
  filenames <- geograbi.download.series.files( 
    path = file.path(dir$output, "series.files"),
    gses = gses$accession, 
    gpls = gses$platform)
  
}, "filenames", redo=F)
filenames <- eval.ret("filenames")
gses$filenames <- ifelse(filenames == "", NA, filenames)

## ----get.pmids -------------------------------------------------------------
eval.save({
  ## get pmids from series files
  series <- lapply(na.omit(gses$filenames), geograbi.get.series) # 2mins
  keep <-sapply(series, function(i) any(grepl("pubmed", names(i))))
  series<- series[keep]
}, "series", redo=F)
series <- eval.ret("series")

series <- lapply(series, function(i){
  keep <- c(grep("accession", names(i)), grep("pubmed", names(i)))[1:2]
  i[keep]
})
series <- as.data.frame(do.call(rbind, series), 
                        stringsAsFactors = F)
names(series)[names(series) == "pubmed_id"] <- "pmid"

## add pmids to gses
gses <- merge(gses, series, 
              by.x="accession", 
              by.y="geo_accession",
              all.x =T)

## ----get.data -------------------------------------------------------------
eval.save({
  dnam <- 
    purr::map2(gses$filenames, gses$accession, ~ {
      path <- file.path(dir$output, "dnam")		
      ret <-tryCatch(geograbi.get.data(filename = .x,
                                       gse = .y,
                                       path = path),
                     error=function(e) e)
      
      if (any(grepl("error", class(ret)))){
        return(c(class(ret), dim(ret), length(ret)))
      }
      
      if (length(ret)==1 && ret==-1){
        gse <- tolower(.y)
        data.filename <- file.path(path, 
                                   paste(tolower(gse), "csv.gz", sep="."))
        ret <- as.data.frame(data.table::fread(data.filename))
      }
      c(class(ret)[1], dim(ret), length(ret))
    })
}, "dnam", redo=F)
dnam <- eval.ret("dnam")

# make a nice summry data frame of  the success status of the download
names(dnam) <- gses$accession
dnam <- lapply(dnam, function(i){
  if (any(grepl("Error|numeric", i))) {
    c("error", rep("0",3))
  } else i
})
dnam <- do.call(rbind, dnam)
identical(gses$accession, rownames(dnam))

colnames(dnam) <- c("class", "nrow", "ncol", "length")
dnam <- data.frame(dnam, stringsAsFactors = F)
dnam[grep("[0-9]", dnam[1,])] <- 
  lapply(dnam[grep("[0-9]", 
                   dnam[1,])], as.numeric)

gses <- cbind(gses, dnam)

## ----supplementary-files------------------------------------------------------

## retrieve supplementary lists for each dataset
eval.save({
  mclapply(1:nrow(gses), function(i) {
    ret = NULL
    if (gses$class[i] != "error" && gses$nrow[i] == 0) {
      try({
        ret = geograbi.list.supplementary.files(gses$accession[i])
      })
      ret
    }})
}, "supplementary-files", redo=F)
supp.files <- eval.ret("supplementary-files")

## filename clues that a supplementary file contains normalized DNA methylation data
## or, if none exists, the raw IDAT files
omit.patterns <- paste0("(",paste(c("intensities","unmeth","unprocessed","non[-_]+normalized","unnormalized","unnormalised","raw","xlsx","tumor"),collapse="|"),")")
keep.patterns <- paste0("(",paste(c("normalized","normalised","methylation","processed","beta"),collapse="|"),")")

## uses patterns to guess appropriate supplementary file
## for datasets without series matrix files
supp.files = sapply(supp.files, function(files) {
  keep = files
  if (length(files) > 0) {
    keep = files[!grepl(omit.patterns,files,ignore.case=T) & grepl(keep.patterns,files,ignore.case=T)]
  }
  if (length(keep) == 0)
    keep = files[grepl("raw.tar$",files,ignore.case=T)]
  keep
})

## omit datasets based on tissue from children or that is not blood
for (i in which(gses$accession %in% readLines("omit.txt")))
  supp.files[[i]] = list()

## adhoc corrections
supp.files[gses$accession == "GSE191082"] = "GSE191082_Normalized_data_Blood.txt.gz"
supp.files[gses$accession == "GSE74414"] = "GSE74414_all_beta_values.txt.gz"
supp.files[gses$accession == "GSE198904"] = "GSE198904_DHRC_Matrix_processed.txt.gz;GSE198904_OBSERVEMDD0001_Matrix_processed.txt.gz"
supp.files[gses$accession == "GSE54882"] = "GSE54882_processed_methylation_matrix.txt.gz"
stopifnot(all(sapply(supp.files,length) <= 1))

gses$supp.file = sapply(supp.files,function(txt) if (length(txt) == 0) "" else unlist(txt))

# ----blood.pred -------------------------------------------------------------
## retrieve sample info for each dataset
eval.save({
  mclapply(1:nrow(gses), function(i) {
    ret = NULL
    	try({
        ret = geograbi.get.samples(gses$filename[i])
      })
      ret
    })
}, "samples", redo=F)
samples <- eval.ret("samples")

## extract just the info on sample source
source <- mclapply(samples, function(i){
		i$source_name_ch1
})

omit.patterns <- paste0("(",paste(c("tumour","spot","tumor", "cancer", "placenta","tissue"),collapse="|"),")")
keep.patterns <- paste0("(",paste(c("blood", "PBMC"),collapse="|"),")")
 
## uses patterns to guess if the source is likely to be  peripheral blood
blood.pred = sapply(source, function(samp) {
	keep = NA
  if (length(samp) > 0) {
    keep = !grepl(omit.patterns,samp,ignore.case=T) & grepl(keep.patterns,samp,ignore.case=T)
  }
  keep
})
gses$blood.pred = sapply(blood.pred, all)

## check on count
#sum(gses[gses$blood.pred == TRUE & !is.na(gses$blood.pred), "ncol"])
# [1] 48357

## ----pubmed -------------------------------------------------------------
eval.save({
  pubmed <- retrieve.papers(gses$pmid) 
}, "pubmed", redo=F)
pubmed <- eval.ret("pubmed")

gses <- pubmed |>
  rename(pub.title = title) |>
  select(all_of(c("pmid", 
                  "pub.title", "abstract")))  %>%
  left_join(gses, ., by = c("pmid"))


## ----write.gses -------------------------------------------------------------
omit.list <- readLines("omit.txt")

data.table::fwrite(
  gses[!gses$accession %in% omit.list,], 
  file.path(dir$output, "gses.csv"),
  row.names = F,
  na = "NA")
