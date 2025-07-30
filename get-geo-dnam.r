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
c("series.files", "dnam", "samples", "supp.files") |>
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
    purrr::map2(gses$filenames, gses$accession, ~ {
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
gses$dnam.file <- 
		file.path(dir$output, "dnam", 
          paste(tolower(gses$accession), "csv.gz", sep="."))	


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
names(samples) <- gses$accession

# write all the sample info to file
samples.file <- purrr::imap(samples, ~ {
		if(!is.null(.x)){
			file <- file.path(dir$output, "samples", 
						paste0(tolower(.y), ".csv.gz"))
			if(!file.exists(file)) data.table::fwrite(.x, file)
			file
		} else NA
}) 
gses$samples.file <- unlist(samples.file)

## extract just the info on sample source
source <- mclapply(samples, function(i){
		i$source_name_ch1
})

omit.patterns <- paste0("(",paste(c("tumour","spot","tumor", "cancer", "placenta","tissue"),collapse="|"),")")
keep.patterns <- paste0("(",paste(c("blood", "WB", "PBMC"),collapse="|"),")")
 
## uses patterns to guess if the source is likely to be  peripheral blood
blood.pred = sapply(source, function(samp) {
	keep = NA
  if (length(samp) > 0) {
    keep = !grepl(omit.patterns,samp,ignore.case=T) & grepl(keep.patterns,samp,ignore.case=T)
  }
  keep
})
gses$blood.pred = sapply(blood.pred, all)

table(gses$blood.pred, useNA = "ifany")
#FALSE  TRUE  <NA> 
#  355   120   201 

## adhoc additions
gses$blood.pred <- (
  gses$blood.pred
  | gses$accession %in% c("GSE40279","GSE42861","GSE51032","GSE51057","GSE74548"))

## corrections
gses$ncol[gses$accession == "GSE42861"] = 689
gses$ncol[gses$accession == "GSE51032"] = 845

table(gses$blood.pred, useNA = "ifany")
## FALSE  TRUE  <NA> 
##   352   125   199 

#str(source[which(gses$blood.pred==F)])
## check on count
#sum(gses[gses$blood.pred == TRUE & !is.na(gses$blood.pred), "ncol"])
# [1] 51050

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

## ----supp.files------------------------------------------------------

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

## omit datasets known to be based on tissue from children or that is not blood
for (i in which(gses$accession %in% readLines("omit.txt")))
  supp.files[[i]] = list()

## adhoc corrections
supp.files[gses$accession == "GSE191082"] = "GSE191082_Normalized_data_Blood.txt.gz"
supp.files[gses$accession == "GSE74414"] = "GSE74414_all_beta_values.txt.gz"
supp.files[gses$accession == "GSE198904"] = "GSE198904_DHRC_Matrix_processed.txt.gz;GSE198904_OBSERVEMDD0001_Matrix_processed.txt.gz"
supp.files[gses$accession == "GSE54882"] = "GSE54882_processed_methylation_matrix.txt.gz"
stopifnot(all(sapply(supp.files,length) <= 1))

gses$supp.file = sapply(supp.files,function(txt) if (length(txt) == 0) "" else unlist(txt))

## ----write.gses -------------------------------------------------------------

# full query results
data.table::fwrite(
  gses, 
  file.path(dir$output, "gses.csv"),
  row.names = F,
  na = "NA")

# subset good blood query results
blood.gses <- gses |>
  subset(blood.pred == TRUE & !is.na(blood.pred)) |>
  subset(nrow >300000 | supp.file != "")

data.table::fwrite(
  blood.gses,
  file.path(dir$output, "blood.gses.csv"),
  row.names = F,
  na = "NA")

## ----download.supp.files ----------------------------------------------------------

supp.indices <- which(blood.gses$supp.file != "")
supp.dir <- file.path(dir$output, "supp.files")
for (i in supp.indices) {
  filenames = unlist(strsplit(blood.gses$supp.file[i], ";"))
  for (filename in filenames) {  
    if (!file.exists(file.path(supp.dir, filename)))
      geograbi.download.supplementary.file(
        path=supp.dir,
        gse=blood.gses$accession[i],
        filename=filename)
  }
}

## ----check.supp.files -------------------------------------------------

supp.filenames <- with(blood.gses, supp.file[grepl("(csv|tsv|txt).gz$", supp.file)])
supp.filenames <- unlist(strsplit(supp.filenames,";"))

## extract the first few rows of each supplementary file
dir.create(peak.dir <- file.path(supp.dir, "peak"))
for (filename in supp.filenames) {
  if (!file.exists(file.path(peak.dir, filename)))
    system(paste(
      "gunzip -c", file.path(supp.dir, filename),
      "| head -n 5 | gzip -c >", file.path(peak.dir, filename)))
}

## check the first few rows of each supplementary file
for (filename in supp.filenames) {
  cat(date(), " checking supplementary file ", filename, "\n")
  dat <- fread(file.path(peak.dir, filename), data.table=F)
  cg.idx = which(colnames(dat) == "ID_REF")
  cg.idx = ifelse(length(cg.idx) == 0, 1, cg.idx[1])
  rownames(dat) <- dat[[cg.idx]]
  dat[[cg.idx]] = NULL
  dat <- as.matrix(dat[,!grepl("(cpg|pval|id_ref)",colnames(dat),ignore.case=T)])
  stopifnot(ncol(dat) >= 100)
  rows <- rownames(dat)
  stopifnot(!any(is.na(rows)))
  stopifnot(all(substring(rows,1,2) == "cg"))
  dat <- apply(dat,2,as.numeric)
  stopifnot(sum(!is.na(dat)) > ncol(dat))
}

## ----clean.supp.files -------------------------------------------------

## clean each supplementary file
supp.filenames <- with(blood.gses, supp.file[grepl("(csv|tsv|txt).gz$", supp.file)])
supp.filenames <- unlist(strsplit(supp.filenames,";"))
ids = unique(sub("_.*", "", supp.filenames))
dir.create(save.dir <- file.path(supp.dir, "clean"))
for (id in ids) {
  out.filename <- file.path(save.dir, paste0(tolower(id),".csv.gz"))
  if (!file.exists(out.filename)) {
    cat(date(), " cleaning supplementary file and saving to ", out.filename, "\n")
    is_filename = grepl(paste0(id,"_"), supp.filenames)
    dats <- lapply(supp.filenames[is_filename], function(filename) {  
      dat <- fread(file.path(supp.dir, filename), data.table=F)
      cg.idx = which(colnames(dat) == "ID_REF")
      cg.idx = ifelse(length(cg.idx) == 0, 1, cg.idx[1])
      rownames(dat) <- dat[[cg.idx]]
      dat[[cg.idx]] = NULL
      dat <- as.matrix(dat[,!grepl("(cpg|pval|id_ref)",colnames(dat),ignore.case=T)])
      rows <- rownames(dat)
      dat <- apply(dat,2,as.numeric)
      rownames(dat) = rows
      ms <- colMeans(dat,na.rm=T)
      stopifnot(!all(ms < 0.1)) ## if no samples have normal looking DNAm
      dat[,ms > 0.1]            ## omit samples with abnormal DNAm
    })
    if (length(dats) > 1) {
      rows <- unlist(lapply(dats,rownames))
      rows <- names(which(table(rows)==length(dats)))
      dats <- lapply(dats, function(dat) dat[match(rows,rownames(dat)),])
      dat <- do.call(cbind, dats)
    }
    else dat <- dats[[1]]
    fwrite(data.frame(cg=rownames(dat), dat, check.names=F), file=out.filename)
  }
}

## ----idats.extraction --------------------------------------------------------------------

arc.filenames <- with(blood.gses, supp.file[grepl("_RAW.tar$", supp.file)])
dir.create(idats.dir <- file.path(supp.dir, "idats"))
for (filename in arc.filenames) {
  contents.out <- file.path(idats.dir, sub(".tar",".out",filename))
  if (!file.exists(contents.out)) {
    cat(date(), " extracting idats from ", filename, "\n")
    system(paste("tar -tf", file.path(supp.dir, filename), ">", contents.out))
    system(paste("tar -xf", file.path(supp.dir, filename), "-C", idats.dir))
  }
}

## ----final.counts -------------------------------------------------------------

## number datasets (not)derived from blood with a sufficient number probes
table(
  blood=gses$blood.pred,
  probes300K=gses$nrow > 300000 | gses$supp.file != "")
##        probes300K
## blood   FALSE TRUE
##   FALSE    42  310
##   TRUE     10  115

## number samples derived from blood and/or with a sufficient number probes
with(gses, {
  data.frame(
    blood=sum(ncol[blood.pred], na.rm=T),
    probes300K=sum(ncol[nrow > 300000 | supp.file != ""]),
    both=sum(ncol[blood.pred & (nrow > 300000 | supp.file != "")]))})
##   blood probes300K  both
## 1 51050     127685 46068

## number of datasets with data coming from a supplementary file
table(blood.gses$supp.file != "")
## FALSE  TRUE 
##    46    69 

## number of datasets with just idats available
sum(grepl("_RAW.tar", blood.gses$supp.file))
## [1] 6

## number of samples with just idats available
sum(blood.gses$ncol[grepl("_RAW.tar", blood.gses$supp.file)])
## [1] 5636
