## ----globals -------------------------------------------------------------
packages <- c(
  "meffil",
  "geograbi", # https://github.com/yousefi138/geograbi
  "eval.save",
  "purrr", # for map2)
  "easyPubMed",
  "dplyr",
  "parallel",
  "meffonym", ## github.com/perishky/meffonym
  "RCurl",
  "clue"
)
lapply(packages, require, character.only=T)


options(mc.cores=20)

# set dirs  
dir <- paths
if(!dir.exists(dir$cache)) dir.create(dir$cache)
eval.save.dir(dir$cache)

## check or make all output dirs
c("series.files", "dnam", "samples", "supp.files", "matching",
  "cell-counts", "dnam-qn","methylation-scores","methylation-scores-qn") |>
  sapply(function(.x){
    if (!dir.exists(file.path(dir$output,.x)))
      dir.create(file.path(dir$output,.x), recursive = T)
    dir.exists(file.path(dir$output,.x))		
  })

## source py written functions
#source(file.path(dir$scripts, "src/get.characteristics.r"))
source(file.path(dir$scripts, "src/retrieve.papers.r"))

## ----gses-------------------------------------------------------------

valid_platforms = c(
  "GPL13534", # 450k
  "GPL16304", # 450k
  "GPL23976", # 850k
  "GPL21145", # EPIC
  "GPL33022") # EPICv2

eval.save({
  ## get gses and basic info on all methylation array records
  gses <- geograbi.retrieve.datasets(
    gpl = valid_platforms) |>
    subset(organism =="Homo sapiens") |>
    subset(samples >100)  
}, "gses", redo=F)
gses <- eval.ret("gses")

gses$multiple_platforms = grepl(" ", gses$platform)

platforms = lapply(strsplit(gses$platform," "), function(ps) intersect(ps,valid_platforms))

gses = gses[rep(1:nrow(gses), sapply(platforms, length)),]
gses$platform = unlist(platforms)

gses$omit = F
omit_accessions = sub("[ #]+.*", "", readLines("omit.txt"))
gses$omit[gses$accession %in% omit_accessions] = T

gses$include = F
include_accessions = sub("[ #]+.*", "", readLines("include.txt"))
gses$include[gses$accession %in% include_accessions] = T

stopifnot(all(include_accessions %in% gses$accession))
stopifnot(all(!(gses$include & gses$omit)))

## ----series.files --------------------------------------------------

get.gse.matrix.url = function(gse,gpl) {
  prefix = sapply(gse, geograbi:::gse.prefix)
  paste(
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
    prefix, "/", 
    gse, "/matrix/",
    gse, ifelse(gpl!="", "-", ""), gpl,
    "_series_matrix.txt.gz", sep = "")
}

platforms = gses$platform
platforms[!gses$multiple_platforms] = ""
gses$series.url = get.gse.matrix.url(
  gses$accession,
  platforms)

series.dir = file.path(dir$output, "series.files")
dir.create(series.dir)
gses$series.filename = file.path(series.dir, basename(gses$series.url))

needs.download = !gses$omit & !file.exists(gses$series.filename)
for (i in which(needs.download)) {  
  cat(date(), "downloading", basename(gses$series.filename[i]),
      which(which(needs.download) == i), "of", sum(needs.download), 
      "\n")
  tryCatch({
    download.file(gses$series.url[i], gses$series.filename[i])
  }, error=function(e) {
    if (file.exists(gses$series.filename[i]))
      unlink(gses$series.filename[i])
  })
  Sys.sleep(1)
}

stopifnot(all(gses$omit | file.exists(gses$series.filename)))

## ----identify.blood.datasets------------------------------------

gses$samples.file = file.path(dir$output, "samples", paste0(tolower(basename(gses$series.filename))))
gses$samples.file = sub("_series_matrix.txt.gz$", ".csv.gz", gses$samples.file)

needs.extraction = !gses$omit & !file.exists(gses$samples.file)
for (i in which(needs.extraction)) {
  cat(date(), "extracting", basename(gses$series.filename[i]),
      which(which(needs.extraction) == i), "of", sum(needs.extraction), 
      "\n")
  tryCatch({
    ret = geograbi.get.samples(gses$series.filename[i])
    fwrite(ret, file=gses$samples.file[i])
  }, error=function(e) {
    if (file.exists(gses$samples.file[i]))
      unlink(gses$samples.file[i])
  })
}

stopifnot(all(gses$omit | file.exists(gses$samples.file)))

sources = mclapply(1:nrow(gses), function(i) {
  if (gses$omit[i]) return(NA)
  fread(gses$samples.file[i], data.table=F)$source_name_ch1
})

stopifnot(all(gses$omit | sapply(sources, length) > 1))

omit.patterns <- paste0("(",paste(c("saliva","buccal","bcell","cells","tumour","spot","tumor", "cancer", "placenta","tissue","FFPE","sarcoma","brain","metastasis","ovary","meningioma","bone","prostate","carcinoma","lung","lymphoma","epidermis","glioma","cortex","plasma","myometrium","lymph","unknown","kidney","liver","breast","nsclc","hcc","colon","organoid","monocyte","ependymoma","pancreas","blastoma","glss","tcga","dlpfc","nucleus","embryonic","villi","derm","thymoma","dfdna","epithelial","hesc","cultured","skin","vein","uterus","bladder","ureter","heart","resection","fibroblast","anal","cervical","T cell","epithelium","lateralis","cd8t","cbmc","myocyte","neutrophil"),collapse="|"),")")
keep.patterns <- paste0("(",paste(c("blood", "WB", "PBMC","PBL","buffy coat","leukocyte"),collapse="|"),")")


gses$is_blood = sapply(1:nrow(gses), function(i) {
  !gses$omit[i] && (all(grepl(keep.patterns,sources[[i]],ignore.case=T)) || gses$include[i])
})
gses$is_not_blood = sapply(1:nrow(gses), function(i) {
  gses$omit[i] || any(grepl(omit.patterns,sources[[i]],ignore.case=T))
})

stopifnot(all(gses$is_not_blood | gses$is_blood))

## ----get.datasets------------------------------------------------------

gses$dnam.file = file.path(dir$output, "dnam", paste(tolower(gses$accession), "csv.gz", sep="."))

needs.data = gses$is_blood & !file.exists(gses$dnam.file)

for (i in which(needs.data)) {
  cat(date(), "extracting DNAm data for", gses$accession[i], 
      which(which(needs.data) == i), "of", sum(needs.data), "\n")
  tryCatch({
    geograbi.get.data(
      filename = gses$series.filename[i],
      gse = gses$accession[i],
      path = file.path(dir$output, "dnam"))
    stopifnot(file.exists(gses$dnam.file[i]))
  }, error=function(e) {
    if (file.exists(gses$dnam.file[i]))
      unlink(gses$dnam.file[i])
  })
}

stopifnot(all(!gses$is_blood | file.exists(gses$dnam.file)))

## ----get.supp----------------------------------------------------------

supp.files = read.csv("supp-files.csv", stringsAsFactors=F)

platforms = strsplit(supp.files$platform,";")
supp.files = supp.files[rep(1:nrow(supp.files), sapply(platforms,length)),]
supp.files$platform = unlist(platforms)

gses_id = paste(gses$accession, gses$platform)
supp_id = paste(supp.files$accession, supp.files$platform)

needs.supp = gses$is_blood & (gses_id %in% supp_id | file.size(gses$dnam.file) < 100000)

stopifnot(all(gses_id[needs.supp] %in% supp_id))

supp.dir = file.path(dir$output, "supp.files")
dir.create(supp.dir)
supp.download.dir = file.path(supp.dir,"downloads")
dir.create(supp.download.dir)

gses$supp.file = ""
gses$supp.file[needs.supp] = supp.files$supp.file[match(gses_id[needs.supp],supp_id)]
gses$supp.file[needs.supp] = sapply(strsplit(gses$supp.file[needs.supp], ";"), function(supp.file)
  paste0(file.path(supp.download.dir, basename(supp.file)), collapse=";"))

cat(sum(!file.exists(file.path(supp.download.dir, unlist(strsplit(supp.files$supp.file,";"))))),
    "supplementary files to download\n")

for (i in 1:nrow(supp.files)) {
  filenames = basename(unlist(strsplit(supp.files$supp.file[i], ";")))
  if (!all(file.exists(file.path(supp.download.dir, filenames)))) {
    cat(date(), "downloading supplementary file(s) for", supp.files$accession[i], "\n")
    for (filename in filenames) {
      geograbi.download.supplementary.file(
        path=supp.download.dir,
        gse=supp.files$accession[i],
        filename=filename)
    }
  }
 }

## ----clean.supp.files -------------------------------------------------

## convert supplementary files to csv files where rows=cpgs and columns=samples and the first column is the cpg identifier

source("gread-function.r")

supp.clean.dir = file.path(supp.dir, "clean")
dir.create(supp.clean.dir)

gses$dnam.file[needs.supp] = file.path(
  supp.clean.dir,
  paste0(tolower(gses$accession[needs.supp]), "_", tolower(gses$platform[needs.supp]), ".csv.gz"))

needs.cleaning = needs.supp & !file.exists(gses$dnam.file) & !endsWith(gses$supp.file,".tar")

for (i in which(needs.cleaning)) {
  cat(date(), " cleaning supplementary file",
      which(i == which(needs.cleaning)), "(of",
      sum(needs.cleaning), "total) for", gses$accession[i], "\n")
  input.filenames = unlist(strsplit(gses$supp.file[i], ";"))
  dats = lapply(input.filenames, function(filename) {
    tryCatch({
      if (grepl(".RDS$",filename,ignore.case=T)) {
        dat = readRDS(filename)
      }
      else {
        sep = ifelse(endsWith(filename,".csv.gz"), "csv", "auto")
        dat = gread(filename, data.table=F, header=T,tmpdir=".")
      }
      if (ncol(dat) > nrow(dat)) {
        cat("transposing meth matrix so rows=cpgs and columns=samples\n")
        ids = dat[[1]]
        dat = t(dat[,-1])
        colnames(dat) = ids
        dat = data.frame(ID_REF=rownames(dat), dat, check.names=F)
      }
      ## identify cpg site identifiers
      cg.idx = which(colnames(dat) == "ID_REF")
      cg.idx = ifelse(length(cg.idx) == 0, 1, cg.idx[1])
      rows = na.omit(unique(dat[[cg.idx]]))
      if (mean(substring(rows,1,2) == "cg", na.rm=T) < 0.9) {
        cg.idx = 2
        rows = na.omit(unique(dat[[cg.idx]]))
      }
      stopifnot(mean(substring(rows,1,2) == "cg",na.rm=T) > 0.9)
      ## remove cpg identifier suffixes (they identify replicate probes)
      rows = sub("_.*$", "", rows)
      unique_rows = unique(rows)
      ## ensure one row per cpg site
      dat = dat[match(unique_rows,rows),]
      rownames(dat) = unique_rows
      ## remove extract non-data columns
      dat[[cg.idx]] = NULL
      if ("V0" %in% colnames(dat)) dat$V0 = NULL
      if ("V1" %in% colnames(dat)) dat$V1 = NULL
      ## remove any columns not supplying dnam levels (e.g. detection p-values, meth/unmeth signal)
      dat = as.matrix(dat[,!grepl("(cpg|pval|detectionp|methylated|id_ref)",colnames(dat),ignore.case=T)])
      ## ensure that methylation levels are numeric
      rows = rownames(dat)
      dat = apply(dat,2,as.numeric)
      rownames(dat) = rows
      ## remove samples with >80% missing values or no variance
      sample.var = apply(dat,2,var,na.rm=T)      
      pct.missing = apply(dat,2,function(v) mean(is.na(v)))
      is.bad = sample.var < 0.0001 | pct.missing > 0.8
      if (any(is.bad, na.rm=T))
        cat("removing", sum(is.bad,na.rm=T), "bad samples from", ncol(dat), "total\n")
      stopifnot(!all(is.bad)) 
      dat[,!is.bad]
      ## convert m-values to beta values
      vals = na.omit(as.vector(dat))
      if (min(vals) < -1 && max(vals) > 2) {
        cat("converting m-values to beta values\n")
        dat = apply(dat, 2, function(val) 2^val/(2^val + 1))
      }
      dat
    }, error=function(e) {
      browser()
      print(e)
    })
  })
  if (length(dats) > 1) {
    rows = unlist(lapply(dats,rownames))
    rows = names(which(table(rows)==length(dats)))
    dats = lapply(dats, function(dat) dat[match(rows,rownames(dat)),])
    dat = do.call(cbind, dats)
  }
  else dat = dats[[1]]
  fwrite(
    data.frame(cg=rownames(dat), dat, check.names=F),
    file=gses$dnam.file[i])
}

stopifnot(!needs.supp | file.exists(gses$dnam.file) | endsWith(gses$supp.file,".tar"))

## ----normalize.idats----------------------------------------

source("fpeak-function.r")

source("normalize-idat-archive-function.r")

normalization.dir = file.path(supp.dir, "normalization-outputs")
dir.create(normalization.dir)

needs.normalization = needs.supp & !file.exists(gses$dnam.file) & endsWith(gses$supp.file,".tar")

for (i in rev(which(needs.normalization))) {
  cat(date(), " normalizing idats for", gses$accession[i],
      ",", which(i==which(needs.normalization)), "of", sum(needs.normalization), "datasets\n")
  normalize.idat.archive(
    archive = gses$supp.file[i],
    norm.dir = file.path(normalization.dir, gses$accession[i]),
    output = gses$dnam.file[i])
}

for (i in which(needs.supp & endsWith(gses$supp.file,".tar"))) {
  samples = fread(gses$samples.file[i],data.table=F, header=T)
  meth = fpeak(gses$dnam.file[i],nrows=10,data.table=F,header=T,tmpdir=".")
  cat(date(), gses$accession[i], nrow(samples), ncol(meth)-1, "\n") 
  if (nrow(samples) > ncol(meth)-1) {
    is_good = samples$geo_accession %in% colnames(meth)
    cat("Removing", sum(!is_good) ,"poor quality samples from", gses$accession[i], "\n")
    stopifnot(sum(!is_good) < 5)
    samples = samples[is_good,]
    fwrite(samples, file=gses$samples.file[i])
  }
}

stopifnot(!gses$is_blood | file.exists(gses$dnam.file))

## ----check-dnam-levels------------------------------------------------

for (i in which(gses$is_blood)) {
  meth = fpeak(gses$dnam.file[i],nrows=20,data.table=F,header=T,tmpdir=".")
  meth_values = na.omit(as.vector(as.matrix(meth[,-1])))
  cat(date(), gses$accession[i], "min:", min(meth_values), "max:", max(meth_values), "\n")
  stopifnot(min(meth_values) >= -0.1 && max(meth_values) <= 1.1) 
}

## ----match.samples.to.meth---------------------

matching.dir = file.path(dir$output, "matching")
dir.create(matching.dir)
gses$matching.file = file.path(matching.dir, basename(gses$samples.file))

good.matching = function(ids,samples,meth) {
    (length(unique(na.omit(ids))) >= min(nrow(samples),ncol(meth)-1)*0.99
      && length(setdiff(ids,colnames(meth))) <= 0.01*length(ids))
}
save.matching = function(samples, samples_col, meth, ids, filename) {
  if (!good.matching(ids,samples,meth))
    return(F)
  matching = data.frame(
    geo_accession=samples$geo_accession,
    samples_col=samples_col,
    meth_id=ids)
  matching[[samples_col]] = samples[[samples_col]]
  fwrite(matching, file=filename)
  return(T)
}
norm_id = function(x,targets=x) {
  good.chars = paste(sort(unique(unlist(strsplit(tolower(targets),"")))),collapse="")
  gsub(paste0("[^",good.chars,"]+"), " ", tolower(x))
}

for (i in which(gses$is_blood)) {
  if (file.exists(gses$matching.file[i])) next
  samples = fread(gses$samples.file[i],data.table=F, header=T)
  cat(date(), "matching", gses$accession[i],
      "(",which(i==which(gses$is_blood)), "of", sum(gses$is_blood), ")\n")
  meth = fpeak(gses$dnam.file[i],nrows=10,data.table=F,header=T,tmpdir=".")
  ## try geo sample accession
  if (save.matching(samples,"geo_accession",meth,samples$geo_accession, gses$matching.file[i]))
    next
  ## try supplementary file name 
  if ("supplementary_file" %in% colnames(samples)) {
    ids = sub(
      "^(|GSM[0-9]+_)([0-9]+_R[0-9]+C[0-9]+)_(Red|Grn).idat.gz$",
      "\\2",
      basename(samples$supplementary_file))
    if (save.matching(samples,"supplementary_file",meth,ids,gses$matching.file[i]))
      next
  }
  ## try description field
  if ("description" %in% colnames(samples)) {
    ids = sub("(|.*: )", "", samples$description)
    if (save.matching(samples,"description",meth,ids,gses$matching.file[i]))
      next
  }
  ## try description.1 field
  if ("description.1" %in% colnames(samples)) {
    ids = sub("(|.*: )", "", samples$description.1)
    if (save.matching(samples,"description.1",meth,ids,gses$matching.file[i]))
      next
  }
  if ("title" %in% colnames(samples)) {
    ## check if meth colnames in title field
    meth_cols = norm_id(colnames(meth)[-1])
    sam_titles = norm_id(samples$title, colnames(meth)[-1])
    idx = sapply(meth_cols, function(col) {
      idx = grep(
        paste0("(^|.*[ ]+)", col, "([ ]+.*|$)"),
        sam_titles)
      ifelse(length(idx) == 1,idx,NA)
    })
    ids = colnames(meth)[-1][match(sam_titles,sam_titles[idx])]
    if (save.matching(samples,"title",meth,ids,gses$matching.file[i]))
      next
  }
  ## matchings for specific datasets #######################
  ids = ""
  samples_col = "title"
  if (gses$accession[i] == "GSE314180")
    ids = sub(".* ([0-9]+)$", "S\\1", samples$title)
  else if (gses$accession[i] == "GSE318177") {
    group = rep("", nrow(samples))
    group[grepl("healthy", samples$title)] = "HC"
    group[grepl("episodic", samples$title)] = "ESA"
    group[grepl("medication", samples$title)] = "MOH"
    ids = paste0(
      group,
      sub(
        ".* subject ([0-9]+) at timepoint (T[0-9]+)$",
        "\\1_\\2_BetaVal",
        samples$title))
  } else if (gses$accession[i] == "GSE308517")
    ids = colnames(meth)[-1] 
  else if (gses$accession[i] == "GSE273970")
    ids = sub("^sample ([0-9]+)_.*","sample_\\1.Noob_Beta",samples$title)
  else if (gses$accession[i] == "GSE293007") {
    ids = paste0(
      sub(".*: ", "", samples$characteristics_ch1),
      "_",
      sub(".*: ", "", samples$characteristics_ch1.1))
    samples_col="characteristics_ch1"
  } else if (gses$accession[i] == "GSE286347")
    ids = sub(".*\\[(.*)\\]$","\\1_BetaVal", samples$title)
  else if (gses$accession[i] == "GSE128235")
    ids = sub(".* sample ([0-9]+)$", "sample\\1", samples$title)
  else if (gses$accession[i] == "GSE93933")
    ids = sub(".* ([A-Za-z0-9-]+)$", "\\1_BetaValue",samples$title)
  else if (gses$accession[i] == "GSE107080") {
    ids = samples$description.2
    samples_col="description.2"
  } else if (gses$accession[i] == "GSE84727") {
    ids = sub("sentrixids: ", "", samples$characteristics_ch1)
    samples_col="characteristics_ch1"
  } else if (gses$accession[i] == "GSE305174")
    ids = sub("whole blood of subject ([0-9]+)","\\1", samples$title)
  else if (gses$accession[i] == "GSE224339") { 
    ids = sub(".* \\[([0-9A-Z]+)\\]$","\\1", samples$title)
    ids[ids=="CEPH"] = "CEPH_control"
  } else if (gses$accession[i] == "GSE240184") {
    group = rep("", nrow(samples))
    group[grepl("three month", samples$title)] = "T3M"
    group[grepl("baseline", samples$title)] = "BL"
    group[grepl("replicate", samples$title)] = paste0(
      group[grepl("replicate", samples$title)], ".1")
    ids = paste0(sub(".* ([0-9]+)(| replicate)$","\\1",samples$title), group)
  }
  if (save.matching(samples,samples_col,meth,ids,gses$matching.file[i]))
    next
  if ("title" %in% colnames(samples)) {
    meth_cols = norm_id(colnames(meth)[-1])
    sam_titles = norm_id(samples$title,colnames(meth)[-1])
    if (length(meth_cols) >= length(sam_titles)) {
      dist = adist(sam_titles,meth_cols,partial=F)
      matching = solve_LSAP(dist) ## each sample index in order => meth index
      ids = colnames(meth)[-1][matching]
    } else {
      dist = adist(meth_cols,sam_titles,partial=F)
      matching = solve_LSAP(dist) ## each meth index in order => sample index
      rev_matching = match(samples$title, samples$title[matching])
      ids = colnames(meth)[-1][rev_matching]
    }
    if (save.matching(samples,"title",meth,ids,gses$matching.file[i])) {
      next
    }
  }
}

stopifnot(all(file.exists(gses$matching.file[gses$is_blood])))


## ----calculate.sample.counts ---------------------------------------

gses$num.matched = NA
gses$num.meth = NA
gses$num.samples = NA

for (i in which(gses$is_blood)) {
  cat(date(), "counting samples in", gses$accession[i],
      "(", which(which(gses$is_blood)==i), "of", sum(gses$is_blood), ")\n")
  samples = fread(gses$samples.file[i],data.table=F, header=T)
  meth = fpeak(gses$dnam.file[i],nrows=10,data.table=F,header=T,tmpdir=".")
  matching = fread(gses$matching.file[i], data.table=F)
  gses$num.matched[i] = nrow(matching)
  gses$num.meth[i] = ncol(meth)-1
  gses$num.samples[i] = nrow(samples)
}

sum(gses$num.matched,na.rm=T)
## 66019
sum(gses$num.meth,na.rm=T)
## 69068

## ----common.cpgs ------------------------------------------------------------

cpg.filename = file.path(dir$output, "common-cpgs.txt")

if (!file.exists(cpg.filename)) {
  cpgs = eval.save({
    extract.cpgs = function(filename) {
      cat(date(), " extracting cpg sites from ", filename, "\n") 
      system(paste("gunzip -c", filename, "| cut -d ',' -f1"), intern=T)[-1]
    }    
    ret = sapply(gses$dnam.file[gses$is_blood], extract.cpgs, simplify=F)
    names(ret) = sub(".csv.gz$", "", basename(names(ret)))
    ret
  }, "cpgs", redo=F) ## ~1.5 hours
  cpg.freqs = table(unlist(sapply(cpgs, function(cpgs) unique(sub("_.*","", cpgs)), simplify=F)))
  ## beware: cg00000029_TC21
  underscores = sapply(cpgs, function(cpgs) sum(grepl("_", cpgs)))
  table(underscores) ## numbers of cpg ids with underscores per dataset
  ## underscores
  ##   0      1      4 854671 883438 937690 
  ## 131     49      1      1      1      1
  names(which(underscores > 850000))
  ## [1] "gse338167_gpl33022" "gse255870_gpl33022" "gse246337_gpl33022"
  common.cpgs = names(cpg.freqs)[which(cpg.freqs > length(cpgs)*0.6)]
  length(common.cpgs)
  ## [1] 409295
  writeLines(common.cpgs, con=cpg.filename)
} else {
  common.cpgs = readLines(cpg.filename)
}

## -----remove.bad.samples---------------------------------------------

for (i in which(gses$accession %in% c("GSE68194","GSE69138","GSE59592"))) {
  cat(date(), "removing samples with too much missing data from ", gses$accession[i], "\n")
  new.file = sub(".csv.gz$", "fixed-samples.csv.gz", gses$dnam.file[i])
  if (file.exists(new.file)) {
    gses$dnam.file[i] = new.file
    next
  }
  meth = fread(gses$dnam.file[i], data.table=F, header=T, tmpdir=".")
  rownames(meth) = meth[[1]]
  meth[[1]] = NULL
  pct.missing = apply(meth,2,function(v) mean(is.na(v)))
  is.bad = pct.missing > 0.8
  if (any(is.bad)) {
    cat("removing", sum(is.bad), "samples with >80% missing data from", ncol(meth), "total\n")
    meth = meth[,!is.bad]
  }
  fwrite(data.frame(cg=rownames(meth), meth, check.names=F), file=new.file)
  gses$dnam.file[i] = new.file
}

## -----cell.counts -----------------------------------------------------------

dir.create(counts.dir <- file.path(dir$output, "cell-counts"))
gses$counts.file = file.path(counts.dir, basename(gses$dnam.file))
needs.counts = gses$is_blood & !file.exists(gses$counts.file) 

for (i in which(needs.counts)) {
  cat(date(), "estimating cell counts for", gses$accession[i],
      "(", which(i==which(needs.counts)), "of", sum(needs.counts), ")\n")
  meth = fread(gses$dnam.file[i], data.table=F, header=T, tmpdir=".")
  rownames(meth) = meth[[1]]
  meth[[1]] = NULL
  try({
    counts = meffil.estimate.cell.counts.from.betas(as.matrix(meth), "blood gse35069 complete")
    fwrite(data.frame(sample=rownames(counts),counts), file=gses$counts.file[i])
  })
} ## about 3 hours

stopifnot(all(sapply(gses$counts.file[gses$is_blood], file.exists)))


## ----quantile-normalize-dnam-----------------------------------------

gses$id = paste0(gses$accession, "_", gses$platform)
targets = sapply(which(gses$is_blood), function(i) {
  in.file = gses$dnam.file[i]
  eval.save({
    meth = fread(in.file, data.table=F,tmpdir=".",header=T)
    rownames(meth) = meth[[1]]
    meth[[1]] = NULL
    meth = as.matrix(meth)
    meth = meth[match(common.cpgs, rownames(meth)),]
    rownames(meth) = common.cpgs
    normalize.quantiles.determine.target(meth, subset=common.cpgs)
  }, paste0("quantile-normalization-target-", gses$id[i]), redo=F)
}) ## 2 hours
colnames(targets) = gses$id[gses$is_blood]

target = targets %*% matrix(gses$num.meth[gses$is_blood],ncol=1)/sum(gses$num.meth[gses$is_blood])
target = as.vector(target)

qn.dir = file.path(dir$output, "dnam-qn")
dir.create(qn.dir, showWarnings = F)
gses$dnam.qn.file = ""
gses$dnam.qn.file[gses$is_blood] = file.path(qn.dir, basename(gses$dnam.file[gses$is_blood]))

for (i in which(gses$is_blood)) {
  in.file = gses$dnam.file[i]
  out.file= gses$dnam.qn.file[i]
  if (file.exists(out.file)) {
    cat(date(), " skipping ", gses$id[i], " as quantile normalized file already exists\n")
    next
  }
  cat(date(), " quantile normalizing ", gses$id[i], "\n")
  meth = fread(in.file, data.table=F, tmpdir=".", header=T)
  rownames(meth) = meth[[1]]
  meth[[1]] = NULL
  meth = as.matrix(meth)
  norm = normalize.quantiles.use.target(meth, target)
  rownames(norm) = rownames(meth)
  colnames(norm) = colnames(meth)
  fwrite(
    data.frame(cg=rownames(norm), norm, check.names=F),
    file=out.file)
} ## 2.5 hours
  
## ----save.methylation.scores-----------------------------------------

save.methylation.scores = function(in.file, studyid, out.file) { 
  models = meffonym.models()
  models = setdiff(models, c("grimage", "grimagev2")) ## require age and sex
  if (file.exists(out.file)) {
    cat(date(), " skipping", studyid, "as scores already exist\n")
    return()
  }
  cat(date(), " generating scores for", studyid, "\n")
  meth = fread(in.file, data.table=F, tmpdir=".", header=T)
  rownames(meth) = meth[[1]]
  meth = as.matrix(meth[,-1])
  scores = mclapply(
    models,
    function(model) {
      tryCatch({
        meffonym.score(meth, model)$score
      }, error=function(e) {
        rep(NA, ncol(meth))
      })
    })
  scores = do.call(cbind, scores)
  rownames(scores) = paste(studyid, colnames(meth), sep="-")
  colnames(scores) = models
  scores = t(scores)
  scores = data.frame(model=rownames(scores), scores, check.names=F)
  fwrite(scores, file=out.file)
}  

scores.dir = file.path(dir$output, "methylation-scores")
dir.create(scores.dir, showWarnings = F)
gses$scores.file = ""
gses$scores.file[gses$is_blood] = file.path(scores.dir, basename(gses$dnam.file[gses$is_blood]))
for (i in which(gses$is_blood)) {
  save.methylation.scores(
    in.file = gses$dnam.file[i],
    studyid = gses$id[i],
    out.file = gses$scores.file[i])
} ## 4 hours

scores.qn.dir = file.path(dir$output, "methylation-scores-qn")
dir.create(scores.qn.dir, showWarnings = F)
gses$scores.qn.file = ""
gses$scores.qn.file[gses$is_blood] = file.path(scores.qn.dir, basename(gses$dnam.file[gses$is_blood]))
for (i in which(gses$is_blood)) {
  save.methylation.scores(
    in.file = gses$dnam.file[i],
    studyid = gses$id[i],
    out.file = gses$scores.qn.file[i])
} ## 4 hours


## ----write.gses --------------------------------------------

# full query results
data.table::fwrite(
  gses, 
  file.path(dir$output, "gses.csv"),
  row.names = F,
  na = "NA")

## datasets derived from blood
data.table::fwrite(
  gses[gses$is_blood,],
  file.path(dir$output, "blood.gses.csv"),
  row.names = F,
  na = "NA")


## ----testing.samples------------------------------------------

testing.dats = c(
"GSE281691_GPL21145",
"GSE280206_GPL21145",
"GSE280206_GPL13534",
"GSE243529_GPL21145",
"GSE183920_GPL21145",
"GSE51032_GPL13534")

testing.samples = sapply(
  gses$samples.file[match(testing.dats, gses$id)],
  fread, data.table=F,
  simplify=F)

names(testing.samples) = testing.dats

testing.samples$GSE281691_GPL21145$cancer = grepl("diseasestate: Case", testing.samples$GSE281691_GPL21145$characteristics_ch1.2)

testing.samples$GSE280206_GPL21145$cancer = grepl("incident_cancer: yes", testing.samples$GSE280206_GPL21145$characteristics_ch1.2)

testing.samples$GSE280206_GPL13534$cancer = grepl("incident_cancer: yes", testing.samples$GSE280206_GPL13534$characteristics_ch1.2)

testing.samples$GSE243529_GPL21145$cancer = grepl("cancer patient", testing.samples$GSE243529_GPL21145$description)

testing.samples$GSE183920_GPL21145$cancer=grepl("ten year-recurrence or death: Yes", testing.samples$GSE183920$characteristics_ch1.3)

testing.samples$GSE51032_GPL13534$cancer = grepl(": C[0-9]+$", testing.samples$GSE51032_GPL13534$characteristics_ch1.2) | grepl(": C[0-9]+$", testing.samples$GSE51032_GPL13534$characteristics_ch1.3) 

testing.samples = do.call(rbind, lapply(testing.dats, function(id) {
  matching = fread(gses$matching.file[match(id,gses$id)], data.table=F)
  matching = matching[,c("geo_accession","meth_id")]
  samples = testing.samples[[id]][,c("geo_accession", "cancer")]
  matching$cancer = samples$cancer[match(matching$geo_accession, samples$geo_accession)]
  cbind(id=id, matching)
}))

with(testing.samples, table(id,cancer))
##                     cancer
## id                   FALSE TRUE
##   GSE183920_GPL21145   230  373
##   GSE243529_GPL21145   268  256
##   GSE280206_GPL13534   357   53
##   GSE280206_GPL21145   208   30
##   GSE281691_GPL21145   260  221
##   GSE51032_GPL13534    424  421

## datasets derived from blood
data.table::fwrite(
  testing.samples, 
  file.path(dir$output, "testing.samples.csv"),
  row.names = F,
  na = "NA")

