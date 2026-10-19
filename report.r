## ----globals -------------------------------------------------------------
dir <- paths

kable_my_defaults <- function(x, ...){
	kable(x) %>%
  	kable_styling( 
  		bootstrap_options = c("striped", "hover"), 
  		full_width = F,
  		position = "left",
  		fixed_thead = T,
  		...)
}

## ----load.data -------------------------------------------------------------
gses <- fread(
  file.path(dir$output, "gses.csv"),
  data.table=F,
  na.strings = "NA")

blood.gses <- fread(
  file.path(dir$output, "blood.gses.csv"),
  data.table=F,
  na.strings = "NA")

testing.samples = fread(
  file.path(dir$output, "testing.samples.csv"),
  data.table=F, na.strings="NA")

## ----n -------------------------------------------------------------
sum(blood.gses$num.meth)

## ----fig.01 -------------------------------------------------------------
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00" )

ggplot(blood.gses, aes(num.meth)) + 
	geom_histogram(binwidth = 20, color=cbPalette[4], 
		fill=cbPalette[4]) +
	geom_vline(aes(xintercept=100), 
		color=cbPalette[3], linetype="dashed")

## ----table.01 -------------------------------------------------------------
blood.gses %>%
	group_by(platform) %>%
	summarise(n = n()) %>%
	mutate(percent = round(100 * (n / sum(n)), 1)) %>%
	rbind(., data.frame(platform = "Total", 
						n = t(colSums(.[2])),
						percent = t(colSums(.[3])))) %>%
	kable()

## ----table.cancer -------------------------------------------------------------

ids = unique(testing.samples$id)

descriptions = c(
  GSE281691_GPL21145="hepatocellular carcinoma vs benign liver disease controls",
  GSE280206_GPL21145="any cancer 15-years after stroke", 
  GSE280206_GPL13534="any cancer 15-years after stroke",
  GSE243529_GPL21145="breast cancer cases vs controls",
  GSE183920_GPL21145="10-year non-muscle invasive bladder cancer recurrence or death",
  GSE51032_GPL13534="breast or colorectal cancer")

kable(
  data.frame(
    dataset=ids,
    n=sapply(ids, function(id) sum(testing.samples$id==id)),
    n.bad=sapply(ids, function(id) sum(testing.samples$id==id & testing.samples$cancer)),
    design=descriptions[ids]),
  row.names=F)

