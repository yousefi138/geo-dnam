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
gses <- as.data.frame(
            data.table::fread(file.path(dir$output, "gses.csv"), 
                            na.strings = "NA"))

blood.gses <- as.data.frame(
            data.table::fread(file.path(dir$output, "blood.gses.csv"), 
                            na.strings = "NA"))

## ----n -------------------------------------------------------------
sum(blood.gses$samples)

## ----fig.01 -------------------------------------------------------------
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00" )

ggplot(blood.gses, aes(samples)) + 
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

