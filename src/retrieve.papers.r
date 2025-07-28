retrieve.papers <- function(pmids, mc.cores = 4) {
	require(easyPubMed)

	pmids <- na.omit(as.integer(gsub(" ", "", pmids)))

        query <- paste0(paste(pmids, 
                    collapse=" "), "[uid]")
        object <- get_pubmed_ids(query)
		papers <- fetch_pubmed_data(object) 

		require(XML)
		papers <- articles_to_list(papers)

		papers <- lapply(1:length(papers), function(i) {
            cat(date(), "retrieving paper", i, "\n", file=stderr())
            paper <- papers[[i]]
            df <- article_to_df(paper, max_chars=-1, getKeywords=T)
            if (is.null(df)) {
                warning(paste(pmids[[i]], "could not be accessed"))
                return(NULL)
            } 
            df$lastname[1] <- paste(df$lastname, collapse=";")
            df$firstname[1] <- paste(df$firstname, collapse=";")
            df[1,]
        })

        do.call(rbind, papers)
} 
