#######################################################
# Processing CD-Batch
#######################################################

#######################################################
# Search the concise hits of a NCBI CD (Conserved Domain) batch search
# 
# Identify input proteins that don't match any domains of interest
#######################################################

# 

# NCBI Batch Web CD-Search:
# https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
#
# Help page on the same:
# https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#SearchMethodBatchQuery

# 1. Put in a fasta file of proteins (1000 max)
#
# 2. For the results Download, click:
#    * Domain Hits, Data Mode: Concise, click "Superfamily Only"
#		 * Align details: BLAST Text
#    * Click Download
#
# 3. (Download Full versions etc. as you like for visualization)


read_concise_CDbatch_results <- function(search_results_fn)
	{
	#tdf = read.table(search_results_fn, header=TRUE, skip=7, sep="\t")

	# Edits
	lines = readLines(search_results_fn)
	for (i in 1:length(lines))
		{
		lines[i] = gsub(pattern="Query", replacement="Num\tQuery", x=lines[i])
		lines[i] = gsub(pattern="Q\\#", replacement="", x=lines[i])
		lines[i] = gsub(pattern=" - >", replacement="\t", x=lines[i])
		}

	newfn = gsub(pattern=".txt", replacement="_refmt.txt", x=search_results_fn)
	writeLines(text=lines, con=newfn)

	#strsplit(lines[8], split="\t")

	tdf = read.table(newfn, header=TRUE, skip=7, sep="\t")
	return(tdf)
	} # END read_concise_CDbatch_results <- function(search_results_fn)


multidomain_prots_to_single_line <- function(tdf, rm.superfamily=TRUE)
	{
	# Assemble a one-protein-per-line file
	tmprows = NULL	# starter
	numdomains = NULL		# starter
	i = 1						# starter
	j = 0						# starter
	ivals = NULL		# starter
	while (i <= nrow(tdf))
		{
		tmprow = tdf[i,]
		TF = tdf$Query %in% tmprow$Query
		if (sum(TF) > 1)
			{
			domains = tdf$Short.name[TF]
			domains_txt = paste(domains, collapse="; ", sep="")
			tmprow$Short.name = domains_txt
			#tmprows[[(j=j+1)]] = tmprow
			tmprows = rbind(tmprows, tmprow)
			if (length(sum(TF)) > 1)
				{
				stop("hey")
				}
			numdomains = c(numdomains, sum(TF))
			ivals = c(ivals, i)
			i = i + (sum(TF) - 0)
			next()
			} else {
		#	tmprows[[(j=j+1)]] = tmprow
			tmprows = rbind(tmprows, tmprow)
			numdomains = c(numdomains, 1)
			ivals = c(ivals, i)
			i = i+1
			}
		} # END while

	tmprows = cbind(tmprows, numdomains)

	domain_hits_per_protein_df = as.data.frame(tmprows, stringsAsFactors=FALSE)
	
	if (rm.superfamily == TRUE)
		{
		domain_hits_per_protein_df$Short.name = gsub(pattern=" superfamily", replacement="", x=domain_hits_per_protein_df$Short.name)
		domain_hits_per_protein_df$Short.name = gdata::trim(domain_hits_per_protein_df$Short.name)
		}
	
	domain_hits_per_protein_df
	} # END multidomain_prots_to_single_line <- function(tdf)


get_domainsTF_in_CDbatch_results <- function(domain_hits_per_protein_df, domains_of_interest = c("MotA","TolQ","ExbB"))
	{
	junk='
	library(gdata)				# for trim
	library(BioGeoBEARS)	# for sourceall
	sourceall("/GitHub/bioinfRhints/Rsrc/")

	wd = "/GitHub/bioinfRhints/flag/589_MotA_rename/CD-search/"
	setwd(wd)
	search_results_fn = "589seqs_hitdata_concise.txt"
	tdf = read_concise_CDbatch_results(search_results_fn)
	domain_hits_per_protein_df = multidomain_prots_to_single_line(tdf)
	domains_of_interest = c("MotA","TolQ","ExbB")
	TF = get_domainsTF_in_CDbatch_results(domain_hits_per_protein_df, domains_of_interest=domains_of_interest)
	proteins_not_matching_df = domain_hits_per_protein_df[TF==FALSE,]
	'
	
	# Setup
	TF = rep(FALSE, times=nrow(domain_hits_per_protein_df))
	for (i in 1:length(domains_of_interest))
		{
		tmpTF = grepl(pattern=domains_of_interest[i], x=domain_hits_per_protein_df$Short.name, ignore.case=FALSE)
		TF[tmpTF] = TRUE
		}
	
	return(TF)	
	} # END get_domainsTF_in_CDbatch_results <- function(domain_hits_per_protein_df, domains_of_interest = c("MotA","TolQ","ExbB"))



get_first_words <- function(tmpstrs)
	{
	firstwords = NULL
	for (i in 1:length(tmpstrs))
		{
		words = strsplit(tmpstrs[i], split=" +")[[1]]
		firstwords = c(firstwords, words[1])
		}
	return(firstwords)
	}

