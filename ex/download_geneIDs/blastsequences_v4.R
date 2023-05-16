
# ALL, All Fields, All terms from all searchable fields
# UID, UID, Unique number assigned to each sequence
# FILT, Filter, Limits the records
# WORD, Text Word, Free text associated with record
# TITL, Title, Words in definition line
# KYWD, Keyword, Nonstandardized terms provided by submitter
# AUTH, Author, Author(s) of publication
# JOUR, Journal, Journal abbreviation of publication
# VOL, Volume, Volume number of publication
# ISS, Issue, Issue number of publication
# PAGE, Page Number, Page number(s) of publication
# ORGN, Organism, Scientific and common names of organism, and all higher levels of taxonomy
# ACCN, Accession, Accession number of sequence
# PACC, Primary Accession, Does not include retired secondary accessions
# GENE, Gene Name, Name of gene associated with sequence
# PROT, Protein Name, Name of protein associated with sequence
# ECNO, EC/RN Number, EC number for enzyme or CAS registry number
# PDAT, Publication Date, Date sequence added to GenBank
# MDAT, Modification Date, Date of last update
# SUBS, Substance Name, CAS chemical name or MEDLINE Substance Name
# PROP, Properties, Classification by source qualifiers and molecule type
# SQID, SeqID String, String identifier for sequence
# GPRJ, Genome Project, Genome Project
# SLEN, Sequence Length, Length of sequence
# FKEY, Feature key, Feature annotated on sequence
# RTYP, Replicon type, Replicon type
# RNAM, Replicon name, Replicon name
# ORGL, Organelle, Organelle

# E.g. mitochondrion[ORGL]
# E.g. srcdb refseq[Properties]
#      &srcdb%20refseq%5BProperties%5D

# E.g. srcdb refseq[PROP]
# Other databases:
# nr
# refseq
# swissprot


# http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.html

# LIMITS
# (txid9608 [ORGN] OR txid9688 [ORGN]) NOT (txid9615 [ORGN] OR txid9606 [ORGN] OR (srcdb refseq model[prop] AND biomol rna[prop]) OR environmental samples[organism] OR metagenomes[orgn] OR txid32644[orgn])

# CODES
# ASCII Encoding Reference
# http://www.w3schools.com/tags/ref_urlencode.asp

# SPECIAL CHARACTERS IN URLs
# space  = %20
# [ = %5B
# ] = %5D
# \ = %5C
# ( = %28
# ) = %29

# OTHERS
# ! = %21
# ! = %21
# " = %22
# # = %23
# $ = %24
# % = %25
# & = %26

convert_txt_to_NCBI <- function(txt)
	{
	defaults='
	txt = "Canis[orgn] AND (srcdb RefSeqGene[Properties] OR srcdb RefSeq[Properties])"
	convert_txt_to_NCBI(txt)
	'
	
	# Spaces
	txt2 = txt
	txt2 = gsub(pattern=" ", replace="%20", x=txt2)
	# [
	txt2 = gsub(pattern="\\[", replace="%5B", x=txt2)
	# ]
	txt2 = gsub(pattern="\\]", replace="%5D", x=txt2)
	# [
	txt2 = gsub(pattern="\\(", replace="%28", x=txt2)
	# ]
	txt2 = gsub(pattern="\\)", replace="%29", x=txt2)
	
	txt2
	return(txt2)
	}


convert_seqs_to_seqCollapse <- function(x)
	{
	# collapse seq into string
	seqCollapse <- paste(toupper(as.character(x)), collapse = "")
	return(seqCollapse)
	}


# Using the blastSequences function
# 
# The blastSequences function appears to work well in most instances. For slow web 
# connections or for large query sequences, it appears to fail. This seem to be due 
# to either a low number of attempts to retrieve hits from the ncbi server, or the 
# short time the routine waits before it throws an error in R. It seem that if we 
# increase the number of attempts to retrieve hits, generally, we can successfully 
# run our BLAST.
# The below code re-defines the blastSequences function, with an additional argument 
# (attempts), to allow users to control how many times the routine should try to 
# retrieve results from the server.
# 
# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html
#
# See also: http://www.inside-r.org/packages/cran/BoSSA/docs/blast
# 
# function definition
blastSeqKK2 <- function (x, database="nr", hitListSize="10", 
                        filter="L", expect="10", program="blastn",
                        organism=NULL, not_organism=NULL, addl_url=NULL, 
                        txt=txt, justurl0=FALSE, attempts=10, retry_wait=5)
	{
	defaults='
	wd = "/drives/Dropbox/_njm/"
	setwd(wd)
	fasta_fn = "Wayne_etal_1997_SysBio_canid_mtDNA_appendix_wolf_cytB.fasta"
	mySeq = read.fasta(fasta_fn)

	x = convert_seqs_to_seqCollapse(x=mySeq[[1]])
	database="nr"
	hitListSize="10"
	filter="L"
	expect="10"
	program="blastn"
	attempts=10
	organism=NULL
	addl_url=NULL
	' # END defaults
	
	baseUrl = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	# Other, non-blast, stuff
	# http://eutils.ncbi.nlm.nih.gov/blast/Blast.cgi
	
	query = paste("QUERY=", as.character(x), "&DATABASE=", database, 
			 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
			 expect, "&PROGRAM=", program, sep = "")
	
	# Put in organism limitation
	entrez_query = FALSE
	q2 = ""
	if (!is.null(organism) || !is.null(not_organism))
		{
		entrez_query = TRUE
		
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		if (length(c(organism, not_organism)) == 1)
			{
			if (length(organism) == 1)
				{
				q2 = paste0(q2, organism, "[orgn]&")
				}
			if (length(not_organism) == 1)
				{
				q2 = paste0(q2, "NOT ", not_organism, "[orgn]&")
				}
			} else {
			
			if ( (length(organism) > 1) || (length(not_organism) > 1) )
				{
				if (length(organism) > 0)
					{
					orgn_txt = rep("[orgn]", times=length(organism))
					organisms_txt1 = paste(unlist(organism), orgn_txt, sep="")
					organisms_txt1
					organisms_txt2 = paste(organisms_txt1, collapse=" OR ", sep="")
					organisms_txt2
					}
				if (length(not_organism) > 0)
					{
					NOT_txt = rep("NOT ", times=length(not_organism))
					orgn_txt = rep("[orgn]", times=length(not_organism))
					not_organisms_txt1 = paste(NOT_txt, unlist(not_organism), orgn_txt, sep="")
					not_organisms_txt1
					not_organisms_txt2 = paste(not_organisms_txt1, collapse=" ", sep="")
					not_organisms_txt2
					}
				}
			
			txt3 = paste0("(", organisms_txt2, " ", not_organisms_txt2, ")")
			q2 = paste0(q2, txt3, "&")
			}
		} # END if (!is.null(organism) || !is.null(not_organism))



	# Put in additional (manual) limitations on the search
	if (!is.null(addl_url))
		{
		entrez_query = TRUE
		
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		q2 = paste0(q2, addl_url, "&")
		} # END if (!is.null(addl_url))

	# Put in additional (text, in ENTREZ format) limitations on the search
	if (!is.null(txt))
		{
		entrez_query = TRUE
		
		# Convert to URL format...
		#addl_txt = convert_txt_to_NCBI(txt)
		# Add e.g. "human[orgn]" to query
		# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node25.html
		q2 = paste0(q2, txt, "&")
		} # END if (!is.null(addl_url))
	
	if (entrez_query == TRUE)
		{
		q2 = paste0("&ENTREZ_QUERY=", q2)
		query = paste0(query, q2)
		} # END if (entrez_query = TRUE)

	# Convert the URL to HTML format
	query = gsub(pattern="&&", replacement="&", x=query)
	query_orig = query
	query = convert_txt_to_NCBI(query_orig)		
	url0 = sprintf("%s?%s&CMD=Put", baseUrl, query)
	url0_orig = sprintf("%s?%s&CMD=Put", baseUrl, query_orig)
	url0 = gsub(pattern="&&", replacement="&", x=url0)
	url0_orig = gsub(pattern="&&", replacement="&", x=url0_orig)
	
	if (justurl0 == TRUE)
		{
		results = NULL
		results$url0_orig = url0_orig
		results$url0 = url0
		return(results)
		}
	
	# Continue if desired
	results = tempfile()
	
	# Wait 5 seconds
	Sys.sleep(150)
	
	
	require(XML)
	
	# Get the search RID
	post = htmlTreeParse(url0, useInternalNodes = TRUE)
	x = post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
	rid = sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
	
	# Put the search RID into a url for searching/downloading
	# ROTE = Estimated waiting time before the BLAST results will be ready
	rtoe = as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
	url1 = sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
	
	# Wait for the waiting time, plus a bit
	Sys.sleep(rtoe + runif(n=1, min=retry_wait, max=2*retry_wait))
	
	# Try to parse
	result = .tryParseResult(url1, attempts, retry_wait=retry_wait)
	
	results = NULL
	if (justurl0 == FALSE)
		{
		results$result = result
		} # END if (justurl0 == FALSE)

	results$url0_orig = url0_orig
	results$url0 = url0
	results$url1 = url1
	return(results)
	} # END blastSeqKK


# Setup a try function
# (url means searching online)
.tryParseResult <- function(url, attempts=10, retry_wait=5)
	{
	total_pause_time = 0
	
	# Try a number of attempts
	for (i in 1:(attempts+1))
		{
		result = tryCatch(
			{
			xmlTreeParse(url, useInternalNodes=TRUE, error = xmlErrorCumulator(immediate=FALSE))
			}, error=function(err) NULL
			)
	
		# Wait if needed
		if (!is.null(result))
			{
			txt = paste0("\nOnline BLAST search returned after ", i, " tries, totalling ", total_pause_time, " seconds.\nResult XML has summary(result)$numNodes=", summary(result)$numNodes, " nodes, ", summary(result)$nameCounts["Hit"], " hits, and ", summary(result)$nameCounts["Hsp"], " Hsps.\n")
			cat(txt)
			return(result)
			} # END if (!is.null(result))
			
		# Wait longer after each try, just in case it comes through...
		retry_wait_time = i * runif(n=1, min=retry_wait, max=2*retry_wait)
		Sys.sleep(retry_wait_time)
		} # END for (i in 1:(attempts+1))
	
	
	stop(paste("no results after ", attempts, 
		   " attempts; please try again later", sep = ""))
	} # END .tryParseResult <- function(url, attempts)



DNAmult_matches <- function(aln)
{
	defaults='
	aln = blastres$aln[[1]]
	'
	# Get the number of IDENTICAL, NUCLEOTIDE positions
	matches_by_nucleotide = consensusMatrix(aln, as.prob=TRUE)
	matches_by_nucleotide[1:18, 1:10]
	# Just take where the frequency of A, C, G, or T (rows 1-4) equals 1
	nummatches = colSums(matches_by_nucleotide[1:4,] == 1)
	match_TF = nummatches == 1
	num_identical = sum(match_TF)
	
	# Get the number of DISAGREEING, NUCLEOTIDE positions
	nummatches = colSums(matches_by_nucleotide[1:4,] == 0.5)
	match_TF = nummatches == 2
	num_disagree = sum(match_TF)
	
	# Number of resolved positions
	num_resolved = num_identical + num_disagree
	
	# Alignment length
	aln_len = dim(aln)[2]
	
	# percent identical
	pctIdent = round( 100 * (num_identical / num_resolved), digits=1)
	
	match_df = c(pctIdent, num_identical, num_disagree, num_resolved, aln_len)
	match_df
	return(match_df)
	}


make_aln_from_qseq_hseq <- function(qseq, hseq)
	{
	# Make the pairwise alignment objects
	aln = DNAMultipleAlignment(
		c(qseq, hseq), 
		rowmask = as(IRanges(), "NormalIRanges"), 
		colmask = as(IRanges(), "NormalIRanges")
		)
	return(aln)
	} # END make_aln_from_hseq_qseq <- function(qseq, hseq)


# Pretty damn slow!
DNAmult_matches_from_dnadf <- function(dnadf)
	{
	defaults='
	dnadf = blastres$dnadf
	'
	alns = mapply(FUN=make_aln_from_qseq_hseq, qseq=dnadf$qseq, hseq=dnadf$hseq)
	tmp_dfs = sapply(X=alns, FUN=DNAmult_matches)
	names(tmp_dfs) = NULL
	match_dfs = t(tmp_dfs)
	match_dfs = as.data.frame(match_dfs, row.names=FALSE, stringsAsFactors=FALSE)
	names(match_dfs) = c("pctIdent", "num_identical", "num_disagree", "num_resolved", "aln_len")
	head(match_dfs)
	
	return(match_dfs)
	} # END DNAmult_matches_from_dnadf <- function(dnadf)


# Get the gi and the first gb from blastres$dnadf$Hit_id
get_first4_from_split_Hit_id <- function(split_Hit_id)
	{
	gi_gbs = split_Hit_id[1:4]
	return(gi_gbs)
	} # END get_first4_from_split_Hit_id <- function(split_Hit_id)
	
get_gi_gb <- function(Hit_ids)
	{
	defaults='
	Hit_ids = blastres$dnadf$Hit_id
	'
	split_Hit_ids = sapply(X=Hit_ids, FUN=strsplit, "\\|")
	gi_gbs = t(sapply(X=split_Hit_ids, FUN=get_first4_from_split_Hit_id))
	gi_gbs = as.data.frame(gi_gbs, row.names=FALSE, stringsAsFactors=FALSE)
	names(gi_gbs) = c("giTF", "gi", "gbTF", "gb")
	gi_gbs$gi = as.numeric(gi_gbs$gi)
	return(gi_gbs)
	} # END get_gi_gb <- function(Hit_ids)

get_sp_ssp_from_hit_id <- function(hit_id, genus)
	{
	require(gdata) # for trim
	
	# Get the start/end of the first hit
	startpos = regexpr(pattern=genus, text=hit_id, ignore.case=FALSE)
	endpos = startpos + nchar(genus) - 1
	str_length = nchar(hit_id)
	print(str_length)
	subtxt = substr(x=hit_id, start=endpos+1, stop=str_length)
	print(subtxt)
	end_string = gdata::trim(subtxt)
	
	# Error check - if genus is at the end, or genus not found
	if (((endpos+1) > str_length) || (startpos == -1))
		{
		genus = genus
		species = ""
		subspecies = ""
		gsp = matrix(c(genus, species, subspecies), nrow=1)
		gsp = as.data.frame(gsp, row.names=NULL, stringsAsFactors=FALSE)
		names(gsp) = c("genus", "species", "subspecies")
		return(gsp)
		} # END if ((endpos+1) > str_length)
	
	# Extract the species (and perhaps subspecies) after the genus
	words = strsplit(end_string, split=" ")[[1]]
	
	# Another error check
	if (length(words[[1]]) == 0)
		{
		genus = genus
		species = ""
		subspecies = ""
		gsp = matrix(c(genus, species, subspecies), nrow=1)
		gsp = as.data.frame(gsp, stringsAsFactors=FALSE, row.names=NULL)
		names(gsp) = c("genus", "species", "subspecies")
		return(gsp)
		} # END if (length(words[[1]]) == 0)
	
	if (length(words) >= 1)
		{
		species = words[[1]]
		
		# Check for equals
		if (grepl(pattern="=", x=species) == TRUE)
			{
			words2 = strsplit(species, split="=")[[1]]
			species = words2[1]
			}
		
		if (length(words) >= 2)
			{
			subspecies = words[[2]]
			} else {
			subspecies = ""
			} # END if (length(words) >= 2)
		} else {
		species = ""
		subspecies = ""
		} # END if (length(words) >= 2)

	gsp = matrix(c(genus, species, subspecies), nrow=1)
	gsp = as.data.frame(gsp, stringsAsFactors=FALSE, row.names=NULL)
	names(gsp) = c("genus", "species", "subspecies")
	return(gsp)
	} # END get_sp_ssp_from_hit_id

dnadf_enhance <- function(dnadf, genera=NULL, calc_match_percentages=TRUE)
	{
	defaults='
	dnadf=blastres$dnadf
	genera = c("Canis", "Lycaon")
	calc_match_percentages=TRUE
	' # END defaults
	
	# Extract the hit ids
	Hit_ids = dnadf$Hit_id
	gi_gbs = get_gi_gb(Hit_ids)
	
	# Find mentions of "mitochond"
	mt1 = grepl(pattern="mitochond", x=dnadf$Hit_def)
	mt2 = grepl(pattern="mtDNA", x=dnadf$Hit_def)
	mt = (mt1 + mt2) > 0

	# Find mentions of "numt" (mtDNA transfered to the nucleus)
	numt = grepl(pattern="numt", x=dnadf$Hit_def, ignore.case=TRUE)

	# Find mentions of "chloroplast"
	cp1 = grepl(pattern="chloroplast", x=dnadf$Hit_def, ignore.case=TRUE)
	cp2 = grepl(pattern="cpDNA", x=dnadf$Hit_def, ignore.case=TRUE)
	cp3 = grepl(pattern="plastid", x=dnadf$Hit_def, ignore.case=TRUE)
	cp = (cp1 + cp2 + cp3) > 0

	# Find mentions of "complete genome"
	cgenome = grepl(pattern="complete genome", x=dnadf$Hit_def, ignore.case=TRUE)
	# Find mentions of "partial genome"
	pgenome = grepl(pattern="partial genome", x=dnadf$Hit_def, ignore.case=TRUE)
	pCDS = grepl(pattern="partial CDS", x=dnadf$Hit_def, ignore.case=TRUE)
	partial = grepl(pattern="partial", x=dnadf$Hit_def, ignore.case=TRUE)
	partial[pgenome==TRUE] = FALSE
	
	# Percentage matches
	pctIdent = round(100*dnadf$Hsp_identity/dnadf$Hsp_align_len, 1)
	pctSim = round(100*dnadf$Hsp_positive/dnadf$Hsp_align_len, 1)
	
	# Make the new fields into a data.frame
	tmpmat = cbind(pctIdent, pctSim, mt, numt, cp, cgenome, pgenome, pCDS, partial)
	typedf = as.data.frame(tmpmat, row.names=NULL, stringsAsFactors=FALSE)
	names(typedf) = c("pctIdent", "pctSim", "mt", "numt", "cp", "cgenome", "pgenome", "pCDS", "partial")
	
	
	# Calculate genus, species from matching genus names
	gsp_df = get_genera_from_dnadf(dnadf, genera=genera)


	# Calculate match percentages
	matches_df = matrix(data=NA, nrow=nrow(dnadf), ncol=5)
	matches_df = as.data.frame(matches_df, row.names=NULL, stringsAsFactors=FALSE)
	names(matches_df) = c("pctIdent", "num_identical", "num_disagree", "num_resolved", "aln_len")
	
	# Calculate the match percentages etc.
	# Slooooooooow....
	# Pretty damn slow!
	if (calc_match_percentages == TRUE)
		{
		txt = "\n\nCalculating pairwise match percentages...\n\n"
		matches_df = DNAmult_matches_from_dnadf(dnadf)
		} # END if (calc_match_percentages == TRUE)
	
	# Combine
	row.names(dnadf) = NULL
	ncols_dnadf = ncol(dnadf)
	dnadf2 = cbind(dnadf[, 1], gi_gbs, matches_df, gsp_df, typedf, dnadf[,2:ncols_dnadf])
	dnadf2[1:5, 1:38]
	
	names(dnadf2)[1] = names(dnadf)[1]
	dnadf2[1:5, 1:38]
	cls.df(dnadf2)
	
	return(dnadf2)
	} # END dnadf_enhance <- function(dnadf)


count_blastresult <- function(result)
	{
	defaults='
	wd = "/drives/Dropbox/_njm/"
	setwd(wd)
	fasta_fn = "Wayne_etal_1997_SysBio_canid_mtDNA_appendix_wolf_cytB.fasta"
	mySeq = read.fasta(fasta_fn)
	x = convert_seqs_to_seqCollapse(x=mySeq[[1]])
	result = blastSeqKK2(x=x, database="nr", hitListSize="10", 
				filter="L", expect="10", program="blastn",
				attempts=10, organism=NULL, addl_url=NULL)
	' # END defaults
	
	# Extract information from the XML in "result"
	Hit_num = xpathApply(result, "//Hit_num", xmlValue)
	Hsp_num = xpathApply(result, "//Hsp_num", xmlValue)
	
	num_Hsps = length(Hsp_num)
	
	# Summarize
	txt = paste0("\n\nThe XML returned by BLAST has ", length(Hit_num), " hits and ", num_Hsps, " Hsps.\n\n")
	cat(txt)
	
	return(num_Hsps)
	} # END count_blastresult


get_genera_from_dnadf <- function(dnadf, genera=NULL)
	{
	# Pull out genus/species names
	gsp_df = matrix(data=NA, nrow=nrow(dnadf), ncol=5)
	gsp_df = as.data.frame(gsp_df, row.names=NULL, stringsAsFactors=FALSE)
	names(gsp_df) = c("genus", "species", "subspecies", "gnsp", "gnspssp")
	
	if (!is.null(genera))
		{
		for (g in 1:length(genera))
			{
			genus = genera[g]
			TF = grepl(pattern=genus, x=dnadf$Hit_def, ignore.case=FALSE)
			
			tmpdf = t(sapply(X=dnadf$Hit_def[TF], FUN=get_sp_ssp_from_hit_id, genus=genus))
			tmpdf = as.data.frame(tmpdf, stringsAsFactors=FALSE, row.names=NULL)
			for (colnum in 1:3)
				{
				gsp_df[TF,colnum] = unlist(tmpdf[1:nrow(tmpdf),colnum])
				} # END for (colnum in 1:3)
			} # END for (g in 1:length(genera))
		} # END if (!is.null(genera))
	
	# Add genus_species and genus_species_ssp
	# (i.e., gnsp and gnspssp)
	gnsp = rep(NA, times=nrow(gsp_df))
	notNA = is.na(gsp_df$genus) == FALSE
	gnsp[notNA] = paste(gsp_df$genus[notNA], gsp_df$species[notNA], sep="_")

	gnspssp = rep(NA, times=nrow(gsp_df))
	gnspssp[notNA] = paste(gsp_df$genus[notNA], gsp_df$species[notNA], gsp_df$subspecies[notNA], sep="_")
	
	gsp_df[,4] = gnsp
	gsp_df[,5] = gnspssp
	
	gsp_df = as.data.frame(gsp_df, row.names=NULL, stringsAsFactors=FALSE)
	names(gsp_df) = c("genus", "species", "subspecies", "gnsp", "gnspssp")
	
	# Convert to character, not factor
	gsp_df$genus = as.character(gsp_df$genus)
	gsp_df$species = as.character(gsp_df$species)
	gsp_df$subspecies = as.character(gsp_df$subspecies)
	gsp_df$gnsp = as.character(gsp_df$gnsp)
	gsp_df$gnspssp = as.character(gsp_df$gnspssp)
	
	return(gsp_df)
	} # END get_genera_from_dnadf <- function(dnadf, genera=NULL


process_blastresult <- function(results, ivals=NULL, genera=NULL)
	{
	defaults='
	wd = "/drives/Dropbox/_njm/"
	setwd(wd)
	fasta_fn = "Wayne_etal_1997_SysBio_canid_mtDNA_appendix_wolf_cytB.fasta"
	mySeq = read.fasta(fasta_fn)
	x = convert_seqs_to_seqCollapse(x=mySeq[[1]])
	result = blastSeqKK2(x=x, database="nr", hitListSize="10", 
				filter="L", expect="10", program="blastn",
				attempts=10, organism=NULL, addl_url=NULL)
	' # END defaults

	result = results$result
	url0 = results$url0 
	url1  = results$url1
		
	# Extract information from the XML in "result"
	Hit_num = xpathApply(result, "//Hit_num", xmlValue)
	Hit_id = xpathApply(result, "//Hit_id", xmlValue)
	Hit_def = xpathApply(result, "//Hit_def", xmlValue)
	Hit_accession = xpathApply(result, "//Hit_accession", xmlValue)
	Hit_len = xpathApply(result, "//Hit_len", xmlValue)
	Hsp_num = xpathApply(result, "//Hsp_num", xmlValue)
	Hsp_bit_score = xpathApply(result, "//Hsp_bit-score", xmlValue)
	Hsp_score = xpathApply(result, "//Hsp_score", xmlValue)
	Hsp_evalue = xpathApply(result, "//Hsp_evalue", xmlValue)
	Hsp_query_from = xpathApply(result, "//Hsp_query-from", xmlValue)
	Hsp_hit_to = xpathApply(result, "//Hsp_hit-to", xmlValue)
	Hsp_query_frame = xpathApply(result, "//Hsp_query-frame", xmlValue)
	Hsp_hit_frame = xpathApply(result, "//Hsp_hit-frame", xmlValue)
	Hsp_identity = xpathApply(result, "//Hsp_identity", xmlValue)
	Hsp_positive = xpathApply(result, "//Hsp_positive", xmlValue)
	Hsp_gaps = xpathApply(result, "//Hsp_gaps", xmlValue)
	Hsp_align_len = xpathApply(result, "//Hsp_align-len", xmlValue)
	qseq = xpathApply(result, "//Hsp_qseq", xmlValue)
	hseq = xpathApply(result, "//Hsp_hseq", xmlValue)
	
	
	# Summarize
	txt = paste0("\n\nThe XML returned by BLAST has ", length(Hit_num), " hits and ", length(hseq), " Hsps.\n\n")
	cat(txt)
	
	require(Biostrings)
	aln = list()
	dnadf = NULL
	
	# Make alignment for each sequence
	current_Hit_num = 1
	previous_Hsp_num = 0
	previous_hit_done = FALSE
	count = 1
	cat("\nProcessing Hsp #: ")
	for (i in seq_len(length(qseq)))
		{
		if (count < 100)
			{
			count = count + 1
			} else {
			cat(i, " ")
			count = 0
			} # END if (count == 100)
		
		current_Hsp_num = Hsp_num[[i]]
		if (previous_Hsp_num >= current_Hsp_num)
			{
			previous_hit_done = TRUE
			}
		if (previous_hit_done == TRUE)
			{
			current_Hit_num = current_Hit_num + 1
			previous_hit_done = FALSE
			}
		previous_Hsp_num = Hsp_num[[i]]
		
		# Make the pairwise alignment objects
		aln[i] = DNAMultipleAlignment(
			c(hseq[[i]], qseq[[i]]), 
			rowmask = as(IRanges(), "NormalIRanges"), 
			colmask = as(IRanges(), "NormalIRanges")
			)
		
		n = current_Hit_num
		tmprow = c(Hit_num[[n]], Hit_id[[n]], Hit_def[[n]], Hit_accession[[n]], Hit_len[[n]], Hsp_num[[i]], Hsp_bit_score[[i]], Hsp_score[[i]], Hsp_evalue[[i]], Hsp_query_from[[i]], Hsp_hit_to[[i]], Hsp_query_frame[[i]], Hsp_hit_frame[[i]], Hsp_identity[[i]], Hsp_positive[[i]], Hsp_gaps[[i]], Hsp_align_len[[i]], qseq[[i]], hseq[[i]])
		dnadf = rbind(dnadf, tmprow)
		} # END for (i in seq_len(length(qseq)))
	
	dnadf = as.data.frame(dnadf, row.names=NULL, stringsAsFactors=FALSE)
	names(dnadf) = c("Hit_num", "Hit_id", "Hit_def", "Hit_accession", "Hit_len", "Hsp_num", "Hsp_bit_score", "Hsp_score", "Hsp_evalue", "Hsp_query_from", "Hsp_hit_to", "Hsp_query_frame", "Hsp_hit_frame", "Hsp_identity", "Hsp_positive", "Hsp_gaps", "Hsp_align_len", "qseq", "hseq")
	
	
	
	# Make appropriate columns numeric
	# This is hella-slow
	#dnadf = dfnums_to_numeric(dnadf)
	dnadf$Hit_num = as.numeric(dnadf$Hit_num)
	dnadf$Hit_len = as.numeric(dnadf$Hit_len)
	dnadf$Hsp_num = as.numeric(dnadf$Hsp_num)
	dnadf$Hsp_bit_score = as.numeric(dnadf$Hsp_bit_score)
	dnadf$Hsp_score = as.numeric(dnadf$Hsp_score)
	dnadf$Hsp_evalue = as.numeric(dnadf$Hsp_evalue)
	dnadf$Hsp_query_from = as.numeric(dnadf$Hsp_query_from)
	dnadf$Hsp_hit_to = as.numeric(dnadf$Hsp_hit_to)
	dnadf$Hsp_query_frame = as.numeric(dnadf$Hsp_query_frame)
	dnadf$Hsp_hit_frame = as.numeric(dnadf$Hsp_hit_frame)
	dnadf$Hsp_identity = as.numeric(dnadf$Hsp_identity)
	dnadf$Hsp_positive = as.numeric(dnadf$Hsp_positive)
	dnadf$Hsp_gaps = as.numeric(dnadf$Hsp_gaps)
	dnadf$Hsp_align_len = as.numeric(dnadf$Hsp_align_len)
	
	cls.df(dnadf)
	dnadf[1:5, 1:17]

	
	dnadf2 = dnadf_enhance(dnadf=dnadf, genera=genera)
	dim(dnadf2)
	#dnadf2[1:5, 1:25]

	
	# Return BLAST results
	blastres = NULL
	blastres$aln = aln
	blastres$dnadf = dnadf2
	blastres$url0 = url0
	blastres$url1 = url1
	
	# Extract
	extract = '
	dnadf = blastres$dnadf
	'
	
	blastres
	return(blastres)
	} # END process_blastresult <- function(result)






# Split cells of a vector / spreadsheet on separation (eg comma)
split_cells_on_sep <- function(vec, split=", ")
	{
	veclist = lapply(X=vec, FUN=split_cell_on_sep, split=split)
	newvec = unlist(veclist)
	return(newvec)
	} # END split_cells_on_sep <- function(vec, sep=",")

split_cell_on_sep <- function(cell, split=",")
	{
	words = strsplit(cell, split=split)[[1]]
	return(words)
	} # END split_cells_on_sep <- function(vec, sep=",")



# Convert Excel column to a vector of items
excel_column_to_vector <- function(vec, ifcommas="takeall")
	{
	defaults='
	vec=xls$CytB_col1
	'
	
	# Check the column for cells with commas; if so, take first or last
	if (ifcommas != "takeall")
		{
		for (v in 1:length(vec))
			{
			TF = grepl(pattern=",", x=vec[v])
			if (TF == TRUE)
				{
				print(vec[v])
				words = strsplit(x=vec[v], split=",")[[1]]
				
				if (ifcommas == "first")
					{
					vec[v] = words[1]
					} # END if (ifcommas == "first")

				if (ifcommas == "last")
					{
					vec[v] = words[length(words)]
					} # END if (ifcommas == "last")
				} # END if (TF == TRUE)
			} # END for (v in 1:length(vec))
		} # END if (ifcommas != "takeall")
	
	
	# Remove blanks
	TF = isblank_TF(vec)==FALSE
	newvec = vec[TF]
	TF = newvec != ".--"
	newvec = newvec[TF]
	
	# Split on ", " or ","
	newvec = split_cells_on_sep(newvec, split=", ")
	newvec = split_cells_on_sep(newvec, split=",")
	
	return(newvec)
	}






#######################################################
# prflag
#######################################################
#' Utility function to conditionally print intermediate results
#'
#' Just a handy shortcut function, allowing other functions to optionally 
#' print, depending on the value of \code{printflag}.
#' 
#' @param x What to print.
#' @param printflag If TRUE, do the printing
#' @return nothing
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
prflag <- function(x, printflag=TRUE)
	{
	# A standard function to print (or not) certain variables,
	#   based on a master printflag
	# This avoids having to comment in/out various code chunks
	#   while debugging.
	if (printflag == TRUE)
		{
		# CAT instead of PRINT if it's a string or numeric
		if (is.character(x))
			{
			cat(x, "\n", sep="")
			}
		if (is.numeric(x))
			{
			cat(x, "\n", sep="")
			} else {
			print(x)
			}
		}
	else
		{
		pass="BLAH"
		}
	}



get_numlines <- function(fn)
	{
	# Get the number of lines in a file
	
	# Check if the file exists
	TF = file.exists(fn)
	if (TF == FALSE)
		{
		txt = paste0("\nWARNING in get_numlines(): file '", fn, "' does not exist in directory\n'", getwd(), "'. Returning linecount=0.")
		cat(txt)
		linecount = 0
		return(linecount)
		}
	
	# intern = TRUE reports the result to R
	cmdstr = paste("wc -l ", fn, sep="")
	wc_return = system(cmdstr, intern=TRUE) 
	
	# split on whitespace (spaces and tabs)
	wc_list = strsplit_whitespace(wc_return)
	
	linecount = as.numeric(wc_list[1])
	return(linecount)
	} # END get_numlines <- function(fn)


#######################################################
# strsplit_whitespace
#######################################################
#' Split strings on whitespace
#' 
#' This function splits strings on whitespace (spaces and tabs), so you don't have
#' to remember the \code{regexp}/\code{grep} format codes.
#' 
#' @param tmpline A string containing text.
#' @return \code{list_of_strs} 
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' tmpline = "Hello world see	my	tabs."
#' strsplit_whitespace(tmpline)
#' 
strsplit_whitespace <- function(tmpline)
	{
	# split on 1 or more whitespaces
	temp = strsplit(tmpline, "[ \t]+")
	
	# get the list
	list_of_strs = temp[[1]]
	
	# remove any leading/trailing ""
	list_of_strs = list_of_strs[list_of_strs != ""]
	
	return(list_of_strs)
	}


#######################################################
# moref
#######################################################
#' print to screen the header of a file
#' 
#' This does the rough equivalent of the \code{UNIX} function \code{more}, but within R.
#' 
#' @param fn A filename.
#' @param printnotcat If \code{TRUE}, use \code{\link[base]{print}} instead of \code{\link[base]{cat}}. Default \code{FALSE}.
#' @return Nothing returned.
#' @export
#' @seealso \code{\link[base]{scan}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
moref <- function(fn, printnotcat = FALSE)
	{
	lines = scan(file=fn, what="character", sep="\n")
	
	if (printnotcat == TRUE)
		{
		for (i in 1:length(lines))
			{
			print(lines[i])
			}
		}
	else
		{
		for (i in 1:length(lines))
			{
			cat(paste(lines[i], "\n", sep=""))
			}
		}
	}





# Convert blanks etc. to 0
convert_blanks_NAs_to_0 <- function(tipdates, newval=0)
	{
	defaults='
	newval=0
	'
	
	TF = tipdates == ""
	tipdates[TF] = newval

	TF = tipdates == " "
	tipdates[TF] = newval

	TF = tipdates == "\t"
	tipdates[TF] = newval

	TF = is.na(tipdates)
	tipdates[TF] = newval

	TF = is.nan(tipdates)
	tipdates[TF] = newval
	
	return(tipdates)
	}


isblank_TF <- function(items)
	{
	TF1 = items == ""
	TF2 = items == " "
	TF3 = items == "\t"
	TF4 = is.na(items)
	TF5 = is.nan(items)
	
	# Keep only the items where none of the above occur
	TFall = TF1 + TF2 + TF3 + TF4 + TF5
	
	# Correct for NA, NaNs
	TFall[is.na(TFall)] = 1
	TFall[is.nan(TFall)] = 1

	blank_TF = TFall > 0
	return(blank_TF)
	}


# Remove blanks etc. from a list
remove_blanks_NAs_etc <- function(items)
	{
	TF1 = items == ""
	TF2 = items == " "
	TF3 = items == "\t"
	TF4 = is.na(items)
	TF5 = is.nan(items)
	
	# Keep only the items where none of the above occur
	TF = TF1 + TF2 + TF3 + TF4 + TF5
	items = items[TF == 0]
	
	return(items)
	}

get_firstword <- function(OTU, split="_")
	{
	firstword = strsplit(x=OTU, split=split)[[1]][1]
	}

