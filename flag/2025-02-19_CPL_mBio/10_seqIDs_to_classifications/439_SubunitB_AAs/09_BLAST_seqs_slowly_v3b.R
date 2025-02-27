#######################################################
# Blast a bunch of sequences in a FASTA file
# (against genbank, slowly)
#
# Separately launch the searches, and retrieve the results
#
# Save everything in a table you are updating and saving
#######################################################
library(openxlsx)
library(ape)
library(phytools)
library(BioGeoBEARS)
library(varhandle) # for check.numeric
library(rentrez) # Fetching full records: entrez_fetch()
library(XML)
library(Biostrings) # For parsing BLAST searches downloaded

sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

# Search on the 133aas_wo_UniProt.fasta file, s1 = search attempt #1
wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/439_SubunitB_AAs/"
setwd(wd)

#######################################################
# SIMPLE: Load and BLAST the FASTA file
#######################################################
fasta_fn = "439_SubunitB_AAs_raw.fasta"
type = "AA"
# This blanks the file
res = build_blastsearch_table_from_fasta(fasta_fn, type=type, overwrite=TRUE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn

# Read in, write out the blast search results table
search_df = read_search_df(fasta_blast_table_fn)
matdf = check_num_still_running(fasta_blast_table_fn); matdf
write_search_df(search_df, fasta_blast_table_fn)

# Construct the BLAST searches (ignoring taxonomic limits,
# apparently these do not speed up the search, but might
# be useful sometimes)
taxids = NULL
database = "nr" 
hitListSize="20"
filter=NULL  # "L" for low-complexity filter
expect="10"
program="blastp"
organism=taxids
not_organism=NULL
addl_url=NULL
txt=NULL
baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
assume_txid=TRUE
res = construct_blasturls_for_table(fasta_blast_table_fn, database=database, hitListSize=hitListSize, filter=filter, expect=expect, program="blastp",
organism=NULL, not_organism=NULL, addl_url=NULL, txt=NULL, baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", assume_txid=TRUE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn

# While loop to gradually fill out table
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

search_df = read_search_df(fasta_blast_table_fn)


search_df[search_df$eTF == TRUE,]

# Re-do search
matdf = check_num_still_running(fasta_blast_table_fn)
matdf
search_df$sTF = FALSE
search_df$pTF = FALSE
search_df$eTF = FALSE
write_search_df(search_df, fasta_blast_table_fn)
matdf = check_num_still_running(fasta_blast_table_fn)
matdf
matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

search_df = read_search_df(fasta_blast_table_fn)
search_df


# SP9SBRUP016 -- example of where an XML was the "wait" result
search_df[2,]


matdf = check_num_still_running(fasta_blast_table_fn)
matdf


# check a specific running RID: SNJKH5BY013
search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)

matchnum = match(x="SNJKH5BY013", table=search_df$rid)
matchnum
search_df[matchnum,]

search_df$sTF[matchnum] = FALSE
search_df$eTF[matchnum] = FALSE

write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf



# Re-do all un-saved ones
search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
search_df$sTF[search_df$eTF == TRUE] = FALSE
search_df$eTF[search_df$eTF == TRUE] = FALSE

write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf

fasta_fn = "44-88aa_of_133_wo_UniProt.fasta"
fasta_blast_table_fn = "44-88aa_wo_UniProt_blast_table_v1.txt"
type = "AA"
res = build_blastsearch_table_from_fasta(fasta_fn, type=type, overwrite=FALSE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

# Shorten for test
write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


taxids = NULL

database = "nr" 
hitListSize="20"
filter=NULL  # "L"
expect="10"
program="blastp"
organism=taxids
not_organism=NULL
addl_url=NULL
txt=NULL
baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
assume_txid=TRUE
res = construct_blasturls_for_table(fasta_blast_table_fn, database=database, hitListSize=hitListSize, filter=filter, expect=expect, program="blastp",
organism=NULL, not_organism=NULL, addl_url=NULL, txt=NULL, baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", assume_txid=TRUE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn

# While loop to gradually fill out table
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf




head(search_df[,1:7])

res = get_some_RIDs_for_table(fasta_blast_table_fn, num_to_get=5, wait_between=runif(1,1.1,2.5), baseUrl="default")
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn

res = launch_some_blast_searches(fasta_blast_table_fn, num_to_launch=5, wait_between=rnorm(n=1, mean=3.8, sd=0.9))
search_df = res$search_df
head(search_df[,1:7])

Sys.sleep(200)

res = download_finished_blast_searches(fasta_blast_table_fn, num_to_get=1, wait_between=1.0)
search_df = res$search_df
head(search_df[,1:7])



# While loop to gradually fill out table
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf




#######################################################
# Check a run
#######################################################
search_df$sTF[c(1:5)] = FALSE
search_df$eTF[c(1:5)] = FALSE
write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



#######################################################
# Revise for a re-run
#######################################################
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")
head(search_df[,1:7])
search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
search_df$sTF[c(3,5)] = FALSE
search_df$eTF[c(3,5)] = FALSE
head(search_df[,1:7])
write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
res = download_finished_blast_searches(fasta_blast_table_fn, num_to_get=5, wait_between=1.0)
search_df = res$search_df
head(search_df[,1:7])

system(search_df$terminal_cmd[3])
system(search_df$terminal_cmd[5])









#######################################################
# Load the FASTA file
#######################################################
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

fasta_fn = "133aas_wo_UniProt.fasta"
type = "AA"
search_df = build_blastsearch_table_from_fasta(fasta_fn, type=type)
head(search_df)
class(search_df)
dim(search_df)



i = 1

cat("\n(Starting search #", i, "/", length(seqs_to_search), ")\n")
Sys.sleep(time=1) # pause for ESC

searchseq = search_df$seqstring[i]

# Database options:
# UniProtKB/Swiss-Prot (swissprot)
# Non-redundant protein sequences (nr)
# Reference proteins (refseq_protein)
# Reference Select proteins (refseq_select)
# Protein Data Bank proteins (pdb)
# 
# See also: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp

# Construct the BLAST URL
taxids = NULL
database = "nr" 
results = construct_blasturl(searchseq=searchseq, database="swissprot", hitListSize="25", 
												filter=NULL, expect="1", program="blastp",
												organism=taxids, not_organism=NULL, addl_url=NULL, 
												txt=NULL, baseUrl=baseUrl)
search_df$cTF[i] = TRUE

# get the RID and waiting time
results = get_blast_RID_waiting_time(results=results, RID_url=NULL, baseUrl=baseUrl)

search_df$RID_url_orig[i] = results$RID_url_orig
search_df$RID_url[i] = results$RID_url
search_df$search_url[i] = results$search_url
search_df$terminal_cmd[i] = results$terminal_cmd
search_df$RIDrequest_fn[i] = results$RIDrequest_fn
search_df$rid[i] = results$rid
search_df$rtoe[i] = results$rtoe

search_df$rTF[i] = TRUE


# Launch the web search; will have to be retrieved later
blastres_multiline = NULL  # BLAST results, may well have >1 per searchseq

search_df$ltime[i] = Sys.time()
start_BLAST_websearch(results)
search_df$lTF[i] = TRUE

# Retrieve results
search_df$search_url[i]

tryParse_results = download_BLAST_xml(results=NULL, search_url=search_df$search_url[i], attempts=1, retry_wait=1)

if ((length(tryParse_results) == 1) && (grepl(pattern="ERROR in download_BLAST_xml", x=tryParse_results) == FALSE))
	{
	#print(tryParse_results$result)
	search_df$sTF[i] = FALSE
	} else {
	search_df$sTF[i] = TRUE
	search_df$stime[i] = Sys.time()
	}


if (search_df$sTF[i] == TRUE)
	{
	try_result = try(process_blastresult(tryParse_results=tryParse_results))

	if (class(try_result) == "try-error")
		{
		search_df$pTF[i] = FALSE
		
		warning_txt = paste0("WARNING: process_blastresult() gave a try-error on i=", i, ", abbr=", abbr[i], ". Inserting a blank line of NAs in blastres$seqdf.")
		
		cat("\n\n")
		warning(warning_txt)
		cat("\n\n")
		
		headers = c("Hit_num", "downID", "taxon", "gi1", "gi2", "gi3", "gi4", "pctIdent", "num_identical", "num_disagree", "num_resolved", "aln_len", "genus", "species", "subspecies", "gnsp", "gnspssp", "pctIdent", "pctSim", "mt", "numt", "cp", "cgenome", "pgenome", "pCDS", "partial", "Hit_id", "Hit_def", "Hit_accession", "accession2", "Hit_len", "Hsp_num", "Hsp_bit_score", "Hsp_score", "Hsp_evalue", "Hsp_query_from", "Hsp_hit_to", "Hsp_query_frame", "Hsp_hit_frame", "Hsp_identity", "Hsp_positive", "Hsp_gaps", "Hsp_align_len", "qseq", "hseq")
		
		blank_line = rep(NA, times=length(headers))
		blank_line_df = as.data.frame(matrix(blank_line, nrow=1), stringsAsFactors=FALSE)
		names(blank_line_df) = headers
		blank_line_df
		
		blastres = NULL
		blastres$aln = list()
		blastres$seqdf = blank_line_df
		
		try_error_is = c(try_error_is, i)
		} else {
		blastres = try_result
		search_df$pTF[i] = TRUE
		}
	
	}
















##########################################################################
# Get seqIDs for ALL mitocogs (most are easily assigned to UniProt)
##########################################################################
mitocogs = c("NAD1", "NAD2", "NAD3", "NAD4", "NAD4L", "NAD5", "NAD6", "NAD7", "NAD8", "NAD9", "NAD10", "NAD11", "SDH2", "SDH3", "SDH4", "COB", "COX1", "COX2", "COX3", "ATP1", "ATP3", "ATP4", "ATP6", "ATP8", "ATP9", "RPS1", "RPS2", "RPS3", "RPS4", "RPS7", "RPS8", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS16", "RPS19", "RPL1", "RPL2", "RPL5", "RPL6", "RPL10", "RPL11", "RPL14", "RPL16", "RPL18", "RPL19", "RPL20", "RPL27", "RPL31", "RPL32", "RPL34", "RPL35", "RPL36", "EFTU")

ivals = NULL
mitocog_names = NULL
seqnames = NULL
allseqs = NULL
seqstrings = NULL
for (i in 1:length(mitocogs))
	{
	catn()
	cat("i=", i, ", mitocog=", mitocogs[i], sep="")
	tmpwd = slashslash(paste0(wd, "/", mitocogs[i]))
	setwd(tmpwd)
	
	seqs = as.character(ape::read.FASTA("seqs.fasta", type="AA"))
	seqstrings = c(seqstrings, sapply(X=seqs, FUN=paste0, collapse=""))
	allseqs = c(allseqs, seqs)
	cat("\n", length(seqs), " sequences read.", sep="")
	
	tmpnames = names(seqs)
	seqnames = c(seqnames, tmpnames)
	
	ivals = c(ivals, rep(i, length(seqs)))
	mitocog_names = c(mitocog_names, rep(mitocogs[i], length(seqs)))
	
	}
allseqs = as.AAbin(allseqs)
head(seqstrings)

# Return to top
setwd(wd)

length(seqnames)
# 1957
length(allseqs)
# 1957

all_seqids = get_leading_seqids_from_name(strings=seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
head(all_seqids)

all_taxnames = get_info_after_leading_seqid_from_name(seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))

# Make a big table
tmptable = cbind(ivals, mitocog_names, all_seqids, all_taxnames, seqstrings)
tmpdf = as.data.frame(tmptable, stringsAsFactors=FALSE)
names(tmpdf) = c("i","mitocog","seqid","taxname","seqstring")
head(tmpdf)




#######################################################
# Load the FASTA file
#######################################################
fasta_fn = "133aas_wo_UniProt.fasta"

seqs = ape::read.FASTA(fasta_fn, type="AA")
seqs_chars = as.character(seqs)
seqs_strings = sapply(X=seqs_chars, FUN=paste0, collapse="")
length(seqs)
seqs_strings

# Extract the seq ids, even though it is delimited by "_" (special code to identify NP_, WP_, etc.)
seqnames = names(seqs)
seqids = get_leading_seqids_from_name(strings=seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
seqids


# Process seqnames to genus, species
taxnames = get_info_after_leading_seqid_from_name(seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))

catn(names(rev(sort(table(taxnames)))))



#######################################################
# Start slow BLAST searches
#######################################################

col_names = c("i","mitocog","seqid","noVersion","seqname","matchto","From","Entry","Reviewed","Entry.Name","Protein.names","Gene.Names","Organism","Length","3D","taxnames")







library(ape)
library(phytools)
library(BioGeoBEARS)
library(varhandle) # for check.numeric
library(rentrez) # Fetching full records: entrez_fetch()
library(XML)
library(Biostrings) # For parsing BLAST searches downloaded

#sourceall("/GitHub/bioinfRhints/Rsrc/")
sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")



wd = "/GitHub/str2phy/ex/MitoCOGs/_03_alphafolds/"
setwd(wd)


#######################################################
# SIMPLE: Load and BLAST the FASTA file
#######################################################
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

wd = "/GitHub/str2phy/ex/MitoCOGs/_04_IDwBLAST_test1/"
setwd(wd)
#fasta_fn = "133aas_wo_UniProt.fasta"
#fasta_fn = "1-43aa_of_133_wo_UniProt.fasta"
#fasta_fn = "44-88aa_of_133_wo_UniProt.fasta"
fasta_fn = "7mitocogs.fasta"
type = "AA"
res = build_blastsearch_table_from_fasta(fasta_fn, type=type, overwrite=TRUE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

# Shorten for test
write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


taxids = NULL

database = "nr" 
hitListSize="20"
filter=NULL  # "L"
expect="10"
program="blastp"
organism=taxids
not_organism=NULL
addl_url=NULL
txt=NULL
baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
assume_txid=TRUE
res = construct_blasturls_for_table(fasta_blast_table_fn, database=database, hitListSize=hitListSize, filter=filter, expect=expect, program="blastp",
organism=NULL, not_organism=NULL, addl_url=NULL, txt=NULL, baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", assume_txid=TRUE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn

# While loop to gradually fill out table
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf



res = get_some_RIDs_for_table(fasta_blast_table_fn, num_to_get=1, wait_between=runif(1,1.1,2.5), baseUrl="default")
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn




res = launch_some_blast_searches(fasta_blast_table_fn, num_to_launch=1, wait_between=rnorm(n=1, mean=3.8, sd=0.9))
search_df = res$search_df


blastXML_fn = "/GitHub/str2phy/ex/MitoCOGs/_03_alphafolds/2025-01-20_14.31.18_SV9K6DFW016_RID_blastres_download.xml"
TF = check_file_for_automatically(blastXML_fn)


# Checks if finished; run repeatedly
res = download_finished_blast_searches(fasta_blast_table_fn, num_to_get=1, wait_between=1.0)
search_df = res$search_df
write_search_df(search_df, fasta_blast_table_fn=fasta_blast_table_fn)

matdf = check_num_still_running(fasta_blast_table_fn)
matdf


head(search_df[,1:7])

res = download_finished_blast_searches(fasta_blast_table_fn, num_to_get=1, wait_between=1.0)
search_df = res$search_df
head(search_df[,1:7])


matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf


# check a specific running RID: SNJKH5BY013
search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)

matchnum = match(x="SNJKH5BY013", table=search_df$rid)
matchnum
search_df[matchnum,]

search_df$sTF[matchnum] = FALSE
search_df$eTF[matchnum] = FALSE

write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf



# Re-do all un-saved ones
search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
search_df$sTF[search_df$eTF == TRUE] = FALSE
search_df$eTF[search_df$eTF == TRUE] = FALSE

write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf

fasta_fn = "44-88aa_of_133_wo_UniProt.fasta"
fasta_blast_table_fn = "44-88aa_wo_UniProt_blast_table_v1.txt"
type = "AA"
res = build_blastsearch_table_from_fasta(fasta_fn, type=type, overwrite=FALSE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

# Shorten for test
write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


taxids = NULL

database = "nr" 
hitListSize="20"
filter=NULL  # "L"
expect="10"
program="blastp"
organism=taxids
not_organism=NULL
addl_url=NULL
txt=NULL
baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
assume_txid=TRUE
res = construct_blasturls_for_table(fasta_blast_table_fn, database=database, hitListSize=hitListSize, filter=filter, expect=expect, program="blastp",
organism=NULL, not_organism=NULL, addl_url=NULL, txt=NULL, baseUrl="http://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", assume_txid=TRUE)
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn

# While loop to gradually fill out table
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf




head(search_df[,1:7])

res = get_some_RIDs_for_table(fasta_blast_table_fn, num_to_get=5, wait_between=runif(1,1.1,2.5), baseUrl="default")
search_df = res$search_df
fasta_blast_table_fn = res$fasta_blast_table_fn

res = launch_some_blast_searches(fasta_blast_table_fn, num_to_launch=5, wait_between=rnorm(n=1, mean=3.8, sd=0.9))
search_df = res$search_df
head(search_df[,1:7])

Sys.sleep(200)

res = download_finished_blast_searches(fasta_blast_table_fn, num_to_get=1, wait_between=1.0)
search_df = res$search_df
head(search_df[,1:7])



# While loop to gradually fill out table
matdf = check_num_still_running(fasta_blast_table_fn)
matdf

matdf2 = fill_out_table_while_loop(fasta_blast_table_fn, wait_between_loops=30, max_num_searches=5)
matdf2

matdf = check_num_still_running(fasta_blast_table_fn)
matdf




#######################################################
# Check a run
#######################################################
search_df$sTF[c(1:5)] = FALSE
search_df$eTF[c(1:5)] = FALSE
write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



#######################################################
# Revise for a re-run
#######################################################
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")
head(search_df[,1:7])
search_df = read.table(file=fasta_blast_table_fn, header=TRUE, sep="\t", quote="", row.names=NULL, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors = FALSE)
search_df$sTF[c(3,5)] = FALSE
search_df$eTF[c(3,5)] = FALSE
head(search_df[,1:7])
write.table(x=search_df, file=fasta_blast_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
res = download_finished_blast_searches(fasta_blast_table_fn, num_to_get=5, wait_between=1.0)
search_df = res$search_df
head(search_df[,1:7])

system(search_df$terminal_cmd[3])
system(search_df$terminal_cmd[5])









#######################################################
# Load the FASTA file
#######################################################
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

fasta_fn = "133aas_wo_UniProt.fasta"
type = "AA"
search_df = build_blastsearch_table_from_fasta(fasta_fn, type=type)
head(search_df)
class(search_df)
dim(search_df)



i = 1

cat("\n(Starting search #", i, "/", length(seqs_to_search), ")\n")
Sys.sleep(time=1) # pause for ESC

searchseq = search_df$seqstring[i]

# Database options:
# UniProtKB/Swiss-Prot (swissprot)
# Non-redundant protein sequences (nr)
# Reference proteins (refseq_protein)
# Reference Select proteins (refseq_select)
# Protein Data Bank proteins (pdb)
# 
# See also: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp

# Construct the BLAST URL
taxids = NULL
database = "nr" 
results = construct_blasturl(searchseq=searchseq, database="swissprot", hitListSize="25", 
												filter=NULL, expect="1", program="blastp",
												organism=taxids, not_organism=NULL, addl_url=NULL, 
												txt=NULL, baseUrl=baseUrl)
search_df$cTF[i] = TRUE

# get the RID and waiting time
results = get_blast_RID_waiting_time(results=results, RID_url=NULL, baseUrl=baseUrl)

search_df$RID_url_orig[i] = results$RID_url_orig
search_df$RID_url[i] = results$RID_url
search_df$search_url[i] = results$search_url
search_df$terminal_cmd[i] = results$terminal_cmd
search_df$blast_RID_fn[i] = results$blast_RID_fn
search_df$rid[i] = results$rid
search_df$rtoe[i] = results$rtoe

search_df$rTF[i] = TRUE


# Launch the web search; will have to be retrieved later
blastres_multiline = NULL  # BLAST results, may well have >1 per searchseq

search_df$ltime[i] = Sys.time()
start_BLAST_websearch(results)
search_df$lTF[i] = TRUE

# Retrieve results
search_df$search_url[i]

tryParse_results = download_BLAST_websearch(results=NULL, search_url=search_df$search_url[i], attempts=1, retry_wait=1)

if ((length(tryParse_results) == 1) && (grepl(pattern="ERROR in download_BLAST_websearch", x=tryParse_results) == FALSE))
	{
	#print(tryParse_results$result)
	search_df$sTF[i] = FALSE
	} else {
	search_df$sTF[i] = TRUE
	search_df$stime[i] = Sys.time()
	}


if (search_df$sTF[i] == TRUE)
	{
	try_result = try(process_blastresult(tryParse_results=tryParse_results))

	if (class(try_result) == "try-error")
		{
		search_df$pTF[i] = FALSE
		
		warning_txt = paste0("WARNING: process_blastresult() gave a try-error on i=", i, ", abbr=", abbr[i], ". Inserting a blank line of NAs in blastres$seqdf.")
		
		cat("\n\n")
		warning(warning_txt)
		cat("\n\n")
		
		headers = c("Hit_num", "downID", "taxon", "gi1", "gi2", "gi3", "gi4", "pctIdent", "num_identical", "num_disagree", "num_resolved", "aln_len", "genus", "species", "subspecies", "gnsp", "gnspssp", "pctIdent", "pctSim", "mt", "numt", "cp", "cgenome", "pgenome", "pCDS", "partial", "Hit_id", "Hit_def", "Hit_accession", "accession2", "Hit_len", "Hsp_num", "Hsp_bit_score", "Hsp_score", "Hsp_evalue", "Hsp_query_from", "Hsp_hit_to", "Hsp_query_frame", "Hsp_hit_frame", "Hsp_identity", "Hsp_positive", "Hsp_gaps", "Hsp_align_len", "qseq", "hseq")
		
		blank_line = rep(NA, times=length(headers))
		blank_line_df = as.data.frame(matrix(blank_line, nrow=1), stringsAsFactors=FALSE)
		names(blank_line_df) = headers
		blank_line_df
		
		blastres = NULL
		blastres$aln = list()
		blastres$seqdf = blank_line_df
		
		try_error_is = c(try_error_is, i)
		} else {
		blastres = try_result
		search_df$pTF[i] = TRUE
		}
	
	}
















##########################################################################
# Get seqIDs for ALL mitocogs (most are easily assigned to UniProt)
##########################################################################
mitocogs = c("NAD1", "NAD2", "NAD3", "NAD4", "NAD4L", "NAD5", "NAD6", "NAD7", "NAD8", "NAD9", "NAD10", "NAD11", "SDH2", "SDH3", "SDH4", "COB", "COX1", "COX2", "COX3", "ATP1", "ATP3", "ATP4", "ATP6", "ATP8", "ATP9", "RPS1", "RPS2", "RPS3", "RPS4", "RPS7", "RPS8", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS16", "RPS19", "RPL1", "RPL2", "RPL5", "RPL6", "RPL10", "RPL11", "RPL14", "RPL16", "RPL18", "RPL19", "RPL20", "RPL27", "RPL31", "RPL32", "RPL34", "RPL35", "RPL36", "EFTU")

ivals = NULL
mitocog_names = NULL
seqnames = NULL
allseqs = NULL
seqstrings = NULL
for (i in 1:length(mitocogs))
	{
	catn()
	cat("i=", i, ", mitocog=", mitocogs[i], sep="")
	tmpwd = slashslash(paste0(wd, "/", mitocogs[i]))
	setwd(tmpwd)
	
	seqs = as.character(ape::read.FASTA("seqs.fasta", type="AA"))
	seqstrings = c(seqstrings, sapply(X=seqs, FUN=paste0, collapse=""))
	allseqs = c(allseqs, seqs)
	cat("\n", length(seqs), " sequences read.", sep="")
	
	tmpnames = names(seqs)
	seqnames = c(seqnames, tmpnames)
	
	ivals = c(ivals, rep(i, length(seqs)))
	mitocog_names = c(mitocog_names, rep(mitocogs[i], length(seqs)))
	
	}
allseqs = as.AAbin(allseqs)
head(seqstrings)

# Return to top
setwd(wd)

length(seqnames)
# 1957
length(allseqs)
# 1957

all_seqids = get_leading_seqids_from_name(strings=seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
head(all_seqids)

all_taxnames = get_info_after_leading_seqid_from_name(seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))

# Make a big table
tmptable = cbind(ivals, mitocog_names, all_seqids, all_taxnames, seqstrings)
tmpdf = as.data.frame(tmptable, stringsAsFactors=FALSE)
names(tmpdf) = c("i","mitocog","seqid","taxname","seqstring")
head(tmpdf)




#######################################################
# Load the FASTA file
#######################################################
fasta_fn = "133aas_wo_UniProt.fasta"

seqs = ape::read.FASTA(fasta_fn, type="AA")
seqs_chars = as.character(seqs)
seqs_strings = sapply(X=seqs_chars, FUN=paste0, collapse="")
length(seqs)
seqs_strings

# Extract the seq ids, even though it is delimited by "_" (special code to identify NP_, WP_, etc.)
seqnames = names(seqs)
seqids = get_leading_seqids_from_name(strings=seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))
seqids


# Process seqnames to genus, species
taxnames = get_info_after_leading_seqid_from_name(seqnames, split="_", recodes=genbank_prefixes(), removetxt=c(">"))

catn(names(rev(sort(table(taxnames)))))



#######################################################
# Start slow BLAST searches
#######################################################

col_names = c("i","mitocog","seqid","noVersion","seqname","matchto","From","Entry","Reviewed","Entry.Name","Protein.names","Gene.Names","Organism","Length","3D","taxnames")






