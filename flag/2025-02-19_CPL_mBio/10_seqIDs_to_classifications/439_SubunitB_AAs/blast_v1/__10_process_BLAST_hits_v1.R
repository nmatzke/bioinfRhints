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
wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/439_SubunitB_AAs/blast_v1/"
setwd(wd)

#######################################################
# SIMPLE: Load and BLAST the FASTA file
#######################################################
fasta_fn = "439_SubunitB_AAs_raw.fasta"
type = "AA"
# This blanks the file
fasta_blast_table_fn = "439_SubunitB_AAs_raw_blast_table_v1_ARCHIVE.txt"

# Read in, write out the blast search results table
search_df = read_search_df(fasta_blast_table_fn)
matdf = check_num_still_running(fasta_blast_table_fn); matdf


# Hmm, many downloaded but not processed.
fns = list.files(wd, pattern="blast_results_rownum")
fns

hits_df_names = c("rownum","RID","seqid","seqname","taxname","searchseq","Hit_num","downID","taxon","gi1","gi2","gi3","gi4","pctIdent","num_identical","num_disagree","num_resolved","aln_len","genus","species","subspecies","gnsp","gnspssp","pctIdent.1","pctSim","mt","numt","cp","cgenome","pgenome","pCDS","partial","Hit_id","Hit_def","Hit_accession","accession2","Hit_len","Hsp_num","Hsp_bit_score","Hsp_score","Hsp_evalue","Hsp_query_from","Hsp_hit_to","Hsp_query_frame","Hsp_hit_frame","Hsp_identity","Hsp_positive","Hsp_gaps","Hsp_align_len","qseq","hseq")

catn("Processing ", nrow(search_df), " rows of the search_df...")

best_blast_hits_df = NULL
for (i in 1:nrow(search_df))
	{
	cat(i, ",", sep="")
	seqid = search_df$seqid[i]
	fn_matchnums = match_grepl(pattern=seqid, x=fns, return_counts=FALSE)
	fn_matchnums

	if (length(fn_matchnums) > 1)
		{
		stop("STOP ERROR: More than one file matched. Delete one and repeat.")
		}
	if ( (length(fn_matchnums) == 1) && (is.na(fn_matchnums) == TRUE) )
		{
		tdf2 = as.data.frame(matrix(rep(NA, times=length(hits_df_names)), nrow=1), stringsAsFactors=FALSE)
		names(tdf2) = c(hits_df_names)
		keepTF1 = FALSE
		} else {
		fn = fns[fn_matchnums]
		tdf = read.table(file=fn, header=TRUE, sep="\t", quote="", strip.white=TRUE, fill=TRUE)
		keepTF1 = tdf$pctIdent.1 == 100.0
		keepTF1[is.na(keepTF1)] = FALSE  # error catch
		}
	
	if (sum(keepTF1) > 0)
		{
		tdf2 = tdf[keepTF1,]
		} else {
		tdf2 = tdf[1,]
		}
	
	df_preface2 = search_df[rep(i, times=nrow(tdf2)),]
	df_to_add = cbind(df_preface2, tdf2)
	
	best_blast_hits_df = rbind(best_blast_hits_df, df_to_add)
	}

catn("...done.")
dim(best_blast_hits_df)

# Error after 465, seqid 428
TF = is.na(best_blast_hits_df$stime)
best_blast_hits_df[]

match(x="stime", table=names(best_blast_hits_df))
# 26
dim(best_blast_hits_df)
# 476  89

best_blast_hits_df[TF, 27:89] = NA

best_blast_hits_df[464:476,-c(15:18,44,88:89)]


best_blast_hits_df$pctIdent.1

best_blast_hits_df$seqid[is.na(best_blast_hits_df$pctIdent.1)]


