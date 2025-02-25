
#######################################################
# Batch CD CDART (conserved domain)
#######################################################

#######################################################
# 08_get_uniprotID_for_each_seq_v1.R
#
# Get UniProt IDs for each sequence
#  (as these will have pre-calculated monomer alphafolds)
#######################################################

#######################################################
# The plan: organized 3di data collection, in stages, filling
# in a TABLE, recording errors.
#
# 1. Get UniProt IDs to match
# 
# Retrieve/ID mapping
# 
# https://www.uniprot.org/id-mapping
#
# 
# 2. When that fails, try BLAST for UniProt
# 3. Download alphafolds for UniProts (as these are pre-calculated)
# 4. Calculate remaining alphafolds, using colabfold
#######################################################

#######################################################
# Download and process GenPept records for each 
# seqid
# 
# process with biofiles: https://github.com/gschofl/biofiles
#######################################################
# devtools::install_github("gschofl/biofiles")

library(reutils)			 # for content()
library(GenomicRanges) # bioconductor
library(biofiles)

library(openxlsx)
library(ape)
library(phytools)
library(BioGeoBEARS)
library(varhandle) # for check.numeric
library(rentrez) # Fetching full records: entrez_fetch()

#sourceall("/GitHub/bioinfRhints/Rsrc/")
sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")


wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/"
setwd(wd)

# Save seqids to:
xlsfn = "439_subunitA_AAs_find_swissprot_UniProt_v1.xlsx"


#######################################################
# OLD STUFF
#######################################################
# OLD: Load FASTA file
old='
fasta_fn = "439_subunitA_AAs_raw.fasta"
seqs = as.character(read_FASTA_safe(fasta_fn, type="AA"))
seqs

# Convert to upper-case
for (i in 1:length(seqs))
	{
	seqs[[i]] = toupper(seqs[[i]])
	}
seqnames = names(seqs)
seqids = get_leading_seqids_from_name(strings=seqnames, split="_")

#######################################################
# Get single record
#######################################################
gb_file <- reutils::efetch(uid="AAC74960", db="protein", rettype = "gp", retmode = "text")
rec <- biofiles::gbRecord(gb_file)
rec
biofiles::uniqueQualifs(rec)
# "organism"     "strain"       "sub_strain"   "db_xref"      "product"      "region_name"  "note"         "gene"         "locus_tag"    "gene_synonym" "coded_by"     "transl_table"

f_single <- biofiles::getFeatures(rec)
class(f_single)
names(f_single)

feature_table_counts = biofiles::featureTable(rec)
feature_table_counts
names(feature_table_counts)

# db_xref: Database cross-references
# https://www.ncbi.nlm.nih.gov/genbank/collab/db_xref/

# biofiles::dbxref(f_single) # FAILS
biofiles::dbxref(f_single["CDS"])
biofiles::dbxref(f_single["Protein"])
biofiles::dbxref(f_single["Region"])
biofiles::dbxref(f_single["source"])
'
#######################################################
# END OLD STUFF
#######################################################


#######################################################
# Current: Load Caroline's tree
#######################################################
trfn = "A379tree.nexus"
tr = read.nexus(trfn)
seqnames = tr$tip.label
seqnames
seqids = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqids[1:10])


#######################################################
# Get lots of records (slower)
#######################################################
runslow = FALSE
tmpfn = "./gb_files.gp"
rds_local_fn = "~/Downloads/recs.rds"
if (runslow)
	{
	gb_files <- reutils::efetch(uid=seqids, db="protein", rettype = "gp", retmode = "text")
	write(reutils::content(gb_files, "text"), file=tmpfn)
	save(gb_files, file="gb_files.Rdata")
	} else {
	# Loads to: gb_files
	#load(file="gb_files.Rdata")
	gb_files = readLines(tmpfn)
	}
runslow = TRUE
if (runslow)
	{
	# parse the efetch object into a gbRecord instance; seems absurdly slow, ~1 per second
	starttime2 = Sys.time()
	recs <- biofiles::gbRecord(rcd=tmpfn, progress=TRUE)
	endtime2 = Sys.time()
	endtime2 - starttime2 # 7.5 minutes
	
	starttime3  = Sys.time()
	saveRDS(recs, file=rds_local_fn) # ~5 minutes
	endtime3 = Sys.time()
	endtime3 - starttime3 # 5 minutes
	
	# Parallelisation doesn't seem to work well in R.app - try Terminal
	# parse (locally saved)
	attempt_parallel_doesnt_work_fast='
	library(doParallel)
	registerDoParallel(cores = 23)
	starttime1 = Sys.time()
	recs <- gbRecord(rcd=tmpfn)
	endtime1 = Sys.time()
	endtime1 - starttime1 # 3.453269 with 1 core in Terminal; 9.294752 mins with 23 cores
	stopImplicitCluster()
	'
	} else {
	# Loads to: recs
	recs = readRDS(file=rds_local_fn) # 1-2 minutes
	}


recs[1:2]

# db_xref and other features
biofiles::uniqueQualifs(recs[[1]])

# Get the features
f <- biofiles::getFeatures(recs)
class(f)
class(f[[1]])
class(f[[2]])
names(f)

# Names and counts of the features:
feature_table_counts = biofiles::featureTable(recs)
allnames = unique(unlist(sapply(X=feature_table_counts, names)))
allnames
# "CDS"     "Protein" "Region"  "source" "sig_peptide"

# Get the db_xrefs for each CDS
biofiles::dbxref(recs[[1]])
biofiles::dbxref(recs[[1]]["CDS"])
biofiles::dbxref(recs[[1]]["Protein"])
biofiles::dbxref(recs[[1]]["Region"])
biofiles::dbxref(recs[[1]]["source"])

#######################################################
# Assemble all of the db_xrefs across all recs
#######################################################
sourceall("/GitHub/str2phy/Rsrc/")
# parse the efetch object into a gbRecord instance.
dbxrefs_table = get_db_xref_for_seqids(seqids, recs)




#######################################################
# Add the accession info, etc.
#######################################################

accessions = biofiles::getAccession(recs)
#> [1] "CP002806"
geneIDs = biofiles::getGeneID(recs)
#> [1] "NA"
seq_labels = biofiles::getDefinition(recs)
#> [1] "Chlamydia psittaci 02DC15, complete genome."
organisms = biofiles::getOrganism(recs)
#> [1] "Chlamydia psittaci 02DC15"
#biofiles::getSequence(recs)
biofiles::getLength(recs)
biofiles::getLocus(recs)
cbind(accessions, geneIDs, seq_labels, organisms)

biofiles::getGeneID(recs, db = "gi")

dbsource = biofiles::getDBSource(recs)
dbsource_prefix = all_but_suffixes(tmptxt=dbsource, split=" ")
dbsource_prefix
dbsource_accession = get_lastwords(OTUs=dbsource, split=" ")
dbsource_accession

catn(dbsource_accession)

seq_labels2 = seq_labels
for (i in 1:length(seq_labels))
	{
	words = strsplit(seq_labels2[i], split="\\[")[[1]]
	seq_labels2[i] = gdata::trim(words[1])
	}
seq_labels2

tgi5_1 = "QDU25289"
tgi5_2 = "AAC74960"

tgi4_1 = "AAG04849"
tgi4_2 = "CAL34488"


t51 = match_grepl(tgi5_1, tr$tip.label, return_counts=FALSE)
t52 = match_grepl(tgi5_2, tr$tip.label, return_counts=FALSE)
t41 = match_grepl(tgi4_1, tr$tip.label, return_counts=FALSE)
t42 = match_grepl(tgi4_2, tr$tip.label, return_counts=FALSE)
t51
t42


tgi5 = getMRCA(phy=tr, tip=c(t51, t52))
tgi4 = getMRCA(phy=tr, tip=c(t41, t42))
fit = getMRCA(phy=tr, tip=c(t41, t51))
# 652

trtable = prt(tr)
ids_in_tgi5 = get_leading_seqids_from_name(strsplit(trtable$tipnames[tgi5], split=",")[[1]])
ids_in_tgi4 = get_leading_seqids_from_name(strsplit(trtable$tipnames[tgi4], split=",")[[1]])
ids_in_fit = get_leading_seqids_from_name(strsplit(trtable$tipnames[fit], split=",")[[1]])

nums_in_tgi5 = match_grepl(ids_in_tgi5, accessions, return_counts=FALSE)
nums_in_tgi4 = match_grepl(ids_in_tgi4, accessions, return_counts=FALSE)
nums_in_fit = match_grepl(ids_in_fit, accessions, return_counts=FALSE)

git_TF = (1:length(tr$tip.label) %in% nums_in_fit) == FALSE
sum(git_TF)
nums_in_git = (1:length(git_TF))[git_TF]

t(t(rev(sort(table(seq_labels2[nums_in_tgi5])))[1:10]))
t(t(rev(sort(table(seq_labels2[nums_in_tgi4])))[1:10]))
t(t(rev(sort(table(seq_labels2[nums_in_fit])))[1:10]))
t(t(rev(sort(table(seq_labels2[nums_in_git])))[1:10]))

t(t(round(rev(sort(table(seq_labels2[nums_in_tgi5])))[1:10] / sum(table(seq_labels2[nums_in_tgi5])), 4) * 100))
t(t(round(rev(sort(table(seq_labels2[nums_in_tgi4])))[1:10] / sum(table(seq_labels2[nums_in_tgi4])), 4) * 100))
t(t(round(rev(sort(table(seq_labels2[nums_in_fit])))[1:10] / sum(table(seq_labels2[nums_in_fit])), 4) * 100))
t(t(round(rev(sort(table(seq_labels2[nums_in_git])))[1:10] / sum(table(seq_labels2[nums_in_git])), 4) * 100))


# Retrieve/ID mapping
# https://www.uniprot.org/id-mapping
# AAC74960 +438 RefSeq_Protein → UniProtKB
# 15 mapped

# AAC74960 +438 UniProtKB_AC-ID → UniProtKB
# 0 mapped

# AAC74960 +438 EMBL-GenBank-DDBJ → UniProtKB
# 0 mapped
 
# Added .1...
# AAC74960.1 +438 Ensembl_Protein → UniProtKB
# 0 hits

# AAC74960.1 +438 Ensembl → UniProtKB-Swiss-Prot
# 0 hits

# AAC74960.1 +438 Ensembl_Genomes_Protein → UniProtKB
# 0 hits











#######################################################
# Query Ensembl bacteria:
# https://rest.ensembl.org/documentation/info/xref_id
#######################################################

library(httr)
library(jsonlite)
library(xml2)
 
server <- "https://rest.ensembl.org"
ext <- "/xrefs/id/AAC74960.1?db_type=pan_homology;"
ext <- "/xrefs/id/EFTU_ECOLI?db_type=pan_ensembl;"
r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
 
stop_for_status(r)
 
# use this if you get a simple nested list back, otherwise inspect its structure
# head(data.frame(t(sapply(content(r),c))))
head(fromJSON(toJSON(content(r))))










Pantaxonomic compara multi-species



library(biomaRt)

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")

datasets <- listDatasets(ensembl)
head(datasets)

refseqids = c("NM_005359","NM_000546")
ipro = getBM(attributes=c("refseq_mrna","interpro","interpro_description"), 
             filters="refseq_mrna",
             values=refseqids, 
             mart=ensembl)
ipro

