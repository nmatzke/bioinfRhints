# 4.2NCBI identical protein retrieval
# https://bioconductor.riken.jp/packages/3.21/bioc/vignettes/ginmappeR/inst/doc/ginmappeR.html
library(ape)
library(BioGeoBEARS)
library(ginmappeR)
library(queryup)

sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/379_MotABs/"
setwd(wd)

motB_seqs_w_MotA_prefix_fn = "motB_filtered.fasta"
motB_seqs_w_MotA_prefix = read_FASTA_safe(motB_seqs_w_MotA_prefix_fn, type="AA")
names(motB_seqs_w_MotA_prefix)
motB_seqnames = all_but_prefixes(names(motB_seqs_w_MotA_prefix), split="_")
motB_seqnames


# MotA Beast2 MCC tree
trfn = "A379tree.nexus"
tr = read.nexus(trfn)
seqnames = tr$tip.label
head(seqnames)
seqidsA = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqidsA[1:10])

#######################################################
# MotB Beast2 MCC tree
#######################################################
trfn = "B379tree.nexus"
trB = read.nexus(trfn)
trB$tip.label = gsub(pattern="'", replacement="", x=trB$tip.label)
seqnames = trB$tip.label
head(seqnames)
seqidsB = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqidsB[1:10])

#######################################################
# MotAB Beast2 MCC tree
#######################################################
trfn = "AB379tree.nexus"
trA = read.nexus(trfn)
trA$tip.label = gsub(pattern="'", replacement="", x=trA$tip.label)
seqnames = trA$tip.label
head(seqnames)
seqidsA = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqidsA[1:10])



run_bigslow = FALSE
if (run_bigslow == TRUE)
{

######################################################
# MotAs from Caroline's paper
#######################################################

# 5-10 minutes the first time - needs Terminal, not R.app
bigdf_379_AQBs = get_uniprot_data_on_seqids(seqidsA, runslow=TRUE, base_fn="379_AQBs", version="v1")
head(bigdf_379_AQBs)
dim(bigdf_379_AQBs)


# 30 seconds to re-load from saved results
bigdf_379_AQBs = get_uniprot_data_on_seqids(seqidsA, runslow=FALSE, base_fn="379_AQBs", version="v1")
head(bigdf_379_AQBs)
dim(bigdf_379_AQBs)


AQBs_379_df = cbind(seqidsA, bigdf_379_AQBs)

write.table(x=AQBs_379_df, file="379_AQBs_uniProt_info_table_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




######################################################
# MotBs from Caroline's paper
#######################################################

# 5-10 minutes the first time
bigdf_379_BRDs = get_uniprot_data_on_seqids(seqidsB, runslow=TRUE, base_fn="379_BRDs", version="v1")
head(bigdf_379_BRDs)
dim(bigdf_379_BRDs)


# 30 seconds to re-load from saved results
bigdf_379_BRDs = get_uniprot_data_on_seqids(seqidsB, runslow=FALSE, base_fn="379_BRDs", version="v1")
head(bigdf_379_BRDs)
dim(bigdf_379_BRDs)


BRDs_379_df = cbind(seqidsB, bigdf_379_BRDs)

write.table(x=BRDs_379_df, file="379_BRDs_uniProt_info_table_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

} # END if (run_bigslow == TRUE)

#######################################################
# Read in from saved tables
#######################################################
AQBs_379_df = read.table(file="379_AQBs_uniProt_info_table_v1.txt", header=TRUE, sep="\t", quote="", row.names=NULL, fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE)

BRDs_379_df = read.table(file="379_BRDs_uniProt_info_table_v1.txt", header=TRUE, sep="\t", quote="", row.names=NULL, fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE)

names_uniProt_info_table = names(AQBs_379_df)
cols_to_replace = names_uniProt_info_table[-1]
cols_to_replace

# Some IDs have no label info
TF1 = isblank_TF(AQBs_379_df$accession)
TF2 = isblank_TF(AQBs_379_df$seqids_wGenBank_codes_txt) == FALSE
TF = (TF1 + TF2) == 2
sum(TF)
tmpinfo = get_uniprot_data_on_seqids(seqids=AQBs_379_df$seqids_wGenBank_codes_txt[TF], runslow=TRUE, base_fn="seqids", version="v1")
names(tmpinfo) = cols_to_replace
dim(tmpinfo)

AQBs_379_df[TF,cols_to_replace] = tmpinfo[,cols_to_replace]
tail(AQBs_379_df)

tmpTF = AQBs_379_df$seqids_wGenBank_codes_txt[TF] == tmpinfo$seqids_wGenBank_codes_txt
AQBs_379_df$seqids_wGenBank_codes_txt[TF][tmpTF==FALSE]
tmpinfo$seqids_wGenBank_codes_txt[tmpTF==FALSE]
newTF = isblank_TF(AQBs_379_df$accession)
sum(newTF)
AQBs_379_df[newTF,]

names(AQBs_379_df)
names(tmpinfo)

res = get_IDs_identical_proteins(seqids=AQBs_379_df$seqid[newTF], runslow=TRUE, identical_protIDs_fn="tmp_protIDs_identical_to_seqIDs_v1.Rdata")

recodes = genbank_prefixes()
res2_reasonable_list = sapply(X=res, FUN=reduce_identical_proteins_to_reasonable_list, codes=recodes)
sum(sapply(X=res2_reasonable_list, FUN=is.null))
sum(sapply(X=res, FUN=is.null))









TF1 = isblank_TF(BRDs_379_df$accession)
TF2 = isblank_TF(BRDs_379_df$seqids_wGenBank_codes_txt) == FALSE
TF = (TF1 + TF2) == 2
sum(TF)
tmpinfo = get_uniprot_data_on_seqids(seqids=, runslow=TRUE, base_fn="seqids", version="v1")


get_uniprot_data_on_seqids(seqids=, runslow=TRUE, base_fn="seqids", version="v1")
