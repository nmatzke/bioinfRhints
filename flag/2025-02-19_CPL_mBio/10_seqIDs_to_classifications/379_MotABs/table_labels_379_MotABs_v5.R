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


bigdf_outdf_A3 = cbind(seqidsA, bigdf_379_AQBs)

write.table(x=bigdf_outdf_A3, file="379_AQBs_uniProt_info_table_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




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


bigdf_outdf_B3 = cbind(seqidsB, bigdf_379_BRDs)

write.table(x=bigdf_outdf_B3, file="379_BRDs_uniProt_info_table_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
