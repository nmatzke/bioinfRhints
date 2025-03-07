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
seqidsAB = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqidsAB[1:10])



