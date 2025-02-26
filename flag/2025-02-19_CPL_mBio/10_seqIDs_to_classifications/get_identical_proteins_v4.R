

# 4.2NCBI identical protein retrieval
# https://bioconductor.riken.jp/packages/3.21/bioc/vignettes/ginmappeR/inst/doc/ginmappeR.html
library(ape)
library(BioGeoBEARS)
library(ginmappeR)

sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/"
setwd(wd)

trfn = "A379tree.nexus"
tr = read.nexus(trfn)
seqnames = tr$tip.label
head(seqnames)
seqids = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqids[1:10])



# 5-10 minutes the first time
bigdf_outdf = get_uniprot_data_on_seqids(seqids, runslow=TRUE, base_fn="379_AQBs", version="v1")
head(bigdf_outdf)
dim(bigdf_outdf)


# 30 seconds to re-load from saved results
bigdf_outdf = get_uniprot_data_on_seqids(seqids, runslow=FALSE, base_fn="379_AQBs", version="v1")
head(bigdf_outdf)
dim(bigdf_outdf)


bigdf_outdf_3 = cbind(seqnames, bigdf_outdf)

write.table(x=bigdf_outdf_3, file="379_AQBs_uniProt_info_table_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

