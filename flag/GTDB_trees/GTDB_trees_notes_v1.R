#######################################################
# Various GTDB trees:
# https://itol.embl.de/shared/ecogenomics
#######################################################

# Release 95, 23,000 genomes
# https://data.gtdb.ecogenomic.org/releases/release95/95.0/genomic_files_reps/


# Release 95, 23,000 genomes
# https://data.gtdb.ecogenomic.org/releases/release95/95.0/genomic_files_reps/

# Release 207, genomes
#https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/
# FASTA file: bac120_ssu_reps_r207.tar.gz

# Latest release is 220:
# https://data.gtdb.ecogenomic.org/releases/release220/220.0/


library(ape)
library(BioGeoBEARS)
library(ginmappeR)

sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

#wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/"
#setwd(wd)

wd = "/GitHub/bioinfRhints/flag/GTDB_trees/"
setwd(wd)

trfn = "bac120_r220.tree"
tr = read.tree(file=trfn)

tipnames = tr$tip.label
length(tipnames)
head(tipnames)

ex='
cd /Users/nickm/GitHub/bioinfRhints/flag/GTDB_trees/
ls
head bac120_ssu_reps_r220.fna

awk 'sub(/^>/, "")' bac120_ssu_reps_r220.fna > bac120_ssu_reps_r220_headers.txt

wc -l bac120_ssu_reps_r220_headers.txt
# 58102

head bac120_ssu_reps_r220_headers.txt

sed 's/;/\t/g' bac120_ssu_reps_r220_headers.txt > bac120_ssu_reps_r220_headers_noSemis.txt
head bac120_ssu_reps_r220_headers_noSemis.txt
'

# Read the headers file to a table
fn = "bac120_ssu_reps_r220_headers_noSemis.txt"
bac120_r220_df = read.table(file=fn, header=FALSE, sep="", quote="", strip.white=TRUE, fill=TRUE)
names(bac120_r220_df) = c("genomeID", "domain", "phylum", "class", "order", "family", "genus", "species1","species2", "locus_tag", "location", "ssu_len", "contig_len")
head(bac120_r220_df)
dim(bac120_r220_df)
# 58102

TF = tipnames %in% bac120_r220_df$genomeID
sum(TF)
length(TF)
length(tipnames)
# 107235

# Reduce the 107235 tree to 58102 tips matching the table
tips_to_remove = tipnames[TF == FALSE]
head(tips_to_remove)
tr2 = drop.tip(phy=tr, tip=tips_to_remove)
length(tr2$tip.label)
# 58102

tr2fn = "bac120_r220_58102tips.newick"
write.tree(tr2, file=tr2fn)


#######################################################
# Now, cross reference the genomes to our Excel table of genomes we've got; remove those phyla
#######################################################
xlsfn = "/GitHub/bioinfRhints/flag/get_MotBs/groupTax_1282_mafftConstr_2023-08-07_edit_wGeneOrder_BestMotBs.xlsx"

# Try the species list
xlsfn = "/GitHub/bioinfRhints/flag/AQB_classification/species_list_10071623_NJMh.xlsx"

xls = openxlsx::read.xlsx(xlsfn)
head(xls)
dim(xls)
# 423  34

#xls_genome_ids = firstwords(strings=xls$genome_id, split="\\.")
xls_genome_ids = firstwords(strings=xls$GenBank.ID, split="\\.")
xls_genome_ids = xls_genome_ids[!is.na(xls_genome_ids)]
head(xls_genome_ids)
tail(xls_genome_ids)

xls_genome_ids_GCF = gsub(pattern="GCA", replacement="GCF", x=xls_genome_ids)
tipnames_GCF = gsub(pattern="GCA", replacement="GCF", x=bac120_r220_df$genomeID)

xls_rows_in_tree_GCF = match_grepl(patterns=xls_genome_ids_GCF, x=tipnames_GCF, return_counts=FALSE)
sum(is.na(xls_rows_in_tree_GCF))
sum(!is.na(xls_rows_in_tree_GCF)) # 249
length(xls_rows_in_tree_GCF) # 420

tipnames_GCF_prefix = firstwords(strings=tipnames_GCF, split="_")
table(tipnames_GCF_prefix)
# GB    RS  (GenBank and RefSeq)
tipnames_GCF2 = gsub(pattern="GB_", replacement="", x=tipnames_GCF)
tipnames_GCF3 = gsub(pattern="RS_", replacement="", x=tipnames_GCF2)
tipnames_GCF4 = firstwords(strings=tipnames_GCF3, split="\\.")
head(tipnames_GCF4)
xls_rownums_of_gtdb_tips_GCF = match_grepl(patterns=tipnames_GCF4, x=xls_genome_ids_GCF, return_counts=FALSE)
sum(is.na(xls_rownums_of_gtdb_tips_GCF))
sum(!is.na(xls_rownums_of_gtdb_tips_GCF)) # 222
length(xls_rownums_of_gtdb_tips_GCF) # 58102


# Phyla already covered in the sampled genomes
tips_in_xls_species_TF = !is.na(xls_rownums_of_gtdb_tips_GCF)
sum(tips_in_xls_species_TF)

tiprows_in_xls_species = bac120_r220_df[TF,]
tiprows_in_xls_species

# tree phyla already included:
phyla_already_included_in_tree = sort(unique(tiprows_in_xls_species$phylum))
length(phyla_already_included_in_tree)
# 149

# Cut tips that are from phyla already included, except for those in species table
phylum_already_used_TF = bac120_r220_df$phylum %in% phyla_already_included_in_tree
sum(phylum_already_used_TF)
length(phylum_already_used_TF)

tips_with_no_sampling_needed = tr2$tip.label[phylum_already_used_TF]
tree_with_all_unsampled_phyla = drop.tip(phy=tr2, tip=tips_with_no_sampling_needed)
length(tree_with_all_unsampled_phyla$tip.label)

plot(tree_with_all_unsampled_phyla)
bac120_r220_df[phylum_already_used_TF==FALSE,]