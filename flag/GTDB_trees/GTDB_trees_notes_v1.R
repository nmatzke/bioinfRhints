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

# Re-order bac120_r220_df to match tree
convert_rows_of_bac120_r220_genomeID_to_tr2 = match(tr2$tip.label, table=bac120_r220_df$genomeID)
sum(is.na(convert_rows_of_bac120_r220_genomeID_to_tr2))

bac120_r220_df = bac120_r220_df[convert_rows_of_bac120_r220_genomeID_to_tr2,]



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

tiprows_in_xls_species = bac120_r220_df[tips_in_xls_species_TF,]
tiprows_in_xls_species

# tree phyla already included:
phyla_already_included_in_tree = sort(unique(tiprows_in_xls_species$phylum))
length(phyla_already_included_in_tree)
# 35

# Sanity check
bac120_r220_df$genomeID[1:10]
tr2$tip.label[1:10]

sort(bac120_r220_df$genomeID[1:10])
sort(tr2$tip.label[1:10])


# Cut tips that are from phyla already included, except for those in species table
phylum_already_used_TF = bac120_r220_df$phylum %in% phyla_already_included_in_tree
sum(phylum_already_used_TF)
length(phylum_already_used_TF)

tips_with_no_sampling_needed = tr2$tip.label[phylum_already_used_TF]
tree_with_all_unsampled_phyla = drop.tip(phy=tr2, tip=tips_with_no_sampling_needed)
length(tree_with_all_unsampled_phyla$tip.label)
# 3231


#plot(tree_with_all_unsampled_phyla)
bac120_r220_wUnsampledPhyla_tiplabels_table = bac120_r220_df[phylum_already_used_TF==FALSE,]
dim(bac120_r220_wUnsampledPhyla_tiplabels_table)

length(unique(bac120_r220_wUnsampledPhyla_tiplabels_table$phylum)) # 136
sort(unique(bac120_r220_wUnsampledPhyla_tiplabels_table$phylum))
rev(sort(table(bac120_r220_wUnsampledPhyla_tiplabels_table$phylum)))
# p__Bacillota_I p__Marinisomatota p__WOR-3 p__Desulfobacterota_B p__Dependentiae p__Bipolaricaulota p__Zixibacteria 
# 1257 138 113 102 88 82 69 

paste0(bac120_r220_wUnsampledPhyla_tiplabels_table$species1, bac120_r220_wUnsampledPhyla_tiplabels_table$species2, collapse="", sep="_")
bac120_r220_wUnsampledPhyla_tiplabels = mapply(FUN=paste, bac120_r220_wUnsampledPhyla_tiplabels_table$genomeID, bac120_r220_wUnsampledPhyla_tiplabels_table$phylum, bac120_r220_wUnsampledPhyla_tiplabels_table$species1, bac120_r220_wUnsampledPhyla_tiplabels_table$species2, MoreArgs=list(sep="_"))
bac120_r220_wUnsampledPhyla_tiplabels = gsub(pattern="p__", replacement="PHYLUM_", x=bac120_r220_wUnsampledPhyla_tiplabels)
bac120_r220_wUnsampledPhyla_tiplabels = gsub(pattern="s__", replacement="sp_", x=bac120_r220_wUnsampledPhyla_tiplabels)
bac120_r220_wUnsampledPhyla_tiplabels
tree_with_all_unsampled_phyla_newLabels = tree_with_all_unsampled_phyla
matchnums_treetips_in_bac120_r220_wUnsampledPhyla_tiplabels = match_grepl(pattern=tree_with_all_unsampled_phyla$tip.label, x=bac120_r220_wUnsampledPhyla_tiplabels, return_counts=FALSE)
matchnums_treetips_in_bac120_r220_wUnsampledPhyla_tiplabels

tail(matchnums_treetips_in_bac120_r220_wUnsampledPhyla_tiplabels)
# 26 3227 3228 3229 3230 3231
# Everything is in order

bac120_r220_wUnsampledPhyla_tiplabels_table[order(bac120_r220_wUnsampledPhyla_tiplabels_table$phylum),]


# Exclude the "sp" (genomes of undescribed species)
species2_starts_with_sp_TF = startsWith(x=bac120_r220_wUnsampledPhyla_tiplabels_table$species2, prefix="sp")
sort(bac120_r220_wUnsampledPhyla_tiplabels_table$species2[species2_starts_with_sp_TF==FALSE])
sum(species2_starts_with_sp_TF)
sum(species2_starts_with_sp_TF==FALSE) # 388 hits on non-sp starts
length(species2_starts_with_sp_TF)



printall(matrix(sort(bac120_r220_wUnsampledPhyla_tiplabels_table$species2[species2_starts_with_sp_TF==TRUE]),ncol=1))


sp_to_add_back = c("spiroformis", "spiroformisspiroformis_Aspumans", "spiroformis", "spiroformis_A", "spumans")

sp_to_add_back_TF = bac120_r220_wUnsampledPhyla_tiplabels_table$species2 %in% sp_to_add_back
species2_starts_with_sp_TF[sp_to_add_back_TF == TRUE] = FALSE
sum(species2_starts_with_sp_TF==FALSE) # 388 is now 391

rev(sort(table(bac120_r220_wUnsampledPhyla_tiplabels_table$phylum[species2_starts_with_sp_TF==FALSE])))

rev(sort(table(bac120_r220_wUnsampledPhyla_tiplabels_table$phylum[species2_starts_with_sp_TF==TRUE])))


# Write out the table of described genomes in phyla we haven't sampled yet
table_of_described_genomes_in_unsampled_phyla = bac120_r220_wUnsampledPhyla_tiplabels_table[species2_starts_with_sp_TF==FALSE,]
dim(table_of_described_genomes_in_unsampled_phyla)

printall(table_of_described_genomes_in_unsampled_phyla[order(table_of_described_genomes_in_unsampled_phyla$phylum),])

outfn = "table_of_described_genomes_in_unsampled_phyla.txt"
write.table(x=table_of_described_genomes_in_unsampled_phyla[order(table_of_described_genomes_in_unsampled_phyla$phylum),], file=outfn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

