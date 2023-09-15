#######################################################
# TASK 04:
# Retrieve gene order dataset for a list of Genbank IDs,
# gids, e.g. from tree tips
#
# * This one uses some functions (e.g. get_adjacent_row() )
#   that (should!) account for the fact that sometimes, the
#   adjacent gene is e.g. a pseudogene that takes up just 
#   one row in the *_feature_table.txt -- most things take
#   up two.  Rarely, I have seen other, weirder things, but
#   we will ignore those here.
# 
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches
library(openxlsx)			# for openxlsx::read.xlsx

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R
sourceall("/GitHub/bioinfRhints/R/BEASTmasteR/") 
# for remove_equals_from_tips()
# for read.beast.table_original


wd = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/"
setwd(wd)

prot_feature_tables_all_fn = "~/Downloads/2023-06-12_prot_feature_tables_all_v1.txt"
genomes_to_spnames_fn = "species_list_10062023_NJM+group+spname_v1.txt"

# Takes a few seconds
prot_feature_tables_all_df = read.table(prot_feature_tables_all_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
head(prot_feature_tables_all_df)

genomes_to_spnames_df = read.table(genomes_to_spnames_fn, header=TRUE, comment.char="%", quote="", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
genomes_to_spnames_df[1:10,]
dim(genomes_to_spnames_df)


# Order an alignment against a tree
trfn = "motA_hmmalign589_full_tr2_order_domainsOrder_cutNonHom_ed1_MiddleConstrained3.fasta.treefile_FigTree.nexus"
tr = read.nexus(trfn)



#######################################################
# Get the flanking genes, save to a tab-delimited table
#######################################################
sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R

list_of_protIDs = gsub(pattern="'", replacement="", x=tr$tip.label)
gene_neighbors_df = get_adjacent_genes(list_of_protIDs, prot_feature_tables_all_df, genomes_to_spnames_df, printwarnings=FALSE)

head(gene_neighbors_df)

outfn = "2023-06-15_gene_neighbors_from_571seqs_v1.txt"
write.table(gene_neighbors_df, file=outfn, sep="\t")



