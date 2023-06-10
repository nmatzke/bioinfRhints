
#######################################################
# Produce various renamings & reorderings
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches
library(openxlsx)			# for openxlsx::read.xlsx

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R


wd = "~/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/"
setwd(wd)

prot_feature_tables_all_fn = "~/Downloads/2023-06-02_prot_feature_tables_all_v1.txt"
genomes_to_spnames_fn = "species_list_10062023_NJM+group+spname_v1.txt"

# Takes a few seconds
prot_feature_tables_all_df = read.table(prot_feature_tables_all_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
head(prot_feature_tables_all_df)

genomes_to_spnames_df = read.table(genomes_to_spnames_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
head(genomes_to_spnames_df)



# Order an alignment against a tree
trfn = "motA_hmmalign589_full_tr2_order_domainsOrder_cutNonHom_ed1_MiddleConstrained3.fasta.treefile_FigTree.nexus"

alnfn = "motA_hmmalign589_full_tr2_order_domainsOrder_cutNonHom_ed1_MiddleConstrained3.fasta"

tr = read.nexus(trfn)
tip_gids = gsub(pattern="'", replacement="", x=tr$tip.label)
head(tip_gids)

aln = seqinr::read.fasta(alnfn, seqtype="AA")
fullnames = fullnames_from_readFasta(aln)
seq_gids = names(aln)
head(fullnames)
head(seq_gids)

TF = tip_gids %in% seq_gids
sum(TF)
length(TF)

TF = seq_gids %in% tip_gids
sum(TF)
length(TF)

