
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


aln_to_tip_order = match(x=tip_gids, table=seq_gids)
aln2 = rev(aln[aln_to_tip_order])
head(fullnames_from_readFasta(aln2))

head(rev(tr$tip.label))


# Find tip_gids in prot_features_table
tip_gids_in_df_matches = match(tip_gids, table=prot_feature_tables_all_df$product_accession)
tip_gids_in_df_matches

nonmatches = tip_gids[is.na(tip_gids_in_df_matches)]
tmplist = c()
for (i in 1:length(nonmatches))
	{
	num = grep(pattern=nonmatches[i], x=fullnames)
	cat("\n")
	cat(fullnames[num])
	tmplist = c(tmplist, fullnames[num])
	}
cat("\n")

sort(unique(extract_last_brackets(tmplist)))


# Myxococcus xanthus DZ2A
id="GCA_020827275.2"
TF = prot_feature_tables_all_df$assembly == id
prot_feature_tables_all_df[TF,]


