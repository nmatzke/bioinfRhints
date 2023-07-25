
#######################################################
# Setup
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches
library(seqinr) 			# for read.fasta

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R

wd = "~/Downloads/iqtree_genomes_wLitLinks/"
setwd(wd)



#######################################################
# Read translation table
# generated at: ~/Downloads/z_genomes_wLitLinks_processing/_cmds_unzip_cat_HMMER_v1.txt
#######################################################
translation_seqsfn = "AQBs_orig+flag1-5_wLitLinks.fasta"
infn = "AQBs_orig+flag1-5_wLitLinks_table.txt"
translate_df = read.table(file=infn, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
head(translate_df)
tail(translate_df)

# Alignment used for tree
alnfn = "1214_seqs_merged_mafftMiddleConstrained2.fasta"
aln = seqinr::read.fasta(alnfn, seqtype="AA")
fullnames = fullnames_from_readFasta(aln)
aln_gids = names(aln)
length(aln_gids)


# Tree
trfn = "1214_seqs_merged_mafftMiddleConstrained2.fasta.treefile"
tr = read.tree(trfn)
tr2 = phytools::midpoint.root(tr)
tr3 = ladderize(phy=tr2, right=TRUE)
plot(tr3, show.tip.label=FALSE)
tr3


# Parse tip labels
tipnames = tr3$tip.label
head(tipnames)
tipnames3 = tipnames

for (i in 1:length(tipnames))
	{
	words = strsplit(x=tipnames[i], split="\\|")[[1]]
	tipnames3[i] = words[2]
	}

tipnames3


# Match to table
matches1 = match(x=tipnames3, table=translate_df$gids)
matches2 = match(x=translate_df$gids, table=tipnames3)
TF = is.na(matches1)
sum(TF)
tipnames3[TF]

TF = is.na(matches2)
sum(TF)
translate_df[TF,]

dim(translate_df)
length(tipnames3)