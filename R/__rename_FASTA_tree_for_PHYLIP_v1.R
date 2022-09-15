
library(XLConnect)
library(ape)
library(seqinr)


# Set working directory (wd)
wd = "/drives/GDrive/__GDrive_projects/2017-08-08_Matt_Baker_flagellum_speaker/01_quicktree_MotB_n2v1/04_rename_v2/"
setwd(wd)

# Load some utility functions
source("http://ib.berkeley.edu/courses/ib200b/scripts/_R_tree_functions_v1.R")
source("parse_sequences_v1.R"



# Read in the Excel file
# Excel file from Matt Baker
xlsfn = "uniprot-yourlist_M201904086746803381A1F0E0DB47453E0216320D57D3F73_repaired_v2.xlsx"
xls = readWorksheetFromFile(file=xlsfn, sheet=1, startRow=1, startCol=1, endCol=8)

head(xls)
dim(xls)


# Check if all the names match the PHYLIP data
library(seqinr)


# Shortnames for FASTA and phylo

# Read in the quicktree tree
trfn = "MotB_757seqs_aln2_kimura_bootstraps2_reroot_FigTree.newick"
tr = read.tree(trfn)
length(tr$tip.label)


# Read in the FASTA alignment file
fasta_fn = "MotB_757seqs_aln.fasta"
seqs = seqinr::read.fasta(fasta_fn)
seqs3 = seqs

names_table = NULL

# Change the names in both
for (i in 1:length(seqs))
	{
	tmpi = sprintf("%04.0f", i)
	
	orig_name = names(seqs)[i]
	words = strsplit(orig_name, split="\\|")[[1]]
	words
	
	if (length(words) == 3)
		{
		tmp = strsplit(words[3], split="_")[[1]]
		if (length(tmp) >= 2)
			{
			tmp2 = tmp[2]
			} else {
			tmp2 = tmp[1]
			}
		
		newname = paste0(tmpi, "|", tmp2, collapse="|")
		}

	if (length(words) == 2)
		{
		tmp = words[2]
			
		newname = paste0(tmpi, "|", tmp, collapse="|")
		}

	
	names(seqs3)[i] = newname
	#newname = paste0(words[1], words[2], collapse="|")	
	
	
	
	# Find the same name in the tree
	TF = tr$tip.label == orig_name
	if (sum(TF) != 1)
		{
		txt = paste0("Error at i=", i, ", ", sum(TF), " sequences found. orig_name was: ", orig_name)
		print(txt)
		stop(txt)
		}
	
	# Apply the new name to tr3
	tr3$tip.label[TF] = newname
	
	tmpline = c(orig_name, newname)
	names_table = rbind(names_table, tmpline)
	}

names_table_df = as.data.frame(names_table, stringsAsFactors=FALSE)
names(names_table_df) = c("orig_name", "10char_name")
write.table(names_table_df, file="names_table_convert_to_10char_v1.txt")


# Check that everything matches
sum(sort(names(seqs3)) == sort(tr3$tip.label))
cbind(sort(names(seqs3)), sort(tr3$tip.label))

# Write to files
write.fasta(sequences=seqs3, names=names(seqs3), file.out="MotB_757seqs_aln_10char.fasta")
write.tree(tr3, file="MotB_757seqs_aln2_kimura_bootstraps2_reroot_FigTree_10char.newick")
write.nexus(tr3, file="MotB_757seqs_aln2_kimura_bootstraps2_reroot_FigTree_10char.nexus")

# Write to PHYLIP
seq3_strings = collapse_seqs(dataset=seqs3)
seq3_dtf = as.data.frame(seq3_strings, stringsAsFactors=FALSE)

dtf_to_phylip(char_dtf=seq3_dtf, outfn="MotB_757seqs_aln_10char.phylip")


