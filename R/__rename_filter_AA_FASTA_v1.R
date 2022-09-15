
library(XLConnect)
library(ape)
library(seqinr)
library(gdata) # for "trim"

# Set working directory (wd)
wd = "/drives/GDrive/__GDrive_projects/2017-08-08_Matt_Baker_flagellum_speaker/01_data/FliI/"
setwd(wd)

# Load some utility functions
source("http://ib.berkeley.edu/courses/ib200b/scripts/_R_tree_functions_v1.R")
source("/GitHub/BEASTmasteR/R/parse_sequences_v1.R")


# Read in the Excel file
# Excel file from Matt Baker
# xlsfn = "uniprot-yourlist_M201904086746803381A1F0E0DB47453E0216320D57D3F73_repaired_v2.xlsx"
# xls = readWorksheetFromFile(file=xlsfn, sheet=1, startRow=1, startCol=1, endCol=8)
# 
# head(xls)
# dim(xls)


# Shortnames for FASTA and phylo

# Read in the quicktree tree
# trfn = "FliI_944seqs_aln2_kimura_bootstraps2_reroot_FigTree.newick"
# tr = read.tree(trfn)
# length(tr$tip.label)


# Read in the FASTA alignment file
seqname = "FliI"
fasta_fn = "FliI_1001Sequences.fasta"
seqs = seqinr::read.fasta(fasta_fn)
seqs3 = seqs

names_table = NULL


#######################################################
# Get all original names
#######################################################
orig_names = rep("", times=length(seqs))
for (i in 1:length(seqs))
	{
	orig_names[i] = attr(seqs[[i]], which="Annot")
	}
	
TF1 = grepl(pattern="MotA", x=orig_names, ignore.case=FALSE)
TF2 = grepl(pattern="motA", x=orig_names, ignore.case=FALSE)
TF3 = grepl(pattern="Motility protein A", x=orig_names, ignore.case=FALSE)
TF4 = grepl(pattern="Motility protein A", x=orig_names, ignore.case=FALSE)
TF = (TF1+TF2+TF3+TF4) > 0
orig_names[TF==FALSE]

# Which are fragments
fragTF = grepl(pattern="ragment", x=orig_names, ignore.case=FALSE)
orig_names[fragTF==TRUE]
fragSeqs_nums = (1:length(seqs))[fragTF]
seqs[fragSeqs_nums]
lengths_of_fragments = sapply(seqs[fragSeqs_nums], FUN=length)


#######################################################
# Get the distribution of lengths
#######################################################
seqlengths = sapply(seqs, FUN=length)
minx = min(pretty(seqlengths))
maxx = max(pretty(seqlengths))
titletxt = paste0("Distribution of lengths of ", seqname, " sequences")
h1 = hist(seqlengths, breaks=50)
h2 = hist(lengths_of_fragments, breaks=h1$breaks)
plot(h1, col=rgb(0,0,1,1/4), xlim=c(minx,maxx), main=titletxt)
plot(h2, col=rgb(1,0,0,1/4), xlim=c(minx,maxx), add=TRUE)
legend(x=minx, y=max(h1$counts), legend=c("all", "fragments", "both"), fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), rgb(1,0,1,1/4)), cex=0.85)

# See the shortest/longest ones
head(sort(sapply(seqs, FUN=length)), 10)
head(rev(sort(sapply(seqs, FUN=length))), 10)







# Cut this as a duplicate, and too long:
name_to_cut = "tr|A0A3M3KG78|A0A3M3KG78_SERPL Uncharacterized protein OS=Serratia plymuthica OX=82996 GN=ALQ63_01277 PE=4 SV=1"
cutTF = orig_names == name_to_cut
seqs2 = seqs[cutTF==FALSE]

# Cut seqs < 425 length:
cutTF1 = sapply(seqs2, FUN=length) < 425
cutTF2 = sapply(seqs2, FUN=length) > 575
cutTF = (cutTF1+cutTF) > 0
seqs2 = seqs2[cutTF==FALSE]

length(seqs)
length(seqs2)
seqs = seqs2


uncharTF = grepl(pattern="nchar", x=orig_names, ignore.case=FALSE)
orig_names[uncharTF==TRUE]


# Change the names
names_table = matrix(data="", ncol=4, nrow=length(seqs))
i=1


for (i in 1:length(seqs))
	{
	tmpi = sprintf("%04.0f", i)
	
	orig_name = attr(seqs[[i]], which="Annot")
	words = strsplit(orig_name, split="\\|")[[1]]
	words
	
	# Split the 3rd term to get geneName and speciesName
	
	# Get OS= species name
	term3 = words[[3]]
	match_OSeq = c(gregexpr(pattern="OS=", text=term3, ignore.case=FALSE))
	start_OSeq = c(match_OSeq[[1]])
	spName_startlength = attr(match_OSeq[[1]], which="match.length")
	start_spName = start_OSeq + spName_startlength
	
	match_OXeq = c(gregexpr(pattern="OX=", text=term3, ignore.case=FALSE))
	end_OSeq = c(match_OXeq[[1]]) - 1
	
	spName = trim(substr(x=term3, start=start_spName, stop=end_OSeq))
	gnName = trim(strsplit(term3, split=" ")[[1]][1])
		
	spName2 = convert_chars_to_underscores(txt=spName)
	
	shortname = paste0("s", tmpi, collapse="")
	medname = paste0(spName2, "|", gnName, collapse="")
	longname = paste0("s", tmpi, "|", spName2, "|", gnName, collapse="")
	shortname
	medname
	longname
	
	names_table[i,] = c(shortname, medname, longname, orig_name)
	}

names_table_df = as.data.frame(names_table, stringsAsFactors=FALSE)
row.names(names_table_df) = NULL
names(names_table_df) = c("short", "medium", "long", "orig")
head(names_table_df)

names(seqs) = names_table_df$medium

# Write to files
write.fasta(sequences=seqs, names=names(seqs), file.out="FliI_944_medNames.fasta")
write.table(names_table_df, file="names_table_FliI_v1.txt")

# Write to PHYLIP
seq_strings = collapse_seqs(dataset=seqs)
seq_dtf = as.data.frame(seq_strings, stringsAsFactors=FALSE)
dtf_to_phylip(char_dtf=seq_dtf, outfn="FliI_974_medNames.phylip")




names_table_df = as.data.frame(names_table, stringsAsFactors=FALSE)
names(names_table_df) = c("orig_name", "10char_name")
write.table(names_table_df, file="names_table_convert_to_10char_v1.txt")



# Write to files
write.fasta(sequences=seqs3, names=names(seqs3), file.out="FliI_944seqs_aln_10char.fasta")
write.tree(tr3, file="FliI_944seqs_aln2_kimura_bootstraps2_reroot_FigTree_10char.newick")
write.nexus(tr3, file="FliI_944seqs_aln2_kimura_bootstraps2_reroot_FigTree_10char.nexus")

# Write to PHYLIP
seq3_strings = collapse_seqs(dataset=seqs3)
seq3_dtf = as.data.frame(seq3_strings, stringsAsFactors=FALSE)

dtf_to_phylip(char_dtf=seq3_dtf, outfn="FliI_944seqs_aln_10char.phylip")


