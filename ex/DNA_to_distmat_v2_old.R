
#######################################################
# 2022-09-15
#
# 1. Read in a sequence alignment
# 2. Construct a distance matrix
#
#######################################################

# Loading a tree, printing to tree table (trtable)
library(ape)


# Data: 18S, ITS1, 5.8S rRNA sequences from
# X-cells (protist parasites of e.g. Antarctic fishes)
# Provided by Craig Miller
#
# Request: distance matrix

seqs_fn = "/GitHub/bioinfRhints/data/seqs.fasta"

# Same effect -- produces an APE "DNAbin" object (DNA in binary format for compactness/speed)
seqs = read.FASTA(file=seqs_fn, type="DNA")
seqs = read.dna(file=seqs_fn, format="fasta", skip=0, as.character=FALSE)
# 69 DNA sequences in binary format stored in a matrix.
# 
# All sequences of same length: 2033 
# 
# Labels:
# MT299786.1 Parvilucifera catillosa strain Kokar2016a isolate...
# KX519761.1 Parvilucifera corolla strain Lanzarote isolate 3 ...
# KF359483.1 Parvilucifera rostrata strain RCC2800 18S ribosom...
# EU502912.1 Parvilucifera sinerae
# KF359485.1 Parvilucifera infectans strain RCC2816 18S riboso...
# AF133909.1 Parvilucifera infectans 18S ribosomal RNA gene, p...
# ...
# 
# Base composition:
#     a     c     g     t 
# 0.269 0.198 0.259 0.275 
# (Total: 140.28 kb)

# What's inside of 'seqs'?

# Gene/locus names, method #1:
rownames(seqs)

# Gene/locus names, method #2:
seqs_rownames = attr(seqs, "dimnames")[[1]]
seqs_rownames

# What all is in a DNAbin object?
names(attributes(seqs))
# "dim"      "dimnames" "class"   

# Alternate methods:
class(seqs)
attr(seqs, "class")

dim(seqs)
attr(seqs, "dim")
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/696/305/GCF_001696305.1_UCN72.1 my_dir/


# Each base as an individual string
seqs_txt = read.dna(file=seqs_fn, format="fasta", skip=0, as.character=TRUE)
dim(seqs_txt)
# 69 2033
seqs_txt[1,]
# Paste a sequence into a single string
first_seq = paste0(seqs_txt[1,], collapse="")
first_seq
#square brakets to specify row and column of matrix, always rows first then columns
# Subset of text sequence
substr(first_seq, start=100, stop=110)

#to assess site variation by counting number of types of nucleotides in each column 

#to do a loop
for (i in 1:100
{
print(i)
}

for (i in 1:100)
{
print(unique(seqs_txt[,i]))
}
#create a vector
unique_nucleotide_counts = rep(0,times = 2033)
#how many columns are variable and how many are not
for (i in 1:2033)
{
ur = unique(seqs_txt[,i])
TF = ur != "-"
number_of_nucleotides = sum(TF)
print(number_of_nucleotides)
unique_nucleotide_counts[i] = number_of_nucleotides

}
sum(unique_nucleotide_counts>1)
# I have written some functions to process/convert DNAbin objects
source("/GitHub/bioinfRhints/R/_R_tree_functions_v2.R")
seqs_as_characters = DNAbin_to_list_of_strings(tmpdata=seqs)
seqs_as_characters[[1]]

# Collapse all sequences to strings
seqs_as_list_of_strings = lapply(X=seqs_as_characters, FUN=paste0, collapse="")
substr(seqs_as_list_of_strings[[1]], start=100, stop=110)




#######################################################
# OK, make the distance matrix
#######################################################
distmat1 = dist.dna(x=seqs, model="K80", gamma=FALSE, as.matrix=TRUE)
dim(distmat1)
distmat1[1:5,1:5]


library(BioGeoBEARS)
library(seqinr)










trstr="((Human:5,Chimp:5):1,Gorilla:6);"
tr = read.tree(file="", text=trstr)
trtable = prt(tr, printflag=FALSE)
trtable
trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)
trtable



# Setup
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)			# for trim
library(pegas)
library(XML)
library(parallel)
library(S4Vectors)
library(stats4)
library(IRanges)
library(XVector)
library(Biostrings)
library(BiocGenerics)
library(annotate)		# for Bioconductor's 
						#"genbank" (efetch by 
						# gene Accession #s)
						# "blastSequences"
#library(reutils)  # for efetch
# EXPIRED, use reutils Use genomes::efetch to download genes
library(genomes)

# rBLAST
library(rBLAST)


# Source this function:
# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html
source("~/Downloads/bioinf/Rcode/blastR_setup/blastsequences_v3.R")

#wd = "/drives/SkyDrive/NIMBioS_projects/2015-01-16_Laura_Saila/data/doggies_combine/Prevosti_2009_seqs_PLUS_NEW/"
#wd = "~/Downloads/bioinf/01_seqs/01_seqs/"
#wd = "~/Downloads/bioinf/scratch/"

setwd(wd)
wd = "/Downloads/bioinf/01_seqs/"

# Sequences to BLAST against:
#fasta_fn = "Eleginops_maclovinus_seqs_to_BLAST_vs_mtDNA_dbs_v1.fasta"
#fasta_fn = "ND6_complete_Histiodraco_velifer_GU214225.1.fasta"
#seqs = read.fasta(file=fasta_fn)
