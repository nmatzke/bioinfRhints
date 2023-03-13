

#######################################################
# 2023-03-13 Lab meeting
# Scripts: manipulating alignments, sequence names / tip names, etc.
#          Extracting gene order from a list of downloaded genomes
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches

# Handy function
extract_last_brackets <- function(list_of_strings, replace_spaces=TRUE)
	{
	species_names = rep("", length(list_of_strings))
	txt = paste0("\nextract_last_brackets() is processing ", length(list_of_strings), " strings. String #")
	cat(txt)
	for (i in 1:length(list_of_strings))
		{
		cat(i, ",", sep="")
		tmptxt = list_of_strings[i]
		tmpstrs = unlist(regmatches(tmptxt, gregexpr("\\[.+?\\]", tmptxt)))
		tmpstr = tmpstrs[length(tmpstrs)] # take the last bracketed text, if more than 1
		
		# Remove "[", "]"
		tmpstr = gsub(pattern="\\[", replacement="", x=tmpstr)
		tmpstr = gsub(pattern="\\]", replacement="", x=tmpstr)
		
		# Replace spaces with "_"
		if (replace_spaces == TRUE)
			{
			tmpstr = gsub(pattern=" ", replacement="_", x=tmpstr)
			}
		species_names[i] = tmpstr
		}
	return(species_names)
	}




# Set working directory
wd = "/GitHub/bioinfRhints/flag/gene_order_example/"
setwd(wd)

# Get the list of genome directories that have been saved somewhere
genome_dirs = list.files(path="genomes", pattern=NULL, recursive=FALSE)
genome_dirs

# Name & load the Excel file containing metadata (e.g. genome filenames, species name etc.)
xlsfn = "species_list_14102022_1page.xlsx"
xls = read.xls(xlsfn)
head(xls)
tail(xls)

# Name and load an MotA alignment file
alnfn = "motA_alignment_refined.fasta"
aln = read.fasta(alnfn)

# "aln" is an R list, each element is a sequence, with a name and some other attributes:
attributes(aln[[1]])

# Number of sequences
length(aln)

# Length of an individual sequence
length(aln[[1]])
length(aln[[2]])

# Length of all the sequences (should be the same, as we have an alignment rather than
# just a list of sequences)
lengths = sapply(X=aln, FUN=length)
lengths

# Are all sequences the same length?
unique(lengths)
all(lengths == unique(lengths))

# For kicks, calculate the percentage gaps in the alignment
numgaps = rep(0.0, times=length(aln))
numdata = rep(0.0, times=length(aln))

for (i in 1:length(aln))
	{
	# Number of gaps in a single sequence
	numgaps[i] = sum(aln[[i]] == "-")
	# Number of non-gaps in a single sequence
	numdata[i] = sum(aln[[i]] != "\\-")
	} # END for (i in 1:length(aln))

# Total number of "-"
sum(numgaps)
sum(numdata)
sum(numgaps) + sum(numdata)
percent_data = sum(numdata) / (sum(numgaps) + sum(numdata)) * 100
percent_data


# Subset alignment to just the columns for which E. coli K-12 MotA has data
# MotA|Escherichia_coli_str._K-12_substr._MG1655|AAC74960.1
aln_names = sapply(X=aln, FUN=attr, "name")
TF = grepl(pattern="AAC74960.1", x=aln_names)
index_num = (1:length(aln))[TF]
aln_names[index_num]
aln[[index_num]]

# Positions to keep from E. coli K-12 MotA
keepTF = aln[[index_num]] != "-"

# Make a copy of the alignment
aln2 = aln
for (i in 1:length(aln))
	{
	aln2[[i]] = aln[[i]][keepTF]
	}
aln2[[1]]
aln2[[index_num]]

# Write out to new fasta file
aln_annotations = sapply(X=aln, FUN=attr, "Annot")
aln_annotations = gsub(pattern="\\>", replacement="", x=aln_annotations) # remove ">" as write.fasta adds ">"
alnfn2 = gsub(pattern=".fasta", replacement="_subset_to_Ecoli_K12.fasta", x=alnfn)
write.fasta(sequences=aln2, names=aln_annotations, file.out=alnfn2)


# Name and load a IQ tree file
trfn = "trees/530_sequences_Alignment_contree_reRootLadder_gIDs.newick"
tr = read.tree(trfn)
tr

# Number of tips
length(tr$tip.label)

# Tree table
trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=2000)
head(trtable)




#######################################################
# Gene order
#######################################################
# Get the list of genome directories that have been saved somewhere
genome_dirs = list.files(path="genomes", pattern=NULL, recursive=FALSE)
genome_dirs

# Get the genome taxon
xlstaxa = rep("", times=nrow(xls))
for (i in 1:nrow(xls))
	{
	tmp = paste(xls$Genus[i], xls$Species[i], xls$Strain[i], sep=" ", collapse="")
	xlstaxa[i] = trim(gsub(pattern="  ", replacement=" ", x=tmp))
	}
xlstaxa

# Name & load the Excel file containing metadata (e.g. genome filenames, species name etc.)
xlsfn = "species_list_14102022_1page.xlsx"
xls = read.xls(xlsfn)
head(xls)
tail(xls)


# find a protein in a genome

genbank_ID = "AAC74960.1"
TF = grepl(pattern=genbank_ID, x=aln_annotations)
match_index = (1:length(aln_annotations))[TF]
match_index

seqname = aln_annotations[match_index]

spname = extract_last_brackets(seqname, replace_spaces=FALSE)
spname  

# Find the species name in the genomes
#match_in_genomes_list = match(x=spname, xlstaxa)

# Approximate String Matching (Fuzzy Matching)
# https://astrostatistics.psu.edu/su07/R/html/base/html/agrep.html
match_in_genomes_list = agrep(pattern=spname, x=xlstaxa, ignore.case=TRUE)

# Genome ID
genomeID_in_xls = xls$GenBank.ID[match_in_genomes_list]
genome_dir_index = agrep(pattern=genomeID_in_xls, x=genome_dirs, ignore.case=TRUE)
genome_dir_index

genome_dir = genome_dirs[genome_dir_index]
genome_dir


# Access that directory, look for protein in list
gene_order_table_zipfn = slashslash(paste0(wd, "/", "genomes/", genome_dir, "/", genome_dir, "_feature_table.txt.gz"))
cmdstr = paste0("gunzip ", gene_order_table_zipfn)
system(cmdstr)

gene_order_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_feature_table.txt")

# Remove "#" from first line
tmplines = readLines(gene_order_table_fn)
tmplines[1] = gsub(pattern="# ", replacement="", x=tmplines[1])
writeLines(tmplines, con=gene_order_table_fn)
gene_order_tmp = read.table(gene_order_table_fn, header=TRUE, comment.char="%", sep="\t", fill=TRUE, stringsAsFactors=FALSE)

# Take every 2nd line
even_lines = seq(2, nrow(gene_order_tmp), by=2)
gene_order_df = gene_order_tmp[even_lines,]
head(gene_order_df)

protein_in_df = match(x=genbank_ID, table=gene_order_df$product_accession)
protein_in_df

x = agrep(pattern=genbank_ID, x=gene_order_df$product_accession, ignore.case=TRUE)
gene_order_df[x,]


x = agrep(pattern=genbank_ID, x=gene_order_tmp$product_accession, ignore.case=TRUE)
gene_order_tmp[x,]
