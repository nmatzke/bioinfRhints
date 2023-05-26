

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



# Set working directory
# wd = "/GitHub/bioinfRhints/flag/gene_order_example/" # example
wd = "~/Downloads/Full_genomes/"	 # full dataset
setwd(wd)

# Get the list of genome directories that have been saved somewhere
genome_dirs = list.files(path="genomes", pattern=NULL, recursive=FALSE)
genome_dirs

# Name & load the Excel file containing metadata (e.g. genome filenames, species name etc.)
xlsfn = "/GitHub/bioinfRhints/flag/gene_order_example/species_list_14102022_1page.xlsx"
xls = read.xls(xlsfn)
head(xls)
tail(xls)

# Name and load an MotA alignment file
alnfn = "/GitHub/bioinfRhints/flag/gene_order_example/motA_alignment_refined.fasta"
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
aln_names = unname(sapply(X=aln, FUN=attr, "name"))
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
trfn = "/GitHub/bioinfRhints/flag/gene_order_example/trees/530_sequences_Alignment_contree_reRootLadder_gIDs.newick"
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
xlsfn = "/GitHub/bioinfRhints/flag/gene_order_example/species_list_14102022_1page.xlsx"
xls = read.xls(xlsfn)
head(xls)
tail(xls)


# find a protein in a genome

protID = "AAC74960.1"
genome_dir = get_genome_name_from_protID(protID, genome_dirs, aln_annotations, xls, xlstaxa)



# Access that directory, look for protein in list
gene_order_table_zipfn = slashslash(paste0(wd, "/", "genomes/", genome_dir, "/", genome_dir, "_feature_table.txt.gz"))
cmdstr = paste0("gunzip ", gene_order_table_zipfn)
system(cmdstr)
gene_order_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_feature_table.txt")

# Get the gene order table 

read_gene_order_table <- function(gene_order_table_fn)
	{
	# Remove "#" from first line
	tmplines = readLines(gene_order_table_fn)
	tmplines[1] = gsub(pattern="# ", replacement="", x=tmplines[1])
	writeLines(tmplines, con=gene_order_table_fn)

	# NOTE: quote="\"" is necessary to avoid "EOF within quoted string"
	#       ...which causes issues with e.g. "2',3'-cyclic phosphodiesterase"
	gene_order_tmp = read.table(gene_order_table_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
	
	# ERROR CHECK
	if (length(tmplines)-1 != nrow(gene_order_tmp))
		{
		txt = paste0("STOP ERROR in read_gene_order_table(", gene_order_table_fn, "): the file has ", length(tmplines), " lines, but read.table() gave only ", nrow(gene_order_tmp), " lines.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}

	# Take every 2nd line
	odd_lines = seq(1, nrow(gene_order_tmp), by=2)
	even_lines = seq(2, nrow(gene_order_tmp), by=2)
	gene_order_df1 = gene_order_tmp[odd_lines,]
	gene_order_df2 = gene_order_tmp[even_lines,]

	gene_order_df1[1987:1992,]
	gene_order_df2[1987:1992,]

	TF = gene_order_df2$symbol == "motB"
	gene_order_df1[TF,]
	gene_order_df2[TF,]

	rbind(gene_order_df1[TF,], gene_order_df2[TF,])


	# This seems to have everything...
	gene_order_df = gene_order_df2
	head(gene_order_df)
	return(gene_order_df)
	}

# Find the MotA protein
gene_order_df = read_gene_order_table(gene_order_table_fn)

gene_num = grep(pattern=protID, x=gene_order_df$product_accession, ignore.case=TRUE)
gene_order_df[gene_num,]

gene_order_df[gene_num+1,]
gene_order_df[gene_num-1,]


list_of_protIDs


nones_to_NA <- function(x)
	{
	if (length(x) == 0)
		return(NA)
	end
	if (is.na(x))
		return(NA)
	else
		return(x)
	end
	return(x)
	}

genome_names_not_found = NULL
protIDs_not_found = NULL
list_of_protIDs = aln_names 
gene_neighbors = NULL
i=1
for (i in 1:length(list_of_protIDs))
	{
	cat("\ni=", i, ",", sep="")
	protID = list_of_protIDs[i]
	both = get_genome_name_from_protID(protID, genome_dirs, aln_annotations, xls, xlstaxa, returnwhat="both")
	if (is.null(both$genome_dir)) # nothing found
		{
		protIDs_not_found = c(protIDs_not_found, protID)
		genome_names_not_found = c(genome_names_not_found, NA)
		next()
		}
	
	spname = both$spname
	genome_dir = both$genome_dir

	# Access that directory, look for protein in list
	gene_order_table_zipfn = slashslash(paste0(wd, "/", "genomes/", genome_dir, "/", genome_dir, "_feature_table.txt.gz"))
	cmdstr = paste0("gunzip ", gene_order_table_zipfn)
	system(cmdstr)
	gene_order_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_feature_table.txt")

	gene_order_df = read_gene_order_table(gene_order_table_fn)

	gene_num = grep(pattern=protID, x=gene_order_df$product_accession, ignore.case=TRUE)
	if (length(gene_num) == 0)
		{
		protIDs_not_found = c(protIDs_not_found, protID)
		genome_names_not_found = c(genome_names_not_found, genome_dir)

		next()
		}
	
	gene_order_df[gene_num,]

	gene_order_df[gene_num+1,]
	gene_order_df[gene_num-1,]
	
	
	# Assemble 2 lines: symbols and accessions
	# symbol line
	strand0 = nones_to_NA(gene_order_df$strand[gene_num])
	sym0 = nones_to_NA(gene_order_df$symbol[gene_num])
	acc0 = nones_to_NA(gene_order_df$product_accession[gene_num])
	name0 = nones_to_NA(gene_order_df$name[gene_num])
	
	strand1 = nones_to_NA(gene_order_df$strand[gene_num+1])
	sym1 = nones_to_NA(gene_order_df$symbol[gene_num+1])
	acc1 = nones_to_NA(gene_order_df$product_accession[gene_num+1])
	name1 = nones_to_NA(gene_order_df$name[gene_num+1])

	strand2 = nones_to_NA(gene_order_df$strand[gene_num+2])
	sym2 = nones_to_NA(gene_order_df$symbol[gene_num+2])
	acc2 = nones_to_NA(gene_order_df$product_accession[gene_num+2])
	name2 = nones_to_NA(gene_order_df$name[gene_num+2])

	strand3 = nones_to_NA(gene_order_df$strand[gene_num+3])
	sym3 = nones_to_NA(gene_order_df$symbol[gene_num+3])
	acc3 = nones_to_NA(gene_order_df$product_accession[gene_num+3])
	name3 = nones_to_NA(gene_order_df$name[gene_num+3])

	strandM1 = nones_to_NA(gene_order_df$strand[gene_num-1])
	symM1 = nones_to_NA(gene_order_df$symbol[gene_num-1])
	accM1 = nones_to_NA(gene_order_df$product_accession[gene_num-1])
	nameM1 = nones_to_NA(gene_order_df$name[gene_num-1])

	strandM2 = nones_to_NA(gene_order_df$strand[gene_num-2])
	symM2 = nones_to_NA(gene_order_df$symbol[gene_num-2])
	accM2 = nones_to_NA(gene_order_df$product_accession[gene_num-2])
	nameM2 = nones_to_NA(gene_order_df$name[gene_num-2])

	strandM3 = nones_to_NA(gene_order_df$strand[gene_num-3])
	symM3 = nones_to_NA(gene_order_df$symbol[gene_num-3])
	accM3 = nones_to_NA(gene_order_df$product_accession[gene_num-3])
	nameM3 = nones_to_NA(gene_order_df$name[gene_num-3])
	
	if (strand0 == "+")
		{
		tmprow = c(i, protID, spname, symM3, symM2, symM1, sym0, sym1, sym2, sym3, strandM3, strandM2, strandM1, strand0, strand1, strand2, strand3, accM3, accM2, accM1, acc0, acc1, acc2, acc3, nameM3, nameM2, nameM1, name0, name1, name2, name3, genome_dir, gene_order_table_fn)
		} else if (strand0 == "-") {
		tmprow = c(i, protID, spname, sym3, sym2, sym1, sym0, symM1, symM2, symM3, strand3, strand2, strand1, strand0, strandM1, strandM2, strandM3, acc3, acc2, acc1, acc0, accM1, accM2, accM3, name3, name2, name1, name0, nameM1, nameM2, nameM3, genome_dir, gene_order_table_fn)		
		}
	
	gene_neighbors = rbind(gene_neighbors, tmprow)
	}


gene_neighbors_df = as.data.frame(gene_neighbors, stringsAsFactors=FALSE)
names(gene_neighbors_df) = c("i", "protID", "spname", "symM3", "symM2", "symM1", "sym0", "sym1", "sym2", "sym3", "strandM3", "strandM2", "strandM1", "strand0", "strand1", "strand2", "strand3", "accM3", "accM2", "accM1", "acc0", "acc1", "acc2", "acc3", "nameM3", "nameM2", "nameM1", "name0", "name1", "name2", "name3", "genome_dir", "gene_order_table_fn")
row.names(gene_neighbors_df) = NULL

gene_neighbors_df

outfn = "gene_neighbors_v1.txt"
write.table(gene_neighbors_df, file=outfn, sep="\t")


cbind(protIDs_not_found, genome_names_not_found)




