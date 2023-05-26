
#######################################################
# Source code for bioinformatics tasks
#######################################################


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




##############
# JUNK: Example code inside 'junk'
##############
junk='
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(xlsx)					# for read.xlsx (replaces gdata::read.xls)
library(stringr) 			# for regmatches

sourceall("/GitHub/bioinfRhints/Rsrc/")

# Set working directory
# wd = "/GitHub/bioinfRhints/flag/gene_order_example/" # example
wd = "~/Downloads/Full_genomes/"	 # full dataset
setwd(wd)

protID = "AAC74960.1"

genome_dirs = list.files(path="genomes", pattern=NULL, recursive=FALSE)
genome_dirs

# Name and load an MotA alignment file
alnfn = "/GitHub/bioinfRhints/flag/gene_order_example/motA_alignment_refined.fasta"
aln = read.fasta(alnfn)
aln_annotations = sapply(X=aln, FUN=attr, "Annot")
aln_annotations = gsub(pattern=">>", replacement=">", x=aln_annotations) # remove ">" as write.fasta adds ">"
aln_annotations = gsub(pattern=">", replacement="", x=aln_annotations) # remove ">" as write.fasta adds ">"

# Name & load the Excel file containing metadata (e.g. genome filenames, species name etc.)
xlsfn = "/GitHub/bioinfRhints/flag/gene_order_example/species_list_14102022_1page.xlsx"
xls = xlsx::read.xlsx(file=xlsfn)

# Get the genome taxon
xlstaxa = rep("", times=nrow(xls))
for (i in 1:nrow(xls))
	{
	tmp = paste(xls$Genus[i], xls$Species[i], xls$Strain[i], sep=" ", collapse="")
	xlstaxa[i] = trim(gsub(pattern="  ", replacement=" ", x=tmp))
	}

protID = "AAC74960.1"
genome_dir = get_genome_name_from_protID(protID, genome_dirs, aln_annotations, xls, xlstaxa)

' # END JUNK

# find a protein in a genome
get_genome_name_from_protID <- function(protID, genome_dirs, aln_annotations, xls, xlstaxa, returnwhat="genome_dir")
	{
	TF = grepl(pattern=protID, x=aln_annotations)
	match_index = (1:length(aln_annotations))[TF]
	match_index

	seqname = aln_annotations[match_index]

	spname = extract_last_brackets(seqname, replace_spaces=FALSE)
	spname  

	# Find the species name in the genomes
	#match_in_genomes_list = match(x=spname, xlstaxa)

	# Approximate String Matching (Fuzzy Matching)
	# https://astrostatistics.psu.edu/su07/R/html/base/html/agrep.html
	
	# Error trap for "=" in name
	if (grepl(pattern="\\=", x=spname) == TRUE)
		{
		words = strsplit(spname, split="\\=")[[1]]
		spname = words[1]
		}
	#tmpname = gsub(pattern="\\=", replacement="EQ", x=spname)
	
	match_in_genomes_list = agrep(pattern=spname, x=xlstaxa, ignore.case=TRUE)
	
	if (length(match_in_genomes_list) > 1)
		{
		txt = paste0("Warning in get_genome_name_from_protID(protID='", protID, "'...): >1 genome matched spname='", spname, "'. Printing them below, taking the first one.")
		cat("\n\n")
		cat(txt)
		cat("\nxlstaxa[match_in_genomes_list]:\n")
		print(xlstaxa[match_in_genomes_list])
		cat("\n")
		warning(txt)
		match_in_genomes_list = match_in_genomes_list[1]
		}
	
	if (length(match_in_genomes_list) == 0)
		{
		return(NULL)
		}
	
	
	# Genome ID
	genomeID_in_xls = xls$GenBank.ID[match_in_genomes_list]
	cat("\nGenome ID found in directories: ", genomeID_in_xls)
	genome_dir_index = grep(pattern=genomeID_in_xls, x=genome_dirs, ignore.case=TRUE)
	genome_dir_index
	
	if (length(genome_dir_index) == 0)
		{
		return(NULL)
		}
	
	genome_dir = genome_dirs[genome_dir_index]
	if (returnwhat == "genome_dir")
		{
		return(genome_dir)
		}
	if (returnwhat == "spname")
		{
		return(spname)
		}
	if (returnwhat == "both")
		{
		extractwith='
		both = get_genome_name_from_protID(protID, genome_dirs, aln_annotations, xls, xlstaxa, returnwhat="both")
		spname = both$spname
		genome_dir = both$genome_dir
		'
		both = list()
		both$spname = spname
		both$genome_dir = genome_dir		
		return(both)
		}		
	return(stop("Shouldn't get here"))
	}
