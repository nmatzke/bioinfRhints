
#######################################################
# Source code for bioinformatics tasks
#######################################################
fullnames_from_readFasta <- function(aln)
	{
	fullnames = sapply(X=aln, FUN=attr, which="Annot")
	fullnames = gsub(pattern=">", replacement="", x=fullnames)
	return(fullnames)
	}


# Assumes: tree tip labels start with a GID
# "fullnames" can be matched by grepl
# "fullnames" has no huge problems (e.g. "(" characters)
relabel_tree_tips_wGID_first <- function(tr, fullnames, split="\\|")
	{
	trgids = rep("", length(tr$tip.label))
	nums = 1:length(fullnames)
	for (i in 1:length(tr$tip.label))
		{
		words = strsplit(tr$tip.label[i], split=split)[[1]]
		trgids[i] = gdata::trim(stringr::str_squish(words[1]))
		TF = grepl(pattern=trgids[i], x=fullnames)
		if (sum(TF) == 1)
			{
			num_in_fullnames = nums[TF]
			tr$tip.label[i] = fullnames[num_in_fullnames]
			} else if (sum(TF) == 0) {
			pass=1
			} else if (sum(TF) > 1) {
			txt = paste0("\nSTOP ERROR in relabel_tree_tips_wGID_first(): gid '", trgids[i], "' had ", sum(TF), " matches in fullnames. It should only have 1 match. Printing the fullnames matches...\n")
			cat(txt)
			print(fullnames[TF])
			stop(txt)
			} # END if (sum(TF) == 1)
		} 
	return(tr)
	}




# Assumes: tree tip labels start with a GID
# "fullnames" can be matched by grepl
# "fullnames" has no huge problems (e.g. "(" characters)
relabel_fastaAln_wGID_first <- function(aln_to_rename, fullnames, split="\\|")
	{
	alnFullnames = fullnames_from_readFasta(aln_to_rename)
	nums = 1:length(fullnames)
	aln2 = list()
	j = 0
	for (i in 1:length(aln_to_rename))
		{
		words = strsplit(alnFullnames[i], split=split)[[1]]
		tmp_gid = gdata::trim(stringr::str_squish(words[1]))
		TF = grepl(pattern=tmp_gid, x=fullnames)
		if (sum(TF) == 1)
			{
			num_in_fullnames = nums[TF]
			aln2[[(j=j+1)]] = aln_to_rename[[i]]
			attr(aln2[[j]], which="Annot") = fullnames[num_in_fullnames]
			} else if (sum(TF) == 0) {
			pass=1
			} else if (sum(TF) > 1) {
			txt = paste0("\nSTOP ERROR in relabel_fastaAln_wGID_first(): gid '", tmp_gid, "' had ", sum(TF), " matches in fullnames. It should only have 1 match. Printing the fullnames matches...\n")
			cat(txt)
			print(fullnames[TF])
			stop(txt)
			} # END if (sum(TF) == 1)
		} 
	return(aln2)
	}





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
library(openxlsx)				# for read.xlsx (replaces gdata::read.xls and xlsx::read.xlsx)
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
xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet=1)

# Get the genome taxon
xlstaxa = rep("", times=nrow(xls))
for (i in 1:nrow(xls))
	{
	tmp = paste(xls$Genus[i], xls$Species[i], xls$Strain[i], sep=" ", collapse="")
	xlstaxa[i] = stringr::str_squish(tmp)
	}

protID = "AAC74960.1"
genome_dir = get_genome_name_from_protID(protID, genome_dirs, aln_annotations, xls, xlstaxa)

' # END JUNK

# find a protein in a genome
# protID -> fasta label in aln_annotations
# -> extract the strain name in [brackets] (at least anything before "=")
# -> search the strain name in xlstaxa (which comes from xls by merging Genus, species, and strain)
# -> from that hit, get the genomeID
# -> return the genome_dir and/or the spname
#
get_genome_name_from_protID2 <- function(protID, genome_dirs, aln_annotations, xls, xlstaxa, returnwhat="genome_dir", max.distance=0.13)
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
	
	# agrep = approximate grep
	match_in_genomes_list = agrep(pattern=spname, x=xlstaxa, max.distance=max.distance, ignore.case=TRUE)
	
	if (length(match_in_genomes_list) > 1)
		{
		tmpdists = rep(0, times=length(match_in_genomes_list))
		
		for (i in 1:length(match_in_genomes_list))
			{
			tmpdists[i] = adist(x=spname, y=xlstaxa[match_in_genomes_list[i]])
			}
		
		TF = order(tmpdists) == 1
		closest_match = xlstaxa[TF]
		
# 		txt = paste0("Warning in get_genome_name_from_protID(protID='", protID, "'...): >1 genome matched spname='", spname, "'. Printing them below, taking the first one.")
# 		cat("\n\n")
# 		cat(txt)
# 		cat("\nxlstaxa[match_in_genomes_list]:\n")
# 		print(xlstaxa[match_in_genomes_list])
# 		cat("\n")
# 		warning(txt)
		match_in_genomes_list = closest_match[1]
		}
	
	if (length(match_in_genomes_list) == 0)
		{
		txt = paste0("WARNING from get_genome_name_from_protID(): no approximate grep (agrep) hit found for protID='", protID, "' at max.distance=", max.distance, ". Returning NA.")
		cat("\n")
		cat(txt)
		cat("\n")
		warning(txt)
		return(NA)
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
