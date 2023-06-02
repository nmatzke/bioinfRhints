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


get_spname_from_assembly_report <- function(assembly_report_table_fn)
	{
	# Also get the species name from the 
	# GCA_000005845.2_ASM584v2_assembly_report.txt
	assembly_report_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_assembly_report.txt")
	assembly_txt = readLines(assembly_report_table_fn)
	
	spname_line = assembly_txt[2]
	spname_line = gsub(pattern="# Organism name:", replacement="", x=spname_line)
	
	# Find & remove anything like (E. coli)
	tmpstrs = unlist(regmatches(spname_line, gregexpr("\\(.+?\\)", spname_line)))
	tmpstr = tmpstrs[length(tmpstrs)] # take the last bracketed text, if more than 1
	spname_line = gsub(pattern=tmpstr, replacement="", x=spname_line)
	spname = trim(gsub(pattern="\\(\\)", replacement="", x=spname_line))
	spname
	return(spname)
	}




# Handy function
extract_last_brackets <- function(list_of_strings, replace_spaces=TRUE)
	{
	species_names = rep("", length(list_of_strings))
	txt = paste0("\nextract_last_brackets() is processing ", length(list_of_strings), " strings. String #")
	cat(txt)
	for (i in 1:length(list_of_strings))
		{
		if (i == 1)
			{
			cat("\nlist_of_strings #", i, ",", sep="")
			} else {
			cat(i, ",", sep="")
			}
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
	} # END get_genome_name_from_protID <- function(protID, genome_dirs, aln_annotations, xls, xlstaxa, returnwhat="genome_dir")



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
wd = "/GitHub/bioinfRhints/flag/gene_order_example/" # example
#wd = "~/Downloads/Full_genomes/"	 # full dataset
setwd(wd)

# Get the list of genome directories that have been saved somewhere
genome_dirs = list.files(path="genomes", pattern=NULL, recursive=FALSE)
genome_dirs

gene_orders_all_fn = "gene_orders_all_v1.txt"

list_of_protIDs_in_genome = list()

unzipTF = FALSE
for (i in 1:length(genome_dirs))
	{
	if (i == 1)
		{
		cat("\ngenome_dirs #", i, ",", sep="")
		} else {
		cat(i, ",", sep="")
		}
	genome_dir = genome_dirs[i]
	# Access that directory, look for protein in list
	gene_order_table_zipfn = slashslash(paste0(wd, "/", "genomes/", genome_dir, "/", genome_dir, "_feature_table.txt.gz"))
	cmdstr = paste0("gunzip ", gene_order_table_zipfn)
	if (unzipTF == TRUE)
		{
		system(cmdstr)
		}
	gene_order_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_feature_table.txt")

	gene_order_df = read_gene_order_table(gene_order_table_fn)
	dim(gene_order_df)

	# Also get the species name from the 
	# GCA_000005845.2_ASM584v2_assembly_report.txt
	assembly_report_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_assembly_report.txt")
	spname = get_spname_from_assembly_report(assembly_report_table_fn)

	# Add some stuff to the front of the line
	id = genome_dir
	gene_order_df = cbind(id, gene_order_df)
	gene_order_df = cbind(spname, gene_order_df)
	
	#gene_orders_all_df = rbind(gene_orders_all_df, gene_order_df)

	if (i == 1)
		{
		write.table(gene_order_df, file=gene_orders_all_fn, sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE)
		} else {
		write.table(gene_order_df, file=gene_orders_all_fn, sep="\t", row.names=FALSE, append=TRUE, quote=FALSE, col.names=FALSE)
		}
	} # END for (i in 1:length(genome_dirs))

gene_orders_all_df = read.table(gene_orders_all_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
dim(gene_orders_all_df)

head(gene_orders_all_df)

gene_orders_all_df[3832:3836,]
dim(gene_orders_all_df)



#######################################################
# Gene order
#######################################################
# Name and load an MotA alignment file
alnfn = "/GitHub/bioinfRhints/flag/gene_order_example/motA_alignment_refined.fasta"
aln = read.fasta(alnfn)
aln_names = unname(sapply(X=aln, FUN=attr, "name"))

genome_names_not_found = NULL
protIDs_not_found = NULL
list_of_protIDs = aln_names 
gene_neighbors = NULL
i=1
for (i in 1:length(list_of_protIDs))
	{
	if (i == 1)
		{
		cat("\nlist_of_protIDs #", i, ",", sep="")
		} else {
		cat(i, ",", sep="")
		}
	protID = list_of_protIDs[i]

	# Access that directory, look for protein in list
	gene_num = grep(pattern=protID, x=gene_orders_all_df$product_accession, ignore.case=TRUE)
	if (length(gene_num) == 0)
		{
		protIDs_not_found = c(protIDs_not_found, protID)
		genome_names_not_found = c(genome_names_not_found, genome_dir)

		next()
		}

	spname = gene_orders_all_df$spname[gene_num]
	id = gene_orders_all_df$id[gene_num]
	# Assemble 2 lines: symbols and accessions
	# symbol line
	strand0 = nones_to_NA(gene_orders_all_df$strand[gene_num])
	sym0 = nones_to_NA(gene_orders_all_df$symbol[gene_num])
	acc0 = nones_to_NA(gene_orders_all_df$product_accession[gene_num])
	name0 = nones_to_NA(gene_orders_all_df$name[gene_num])
	
	strand1 = nones_to_NA(gene_orders_all_df$strand[gene_num+1])
	sym1 = nones_to_NA(gene_orders_all_df$symbol[gene_num+1])
	acc1 = nones_to_NA(gene_orders_all_df$product_accession[gene_num+1])
	name1 = nones_to_NA(gene_orders_all_df$name[gene_num+1])

	strand2 = nones_to_NA(gene_orders_all_df$strand[gene_num+2])
	sym2 = nones_to_NA(gene_orders_all_df$symbol[gene_num+2])
	acc2 = nones_to_NA(gene_orders_all_df$product_accession[gene_num+2])
	name2 = nones_to_NA(gene_orders_all_df$name[gene_num+2])

	strand3 = nones_to_NA(gene_orders_all_df$strand[gene_num+3])
	sym3 = nones_to_NA(gene_orders_all_df$symbol[gene_num+3])
	acc3 = nones_to_NA(gene_orders_all_df$product_accession[gene_num+3])
	name3 = nones_to_NA(gene_orders_all_df$name[gene_num+3])

	strandM1 = nones_to_NA(gene_orders_all_df$strand[gene_num-1])
	symM1 = nones_to_NA(gene_orders_all_df$symbol[gene_num-1])
	accM1 = nones_to_NA(gene_orders_all_df$product_accession[gene_num-1])
	nameM1 = nones_to_NA(gene_orders_all_df$name[gene_num-1])

	strandM2 = nones_to_NA(gene_orders_all_df$strand[gene_num-2])
	symM2 = nones_to_NA(gene_orders_all_df$symbol[gene_num-2])
	accM2 = nones_to_NA(gene_orders_all_df$product_accession[gene_num-2])
	nameM2 = nones_to_NA(gene_orders_all_df$name[gene_num-2])

	strandM3 = nones_to_NA(gene_orders_all_df$strand[gene_num-3])
	symM3 = nones_to_NA(gene_orders_all_df$symbol[gene_num-3])
	accM3 = nones_to_NA(gene_orders_all_df$product_accession[gene_num-3])
	nameM3 = nones_to_NA(gene_orders_all_df$name[gene_num-3])
	
	
	if (strand0 == "+")
		{
		tmprow = c(i, protID, spname, symM3, symM2, symM1, sym0, sym1, sym2, sym3, strandM3, strandM2, strandM1, strand0, strand1, strand2, strand3, accM3, accM2, accM1, acc0, acc1, acc2, acc3, nameM3, nameM2, nameM1, name0, name1, name2, name3, id)
		} else if (strand0 == "-") {
		tmprow = c(i, protID, spname, sym3, sym2, sym1, sym0, symM1, symM2, symM3, strand3, strand2, strand1, strand0, strandM1, strandM2, strandM3, acc3, acc2, acc1, acc0, accM1, accM2, accM3, name3, name2, name1, name0, nameM1, nameM2, nameM3, id)		
		}
	
	gene_neighbors = rbind(gene_neighbors, tmprow)
	}


gene_neighbors_df = as.data.frame(gene_neighbors, stringsAsFactors=FALSE)
names(gene_neighbors_df) = c("i", "protID", "spname", "symM3", "symM2", "symM1", "sym0", "sym1", "sym2", "sym3", "strandM3", "strandM2", "strandM1", "strand0", "strand1", "strand2", "strand3", "accM3", "accM2", "accM1", "acc0", "acc1", "acc2", "acc3", "nameM3", "nameM2", "nameM1", "name0", "name1", "name2", "name3", "id")
row.names(gene_neighbors_df) = NULL

gene_neighbors_df

outfn = "gene_neighbors_v1.txt"
write.table(gene_neighbors_df, file=outfn, sep="\t")


IDs_not_found = cbind(protIDs_not_found, genome_names_not_found)
write.table(IDs_not_found, file="IDs_not_found.txt", sep="\t")


gene_neighbors_df_orig = gene_neighbors_df


#######################################################
# Rename ExbDs etc.
#######################################################

# MotB homolog names - print to screen
sort(table(gene_neighbors_df$name1))

# Translate names to symbols (Excel file)
translate_motBs_fn = "~/Downloads/Full_genomes/MotB_translation_v1.xlsx"

# Works best if you convert the longest strings first
translate_motBs = read.xls(translate_motBs_fn, header=TRUE)
translate_motBs$orig = trim(translate_motBs$orig)
translate_motBs$new = trim(translate_motBs$new)
biggest_first = rev(order(str_length(translate_motBs$orig)))
translate_motBs = translate_motBs[biggest_first,]
head(translate_motBs)


TF = translate_motBs$orig == "flagellar motor protein MotB"
translate_motBs[TF,]

i=70
for (i in 1:nrow(translate_motBs))
	{
	# For the "blank" original
	if (translate_motBs$orig[i] == "")
		{
		gene_neighbors_df$name1[gene_neighbors_df$name1 == translate_motBs$orig[i]] = translate_motBs$new[i]
		gene_neighbors_df$name2[gene_neighbors_df$name2 == translate_motBs$orig[i]] = translate_motBs$new[i]
		gene_neighbors_df$name3[gene_neighbors_df$name3 == translate_motBs$orig[i]] = translate_motBs$new[i]
		gene_neighbors_df$name0[gene_neighbors_df$name0 == translate_motBs$orig[i]] = translate_motBs$new[i]
		gene_neighbors_df$nameM1[gene_neighbors_df$nameM1 == translate_motBs$orig[i]] = translate_motBs$new[i]
		gene_neighbors_df$nameM2[gene_neighbors_df$nameM2 == translate_motBs$orig[i]] = translate_motBs$new[i]
		gene_neighbors_df$nameM3[gene_neighbors_df$nameM3 == translate_motBs$orig[i]] = translate_motBs$new[i]
		next()
		}

	gene_neighbors_df$name1 = gsub(pattern=translate_motBs$orig[i], replacement=translate_motBs$new[i], x=gene_neighbors_df$name1)
	gene_neighbors_df$name2 = gsub(pattern=translate_motBs$orig[i], replacement=translate_motBs$new[i], x=gene_neighbors_df$name2)
	gene_neighbors_df$name3 = gsub(pattern=translate_motBs$orig[i], replacement=translate_motBs$new[i], x=gene_neighbors_df$name3)
	gene_neighbors_df$name0 = gsub(pattern=translate_motBs$orig[i], replacement=translate_motBs$new[i], x=gene_neighbors_df$name0)
	gene_neighbors_df$nameM1 = gsub(pattern=translate_motBs$orig[i], replacement=translate_motBs$new[i], x=gene_neighbors_df$nameM1)
	gene_neighbors_df$nameM2 = gsub(pattern=translate_motBs$orig[i], replacement=translate_motBs$new[i], x=gene_neighbors_df$nameM2)
	gene_neighbors_df$nameM3 = gsub(pattern=translate_motBs$orig[i], replacement=translate_motBs$new[i], x=gene_neighbors_df$nameM3)
	
	}

# Insert the symbols in the "sym" columns (when blank)
gene_neighbors_df$sym1[gene_neighbors_df$sym1[i] == ""] = gene_neighbors_df$name1[gene_neighbors_df$sym1[i] == ""]
gene_neighbors_df$sym2[gene_neighbors_df$sym1[i] == ""] = gene_neighbors_df$name2[gene_neighbors_df$sym1[i] == ""]
gene_neighbors_df$sym3[gene_neighbors_df$sym1[i] == ""] = gene_neighbors_df$name3[gene_neighbors_df$sym1[i] == ""]
gene_neighbors_df$sym0[gene_neighbors_df$sym1[i] == ""] = gene_neighbors_df$name0[gene_neighbors_df$sym1[i] == ""]
gene_neighbors_df$symM1[gene_neighbors_df$sym1[i] == ""] = gene_neighbors_df$nameM1[gene_neighbors_df$sym1[i] == ""]
gene_neighbors_df$symM2[gene_neighbors_df$sym1[i] == ""] = gene_neighbors_df$nameM2[gene_neighbors_df$sym1[i] == ""]
gene_neighbors_df$symM3[gene_neighbors_df$sym1[i] == ""] = gene_neighbors_df$nameM3[gene_neighbors_df$sym1[i] == ""]


outfn = "gene_neighbors_v1_translated.txt"

write.table(gene_neighbors_df, file=outfn, sep="\t")







