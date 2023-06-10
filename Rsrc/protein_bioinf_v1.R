
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

lastword <- function(string, split="/")
	{
	words = strsplit(gdata::trim(string), split=split)[[1]]
	return(gdata::trim(words[length(words)]))
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
get_genome_name_from_protID2 <- function(protID, genome_dirs, aln_annotations, xls, xlstaxa, returnwhat="genome_dir", genomes_path="~/Downloads/Full_genomes/", max.distance=0.2)
	{
	# Remove ".1" from protID
	words = strsplit(protID, split="\\.")[[1]]
	protID = words[1]
	
	TF = grepl(pattern=protID, x=aln_annotations)
	match_index = (1:length(aln_annotations))[TF]
	match_index

	seqname = aln_annotations[match_index]
	print(seqname)
	
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
	
	# There is more than one match 
	# (e.g. "Vibrio cholerae" matches:
	# Vibrio cholerae RFB16 chromosome 1
	# Vibrio cholerae RFB16 chromosome 2
	# Vibrio cholerae RFB16 plasmid
	# ...then search through the protein table for each one
	if (length(match_in_genomes_list) > 1)
		{
		genomeIDs_in_xls = xls$GenBank.ID[match_in_genomes_list]
		cat("\nGenome IDs found in directories: '", paste0(genomeIDs_in_xls, sep="', '"), "\n", sep="")
		
		genome_dir_indices = rep(0, times=length(match_in_genomes_list))
		for (i in 1:length(genome_dir_indices))
			{
			genome_dir_indices[i] = grep(pattern=genomeIDs_in_xls[i], x=genome_dirs, ignore.case=TRUE)
			
			genome_dir = genome_dirs[genome_dir_indices[i]]
			gene_order_table_fn = slashslash(paste0(genomes_path, "/genomes/", genome_dir, "/", genome_dir, "_feature_table.txt"))
			
			# Unzip if needed
			if (file.exists(gene_order_table_fn) == FALSE)
				{
				gene_order_table_zipfn = slashslash(paste0(genomes_path, "genomes/", genome_dir, "/", genome_dir, "_feature_table.txt.gz"))
				
				if (file.exists(gene_order_table_fn) == TRUE)
					{
					cmdstr = paste0("gunzip ", gene_order_table_zipfn)
					system(cmdstr)
					} else {
					txt = paste0("WARNING ERROR in get_genome_name_from_protID2: for protID='", protID, "', seqname='", seqname, "'. Neither of these files was found:\n'", gene_order_table_fn, "'\n'", gene_order_table_zipfn, "'\n Returning NULL.\n")
					cat("\n")
					cat(txt)
					warning(txt)
					return(NULL)
					}
				} # END if (file.exists(gene_order_table_fn) == FALSE)
			
			# Find the MotA protein
			gene_order_df = read_gene_order_table(gene_order_table_fn)
			gene_num = grep(pattern=protID, x=gene_order_df$product_accession2, ignore.case=TRUE)
			
			if (length(gene_num) == 1)
				{
				match_in_genomes_list = match_in_genomes_list[i]
				break()
				}
			} # END for (i in 1:length(genome_dir_indices))
		} # END if (length(match_in_genomes_list) > 1)

	# Multiple genomes searched, but none found
	if (length(match_in_genomes_list) > 1)
		{
		genomeIDs_in_xls = xls$GenBank.ID[match_in_genomes_list]
		missing_genomes_txt = paste0(c(xlstaxa[match_in_genomes_list]), sep="', '")
		txt= paste0("WARNING ERROR in get_genome_name_from_protID2: for protID='", protID, "', seqname='", seqname, "', no matches to the protID were found in any of the genomes matching spname='", spname, "'. Genomes searched:\n", missing_genomes_txt, "\nReturning NA.\n")
		cat("\n")
		cat(txt)
		warning(txt)
		return(NA)
		}
	
	# None found
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
	cat("\nGenome ID found in directories: ", genomeID_in_xls, sep="")
	genome_dir_index = grep(pattern=genomeID_in_xls, x=genome_dirs, ignore.case=TRUE)
	genome_dir_index
	
	if (length(genome_dir_index) == 0)
		{
		warning("WARNING ERROR: no genome_dir_index found. Returning NULL.")
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
	} # END get_genome_name_from_protID2 <- function(protID, genome_dirs, aln_annotations, xls, xlstaxa, returnwhat="genome_dir", max.distance=0.2)





# Get the gene order table 
read_gene_order_table_ATTEMPT <- function(gene_order_table_fn)
	{
	# Remove "#" from first line
	tmplines = readLines(gene_order_table_fn)
	tmplines[1] = gsub(pattern="# ", replacement="", x=tmplines[1])
	writeLines(tmplines, con=gene_order_table_fn)

	# NOTE: quote="\"" is necessary to avoid "EOF within quoted string"
	#       ...which causes issues with e.g. "2',3'-cyclic phosphodiesterase"
	gene_order_tmp = read.table(gene_order_table_fn, header=TRUE, comment.char="%", quote="", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
	
	# ERROR CHECK
	if (length(tmplines)-1 != nrow(gene_order_tmp))
		{
		txt = paste0("STOP ERROR in read_gene_order_table(", gene_order_table_fn, "): the file has ", length(tmplines), " lines, but read.table() gave only ", nrow(gene_order_tmp), " lines.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}

	# Take the genes, add protein name
	# Lines can include
# 	gene
# 	CDS with protein
# 	
# 	gene
# 	tRNA, rRNA
# 	
#   $class: protein_coding" "with_protein"   "tRNA"           ""               "pseudogene"   "rRNA" 
	
	# Find the pseudogenes, convert "gene" to "pseudogene" in the $feature category:
	all_features_TF = gene_order_tmp$feature == "gene"
	gene_order_df = gene_order_tmp[all_features_TF,]
	gene_order_df$feature[gene_order_df$class == "pseudogene"] = "pseudogene"
	
	# Relabel the pseudogenes in the raw table
	TF = gene_order_tmp$class == "pseudogene"
	gene_order_tmp$feature[TF] = "pseudogene"

#	TF = gene_order_tmp$class == "other"
#	gene_order_tmp$feature[TF] = "other"

	
	# Find the genes:
	nums = 1:nrow(gene_order_tmp)
	TF = gene_order_tmp$feature == "gene"
	gene_nums = nums[TF]
	
	# Find the products:
	product_nums = gene_nums + 1
	gene_order_df1 = gene_order_tmp[gene_nums,]
	gene_order_df2 = gene_order_tmp[product_nums,]
	
	# Error check
	if (nrow(gene_order_df1) != nrow(gene_order_df2))
		{
		txt = paste0("STOP ERROR in read_gene_order_table(", gene_order_table_fn, "): nrow(gene_order_df1) != nrow(gene_order_df2)")
		cat("\n")
		cat(txt)
		cat("\n")
		stop(txt)
		}
	
	# columns from protein
	cols_from_protein = c("start", "end", "product_accession", "non.redundant_refseq", "related_accession", "name", "symbol", "GeneID", "locus_tag", "feature_interval_length", "product_length", "attributes")
	tmpdf = gene_order_df2[,cols_from_protein]
	newcol_names = paste0(cols_from_protein, "2")
	names(tmpdf) = newcol_names
	gene_order_df12 = cbind(gene_order_df1, tmpdf)
	
	# Error check
	if (all(gene_order_df12$locus_tag == gene_order_df12$locus_tag2, na.rm=TRUE) == FALSE)
		{
		txt = paste0("STOP ERROR in read_gene_order_table(", gene_order_table_fn, "): not all gene_order_df12$locus_tag == gene_order_df12$locus_tag2)")
		cat("\n")
		cat(txt)
		cat("\n")
		stop(txt)
		}
	
	# Merge back into the main dataframe
	tmpmat2 = matrix(data="", nrow=nrow(gene_order_df), ncol=length(cols_from_protein))
	tmpdf2 = as.data.frame(tmpmat2, stringsAsFactors=FALSE)
	names(tmpdf2) = newcol_names
	gene_order_df = cbind(gene_order_df, tmpdf2)
	
	# Unique matches (locus_tag not always filled out)
	matchvals1 = c(apply(X=cbind(gene_order_df$symbol, gene_order_df$GeneID, gene_order_df$locus_tag), MARGIN=1, FUN=paste0, sep="_"))
	matchvals2 = c(apply(X=cbind(gene_order_df12$symbol, gene_order_df12$GeneID, gene_order_df12$locus_tag), MARGIN=1, FUN=paste0, sep="_"))
	
	matches = match(x=matchvals1, table=matchvals2)
	matches = matches[is.na(matches) == FALSE]
	
	gene_order_df[matches,] = gene_order_df12
	
	# column name to force to character
	colnames_for_char = c("feature","class","assembly","assembly_unit","seq_type","chromosome","genomic_accession","strand","product_accession","non.redundant_refseq","related_accession","name","symbol","GeneID","locus_tag","attributes","product_accession2","non.redundant_refseq2","related_accession2","name2","symbol2","GeneID2","locus_tag2","attributes2")
	
	cls.df(gene_order_df)
	for (i in 1:length(colnames_for_char))
		{
		gene_order_df[,colnames_for_char[i]] = as.character(gene_order_df[,colnames_for_char[i]])
		}
	gene_order_df$start = as.integer(gene_order_df$start)
	gene_order_df$end = as.integer(gene_order_df$end)
	gene_order_df$start2 = as.integer(gene_order_df$start2)
	gene_order_df$end2 = as.integer(gene_order_df$end2)
	gene_order_df$feature_interval_length = as.integer(gene_order_df$feature_interval_length)
	gene_order_df$product_length = as.integer(gene_order_df$product_length)
	gene_order_df$feature_interval_length2 = as.integer(gene_order_df$feature_interval_length2)
	gene_order_df$product_length2 = as.integer(gene_order_df$product_length2)
	cls.df(gene_order_df)
	
	# This seems to have everything...
	return(gene_order_df)
	}



# Read a gene order table (aka a protein feature table)
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

	# This seems to have everything...
	gene_order_df = gene_order_tmp
	head(gene_order_df)
	return(gene_order_df)
	}


# Get the gene order table 
read_gene_order_table_HALF <- function(gene_order_table_fn)
	{
	# Remove "#" from first line
	tmplines = readLines(gene_order_table_fn)
	tmplines[1] = gsub(pattern="# ", replacement="", x=tmplines[1])
	writeLines(tmplines, con=gene_order_table_fn)

	# NOTE: quote="\"" is necessary to avoid "EOF within quoted string"
	#       ...which causes issues with e.g. "2',3'-cyclic phosphodiesterase"
	gene_order_tmp = read.table(gene_order_table_fn, header=TRUE, comment.char="%", quote="", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
	
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



	# This seems to have everything...
	gene_order_df = gene_order_df2
	head(gene_order_df)
	return(gene_order_df)
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










# prokka outputs have just e.g.:
# 
junk='
locus_tag	ftype	length_bp	gene	EC_number	COG	product
DGBBEKCF_00001	CDS	2412	mshA_1	2.4.1.250		D-inositol-3-phosphate glycosyltransferase
DGBBEKCF_00002	CDS	762				hypothetical protein
'

# For merging into a huge protein features table,
# we want to have a _feature_table.txt file with:
junk='
# feature	class	assembly	assembly_unit	seq_type	chromosome	genomic_accession	start	end	strand	product_accession	non-redundant_refseq	related_accession	name	symbol	GeneID	locus_tag	feature_interval_length	product_length	attributes
gene	protein_coding	GCA_900093645.1	Primary Assembly	unplaced scaffold		FLYF01000089.1	3	122	+							AB751O23_DK_00010	120		partial
CDS	with_protein	GCA_900093645.1	Primary Assembly	unplaced scaffold		FLYF01000089.1	3	122	+	SCA59112.1			hypothetical protein			AB751O23_DK_00010	120	39	partial
'

convert_prokka_tsv_to_prot_feature_table <- function(tsv_fn, write_to_file=TRUE)
	{
	junk='
	tsv_fn = "~/Downloads/Full_genomes/genomes/PRJEB38681_Verrucomicrobium_sp_CAISZB01/PRJEB38681_Verrucomicrobium_sp_CAISZB01.tsv"
	
	prot_feature_table_df = convert_prokka_tsv_to_prot_feature_table(tsv_fn)
	
	prot_feature_table_fn = gsub(pattern=".tsv", replacement="_feature_table.txt", x=tsv_fn)
	write.table(prot_feature_table_df, file=prot_feature_table_fn, sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE)

	tsv_fn = "~/Downloads/Full_genomes/genomes/PRJEB3868_Verrucomicrobium_sp_CAIZXV01/PRJEB3868_Verrucomicrobium_sp_CAIZXV01.tsv"
	
	prot_feature_table_df = convert_prokka_tsv_to_prot_feature_table(tsv_fn)


	'
	
	headers = c("feature", "class", "assembly", "assembly_unit", "seq_type", "chromosome", "genomic_accession", "start", "end", "strand", "product_accession", "non-redundant_refseq", "related_accession", "name", "symbol", "GeneID", "locus_tag", "feature_interval_length", "product_length", "attributes")
	
	tsvdf = read.table(tsv_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
	unique(tsvdf$ftype)
	
	# Fill in the fields
	feature = tsvdf$ftype
	class = rep("with_protein", times=nrow(tsvdf))
	class[feature != "CDS"] = ""
	
	filename_only_tsv_fn = lastword(tsv_fn)
	assembly_name = gsub(pattern=".tsv", replacement="", x=filename_only_tsv_fn)
	assembly = rep(assembly_name, times=nrow(tsvdf))
	
	assembly_unit = rep("prokka assembly", times=nrow(tsvdf))
	seq_type = rep("prokka unplaced scaffold", times=nrow(tsvdf))
	chromosome = rep("", times=nrow(tsvdf))
	genomic_accession = tsvdf$locus_tag
	
	# Make fake start and end, based on order in file
	nums = 1:nrow(tsvdf)
	startvals = rep(0, nrow(tsvdf))
	endvals = rep(0, nrow(tsvdf))
	
	for (i in 1:nrow(tsvdf))
		{
		if (i == 1)
			{
			startvals[i] = 1
			endvals[i] = startvals[i] + tsvdf$length_bp[i] - 1
			} else {
			startvals[i] = endvals[i-1] + 1
			endvals[i] = startvals[i] + tsvdf$length_bp[i] - 1			
			}
		}
	
	start = startvals
	end = endvals
	strand = rep("+", times=nrow(tsvdf))
	product_accession = tsvdf$locus_tag
	non.redundant_refseq = rep("", times=nrow(tsvdf))
	related_accession = rep("", times=nrow(tsvdf))
	name = tsvdf$product
	symbol = tsvdf$gene
	GeneID = rep("", times=nrow(tsvdf))
	locus_tag = tsvdf$locus_tag
	feature_interval_length = tsvdf$length_bp
	product_length = tsvdf$length_bp / 3
	
	# enzyme pathway
	EC_numbers = paste0("EC_number:", tsvdf$EC_number, sep="")
	EC_numbers[EC_numbers == "EC_number:"] = ""
	EC_numbers
	
	# COGs etc.
	COGs = paste0("COG:", tsvdf$COG, sep="")
	COGs[COGs == "COG:"] = ""
	COGs
	
	# merge
	tmptxt = paste0(EC_numbers, " | ", COGs, sep="")
	tmptxt[tmptxt == " | "] = ""
	attributes = tmptxt
	
	feature_table_df = cbind(feature, class, assembly, assembly_unit, seq_type, chromosome, genomic_accession, start, end, strand, product_accession, non.redundant_refseq, related_accession, name, symbol, GeneID, locus_tag, feature_interval_length, product_length, attributes)
	prot_feature_table_df = as.data.frame(feature_table_df, stringsAsFactors=FALSE)
	names(prot_feature_table_df) = headers
	
	prot_feature_table_df$start = as.numeric(prot_feature_table_df$start)
	prot_feature_table_df$end = as.numeric(prot_feature_table_df$end)
	prot_feature_table_df$feature_interval_length = as.numeric(prot_feature_table_df$feature_interval_length)
	prot_feature_table_df$product_length = as.numeric(prot_feature_table_df$product_length)
	
	# Also, write to file
	if (write_to_file == TRUE)
		{
		prot_feature_table_fn = gsub(pattern=".tsv", replacement="_feature_table.txt", x=tsv_fn)
		write.table(prot_feature_table_df, file=prot_feature_table_fn, sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE)
		}

	return(prot_feature_table_df)
	}
