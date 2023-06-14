
#######################################################
# Source code for bioinformatics tasks
#######################################################

# Run "trim" on each column of a data.frame or similar
trim_each_col_in_df <- function(xlsx)
	{
	for (i in 1:ncol(xlsx))
		{
		xlsx[,i] = stringr::str_trim(xlsx[,i])
		}
	return(xlsx)
	}

# In a subsection of a feature_table, copy the info from a nonblank cell 
# to a blank cell
copy_over_blankcell <- function(tmpdf2)
	{
	tmpnums = 1:nrow(tmpdf2)
	
	# Names
	unique_names = unique(tmpdf2$name[tmpnums])
	if (sum(isblank_TF(unique_names)) > 0)
		{
		unique_names = unique_names[isblank_TF(unique_names) == FALSE]
		}
	if (length(unique_names) == 0)
		{
		tmpdf2$name[tmpnums] = ""
		}
	if (length(unique_names) == 1)
		{
		tmpdf2$name[tmpnums] = unique_names
		}
	if (length(unique_names) > 1)
		{
		tmpdf2$name[tmpnums] = paste0(unique_names, collapse="_")
		}


	# Product accession
	unique_product_accessions = unique(tmpdf2$product_accession[tmpnums])
	if (sum(isblank_TF(unique_product_accessions)) > 0)
		{
		unique_product_accessions = unique_product_accessions[isblank_TF(unique_product_accessions) == FALSE]
		}
	if (length(unique_product_accessions) == 0)
		{
		tmpdf2$product_accession[tmpnums] = ""
		}
	if (length(unique_product_accessions) == 1)
		{
		tmpdf2$product_accession[tmpnums] = unique_product_accessions
		}
	if (length(unique_product_accessions) > 1)
		{
		tmpdf2$product_accession[tmpnums] = paste0(unique_product_accessions, collapse="_")
		}

	
	# Symbol
	unique_symbols = unique(tmpdf2$symbol[tmpnums])
	if (sum(isblank_TF(unique_symbols)) > 0)
		{
		unique_symbols = unique_symbols[isblank_TF(unique_symbols) == FALSE]
		}
	if (length(unique_symbols) == 0)
		{
		tmpdf2$symbol[tmpnums] = ""
		}
	if (length(unique_symbols) == 1)
		{
		tmpdf2$symbol[tmpnums] = unique_symbols
		}
	if (length(unique_symbols) > 1)
		{
		tmpdf2$symbol[tmpnums] = paste0(unique_symbols, collapse="_")
		}

	# Strand
	unique_strands = unique(tmpdf2$strand[tmpnums])
	if (sum(isblank_TF(unique_strands)) > 0)
		{
		unique_strands = unique_strands[isblank_TF(unique_strands) == FALSE]
		}
	if (length(unique_strands) == 0)
		{
		tmpdf2$strand[tmpnums] = ""
		}
	if (length(unique_strands) == 1)
		{
		tmpdf2$strand[tmpnums] = unique_strands
		}
	if (length(unique_strands) > 1)
		{
		tmpdf2$strand[tmpnums] = paste0(unique_strands, collapse="_")
		}
	
	return(tmpdf2)
	} # END copy_over_blankcell


#######################################################
# Get an adjacent item from a full *_feature_table.txt file
# 
# Usually you go up (or down) two rows, but this can
# break down for e.g. pseudogenes
#######################################################

get_adjacent_row <- function(rownum, shift, prot_feature_tables_all_df)
	{
	countshift = shift
	nums = rownum:(rownum+3*shift)

	# symbol + locus tag should give something unique to count
	symbol_locus = paste0(prot_feature_tables_all_df$symbol[nums], "|", prot_feature_tables_all_df$locus_tag[nums])
	core_gene_symbol_locus = symbol_locus[1]
	unique_symbol_locus = unique(symbol_locus)
	keepTF = unique_symbol_locus != core_gene_symbol_locus
	unique_symbol_locus = unique_symbol_locus[keepTF]
	
	TF = symbol_locus %in% unique_symbol_locus
	tmpnums = (1:length(TF))[TF]
	index_to_rownum_shifts = tmpnums
	
	# Copy info over blanks ""
	if (shift < 1)
		{
		countshift = -1 * shift # for counting backwards if needed
		}
	tmpdf = prot_feature_tables_all_df[nums,]
	usl = unique_symbol_locus[countshift]
	TF = symbol_locus %in% usl
	usl_nums = (1:length(TF))[TF]
	tmpdf2 = tmpdf[usl_nums,]
	tmpdf2 = copy_over_blankcell(tmpdf2) # fix the issue of blanks
	
	# Extract the info
	strand = tmpdf2$strand[1]
	symbol = tmpdf2$symbol[1]
	product_accession = tmpdf2$product_accession[1]
	name = tmpdf2$name[1]
	info = c(strand, symbol, product_accession, name)
	names(info) = c("strand", "symbol", "product_accession", "name")
	
	return(info)
	}	# END get_adjacent_row






get_adjacent_genes <- function(list_of_protIDs, prot_feature_tables_all_df, genomes_to_spnames_df, printwarnings=FALSE)
	{
	genome_names_not_found = NULL
	protIDs_not_found = NULL
	gene_neighbors = NULL
	i=1
	txt = paste0("\nget_adjacent_genes() is extracting adjacent gene information for ", length(list_of_protIDs), " input protein IDs (GenBank IDs, gids).\ni=")
	cat(txt)
	for (i in 1:length(list_of_protIDs))
		{
		cat(i, ",", sep="")
		protID = list_of_protIDs[i]
	
		# Slow
		#gene_num = grep(pattern=protID, x=prot_feature_tables_all_df$product_accession, ignore.case=TRUE)
		TF = prot_feature_tables_all_df$product_accession == protID
		gene_num = (1:length(TF))[TF]
		if (length(gene_num) == 0)
			{
			protIDs_not_found = c(protIDs_not_found, protID)
			genome_names_not_found = c(genome_names_not_found, genome_dir)

			next()
			}
	
		# Get basic info
		assembly = prot_feature_tables_all_df$assembly[gene_num]
		
		if (printwarnings == FALSE)
			{
			group = suppressWarnings(get_group_from_genome_ID(assembly=assembly, genomes_to_spnames_df=genomes_to_spnames_df, printwarnings=FALSE))
			spname = suppressWarnings(get_spname_from_genome_ID(assembly=assembly, genomes_to_spnames_df=genomes_to_spnames_df, printwarnings=FALSE))	
			} else {
			group = get_group_from_genome_ID(assembly=assembly, genomes_to_spnames_df=genomes_to_spnames_df, printwarnings=TRUE)
			spname = get_spname_from_genome_ID(assembly=assembly, genomes_to_spnames_df=genomes_to_spnames_df, printwarnings=TRUE)	
			}

		# Look at near-neighbors
		prot_feature_tables_all_df[gene_num,]
		prot_feature_tables_all_df[gene_num+1,]
		prot_feature_tables_all_df[gene_num-1,]

		rownum = gene_num
		info1 = get_adjacent_row(rownum, shift=1, prot_feature_tables_all_df)
		info2 = get_adjacent_row(rownum, shift=2, prot_feature_tables_all_df)
		info3 = get_adjacent_row(rownum, shift=3, prot_feature_tables_all_df)
	
		infoM1 = get_adjacent_row(rownum, shift=-1, prot_feature_tables_all_df)
		infoM2 = get_adjacent_row(rownum, shift=-2, prot_feature_tables_all_df)
		infoM3 = get_adjacent_row(rownum, shift=-3, prot_feature_tables_all_df)

		# Assemble 2 lines: symbols and accessions
		# symbol line
		strand0 = nones_to_NA(prot_feature_tables_all_df$strand[gene_num])
		sym0 = nones_to_NA(prot_feature_tables_all_df$symbol[gene_num])
		acc0 = nones_to_NA(prot_feature_tables_all_df$product_accession[gene_num])
		name0 = nones_to_NA(prot_feature_tables_all_df$name[gene_num])
	
		strand1 = nones_to_NA(info1["strand"])
		sym1 = nones_to_NA(info1["symbol"])
		acc1 = nones_to_NA(info1["product_accession"])
		name1 = nones_to_NA(info1["name"])

		strand2 = nones_to_NA(info2["strand"])
		sym2 = nones_to_NA(info2["symbol"])
		acc2 = nones_to_NA(info2["product_accession"])
		name2 = nones_to_NA(info2["name"])

		strand3 = nones_to_NA(info3["strand"])
		sym3 = nones_to_NA(info3["symbol"])
		acc3 = nones_to_NA(info3["product_accession"])
		name3 = nones_to_NA(info3["name"])

		strandM1 = nones_to_NA(infoM1["strand"])
		symM1 = nones_to_NA(infoM1["symbol"])
		accM1 = nones_to_NA(infoM1["product_accession"])
		nameM1 = nones_to_NA(infoM1["name"])

		strandM2 = nones_to_NA(infoM2["strand"])
		symM2 = nones_to_NA(infoM2["symbol"])
		accM2 = nones_to_NA(infoM2["product_accession"])
		nameM2 = nones_to_NA(infoM2["name"])

		strandM3 = nones_to_NA(infoM3["strand"])
		symM3 = nones_to_NA(infoM3["symbol"])
		accM3 = nones_to_NA(infoM3["product_accession"])
		nameM3 = nones_to_NA(infoM3["name"])
	
		if (strand0 == "+")
			{
			tmprow = c(i, protID, spname, group, assembly, symM3, symM2, symM1, sym0, sym1, sym2, sym3, strandM3, strandM2, strandM1, strand0, strand1, strand2, strand3, accM3, accM2, accM1, acc0, acc1, acc2, acc3, nameM3, nameM2, nameM1, name0, name1, name2, name3, prot_feature_tables_all_fn)
			} else if (strand0 == "-") {
			tmprow = c(i, protID, spname, group, assembly, sym3, sym2, sym1, sym0, symM1, symM2, symM3, strand3, strand2, strand1, strand0, strandM1, strandM2, strandM3, acc3, acc2, acc1, acc0, accM1, accM2, accM3, name3, name2, name1, name0, nameM1, nameM2, nameM3, prot_feature_tables_all_fn)		
			}
	
		gene_neighbors = rbind(gene_neighbors, tmprow)
		}
	cat("...done\n")

	gene_neighbors_df = as.data.frame(gene_neighbors, stringsAsFactors=FALSE)
	names(gene_neighbors_df) = c("i", "protID", "spname", "group", "assembly", "symM3", "symM2", "symM1", "sym0", "sym1", "sym2", "sym3", "strandM3", "strandM2", "strandM1", "strand0", "strand1", "strand2", "strand3", "accM3", "accM2", "accM1", "acc0", "acc1", "acc2", "acc3", "nameM3", "nameM2", "nameM1", "name0", "name1", "name2", "name3", "prot_feature_tables_all_fn")
	row.names(gene_neighbors_df) = NULL
	return(gene_neighbors_df)
	}








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


# Create "spname", a short & program-universal species/strain name
# ESPECIALLY, REMOVE THESE CHARACTERS, WHICH MESS UP NEWICK FILES:
# ( ) ; , \ / [ ]
clean_gb_spname <- function(spname)
	{
	spname = gsub(pattern="serovar", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="unclassified", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="uncultured", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Candidatus", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="chromosome ", replacement="chrom", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="sp.", replacement="sp", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="WGS isolate:", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="plasmid ", replacement="pl_", x=spname, ignore.case=TRUE)

	spname = stringr::str_trim(stringr::str_squish(spname))
	sum(grepl("Wolfebacteria Wolfebacteria", x=spname))
	sum(grepl("Wolfebacteria_Wolfebacteria", x=spname))

	spname = gsub(pattern="Korarchaeota Candidatus", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Chlamydiales bacterium Chlamydiales bacterium", replacement="Chlamydiales bacterium", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Candidatus Prometheoarchaeum Candidatus Prometheoarchaeum", replacement="Prometheoarchaeum", x=spname, ignore.case=TRUE)

	spname = gsub(pattern="of Drosophila melanogaster isolate wMel", replacement="Dmelanogaster", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Beckwithbacteria Beckwithbacteria", replacement="Beckwithbacteria", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Campbellbacteria Campbellbacteria", replacement="Campbellbacteria", x=spname, ignore.case=TRUE)

	spname = gsub(pattern="Moranbacteria Moranbacteria", replacement="Moranbacteria", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Uhrbacteria Uhrbacteria", replacement="Uhrbacteria", x=spname, ignore.case=TRUE)

	spname = gsub(pattern="Wolfebacteria Wolfebacteria", replacement="Wolfebacteria", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Woesebacteria Woesebacteria", replacement="Woesebacteria", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Wolfebacteria_Wolfebacteria", replacement="Wolfebacteria", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="Woesebacteria_Woesebacteria", replacement="Woesebacteria", x=spname, ignore.case=TRUE)

	sum(grepl("Wolfebacteria Wolfebacteria", x=spname))
	sum(grepl("Wolfebacteria_Wolfebacteria", x=spname))

# 	spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)

	spname = gsub(pattern=";", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="\\(", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="\\)", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="\\[", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="\\]", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="=", replacement="EQ", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="\\\\", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="str. ", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="\\.", replacement="", x=spname, ignore.case=TRUE)
	spname = gsub(pattern="/", replacement="-", x=spname, ignore.case=TRUE)
	spname = gsub(pattern=":", replacement="", x=spname, ignore.case=TRUE)

	# Remove extra spaces
	spname = stringr::str_trim(stringr::str_squish(spname))

	# Add underscores
	spname = gsub(pattern=" ", replacement="_", x=spname, ignore.case=TRUE)
	#sum(grepl("Wolfebacteria_Wolfebacteria", x=spname))
	
	return(spname)
	}


#######################################################
# Parse the names for protein name info
#######################################################
extract_protein_name_info <- function(fullnames)
	{
	example_run='
	library(ape)
	library(phytools)
	library(seqinr)				# for read.fasta
	library(BioGeoBEARS)	# for cls.df
	library(gdata)				# for trim
	library(stringr) 			# for regmatches
	library(openxlsx)			# for openxlsx::read.xlsx

	sourceall("/GitHub/bioinfRhints/Rsrc/") # for 
	fullnames = ">ACF14374.1 MotA/TolQ/ExbB proton channel [Chloroherpeton thalassium ATCC 35110]"
	protein_name = extract_protein_name_info(fullnames)
	protein_name
	'
	
	protein_name = rep("", times=length(fullnames))
	species_names_wSpaces = extract_last_brackets(list_of_strings=fullnames, replace_spaces=FALSE)

	for (i in 1:length(fullnames))
		{
		orig_name = fullnames[i]

		# Remove ">", ], [
		subname = gsub(pattern=">", replacement="", x=orig_name)
		subname = gsub(pattern="\\[", replacement="", x=subname)
		subname = gsub(pattern="\\]", replacement="", x=subname)

		# Remove GenBank ID:
		gid = firstword(subname)
		subname = gsub(pattern=gid, replacement="", x=subname)
	
		# Remove species name:
		subname = gsub(pattern=species_names_wSpaces[i], replacement="", x=subname)

		# REMOVE SEMICOLONS, GD IT
		subname = gsub(pattern="\\;", replacement="", x=subname, ignore.case=TRUE)
		# REMOVE EQUALS, GD IT
		subname = gsub(pattern="\\=", replacement="EQ", x=subname, ignore.case=TRUE)

	
		# fix spaces
		subname = gsub(pattern="  ", replacement=" ", x=subname)
		subname = gsub(pattern="  ", replacement=" ", x=subname)
		subname = gsub(pattern="  ", replacement=" ", x=subname)
		subname = gsub(pattern="  ", replacement=" ", x=subname)
		subname = gdata::trim(subname)
		
		protein_name[i] = subname
		}
	return(protein_name)
	}


classify_MotAfam_labels <- function(list_of_strings)
	{
	junk='
	
	list_of_strings = protein_name
	short_protname = classify_MotAfam_labels(list_of_strings=protein_name)
	short_protname
	'
	
	
	short_protname = rep("", length(list_of_strings))
	for (i in 1:length(list_of_strings))
		{
		tmpstr = list_of_strings[i]

		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)

		
		tmpstr = gsub(pattern="probable biopolymer transport protein", replacement="biopoly_transp", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Biopolymer transport proteins", replacement="biopoly_transp", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbB protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="related to flagellar apparatus MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)
		
		# Add semicolon on the end
		
		# Many are just MotA/TolQ/ExbB
		if (grepl(pattern="MotA/TolQ/ExbB", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "AQB"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="gliding", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "gliding"
			short_protname[i] = tmpstr
			next()
			}
		
		# Remove obvious words - flagellum
		tmpstr = gsub(pattern="transport membrane proton channel, TolQ-related protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport membrane proton channel, TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)


		tmpstr = gsub(pattern="flagellar motor rotation protein ", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar motor protein ", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MAG: flagellar motor protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="chemotaxis protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="chemotaxis", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="proton channel family", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar motor component", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MotA protein", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="motility protein A", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar basal body stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="endoflagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="flagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="MAG\\: flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="MotA\\; MotA component of the H\\+\\-coupled stator flagellum complex", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="Mot family proton \\(H\\+\\) or sodium \\(Na\\+\\) transporter MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="sodium channel stator-force generator subunit of flagellar rotation", replacement="sodium_MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="signal recognition particle-docking protein FtsY", replacement="SRP_FtsY", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="sodium-driven polar flag MotA", replacement="sodiumPolarMotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="motA\\, \\(MotA\\)", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="MAG: flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="pomA protein", replacement="PomA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="protein of unassigned function", replacement="unk", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="multi-sensor hybrid histidine kinase", replacement="HistKin", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MotA\\; MotA component of the H\\+\\-coupled stator flagellum complex", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Mot family proton \\(H+\\) or sodium \\(Na+\\) transporter MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="membrane spanning protein in TonB-ExbB-ExbD complex", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		#tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)
		#tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="flagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="putative", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="motor", replacement="", x=tmpstr, ignore.case=TRUE)



		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)
		
		# Any blank here is just a motor protein
		if ((tmpstr == "") || (tmpstr == "protein"))
			{
			tmpstr = "motor"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="tonB system transport protein ExbB/TolQ", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "TolQ_ExbB_wTonB"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="tonB system transport protein ExbB/TolQ", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "TolQ_ExbB_wTonB"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="ExbB/TolQ", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "ExbB_TolQ"
			short_protname[i] = tmpstr
			next()
			}
		tmpstr = gsub(pattern="system transport protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="accessory protein", replacement="", x=tmpstr, ignore.case=TRUE)


		if (tolower(tmpstr) == "probable biopolymer transport protein")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}

		if (tolower(tmpstr) == "biopolymer transporter")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}
		if (tolower(tmpstr) == "biopolymer transport proteins")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}
		if (tolower(tmpstr) == "biopolymer transport protein")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}


		tmpstr = gsub(pattern="hypothetical protein ", replacement="HYP_", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="outer membrane transport energization protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Biopolymer", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="domain-containing protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="colicin uptake protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbB-related protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ferric siderophore transport system,", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable tolQ-type", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable tolQ protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Cell division and transport-associated protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolQ protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Tol-Pal system protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="(ExbB-like)", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Ton complex subunit ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQ-like protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQ transporter", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Tol-Pal system subunit TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="\\(ExbB\\)", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbB protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tonB\\-system energizer ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TonB exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)


		tmpstr = gsub(pattern="Cell division and transportassociated TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tonBsystem energizer ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolPal system TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQtype", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolPal system subunit TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQlike protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQrelated protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="sodium channel statorforce generator subunit of flagellar rotation", replacement="sodium MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="related to TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="related to \\(TolQ\\)", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolerance to group A colicins, singlestranded DNA filamentous phage, required for OM integrity", replacement="TolA_OM", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="\\(ExbBlike\\)", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBrelated protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBlike", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBrelated protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBlike protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Cell division and transportassociated TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tonBsystem energizer ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolPal system TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolpal systemassociated acylCoA thioesterase", replacement="tolpal thioest", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="sodiumdriven polar flag MotA", replacement="sodium MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MotA MotA component of the H+coupled stator flagellum complex", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical proteintransmembrane prediction", replacement="HYP", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="membrane spanning protein in TonBExbBExbD complex", replacement="ExbB", x=tmpstr, ignore.case=TRUE)




		tmpstr = gsub(pattern="MAG: TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical membrane protein", replacement="HYP", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical protein-transmembrane prediction", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern=", tolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolerance to group A colicins, single-stranded DNA filamentous phage, required for OM integrity", replacement="colicin-tol", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tol-pal system-associated acyl-CoA thioesterase ", replacement="TolPal_thioester", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="two-component sensor histidine kinase YdfI", replacement="2compSens_hisKin_YdfI", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="protein of unassigned function", replacement="unk", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="conserved hypothetical protein", replacement="HYP", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="motA\\, \\(MotA\\)", replacement="MotA", x=tmpstr, ignore.case=TRUE)


		tmpstr = gsub(pattern="flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MAG: MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar proton channel", replacement="flag", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MAG\\: MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)

		tmpstr = gsub(pattern="tol\\-pal system\\-associated acyl\\-CoA thioesterase", replacement="tolPal_thioest", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="probable ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="TonB ExbB2", replacement="ExbB2", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TonB exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TonB ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="LafT protein", replacement="LafT", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Tolq", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transporter ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="LafT protein", replacement="LafT", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="uncharacterized protein", replacement="unk", x=tmpstr)


		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)

		tmpstr = gsub(pattern=" ", replacement="_", x=tmpstr)

		tmpstr = gsub(pattern="MotA_MotA_component_of_the_H\\+coupled_stator_flagellum_complex", replacement="MotA", x=tmpstr)
		tmpstr = gsub(pattern="MotA_MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Mot_family_proton_H\\+_or_sodium_Na\\+_transporter_MotA", replacement="MotA", x=tmpstr)
		tmpstr = gsub(pattern="related_to_flagellar_apparatus_MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="related_to_ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transporter_ExbD", replacement="ExbD", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="conserved_TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MotA_MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		


		if ((tmpstr == "") || (tmpstr == "transporter"))
			{
			tmpstr = "transport"
			short_protname[i] = tmpstr
			next()
			}

		
		short_protname[i] = tmpstr
		}
	return(short_protname)
	} # END classify_MotAfam_labels <- function(list_of_strings)





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
	if (length(tmpstr) == 0)
		{
		spname = gdata::trim(stringr::str_squish(spname_line))
		} else {
		spname_line = gsub(pattern=tmpstr, replacement="", x=spname_line)
		spname_line = stringr::str_squish(spname_line)
		spname = trim(gsub(pattern="\\(\\)", replacement="", x=spname_line))
		spname
		}
	return(spname)
	}


#######################################################
# Get the sequence lengths (ignoring indels "-")
#######################################################
get_seqlengths <- function(aln)
	{
	seqlengths = rep(0, times=length(aln))
	for (i in 1:length(aln))
		{
		# Length of this row of the alignment, minus indels ("-")
		seqlengths[i] = length(aln[[i]]) - sum(aln[[i]] == "-")
		}
	return(seqlengths)
	}


# Handy function
extract_last_brackets <- function(list_of_strings, replace_spaces=TRUE)
	{
	example_code='
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	list_of_strings = tmptxt; replace_spaces=TRUE
	species_names = extract_last_brackets(list_of_strings=list_of_strings, replace_spaces=replace_spaces)
	species_names
	'
	
	
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



firstword <- function(string, split=" ")
	{
	example_code='
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	firstword(tmptxt, split=" ")
	'

	words = strsplit(gdata::trim(string), split=split)[[1]]
	return(gdata::trim(words[1]))
	}


lastword <- function(string, split="/")
	{
	example_code='
	tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
	lastword(tmptxt, split=" ")
	
	# Run on multiple inputs
	tmptxts = c(tmptxt, tmptxt)
	results = sapply(X=tmptxts, FUN=lastword, split=" ")
	results
	unname(results)
	'
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
	} # END convert_prokka_tsv_to_prot_feature_table <- function(tsv_fn, write_to_file=TRUE)








get_phylum_from_genome_ID <- function(assembly, genomes_to_spnames_df, printwarnings=FALSE)
	{
	junk='
	wd = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/"
	setwd(wd)
	genomes_to_spnames_fn = "species_list_10062023_NJM+group+spname_v1.txt"
	genomes_to_spnames_df = read.table(genomes_to_spnames_fn, header=TRUE, comment.char="%", quote="", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
	assembly = "GCA_000012685.1"
	phylum = get_phylum_from_genome_ID(assembly, genomes_to_spnames_df)	
	'

	assembly = firstword(assembly, split="\\.")
	
	# Find assembly in genome GenBank IDs
	genomes_to_spnames_matchnum = grep(pattern=assembly, x=genomes_to_spnames_df$GenBank.ID)

	if (length(genomes_to_spnames_matchnum) == 1)
		{
		phylum = genomes_to_spnames_df$Phylum[genomes_to_spnames_matchnum]
		}

	if (length(genomes_to_spnames_matchnum) > 1)
		{
		txt = paste0("Warning: more than 1 match to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. Probably this is due to multiple chromosomes/plasmids. Taking first hit, this is sufficient to identify the phylum. Taking first hit, but printing matches:")
		
		if (printwarnings == TRUE)
			{
			cat("\n")
			cat(txt)
			cat("\n")
			print(genomes_to_spnames_df[genomes_to_spnames_matchnum,])
			cat("\n")
			}
		genomes_to_spnames_matchnum = genomes_to_spnames_matchnum[1]
		phylum = genomes_to_spnames_df$Phylum[genomes_to_spnames_matchnum]
		warning(txt)
		}

	if (length(genomes_to_spnames_matchnum) == 0)
		{
		# Find assembly in genome RefSeq IDs
		genomes_to_spnames_matchnum = grep(pattern=assembly, x=genomes_to_spnames_df$RefSeq)
		if (length(genomes_to_spnames_matchnum) == 1)
			{
			phylum = genomes_to_spnames_df$Phylum[genomes_to_spnames_matchnum]
			}
		if (length(genomes_to_spnames_matchnum) > 1)
			{
			txt = paste0("Warning: more than 1 match to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. Probably this is due to multiple chromosomes/plasmids. Taking first hit, this is sufficient to identify the phylum, but printing matches:")

			if (printwarnings == TRUE)
				{
				cat("\n")
				cat(txt)
				cat("\n")
				print(genomes_to_spnames_df[genomes_to_spnames_matchnum,])
				cat("\n")
				}
			genomes_to_spnames_matchnum = genomes_to_spnames_matchnum[1]
			phylum = genomes_to_spnames_df$Phylum[genomes_to_spnames_matchnum]
			warning(txt)
			}
	
		if (length(genomes_to_spnames_matchnum) == 0)
			{
			txt = paste0("Warning: 0 matches to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. The 'Phylum' part of the label will be blank.")

			if (printwarnings == TRUE)
				{
				cat("\n")
				cat(txt)
				cat("\n")
				}
			warning(txt)
			phylum = ""
			}
		} # END if (length(genomes_to_spnames_matchnum) == 0)
	return(phylum)
	} # END get_phylum_from_genome_ID <- function(assembly, genomes_to_spnames_df)



get_group_from_genome_ID <- function(assembly, genomes_to_spnames_df, printwarnings=FALSE)
	{
	junk='
	wd = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/"
	setwd(wd)
	genomes_to_spnames_fn = "species_list_10062023_NJM+group+spname_v1.txt"
	genomes_to_spnames_df = read.table(genomes_to_spnames_fn, header=TRUE, comment.char="%", quote="", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
	assembly = "GCA_000012685.1"
	group = get_group_from_genome_ID(assembly, genomes_to_spnames_df)	
	'

	assembly = firstword(assembly, split="\\.")
	
	# Find assembly in genome GenBank IDs
	genomes_to_spnames_matchnum = grep(pattern=assembly, x=genomes_to_spnames_df$GenBank.ID)

	if (length(genomes_to_spnames_matchnum) == 1)
		{
		group = genomes_to_spnames_df$group[genomes_to_spnames_matchnum]
		}

	if (length(genomes_to_spnames_matchnum) > 1)
		{
		txt = paste0("Warning: more than 1 match to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. Probably this is due to multiple chromosomes/plasmids. Taking first hit, this is sufficient to identify the group. Taking first hit, but printing matches:")
		
		if (printwarnings == TRUE)
			{
			cat("\n")
			cat(txt)
			cat("\n")
			print(genomes_to_spnames_df[genomes_to_spnames_matchnum,])
			cat("\n")
			}
		genomes_to_spnames_matchnum = genomes_to_spnames_matchnum[1]
		group = genomes_to_spnames_df$group[genomes_to_spnames_matchnum]
		warning(txt)
		}

	if (length(genomes_to_spnames_matchnum) == 0)
		{
		# Find assembly in genome RefSeq IDs
		genomes_to_spnames_matchnum = grep(pattern=assembly, x=genomes_to_spnames_df$RefSeq)
		if (length(genomes_to_spnames_matchnum) == 1)
			{
			group = genomes_to_spnames_df$group[genomes_to_spnames_matchnum]
			}
		if (length(genomes_to_spnames_matchnum) > 1)
			{
			txt = paste0("Warning: more than 1 match to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. Probably this is due to multiple chromosomes/plasmids. Taking first hit, this is sufficient to identify the group, but printing matches:")

			if (printwarnings == TRUE)
				{
				cat("\n")
				cat(txt)
				cat("\n")
				print(genomes_to_spnames_df[genomes_to_spnames_matchnum,])
				cat("\n")
				}
			genomes_to_spnames_matchnum = genomes_to_spnames_matchnum[1]
			group = genomes_to_spnames_df$group[genomes_to_spnames_matchnum]
			warning(txt)
			}
	
		if (length(genomes_to_spnames_matchnum) == 0)
			{
			txt = paste0("Warning: 0 matches to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. The 'group' part of the label will be blank.")

			if (printwarnings == TRUE)
				{
				cat("\n")
				cat(txt)
				cat("\n")
				}
			warning(txt)
			group = ""
			}
		} # END if (length(genomes_to_spnames_matchnum) == 0)
	return(group)
	} # END get_group_from_genome_ID <- function(assembly, genomes_to_spnames_df)



get_spname_from_genome_ID <- function(assembly, genomes_to_spnames_df, printwarnings=FALSE)
	{
	junk='
	wd = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/"
	setwd(wd)
	genomes_to_spnames_fn = "species_list_10062023_NJM+spname+spname_v1.txt"
	genomes_to_spnames_df = read.table(genomes_to_spnames_fn, header=TRUE, comment.char="%", quote="", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
	assembly = "GCA_000012685.1"
	spname = get_spname_from_genome_ID(assembly, genomes_to_spnames_df)	
	'

	assembly = firstword(assembly, split="\\.")
	
	# Find assembly in genome GenBank IDs
	genomes_to_spnames_matchnum = grep(pattern=assembly, x=genomes_to_spnames_df$GenBank.ID)

	if (length(genomes_to_spnames_matchnum) == 1)
		{
		spname = genomes_to_spnames_df$spname[genomes_to_spnames_matchnum]
		}

	if (length(genomes_to_spnames_matchnum) > 1)
		{
		txt = paste0("Warning: more than 1 match to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. Probably this is due to multiple chromosomes/plasmids. Taking first hit, this is sufficient to identify the spname. Taking first hit, but printing matches:")
		
		if (printwarnings == TRUE)
			{
			cat("\n")
			cat(txt)
			cat("\n")
			print(genomes_to_spnames_df[genomes_to_spnames_matchnum,])
			cat("\n")
			}
		genomes_to_spnames_matchnum = genomes_to_spnames_matchnum[1]
		spname = genomes_to_spnames_df$spname[genomes_to_spnames_matchnum]
		warning(txt)
		}

	if (length(genomes_to_spnames_matchnum) == 0)
		{
		# Find assembly in genome RefSeq IDs
		genomes_to_spnames_matchnum = grep(pattern=assembly, x=genomes_to_spnames_df$RefSeq)
		if (length(genomes_to_spnames_matchnum) == 1)
			{
			spname = genomes_to_spnames_df$spname[genomes_to_spnames_matchnum]
			}
		if (length(genomes_to_spnames_matchnum) > 1)
			{
			txt = paste0("Warning: more than 1 match to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. Probably this is due to multiple chromosomes/plasmids. Taking first hit, this is sufficient to identify the spname, but printing matches:")

			if (printwarnings == TRUE)
				{
				cat("\n")
				cat(txt)
				cat("\n")
				print(genomes_to_spnames_df[genomes_to_spnames_matchnum,])
				cat("\n")
				}
			genomes_to_spnames_matchnum = genomes_to_spnames_matchnum[1]
			spname = genomes_to_spnames_df$spname[genomes_to_spnames_matchnum]
			warning(txt)
			}
	
		if (length(genomes_to_spnames_matchnum) == 0)
			{
			txt = paste0("Warning: 0 matches to assembly='", assembly, "' in '", genomes_to_spnames_fn, "'. The 'spname' part of the label will be blank.")

			if (printwarnings == TRUE)
				{
				cat("\n")
				cat(txt)
				cat("\n")
				}
			warning(txt)
			spname = ""
			}
		} # END if (length(genomes_to_spnames_matchnum) == 0)
	return(spname)
	} # END get_spname_from_genome_ID <- function(assembly, genomes_to_spnames_df)




# Trim each column of a data.frame
# # E.g. openxls leaves some annoying \t characters
trim_each_column <- function(xlsx)
	{
	example_code='
	xlsfn = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/species_list_10062023_NJM.xlsx"
	xlsx[4:8,1:5]	
	xlsx = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet=1)
	xlsx[4:8,1:5]	
	'

	for (i in 1:ncol(xlsx))
		{
		xlsx[,i] = stringr::str_trim(xlsx[,i])
		}
	return(xlsx)
	} # END trim_each_column <- function(xlsx)



all_but_suffix <- function(fn, split="\\.")
	{
	words = strsplit(fn, split=split)[[1]]
	word_to_remove = paste0(split, words[length(words)])
	subfn = gsub(pattern=word_to_remove, replacement="", x=fn)
	return(subfn)
	}

