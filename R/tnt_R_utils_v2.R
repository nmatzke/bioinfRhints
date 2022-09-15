# TNT-manipulation functions for R
# Load with:
#sourcedir = '/drives/Dropbox/_njm/'
#source4 = 'tnt_R_utils_v1.R'
#source(paste(sourcedir, source4, sep=""))


#######################################################
# Generic utilities
#######################################################

# Run the command-line sed on a file to find/replace
rsed <- function(oldtxt, newtxt, infn, outfn="", print_warnings=TRUE)
	{
	if (print_warnings == TRUE)
		{
		txt = "Warning: running system command 'sed' from R. This is a UNIX command, so\nthis should work on Linux and Mac OS X. On windows, you will\nprobably have to manually download/install a 'sed.exe' program and then make sure its PATH is accessible."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		} # END if (print_warnings == TRUE)
	
	# Assemble the command
	
	# Convert " " to "\w"
	#oldtxt = gsub(pattern=" ", replacement="\\\\w", x=" ;", perl=TRUE)
	
	# Redirect so if sed fails it doesn't wipe out file:
	# http://stackoverflow.com/questions/2585438/redirect-output-from-sed-s-c-d-myfile-to-myfile
	
	if (outfn == "")
		{
		cmdtxt = paste0("sed 's/", oldtxt, "/", newtxt, "/g' ", infn, " > tmp && mv tmp ", infn)
		} else {
		cmdtxt = paste0("sed 's/", oldtxt, "/", newtxt, "/g' ", infn, " > tmp && mv tmp ", outfn)
		}
	

	if (print_warnings == TRUE)
		{
		cat("\n")
		cat("Running this system command at the terminal:\n")
		cat(cmdtxt)
		cat("\n\n")
		}
	
	# Run the command
	system(cmdtxt)
	} # END rsed <- function(oldtxt, newtxt, infn, outfn="", print_warnings=TRUE)


# Read a TNT-output NEXUS tree file
# Using 'sed' (in UNIX terminal) or gsub (inside of R) to correct e.g. "begin ;" -> "begin;"
# correction="sed_permafix" means: use sed, overwrite the input NEXUS file with the new one
# correction="sed_permafix" means: use sed, save corrected file in outfn
# correction="gsub_permafix" means: use gsub in R, overwrite the input NEXUS file with the new one
# correction="gsub_nofix" means: use gsub in R,, save corrected file in outfn
read_TNT_nexus <- function(file, tree.names=NULL, correction="sed_nofix", outfn="tmp.nex", print_warnings=FALSE) 
	{
	if (correction == "sed_permafix")
		{
		# Run UNIX 'sed' command to fix
		rsed(oldtxt=" ;", newtxt=";", infn=file, outfn="", print_warnings=print_warnings)
		outfn = file
		}

	if (correction == "sed_nofix")
		{
		# Run UNIX 'sed' command to fix
		rsed(oldtxt=" ;", newtxt=";", infn=file, outfn=outfn, print_warnings=print_warnings)
		}

	if (correction == "gsub_permafix")
		{
		lines = readLines(file, warn=print_warnings)
		lines = gsub(pattern=" ;", replacement=";", x=lines)
		writeLines(text=lines, con=file)
		}

	if (correction == "gsub_nofix")
		{
		lines = readLines(file, warn=print_warnings)
		lines = gsub(pattern=" ;", replacement=";", x=lines)
		if (outfn == "")
			{
			outfn = "tmp.nex"
			}
		writeLines(text=lines, con=outfn)
		}

	
	# Otherwise, do nothing...
	
	# Then read the file like usual:
	# Read NEXUS like usual
	trs = read.nexus(file=outfn, tree.names=tree.names)
	
	# Remove any final "/" tags
	if (class(trs) == "multiPhylo")
		{
		for (i in 1:length(trs))
			{
			trs[[i]]$tip.label = sapply(X=trs[[i]]$tip.label, FUN=remove_endchar, pattern="/")
			trs[[i]]$node.label = sapply(X=trs[[i]]$node.label, FUN=remove_endchar, pattern="/")
			} # END for (i in 1:length(trs))
		
		
		} else {
		trs$tip.label = sapply(X=trs$tip.label, FUN=remove_endchar, pattern="/")
		trs$node.label = sapply(X=trs$node.label, FUN=remove_endchar, pattern="/")
		} # END if (class(trs) == "multiPhylo")
		
	return(trs)
	} # END read_TNT_nexus <- function(file, tree.names=NULL, correction="sed_permafix", outfn="tmp.nex", print_warnings=FALSE) 


# Remove the first character, if it matches the pattern
remove_startchar <- function(txt, pattern="=")
	{	
	startval = 1
	stopval = 1
	firstchar = substr(x=txt, start=startval, stop=stopval)
	if (firstchar == pattern)
		{
		newtxt = substr(x=txt, start=2, stop=nchar(txt))
		} else {
		newtxt = txt
		}
	return(newtxt)	
	} # END remove_endchar <- function(txt, pattern="/")



# Remove the last character, if it matches the pattern
remove_endchar <- function(txt, pattern="/")
	{	
	startval = nchar(txt)
	stopval = nchar(txt)
	lastchar = substr(x=txt, start=startval, stop=stopval)
	if (lastchar == pattern)
		{
		newtxt = substr(x=txt, start=1, stop=stopval-1)
		} else {
		newtxt = txt
		}
	return(newtxt)	
	} # END remove_endchar <- function(txt, pattern="/")


# Does the line start with pattern?
string_startchar <- function(txt, pattern="(")
	{
	firstchar = substr(x=txt, start=1, stop=1)
	if (firstchar == pattern)
		{
		return(TRUE)
		} else {
		return(FALSE)
		}
	} # END string_startchar <- function(txt, pattern="(")

# Does the end with pattern?
string_endchar <- function(txt, pattern=")")
	{
	lastchar = substr(x=txt, start=nchar(txt), stop=nchar(txt))
	if (lastchar == pattern)
		{
		return(TRUE)
		} else {
		return(FALSE)
		}
	} # END string_startchar <- function(txt, pattern="(")





#######################################################
# Write a dataframe to TNT
#######################################################
# data_df is assumed to have *COLS* as OTUs, *ROWS* as characters
# 
# nstates can be e.g. "8", "dna", etc.
# 
# NOTE: YOU MAY NEED TO CONVERT ( AND ) TO [ AND ] TO MAKE 
# THE FILE TNT-READABLE
#
df2tnt <- function(data_df, nstates="8", tntfn="data_df.tnt", truncate_OTUnames=FALSE, maxlength=32, force_unique=TRUE)
	{
	data_df = t(data_df)
	
	numrows = nrow(data_df)
	numcols = ncol(data_df)
	OTU_names = row.names(data_df)
	OTU_names = check_OTU_name_length(OTU_names, truncate_OTUnames=truncate_OTUnames, maxlength=maxlength, force_unique=force_unique)
	row.names(data_df) = OTU_names
	
	
	# Write the header
	txt = paste0("nstates ", nstates, " ;")
	write(x=txt, file=tntfn, append=FALSE)
	txt = paste0("xread 'Data converted from R data frame to TNT format by df2tnt()'")
	write(x=txt, file=tntfn, append=TRUE)
	txt = paste0(numcols, " ", numrows)
	write(x=txt, file=tntfn, append=TRUE)
	
	# Collapse the data
	chars_txt = apply(X=data_df, MARGIN=1, FUN=paste0, collapse="")
	OTU_names_plus_chars = cbind(OTU_names, chars_txt)
	write.table(x=OTU_names_plus_chars, file=tntfn, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	
	txt = ";"
	write(x=txt, file=tntfn, append=TRUE)
	txt = ""
	write(x=txt, file=tntfn, append=TRUE)
	txt = "proc /;"
	write(x=txt, file=tntfn, append=TRUE)
	txt = "comments 0"
	write(x=txt, file=tntfn, append=TRUE)
	txt = ";"
	write(x=txt, file=tntfn, append=TRUE)
	
	return(tntfn)
	} # END df2tnt <- function(data_df, nstates="8", tntfn="data_df.tnt", truncate_OTUnames=FALSE)


check_OTU_name_length <- function(OTU_names, truncate_OTUnames=FALSE, maxlength=32, force_unique=TRUE)
	{
	# Check OTU_names length
	lengths = sapply(X=OTU_names, FUN=nchar)
	if (any(lengths > maxlength))
		{
		TF = lengths > maxlength
		OTUnames_to_chop = OTU_names[TF]
		
		# Truncate the names if desired
		if (truncate_OTUnames == TRUE)
			{
			OTUnames_to_chop = chop_OTUnames(OTUnames_to_chop=OTUnames_to_chop, maxlength=maxlength, force_unique=force_unique)
			
			numchanged = length(OTUnames_to_chop)
			
			txt = paste0("WARNING: ", numchanged, " OTU names were greater than maxlength=", maxlength, ". These were changed:\n\n")
			cat("\n\n")
			cat(txt)
			
			for (i in 1:numchanged)
				{
				txt = paste0(OTU_names[TF][i], "\t\t->\t\t", OTUnames_to_chop[i])
				cat(txt)
				cat("\n")
				} # END for (i in 1:numchanged)
			} # END if (truncate_OTUnames == TRUE)
		OTU_names[TF] = OTUnames_to_chop
		}
	return(OTU_names)
	} # END check_OTU_name_length <- function(OTU_names, truncate_OTUnames=FALSE, maxlength=32, force_unique=TRUE)


chop_OTUnames <- function(OTUnames_to_chop, maxlength=32, force_unique=TRUE)
	{
	OTUnames_chopped = sapply(X=OTUnames_to_chop, FUN=chop_OTUname, maxlength=maxlength)
	
	# Force them to be unique, if desired
	if (force_unique == TRUE)
		{
		uniq_OTUs = unique(OTUnames_chopped)
		if (length(uniq_OTUs) != length(OTUnames_chopped))
			{
			# Get counts of each OTU
			counts = c(unlist(uniq_OTUs))
			TF = counts > 1
			OTUs_to_shorten = uniq_OTUs[TF]
			counts = counts[TF]
			
			for (i in 1:length(OTUs_to_shorten))
				{
				num_matches = counts[i]
				
				# Get the length of the nums to add
				length_nums = 1+nchar(as.character(num_matches))
				
				# Go through the identical_OTUs
				identical_OTUs_TF = OTUnames_chopped == OTUs_to_shorten[i]
				
				for (j in 1:num_matches)
					{
					tmp_OTUname = OTUnames_chopped[identical_OTUs_TF][j]
					tmp_nchars = nchar(tmp_OTUname)
					
					length_final = tmp_nchars + length_nums
					
					# Make the number suffix
					fmt_string = paste0("%", (length_nums-1), ".0f")
					newnum1 = sprintf(fmt=fmt_string, j)
					newnum2 = paste0("_", newnum1)
					
					# Do one thing if adding numbers makes it too long
					if (length_final > maxlength)
						{
						tmp_OTUname2 = substr(x=tmp_OTUname, start=1, stop=maxlength-length_nums)
						
						new_OTUname = paste0(tmp_OTUname2, newnum2)
						} else {
						# Just paste the numbers on
						new_OTUname = paste0(tmp_OTUname, newnum2)
						} # END if (length_final > maxlength)
					
					# Save the new name
					OTUnames_chopped[identical_OTUs_TF][j] = new_OTUname
					} # END for (j in 1:num_matches)
				} # END for (i in 1:length(OTUs_to_shorten))
			} # END if (length(uniq_OTUs) == length(OTUnames_chopped))
		} # END if (force_unique == TRUE)
	
	# Return the chopped names
	return(OTUnames_chopped)
	} # END chop_OTUnames <- function(OTUnames_to_chop, maxlength=32, force_unique=TRUE)

chop_OTUname <- function(OTUname_to_chop, maxlength=32)
	{
	OTUname_chopped = substr(x=OTUname_to_chop, start=1, stop=maxlength)
	return(OTUname_chopped)
	} # END chop_OTUname <- function(OTUname_to_chop, maxlength=32)

###########################################
# TNT functions
###########################################
# options = text_convert, via_Rphylo, nexus

setup = '
tntfn = "majority_bremer.tnttree"
'
tntfile2newick <- function(tntfn, brlens=FALSE, options="text_convert", translate=FALSE, branchlabels="=")
	{
	# Convert a TNT output tree to Newick format.
	# Following the directions online here:
	# http://tnt.insectmuseum.org/index.php/FAQ#Can_TNT_export_Newick_formatted_trees.3F
	#
	
	# Read a text file into a list of strings
	lines = scan(tntfn, what="character", sep="\n")
	
	#tntstr = lines[2]
	
	TF = sapply(X=lines, FUN=string_startchar, pattern="(")
	lines = lines[TF]
	numlines = length(lines)

	if (numlines < 1)
		{
		txt = paste0("No lines found that start with '(', so no TNT trees found in '", tntfn, "'. Returning NULL.")
		cat("\n")
		cat("txt")
		cat("\n")
		return(NULL)
		}
	

	if (numlines == 1)
		{
		tntstr = lines
		newickstr = tntstr2newick(tntstr, brlens=brlens, branchlabels=branchlabels)
		#tnttr = read.tree(file="", text=newickstr)
		
		if (options == "text_convert")
			{
			newick_outfn = paste(tntfn, ".newick", sep="")

			# write straight to text
			write(newickstr, file=newick_outfn, sep="\n")
			cat("\nWrote ", tntfn, " to ", newick_outfn, " (via text conversion)\n\n", sep="")
			return(newick_outfn)
			}

		if (options == "via_Rphylo")
			{
			newick_outfn = paste(tntfn, ".newick", sep="")

			# write to newick text via R
			tr = read.tree(file="", text=newickstr)
			
			write.tree(tr, file=newick_outfn)
			cat("\nWrote ", tntfn, " to ", newick_outfn, " (via R phylo object)\n\n", sep="")
			return(newick_outfn)
			}
			
		if (options == "nexus")
			{
			nexus_outfn = paste(tntfn, ".nex", sep="")

			# write to newick text via R
			tr = read.tree(file="", text=newickstr)
			
			write.nexus(tr, file=nexus_outfn, translate=translate)
			cat("\nWrote ", tntfn, " to ", nexus_outfn, " (via R phylo object to NEXUS)\n\n", sep="")
			return(nexus_outfn)
			} # END if (nexus == FALSE)
		} 
	
	if (numlines > 1)
		{
		# Otherwise...
		cat(paste0("\n\ntntfile2R() detected ", numlines, " TNT trees...\n\n"))
		trs = NULL
		for (i in 1:numlines)
			{
			txt = paste0("Reading TNT tree #", i, "/", numlines, "...\n")
			tntstr = lines[i]
			newickstr = tntstr2newick(tntstr, brlens=brlens, branchlabels=branchlabels)
			tnttr = read.tree(file="", text=newickstr)
			trs[[i]] = tnttr
			} # END for (i in 1:numlines)
		class(trs) = "multiPhylo"
		
		if (options != "nexus")
			{
			newick_outfn = paste(tntfn, ".newick", sep="")
			write.tree(trs, file=newick_outfn)
			cat("\nWrote ", tntfn, " to ", newick_outfn, " (via R phylo object)\n\n", sep="")
			return(newick_outfn)
			} else {
			nexus_outfn = paste(tntfn, ".nex", sep="")
			write.nexus(trs, file=nexus_outfn, translate=translate)
			cat("\nWrote ", tntfn, " to ", nexus_outfn, " (via R phylo object to NEXUS)\n\n", sep="")
			return(nexus_outfn)
			} # END if (options != "nexus")
		} # END if (numlines > 1)
	
	stop("\n\nSTOP ERROR: tntfile2newick shouldn't get here!\n\n")
	} # END tntfile2newick <- function(tntfn, brlens=FALSE, options="text_convert", translate=FALSE)


tntfile2R <- function(tntfn, brlens=FALSE, branchlabels="=", commas_to="|")
	{
	# Convert a TNT output tree to Newick format & read into R
	# Following the directions online here:
	# http://tnt.insectmuseum.org/index.php/FAQ#Can_TNT_export_Newick_formatted_trees.3F
	#
	
	# Read a text file into a list of strings
	lines = scan(tntfn, what="character", sep="\n")
	
	# If you have e.g. synapomorphies, there could be commas.
	# These need to be removed.
	lines = gsub(pattern=",", replacement="|", x=lines)
	
	#tntstr = lines[2]
	# Lines that start with "("
	TF = sapply(X=lines, FUN=string_startchar, pattern="(")
	lines = lines[TF]
	numlines = length(lines)

	if (numlines < 1)
		{
		txt = paste0("No lines found that start with '(', so no TNT trees found in '", tntfn, "'. Returning NULL.")
		cat("\n")
		cat("txt")
		cat("\n")
		return(NULL)
		}
	

	if (numlines == 1)
		{
		tntstr = lines
		newickstr = tntstr2newick(tntstr, brlens=brlens, branchlabels=branchlabels)
		
		tnttr = read.tree(file="", text=newickstr)
		
		if (is.null(tnttr$node.label))
			{
			ntips = length(tnttr$tip.label)
			nodenums = (ntips+1) : (ntips+tnttr$Nnode)
			tnttr$node.label = nodenums
			} # END if (is.null(tnttr$node.label))
		
		# Remove any "=" tags (for brlens==FALSE
		if (brlens == FALSE)
			{
			# Node labels: easy
			tnttr$node.label = sapply(X=tnttr$node.label, FUN=remove_startchar, pattern=branchlabels)
			
			# Probably pointless; tip labels have to be post-processed later
			tnttr$tip.label = sapply(X=tnttr$tip.label, FUN=remove_startchar, pattern=branchlabels)
			
			# Tip labels: split on "="
#			tip_labels_wEquals_TF = grepl(pattern=branchlabels, x=tnttr$tip.label)
#			numhits = sum(tip_labels_wEquals_TF)
#			if (numhits > 0)
#				{
#				data_below_tip = rep("", times=numhits)
#				for (i in 1:numhits)
#					{
#					# Convert the labels by removing "="
#					words = strsplit(x=tnttr$tip.label[tip_labels_wEquals_TF][i], split=branchlabels)[[1]]
#					
#					# Store the OTU name
#					tnttr$tip.label[tip_labels_wEquals_TF][i]
#					
#					if (length(words) > 1)
#						{
#						tmp = paste0(words[2:length(words)], collapse="_", sep="")
#						data_below_tip[i] = tmp
#						}
#					}
#				}
			
			} # END if (brlens == FALSE)
		return(tnttr)
		}
		
	# Otherwise...
	cat(paste0("\n\ntntfile2R() detected ", numlines, " TNT trees...\n\n"))
	trs = NULL
	for (i in 1:numlines)
		{
		txt = paste0("Reading TNT tree #", i, "/", numlines, "...\n")
		tntstr = lines[i]
		newickstr = tntstr2newick(tntstr, brlens=brlens)
		tnttr = read.tree(file="", text=newickstr)

		if (is.null(tnttr$node.label))
			{
			ntips = length(tnttr$tip.label)
			nodenums = (ntips+1) : (ntips+tnttr$Nnode)
			tnttr$node.label = nodenums
			} # END if (is.null(tnttr$node.label))


		# Convert the labels by removing "="
		tnttr$node.label = sapply(X=tnttr$node.label, FUN=remove_startchar, pattern=branchlabels)
		# Probably pointless; tip labels have to be post-processed later
		tnttr$tip.label = sapply(X=tnttr$tip.label, FUN=remove_startchar, pattern=branchlabels)
		
		# Save
		trs[[i]] = tnttr
		} # END for (i in 1:numlines)
	
	class(trs) = "multiPhylo"
	return(trs)
	} # END tntfile2R <- function(tntfn, brlens=FALSE)


# Convert a Newick string to TNT string
# (make sure no branchlengths, not sure they would be read anyway)
newickstr2tntstr <- function(newickstr)
	{
	# 1. Change any ; to *
	tntstr1 = gsub(pattern=";", replacement="*", x=newickstr)
	
	# 2. Change any commas to whitespace
	tntstr2 = gsub(pattern=",", replacement=" ", x=tntstr1)

	# 3. # replace "),(" with ")("
	tntstr3 = gsub(pattern="\\),\\(", replacement="\\)\\(", x=tntstr2)

	# 3. # replace ") (" with ")("
	tntstr3a = gsub(pattern="\\) \\(", replacement="\\)\\(", x=tntstr3)

	# 4. replace ")" with " )"
	tntstr4 = gsub("\\)", " \\)", tntstr3a)
		
	# 4. replace ") )" with ")"
	tntstr5 = gsub("\\) \\)", "\\)\\)", tntstr4)
	tntstr6 = gsub("\\) \\)", "\\)\\)", tntstr5)
	
	tntstr = tntstr6
	
	return(tntstr6)
	}


# Write TNT strings to file
tntstr2file <- function(tntstr, file="tree.tnt")
	{
	numlines = length(tntstr)
	numlines = numlines + 4
	
	startline = "tread 'tree(s) from TNTR function tntstr2file() by Nick Matzke'"
	treelines = tntstr
	endlines = c(";", "proc-;", "")
	
	lines_to_write = c(startline, treelines, endlines)
	writeLines(text=lines_to_write, con=file)
	return(file)
	} # END tntstr2file <- function(tntstr, file="tree.tnt")


# Write tree(s) to TNT file or string
write_tree_TNT <- function(tr, file)
	{
	# Blank out the branch lengths, if any
	if (class(tr) == "multiPhylo")
		{
		for (i in 1:length(tr))
			{
			# Blank out the branch lengths, if any
			tr[[i]]$edge.length = NULL
			tr[[i]]$node.label = NULL
			} # END for (i in 1:length(tr))
		} else {
		# Blank out the branch lengths, if any
		tr$edge.length = NULL
		tr$node.label = NULL
		} # END if (class(tr) == "multiPhylo")
	
	# Write tree to Newick
	newickstr = write.tree(tr, file="")
	
	# Convert the Newick string to a TNT string
	tntstr = newickstr2tntstr(newickstr)
	
	# Return string if desired
	if (file == "")
		{
		return(tntstr)
		}
	
	# Otherwise...
	file = tntstr2file(tntstr=tntstr, file=file)
	
	return(file)
	}


sep_tipnames_branchlabels <- function(tip_labels, branchlabels="=")
	{
	# Check that everything has branchlabels
	# e.g., one "=" in each label
	hit_locations_on_each_tiplabel = gregexpr(pattern=branchlabels, text=tip_labels)
	numhits_on_branchlabels = sapply(X=hit_locations_on_each_tiplabel, FUN=length)
	
	# Also -- misses have length 1, but they have value -1. Set numhits to 0.
	numhits_on_branchlabels[sapply(X=hit_locations_on_each_tiplabel, FUN=isNegOne)] = 0
	TF = (numhits_on_branchlabels == 1)

	# Process the hits
	if (all(TF == TRUE))
		{
		# Put the data in a matrix format
		tmpmat = matrix(data=unlist(strsplit(x=tip_labels, split=branchlabels)), ncol=2, byrow=TRUE)
		} else {
		tmpmat = NULL
		for (i in 1:length(tip_labels))
			{
			words = strsplit(tip_labels[i], split=branchlabels)[[1]]
			if (length(words) == 1)
				{
				tmpline = c(words, "")
				}
			if (length(words) == 2)
				{
				tmpline = words
				}
			if (length(words) > 2)
				{
				tmpline = c(words[1], paste0(words[2:length(words)], collapse=branchlabels, sep=""))
				}
			tmpmat = rbind(tmpmat, tmpline)
			}
		}
	
	tip_branchlabels_df = as.data.frame(tmpmat, stringsAsFactors=FALSE, row.names=NULL)
	names(tip_branchlabels_df) = c("tipnames", "branch_data")
	row.names(tip_branchlabels_df) = NULL
	
	return(tip_branchlabels_df)
	} # END sep_tipnames_branchlabels <- function(tip_labels, branchlabels="=")



# Split a series of ttags on slash
# deals with inconsistencies, e.g.
# /
# /
# 100/90
# 100/90/80
split_ttags_on_slash <- function(ttags, slash="/")
	{
	defaults='
	ttags=tmp_labels
	slash="/"
	'
	
	hit_locations_on_each_tiplabel = gregexpr(pattern=slash, text=ttags)
	numhits_on_branchlabels = sapply(X=hit_locations_on_each_tiplabel, FUN=length)
	numhits_on_branchlabels[sapply(X=hit_locations_on_each_tiplabel, FUN=isNegOne)] = 0
	max_numhits = max(numhits_on_branchlabels)
	
	# Blank matrix
	ttags_mat = matrix(data="", ncol=(max_numhits+1), nrow=length(ttags))
	
	# Cells with 0 slashes just go into the first column
	TF0 = numhits_on_branchlabels == 0
	ttags_mat[TF0, ][,1] = ttags[TF0]
	
	# Fix various weirdness
	ttags = gsub(pattern="//", replacement="/_/", x=ttags)
	ttags = gsub(pattern="//", replacement="/_/", x=ttags)
	
	TF = sapply(X=ttags, FUN=string_endchar, pattern="/")
	TF[ttags="/"] = FALSE
	ttags[TF] = paste0(ttags[TF], "_")

	TF = sapply(X=ttags, FUN=string_startchar, pattern="/")
	TF[ttags="/"] = FALSE
	ttags[TF] = paste0("_", ttags[TF])

	
	for (i in 1:max_numhits)
		{
		TF1 = numhits_on_branchlabels == i
		TF2 = ttags != paste0(rep(slash, times=i), collapse="")
		TF = (TF1 + TF2) == 2
		
		# Put the data in a matrix format
		tmpmat = matrix(data=unlist(strsplit(x=ttags[TF], split=slash)), ncol=i+1, nrow=sum(TF), byrow=TRUE)
		# Tricky bit to get around single rows that are mis-cast as vectors rather than matrices
		ttags_mat2 = matrix(ttags_mat[TF,], ncol=ncol(ttags_mat), byrow=TRUE)
		ttags_mat2[, 1:(i+1)] = tmpmat
		ttags_mat[TF,] = ttags_mat2
		}
	
	ttags_df = as.data.frame(ttags_mat, stringsAsFactors=FALSE, row.names=NULL)
	names(ttags_df) = paste0("c", 1:(max_numhits+1), sep="")
	
	# Take out the "_"
	ttags_df[ttags_df=="_"] = ""
	
	return(ttags_df)
	} # END split_ttags_on_slash <- function(ttags, slash="/")


tntstr2newick <- function(tntstr, brlens=FALSE, branchlabels="=")
	{
	# Convert a TNT output tree (as a string) to Newick format.
	# Following the directions online here:
	# http://tnt.insectmuseum.org/index.php/FAQ#Can_TNT_export_Newick_formatted_trees.3F
	#
	
	# Remove " ;" from the string
	tntstr = gsub(pattern=" ;", replacement=";", x=tntstr)

	# Remove " *" from the string
	tntstr = gsub(pattern=" \\*", replacement="*", x=tntstr)

	
	if (brlens == FALSE)
		{
		# 1. replace any * at the end of the lines with ";" 
		tntchars = strsplit(tntstr, "")[[1]]
		if (tntchars[length(tntchars)] == "*")
			{
			tntchars[length(tntchars)] = ";"
			}
		tntstr2 = list2str(tntchars, spacer="")
		
		# 2. replace all whitespace with commas
		tntstr3 = gsub("[: :]", ",", tntstr2)
		
		# 3. # replace ")(" with "),("
		tntstr4 = gsub("\\)\\(", "\\),\\(", tntstr3)
		
		# 4. replace ",)" with ")"
		tntstr5 = gsub(",\\)", "\\)", tntstr4)
		
		# 5. Branch lengths on tip branches might get mis-translated as separate tips
		branchlabels_pattern = paste0(",", branchlabels)
		newickstr = gsub(pattern=branchlabels_pattern, replacement=branchlabels, x=tntstr5)
		}

	# Trees WITH branch lengths are slightly different
	if (brlens == TRUE)
		{
		# 1. replace any * at the end of the lines with ";" 
		tntchars = strsplit(tntstr, "")[[1]]
		if (tntchars[length(tntchars)] == "*")
			{
			tntchars[length(tntchars)] = ";"
			}
		tntstr2 = list2str(tntchars, spacer="")
		
		# 2. replace all whitespace with commas
		tntstr3 = gsub("[: :]", ",", tntstr2)
		
		# 3. # replace ")(" with "),("
		tntstr4 = gsub("\\)\\(", "\\),\\(", tntstr3)
		
		# 4. replace ",)" with ")"
		tntstr5 = gsub(",\\)", "\\)", tntstr4)

		# 5. replace "=" with ":"
		tntstr6 = gsub(branchlabels, ":", tntstr5)
		
		# 6. replace ",:" with ":"
		newickstr = gsub(",:", ":", tntstr6)
		newickstr
		}

		
	#plot(read.tree(text=tntstr5))
	
	return(newickstr)
	}









#######################################################
# Read the results of "auto.run"
#######################################################
read_auto_results <- function(settings=NULL)
	{
	defaults='
	wd = "/drives/GDrive/__GDrive_projects/2015-08-01_AFA_cladogram/_02_TNT/2015-08-25_runs/allchars_MP/"
	# Working directory
	settings$wd = "."
	# Filenames
	settings$script_fn = "auto.run"
	settings$MP_topologies_fn = "MP_topologies.tnt"
	settings$indata_fn = "tmpdata.tnt"
	#settings$intrees_fn = "tmptrees.tnt"
	settings$nodenums_fn = "auto_nodenums.tnt"
	settings$branchlengths_fn = "auto_branchlengths.tnt"
	settings$strictcon_fn = "auto_strict_consensus.tnt"
	settings$bremer_abs_fn = "auto_BremerAbsolute.tnt"
	settings$bremer_rel_fn = "auto_BremerRelative.tnt"
	settings$bootstraps_fn = "auto_Bootstraps.tnt"
	'
	
	if (is.null(settings))
		{
		# Working directory
		settings$wd = "."
		# Filenames
		settings$script_fn = "auto.run"
		settings$MP_topologies_fn = "MP_topologies.tnt"
		settings$indata_fn = "tmpdata.tnt"
		#settings$intrees_fn = "tmptrees.tnt"
		settings$nodenums_fn = "auto_nodenums.tnt"
		settings$synapos_fn = "auto_synapos.tnt"
		settings$branchlengths_fn = "auto_branchlengths.tnt"
		settings$strictcon_fn = "auto_strict_consensus.tnt"
		settings$bremer_abs_fn = "auto_BremerAbsolute.tnt"
		settings$bremer_rel_fn = "auto_BremerRelative.tnt"
		settings$bootstraps_fn = "auto_Bootstraps.tnt"
		} # END if (settings = NULL)
	
	# Abbreviate settings
	s = settings
	
	# Get original working directory
	orig_wd = getwd()
	# Change wd
	setwd(s$wd)
	
	cat("\nSearching for files in settings$wd=\n")
	cat(s$wd)
	cat("\n")
	
	if (file.exists(s$branchlengths_fn))
		{
		cat("\nSearching for '", s$branchlengths_fn, "', a TNT tree with branchlengths...found.")
	
		# Read the branchlengths tree first
		tr_branchlengths = tntfile2R(tntfn=s$branchlengths_fn, brlens=TRUE, branchlabels="=")
		tr_branchlengths
		trtable = prt(tr_branchlengths, printflag=FALSE, get_tipnames=TRUE)
		trtable
		} else {
		cat("\nSearching for '", s$branchlengths_fn, "', a TNT tree with branchlengths...not found. Instead...")
		if (file.exists(s$strictcon_fn))
			{
			cat("\nsearching for '", s$strictcon_fn, "', a TNT strict consensus tree without branchlengths...found.")
	
			# Read the branchlengths tree first
			tr_branchlengths = tntfile2R(tntfn=s$strictcon_fn, brlens=TRUE, branchlabels="=")
			tr_branchlengths
			trtable = prt(tr_branchlengths, printflag=FALSE, silence_warnings=FALSE)
			trtable
			} else {
			stoptxt = "STOP ERROR in read_auto_results().\nNeither settings$branchlengths_fn nor settings$strictcon_fn can be found."
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			} # END if (file.exists(s$strictcon_fn))
		} # END if (file.exists(s$branchlengths_fn))
	
	
	if (file.exists(s$nodenums_fn))
		{
		cat("\nSearching for '", s$nodenums_fn, "', a TNT tree with TNT node numbers...found.")
	
		# Read the branchlengths tree first
		tr_nodenums = tntfile2R(tntfn=s$nodenums_fn, brlens=FALSE, branchlabels="=")
		tr_nodenums
		tmptable = prt(tr_nodenums, printflag=FALSE, silence_warnings=TRUE)
	
		ntips = length(tr_nodenums$tip.label)
		internal_Rnodenums = (ntips+1):(ntips+tr_nodenums$Nnode)
		tip_Rnodenums = 1:ntips
	
		# Extract the TNT node numbers
		tnt_nodenums_internal = tmptable$label[internal_Rnodenums]
		tnt_nodenums_tips = sep_tipnames_branchlabels(tip_labels=tmptable$label[tip_Rnodenums], branchlabels="=")
		tnt_nodenums_tips
		TNT_nodenums = c(tnt_nodenums_tips$branch_data, tnt_nodenums_internal)
		
		# Add to trtable
		trtable = cbind(trtable, TNT_nodenums)
		} else {
		cat("\nsearching for '", s$nodenums_fn, "', a TNT tree with TNT node numbers...not found.")
		TNT_nodenums = NULL
		} # END if (file.exists(s$nodenums_fn))




	if (file.exists(s$bremer_abs_fn))
		{
		cat("\nSearching for '", s$bremer_abs_fn, "', a TNT tree with absolute Bremer supports...found.")
	
		# Read the branchlengths tree first
		tr_bremer_abs = tntfile2R(tntfn=s$bremer_abs_fn, brlens=FALSE, branchlabels="=")
		tr_bremer_abs
		tmptable = prt(tr_bremer_abs, printflag=FALSE, silence_warnings=TRUE)
	
		# Extract the TNT Bremer supports
		ntips = length(tr_bremer_abs$tip.label)
		internal_Rnodenums = (ntips+1):(ntips+tr_bremer_abs$Nnode)
		tnt_bremers_internal = tmptable$label[internal_Rnodenums]
		tnt_bremers_tips = rep("", times=ntips)
		absBremer = c(tnt_bremers_tips, tnt_bremers_internal)

		# In TNT Bremer support, something like "9?" means "9 or higher", i.e. 9+
		# So, we will change "?" to "+"
		absBremer = gsub(pattern="\\?", replacement="+", x=absBremer)
		
		# Add to trtable
		trtable = cbind(trtable, absBremer)
		} else {
		cat("\nSearching for '", s$bremer_abs_fn, "', a TNT tree with absolute Bremer supports...not found.")
		absBremer = NULL
		} # END if (file.exists(s$nodenums_fn))

	if (file.exists(s$bremer_rel_fn))
		{
		cat("\nSearching for '", s$bremer_rel_fn, "', a TNT tree with relative Bremer supports...found.")
	
		# Read the branchlengths tree first
		tr_bremer_rel = tntfile2R(tntfn=s$bremer_rel_fn, brlens=FALSE, branchlabels="=")
		tr_bremer_rel
		tmptable = prt(tr_bremer_rel, printflag=FALSE, silence_warnings=TRUE)
	
		# Extract the TNT Bremer supports
		ntips = length(tr_bremer_rel$tip.label)
		internal_Rnodenums = (ntips+1):(ntips+tr_bremer_rel$Nnode)
		tnt_bremers_internal = tmptable$label[internal_Rnodenums]
		tnt_bremers_tips = rep("", times=ntips)
		relBremer = c(tnt_bremers_tips, tnt_bremers_internal)
		
		# In TNT Bremer support, something like "9?" means "9 or higher", i.e. 9+
		# So, we will change "?" to "+"
		relBremer = gsub(pattern="\\?", replacement="+", x=relBremer)
		
		# Add to trtable
		trtable = cbind(trtable, relBremer)
		} else {
		cat("\nSearching for '", s$bremer_abs_fn, "', a TNT tree with relative Bremer supports...not found.")
		relBremer = NULL
		} # END if (file.exists(s$nodenums_fn))

	if (file.exists(s$bootstraps_fn))
		{
		cat("\nSearching for '", s$bootstraps_fn, "', a TNT tree with Bootstrap supports...found.")
		
		# Read the branchlengths tree first
		tr_bootstraps_abs = tntfile2R(tntfn=s$bootstraps_fn, brlens=FALSE, branchlabels="=")
		tr_bootstraps_abs
		tmptable = prt(tr_bootstraps_abs, printflag=FALSE, silence_warnings=TRUE)
	
		# Extract the TNT Bremer supports
		ntips = length(tr_bootstraps_abs$tip.label)
		internal_Rnodenums = (ntips+1):(ntips+tr_bootstraps_abs$Nnode)
		tip_Rnodenums = 1:ntips

		tnt_bootstraps_internal = tmptable$label[internal_Rnodenums]
		tnt_bootstraps_tips = sep_tipnames_branchlabels(tip_labels=tmptable$label[tip_Rnodenums], branchlabels="=")
		
		ttags = c(tnt_bootstraps_tips$branch_data, tnt_bootstraps_internal)
		ttags
		
		ttags_df = split_ttags_on_slash(ttags=ttags, slash="/")
		bootnames = c("freqBoot", "gcBoot", "otherBoot", "other1", "other2", "other3", "other4", "other5")
		
		names(ttags_df) = bootnames[1:ncol(ttags_df)]
		ttags_df
		
		# Add to trtable
		trtable = cbind(trtable, ttags_df)
		} else {
		cat("\nSearching for '", s$bootstraps_fn, "', a TNT tree with Bootstrap supports...not found.")
		bootstraps = NULL
		} # END if (file.exists(s$nodenums_fn))



	if (file.exists(s$synapos_fn))
		{
		cat("\nSearching for '", s$synapos_fn, "', a TNT tree with synapomorphies plotted...found.")
	
		# Read the branchlengths tree first
		tr_synapos = tntfile2R(tntfn=s$synapos_fn, brlens=FALSE, branchlabels="=", commas_to="|")
		tr_synapos
		tmptable = prt(tr_synapos, printflag=FALSE, silence_warnings=TRUE, get_tipnames=TRUE)
		

	
		ntips = length(tr_synapos$tip.label)
		internal_Rnodenums = (ntips+1):(nrow(tmptable))
		tip_Rnodenums = 1:ntips
	
		# It may be that some tips have synapomorphies attached, and others don't
		# Fix this here
		tnt_synapos_tips = sep_tipnames_branchlabels(tip_labels=tmptable$label[tip_Rnodenums], branchlabels="=")
		
		# Extract the TNT synapomorphies
		tnt_synapos_internal = tmptable$label[internal_Rnodenums]
		TNT_synapos = c(tnt_synapos_tips$branch_data, tnt_synapos_internal)
		
		# Add to trtable
		# THIS CRASHES SOMETIMES DUE TO TNT INCONSISTENCY
		
		# Make new trtable
		tr_synapos$tip.label = tnt_synapos_tips$tipnames
		tmptable = prt(tr_synapos, printflag=FALSE, silence_warnings=TRUE, get_tipnames=TRUE)
		
		#trtable = cbind(trtable, TNT_synapos)
		rows_in_trtable = match(x=tmptable$tipnames, table=trtable$tipnames)
		rows_in_trtable = rows_in_trtable[!is.na(rows_in_trtable)]
		TNT_synapos = TNT_synapos[rows_in_trtable]
		
		# Add to trtable
		trtable = cbind(trtable, TNT_synapos)
		} else {
		cat("\nsearching for '", s$synapos_fn, "', a TNT tree with synapomorphies plotted...not found.")
		TNT_synapos = NULL
		} # END if (file.exists(s$synapos_fn))

	cat("\n\n")

	# Change back to original working directory
	setwd(orig_wd)
	
	# Change anything that is a factor, back into a character
	classes = cls.df(trtable)
	TF = classes$cls_col_list == "factor"
	if (sum(TF) > 0)
		{
		if (sum(TF) > 1)
			{
			for (i in 1:sum(TF))
				{
				trtable[,TF][,i] = as.character(trtable[,TF][,i])
				}
			} else {
			trtable[,TF] = as.character(trtable[,TF])
			}
		} # END if (sum(TF) > 0)
	cls.df(trtable)
	
	# Return the table
	return(trtable)
	} # END read_auto_results <- function(settings=NULL)


#######################################################
# Check for -1, in a list of gregexpr() matches
# Best used with:
# sapply(X=matches, FUN=isNegOne)
#######################################################
isNegOne <- function(matches)
	{
	val = attr(x=matches, which="match.length")
	if (length(val) == 1)
		{
		if (val == -1)
			{
			return(TRUE)
			} else {
			return(FALSE)
			}
		} else {
		return(FALSE)
		}
	}





# This will only extract 1 group between 1st startflag & 1st doneflag
extract_lines_startstr_to_endstr <- function(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	{
	startflag = 0
	doneflag = 0
	lines_to_keep = c()
	lines_to_keep_i = 0
	instance = 0
	for (i in 1:length(lines))
		{
		line = lines[i]
		
		#print(paste(string_to_start_at, line))
		if (grepl(string_to_start_at, line))
			{
			instance = instance+1
			if (instance == instance_to_find)
				{
				startflag = 1
				}
			}
		
		# Only check for the doneflag, *IF YOU HAVE ALREADY FOUND THE STARTFLAG*
		if (startflag == 1)
			{
			if (grepl(string_to_end_at, line))
				{
				doneflag = 1
				}
			}
		
		
		# end if done; if not, check startflag
		if (doneflag == 1)
			{
			blah = 0
			if(printflag)
				{			
				print("Found string_to_end_at, exiting extract_lines_startstr_to_endstr()...")
				}
			
			# Include the last line, if desired
			if (include_endstring == TRUE)
				{
				lines_to_keep_i = lines_to_keep_i + 1
				if(printflag)
					{			
					print(paste("Storing ", lines_to_keep_i, "th line, line #", i, sep=""))
					}
				lines_to_keep[[lines_to_keep_i]] = line
				
				return(lines_to_keep)
				}
			}
		else
			{
			# MAIN inclusion of lines
			if (startflag == 1)
				{
				lines_to_keep_i = lines_to_keep_i + 1
				#print(lines_to_keep_i)
				if(printflag)
					{			
					print(paste("Storing ", lines_to_keep_i, "th line, line #", i, sep=""))
					}
				lines_to_keep[[lines_to_keep_i]] = line				
				}
			else
				{
				blah = 0
				}
			}
		
		}
	
	# Check if you've got the correct instance
	
	return(lines_to_keep)
	}




#######################################################
# Read the saved list of synapomorphies from the logfile
#######################################################
extract_synapos <- function(logfn, trtable=NULL, string_to_start_at="START FULL LIST OF INFERRED SYNAPOMORPHIES", string_to_end_at="END FULL LIST OF INFERRED SYNAPOMORPHIES")
	{
	defaults='
	logfn = "auto_logfile.txt"
	trtable=NULL
	string_to_start_at = "START FULL LIST OF INFERRED SYNAPOMORPHIES"
	string_to_end_at = "END FULL LIST OF INFERRED SYNAPOMORPHIES"
	'
	
	# Read the lines
	lines = readLines(con=logfn)
	lines_to_keep = extract_lines_startstr_to_endstr(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	lines_to_keep = gdata::trim(lines_to_keep)
	lines_to_keep

	# Remove useless stuff
	TF1 = lines_to_keep == ""
	#TF2 = lines_to_keep == "All trees:"
	TF3 = lines_to_keep == "START FULL LIST OF INFERRED SYNAPOMORPHIES"
	TF4 = grepl(pattern="\\*\\*\\*\\*\\*", x=lines_to_keep)
	TF5a = grepl(pattern="apo ", x=lines_to_keep)
	TF5b = grepl(pattern=" secs.", x=lines_to_keep)
	TF5 = (TF5a + TF5b) == 2
	TF6 = grepl(pattern="Synapomorphies common to ", x=lines_to_keep)
	TF7 = grepl(pattern="Node numbers refer to nodes in consensus", x=lines_to_keep)
	TF = (TF1 + TF3 + TF4 + TF5 + TF6 + TF7) > 0
	lines_to_keep = lines_to_keep[TF == FALSE]
	lines_to_keep
	
	# Exit if none found
	if (is.null(lines_to_keep))
		{
		return(NULL)
		}
	
	# Identify types of information
	no_autapos_TF = lines_to_keep == "No autapomorphies:"
	no_synapos_TF = lines_to_keep == "No synapomorphies"
	char_TF = grepl(pattern="Char. ", x=lines_to_keep)
	node_TF = grepl(pattern="Node ", x=lines_to_keep)
	all_TF = grepl(pattern="All trees:", x=lines_to_keep)
	some_TF = grepl(pattern="Some trees:", x=lines_to_keep)
	tipnames_TF = (no_autapos_TF + no_synapos_TF + char_TF + node_TF + all_TF + some_TF) == 0

	# Remove " :"
	lines_to_keep[tipnames_TF] = gsub(pattern=" :", replacement="", x=lines_to_keep[tipnames_TF])
	tipnames = lines_to_keep[tipnames_TF]

	# Parse the list
	tmpmat = matrix(data="", nrow=length(lines_to_keep), ncol=9)
	nodecount = 0
	apocount = 0
	for (i in 1:length(lines_to_keep))
		{
		# If tip or node, start new row
		if (tipnames_TF[i] == TRUE)
			{
			type = "tip"
			nodecount = nodecount + 1
			node = lines_to_keep[i]
			}
		if (node_TF[i] == TRUE)
			{
			type = "node"
			nodecount = nodecount + 1
			node = gsub(pattern="Node ", replacement="", lines_to_keep[i])
			node = gsub(pattern=" :", replacement="", node)
			}
		if (some_TF[i] == TRUE)
			{
			someall = "some"
			}
		if (all_TF[i] == TRUE)
			{
			someall = "all"
			}
		if (char_TF[i] == TRUE)
			{
			# parse character
			words = strsplit(lines_to_keep[i], split=":")[[1]]
			words = trim(words)
			# get character num
			charnum = gsub(pattern="Char. ", replacement="", x=words[1])
			# get character change
			words3 = strsplit(words[2], split=" ")[[1]]
			ancstate = words3[1]
			decstate = words3[3]
			apocount = apocount + 1
			tmprow = c(i, nodecount, apocount, type, node, someall, charnum, ancstate, decstate)
			tmpmat[apocount,] = tmprow
			} # END if (char_TF[i] == TRUE)
		} # END for (i in 1:length(lines_to_keep))
		# Subset table
	synapos = as.data.frame(tmpmat[1:apocount,], stringsAsFactors=FALSE, row.names=NULL)
	names(synapos) = c("i", "nodecount", "apocount", "type", "node", "someall", "charnum", "ancstate", "decstate")
	
	# If the trtable is available, with TNT node numbers, add to that:
	if (is.null(trtable) == FALSE)
		{
		# Link to trtable / APE phylo object
		tipsTF = synapos$type == "tip"
		tipmatches = match(x=synapos$node[tipsTF], trtable$label)
		nodesTF = synapos$type == "node"
		nodematches = match(x=synapos$node[nodesTF], trtable$TNT_nodenums)

		APEnode = rep(NA, times=nrow(synapos))
		APEnode[tipsTF] = trtable$node[tipmatches]
		APEnode[nodesTF] = trtable$node[nodematches]
		
		synapos = cbind(synapos, APEnode)
		} # END if (is.null(trtable) == FALSE)
	
	return(synapos)
	}




# Consistency Index
calc_CI <- function(cscores, minsteps)
	{
	CIs = minsteps / cscores
	return(CIs)
	}

# Retention Index
calc_RI <- function(cscores, minsteps, maxsteps)
	{
	RIs = (maxsteps - cscores) / (maxsteps - minsteps)
	return(RIs)
	}

# Rescaled Consistency Index (RC or RCI)
calc_RCI <- function(cscores, minsteps, maxsteps)
	{
	CIs = calc_CI(cscores, minsteps)
	RIs = calc_RI(cscores, minsteps, maxsteps)
	RCIs = RIs * CIs
	return(RCIs)
	}

# Homoplasy Index: HI = 1 - CI
calc_HI <- function(cscores, minsteps)
	{
	CIs = calc_CI(cscores, minsteps)
	HIs = 1 - CIs
	return(HIs)
	}

autoget_charstats <- function(logfn)
	{
	defaults='
	logfn = "auto_logfile.txt"
	'
	cscores = extract_cscores(logfn)
	chomo = extract_chomo(logfn)
	maxsteps = extract_maxsteps(logfn)
	minsteps = extract_minsteps(logfn)

	cnum = 1:length(cscores)
	TNTcnum = cnum-1

	charstats = cbind(cnum, TNTcnum, minsteps, maxsteps, cscores, chomo)
	charstats = as.data.frame(charstats, stringsAsFactors=FALSE, row.names=NULL)
	head(charstats)
	tail(charstats)

	CIs = calc_CI(cscores, minsteps)
	RIs = calc_RI(cscores, minsteps, maxsteps)
	RCIs = calc_RCI(cscores, minsteps, maxsteps)

	charstats = cbind(charstats, CIs, RIs, RCIs)
	charstats
	
	return(charstats)
	}


autoget_treestats <- function(logfn)
	{
	defaults='
	logfn = "auto_logfile.txt"
	'
	res = extract_CI_RI(logfn)
	CI = res$CI
	RI = res$RI	
	TL = extract_TL(logfn)
	RCI = CI*RI
	
	treestats = c(TL, CI, RI, RCI)
	tmpmat = matrix(data=treestats, nrow=1, ncol=length(treestats), byrow=TRUE)
	tmpmat
	treestats = as.data.frame(tmpmat, stringsAsFactors=FALSE, row.names=NULL)
	treestats
	names(treestats) = c("TL", "CI", "RI", "RCI")
	
	return(treestats)
	}

extract_CI_RI <- function(logfn, string_to_start_at="Reading from stats.run ", string_to_end_at="Again reading from ")
	{
	defaults='
	logfn = "auto_logfile.txt"
	string_to_start_at = "Reading from stats.run"
	string_to_end_at = "Again reading from "
	'
	
	# Read the lines
	lines = readLines(con=logfn)
	lines_to_keep = extract_lines_startstr_to_endstr(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	lines_to_keep = gdata::trim(lines_to_keep)
	lines_to_keep
	
	# Remove blanks
	lines_to_keep = lines_to_keep[lines_to_keep != ""]
	
	# Subset lines
	TF = lines_to_keep == "Consistency index"
	linenum = (1:length(lines_to_keep))[TF]
	linenums = (linenum+1):length(lines_to_keep)
	lines_to_keep = lines_to_keep[linenums]
	
	TF = lines_to_keep == "Retention index"
	linenum = (1:length(lines_to_keep))[TF]
	CI_lines = 1:(linenum-1)
	RI_lines = (linenum+1):length(lines_to_keep)
	
	CItxt = lines_to_keep[CI_lines]
	RItxt = lines_to_keep[RI_lines]
	
	CIwords = trim(strsplit(paste0(CItxt, collapse=" "), split=" "))[[1]]
	CIwords = CIwords[CIwords != ""]
	CInums = as.numeric(CIwords)
	CI = CInums[CInums==max(CInums)][1]
	
	RIwords = trim(strsplit(paste0(RItxt, collapse=" "), split=" "))[[1]]
	RIwords = RIwords[RIwords != ""]
	RInums = as.numeric(RIwords)
	RI = RInums[RInums==max(RInums)][1]
	
	res = NULL
	res$CI = CI
	res$RI = RI
	
	extract='
	CI = res$CI
	RI = res$RI	
	'
	
	return(res)
	} # END extract_stats
	

# Extract total length
extract_TL <- function(logfn, string_to_start_at="Tree lengths", string_to_end_at="length ")
	{
	defaults='
	logfn = "auto_logfile.txt"
	string_to_start_at = "Tree lengths"
	string_to_end_at = "length "
	'
	
	# Read the lines
	lines = readLines(con=logfn)
	lines_to_keep = extract_lines_startstr_to_endstr(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	lines_to_keep = gdata::trim(lines_to_keep)
	lines_to_keep
	
	# Remove blanks
	lines_to_keep = lines_to_keep[lines_to_keep != ""]
	
	# Subset lines
	TF = grepl(pattern=string_to_start_at, x=lines_to_keep)
	linenum = (1:length(lines_to_keep))[TF]
	linenums = (linenum+1):length(lines_to_keep)
	lines_to_keep = lines_to_keep[linenums]
	TLtxt = lines_to_keep
	
	TLwords = trim(strsplit(paste0(TLtxt, collapse=" "), split=" "))[[1]]
	TLwords = TLwords[TLwords != ""]
	TLnums = as.numeric(TLwords)
	TL = TLnums[TLnums==max(TLnums)][1]
	TL
	
	return(TL)
	} # END extract_stats
	


# Extract cscores (character length for each character)
extract_cscores <- function(logfn, string_to_start_at="Tree 0, total", string_to_end_at="cscores ")
	{
	defaults='
	logfn = "auto_logfile.txt"
	string_to_start_at = "Tree 0, total"
	string_to_end_at = "cscores "
	'
	
	# Read the lines
	lines = readLines(con=logfn)
	lines_to_keep = extract_lines_startstr_to_endstr(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	lines_to_keep = gdata::trim(lines_to_keep)
	lines_to_keep
	
	# Remove blanks
	lines_to_keep = lines_to_keep[lines_to_keep != ""]
	
	# Subset lines
	TF = grepl(pattern=string_to_start_at, x=lines_to_keep)
	linenum = (1:length(lines_to_keep))[TF]
	linenums = (linenum+1):length(lines_to_keep)
	lines_to_keep = lines_to_keep[linenums]
	
	# Remove "0      1      2      3      4      5      6      7      8      9"
	TLtxt = lines_to_keep[2:length(lines_to_keep)]
	TLtxt
	
	cscores = NULL
	for (i in 1:length(TLtxt))
		{
		TLwords = trim(strsplit(paste0(TLtxt[i], collapse=" "), split=" "))[[1]]
		TLwords = TLwords[TLwords != ""]
		TLnums = as.numeric(TLwords[2:length(TLwords)])
		cscores = c(cscores, TLnums)
		} # END for (i in 1:length(TLtxt))
	cscores
	
	return(cscores)
	} # END extract_stats
	



# Extract chomo (homoplasy: extract steps beyond the minimum)
extract_chomo <- function(logfn, string_to_start_at="Tree 0, homoplasy", string_to_end_at="chomo ")
	{
	defaults='
	logfn = "auto_logfile.txt"
	string_to_start_at = "Tree 0, homoplasy"
	string_to_end_at = "chomo "
	'
	
	# Read the lines
	lines = readLines(con=logfn)
	lines_to_keep = extract_lines_startstr_to_endstr(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	lines_to_keep = gdata::trim(lines_to_keep)
	lines_to_keep
	
	# Remove blanks
	lines_to_keep = lines_to_keep[lines_to_keep != ""]
	
	# Subset lines
	TF = grepl(pattern=string_to_start_at, x=lines_to_keep)
	linenum = (1:length(lines_to_keep))[TF]
	linenums = (linenum+1):length(lines_to_keep)
	lines_to_keep = lines_to_keep[linenums]
	
	# Remove "0      1      2      3      4      5      6      7      8      9"
	TLtxt = lines_to_keep[2:length(lines_to_keep)]
	TLtxt
	
	chomo = NULL
	for (i in 1:length(TLtxt))
		{
		TLwords = trim(strsplit(paste0(TLtxt[i], collapse=" "), split=" "))[[1]]
		TLwords = TLwords[TLwords != ""]
		TLnums = as.numeric(TLwords[2:length(TLwords)])
		chomo = c(chomo, TLnums)
		} # END for (i in 1:length(TLtxt))
	chomo
	
	return(chomo)
	} # END extract_chomo
	



# Extract maxsteps (character length for each character)
extract_maxsteps <- function(logfn, string_to_start_at="Maximum possible steps ", string_to_end_at="\\*\\*\\*\\*\\*")
	{
	defaults='
	logfn = "auto_logfile.txt"
	string_to_start_at = "Maximum possible steps "
	string_to_end_at = "\\*\\*\\*\\*\\*"
	'
	
	# Read the lines
	lines = readLines(con=logfn)
	lines_to_keep = extract_lines_startstr_to_endstr(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	lines_to_keep = gdata::trim(lines_to_keep)
	lines_to_keep
	
	# Remove blanks
	lines_to_keep = lines_to_keep[lines_to_keep != ""]
	
	# Subset lines
	TF = grepl(pattern=string_to_start_at, x=lines_to_keep)
	linenum = (1:length(lines_to_keep))[TF]
	linenums = (linenum+1):length(lines_to_keep)
	lines_to_keep = lines_to_keep[linenums]
	
	# Remove "0      1      2      3      4      5      6      7      8      9"
	TLtxt = lines_to_keep[2:length(lines_to_keep)]
	TLtxt
	
	maxsteps = NULL
	for (i in 1:length(TLtxt))
		{
		TLwords = trim(strsplit(paste0(TLtxt[i], collapse=" "), split=" "))[[1]]
		TLwords = TLwords[TLwords != ""]
		TLnums = as.numeric(TLwords[2:length(TLwords)])
		maxsteps = c(maxsteps, TLnums)
		} # END for (i in 1:length(TLtxt))
	maxsteps
	
	return(maxsteps)
	} # END extract_stats
	



# Extract minsteps (character length for each character)
extract_minsteps <- function(logfn, string_to_start_at="Minimum possible steps ", string_to_end_at="Maximum possible steps ")
	{
	defaults='
	logfn = "auto_logfile.txt"
	string_to_start_at = "Minimum possible steps "
	string_to_end_at = "Maximum possible steps "
	'
	
	# Read the lines
	lines = readLines(con=logfn)
	lines_to_keep = extract_lines_startstr_to_endstr(lines, string_to_start_at, string_to_end_at, printflag=FALSE, include_endstring=FALSE, instance_to_find=1)
	lines_to_keep = gdata::trim(lines_to_keep)
	lines_to_keep
	
	# Remove blanks
	lines_to_keep = lines_to_keep[lines_to_keep != ""]
	
	# Subset lines
	TF = grepl(pattern=string_to_start_at, x=lines_to_keep)
	linenum = (1:length(lines_to_keep))[TF]
	linenums = (linenum+1):length(lines_to_keep)
	lines_to_keep = lines_to_keep[linenums]
	
	# Remove "0      1      2      3      4      5      6      7      8      9"
	TLtxt = lines_to_keep[2:length(lines_to_keep)]
	TLtxt
	
	minsteps = NULL
	for (i in 1:length(TLtxt))
		{
		TLwords = trim(strsplit(paste0(TLtxt[i], collapse=" "), split=" "))[[1]]
		TLwords = TLwords[TLwords != ""]
		TLnums = as.numeric(TLwords[2:length(TLwords)])
		minsteps = c(minsteps, TLnums)
		} # END for (i in 1:length(TLtxt))
	minsteps
	
	return(minsteps)
	} # END extract_stats
	











# Write MSM.run with specific input/output files
write_TNT_MSM_runfile <- function(fn_to_copy="MSM.run", msm_in_dates_fn, msm_run_fn, msm_out_tree_fn, msm_out_fn)
	{
	
	# Copy the original MSM.run file, change the text in the relevant lines for input/output files
	
	# Read a text file into a list of strings
	lines = scan(fn_to_copy, what="character", sep="\n")

	# Replace output date-calibrated parsimony trees file
	lines = gsub("datamsm.tnt", msm_in_dates_fn, lines)

	# Replace output date-calibrated parsimony trees file
	lines = gsub("MSM.run", msm_run_fn, lines)
	
	# Replace output date-calibrated parsimony trees file
	lines = gsub("calibrated_tree.tre", msm_out_tree_fn, lines)

	# Replace output date-calibrated parsimony trees file
	lines = gsub("MSM.OUT", msm_out_fn, lines)
	
	
	# Write the list of commands to a file
	write.table(lines, file=paste(msm_run_fn, sep=""), quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	return(msm_run_fn)
	}


# Concatenate a list of TNT tree files into one big file
cat_tree_files <- function(list_of_treefns, outfn)
	{
	new_lines = c()
	list_num_files = c()
	
	for (i in 1:length(list_of_treefns))
		{
		# input file
		infn = list_of_treefns[i]
		
		lines = scan(infn, what="character", sep="\n")
		
		list_num_files = c(list_num_files, length(lines))

		new_lines = c(new_lines, lines)
		}

	cat("# of lines in raw cat of tree file: ", length(new_lines), "\n", sep="")

	# save the first ("tread...") and last ("proc -") commands
	first_line_perm = new_lines[1]
	last_line_perm = new_lines[length(new_lines)]

	# remove semi-colons
	new_lines2 = replace_str_lines(new_lines, ");", ")*")

	# save the first ("tread...") and last ("proc -") commands
	first_line = new_lines2[1]
	last_line = new_lines2[length(new_lines2)]
	
	new_lines3 = null_lines(new_lines2, "tread ")
	new_lines4 = null_lines(new_lines3, last_line)

	cat("# of lines in cat of tree file, minus 'tread...' and 'proc -;' commands: ", length(new_lines4), "\n", sep="")

	
	last_tree = replace_str_lines(new_lines4[length(new_lines4)], ")\\*", ");")
	new_lines5 = c(first_line_perm, new_lines4[1:(length(new_lines4)-1)], last_tree, last_line_perm)
	cat("# of lines in cat of tree file, with start/end lines: ", length(new_lines5), "\n", sep="")

	
	
	write.table(new_lines5, file=outfn, quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	return(list_num_files)
	}


extract_lines_startnum_to_flag <- function(lines, row_to_start, string_to_stop_at)
	{
	doneflag = 0
	lines_to_keep = c()
	lines_to_keep_i = 0
	for (i in 1:length(lines))
		{
		line = lines[i]
	
		# skip the 1st two lines
		if (i < row_to_start)
			{
			blah = 0
			}
		else
			{		
			if (doneflag == 1)
				{
				blah = 0
				}
			else
				{
				# if line starts with a dash, you are done...
				if (grepl(string_to_stop_at, line))
					{
					#print("yo!")
					doneflag = 1
					}
				else
					{
					lines_to_keep_i = lines_to_keep_i + 1
					print(paste("Storing line #", lines_to_keep_i, sep=""))
					lines_to_keep[[lines_to_keep_i]] = line
					}
				}
			}
		
		}
	return(lines_to_keep)
	}



# Parse an MSM.out file that results from MSM.run
parse_MSM_out <- function(fn, outfn, row_to_start=4, string_to_stop_at = "----")
	{

	# Read a text file into a list of strings
	lines = scan(fn, what="character", sep="\n")
	
	# This seems to change a little depending on log file etc...
	#row_to_start=4
	#string_to_stop_at = "----"
	lines_to_keep = extract_lines_startnum_to_flag(lines, row_to_start, string_to_stop_at)
	print(length(lines_to_keep))
	
	# parse the summary stat lines from msm.out:
	newlines = c()
	for (i in 1:length(lines_to_keep))
		{
		line = lines_to_keep[i]
		# split on whitespace
		list_of_strs = strsplit_whitespace(line)
		
		MSM = as.numeric(list_of_strs[2])
		GER = as.numeric(list_of_strs[3])
		LmLo = list_of_strs[7]
		LmLo_str = strsplit(LmLo, "/")			# L = length in terms of years etc. gapped
		Lm = as.numeric(LmLo_str[[1]][1])		# Lm = Lmin = minimum possible tree length
		Lo = as.numeric(LmLo_str[[1]][2])		# Lo = Lobserved = observed length of MP tree
		
		finalcalcs = list_of_strs[12]
		finalcalcs1 = strsplit(finalcalcs, "-")
		finalcalcs2 = finalcalcs1[[1]][1]
		finalcalcs3 = sub("[[]", "", finalcalcs2)		# [[] = regexp for [
		LM = as.numeric(finalcalcs3)			# LM = LMax = maximum possible tree length
		
		treenum = i-1
		
		newline = c(treenum, MSM, GER, Lm, Lo, LM)
		newlines = rbind(newlines, newline)
		}
	
	newlines = data.frame(newlines, row.names = NULL)
	names(newlines) = c("Tree", "MSMs", "GER", "Lm", "Lo", "LM")
	
	(newlines)
	
	write.table(newlines, file=outfn, append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="	")

	return(newlines)
	}


write_TNT_MSM_input_dates_data <- function(fn=NULL, msm_in_dates_fn, taxon_years=NULL, taxon_names=NULL, present_year=2010)
	{	
	# SET UP THE MSM ANALYSIS
	
	# need gdata for read.xls
	require(gdata)
	
	#dir = "/Users/nick/Documents/NCSE/_issues/_issues_history/Darwin_quote/__v3/_3_tntall/"
	#dir = "/Users/nick/Documents/NCSE/_issues/_issues_history/Darwin_quote/"
	#setwd(dir)
	
	# =========================================================
	# Read in the table with the row names (which have to contain the dates, in this script)
	# =========================================================# version 24 removed all "#" signs from links and elsewhere, these apparently scrambled read.table
	#fn = "2010-08-04_Darwin_quote_db_v34.txt"
	#fn = "input_chars.txt"
	#fn = rownames_fn
	#fnxls = "2010-08-04_Darwin_quote_db_v34.xls"
	
	#dirfn = paste(dir, fn, sep="")
	#dirfnxls = paste(dir, fnxls, sep="")
	
	
	# ===================
	# First run the standard Bremer or Bootstrap script to make a list of MP trees.
	# Then, the following...
	# ===================
	
	# ==========================================================================================
	# Comparing the most parsimonious trees using the MSM* statistic, calculated in TNT
	# 
	# Make a list of the dates to write a "dataMSM.tnt" file that contains the 
	#   age character (an irreversible Sankoff character, I think, where the 
	#   number of steps = distance in my between dates
	# ==========================================================================================
	
	
	#agefn = "dataMSM.tnt"
	agefn = msm_in_dates_fn
	
	# a list/matrix of the strings to output to the text file
	s = c()
	
	# set the maximum amount of RAM, and number of states in character
	s[1] = "mxram 50; nstates 32;"
	s[2] = "xread"
	s[3] = "'Age character matrix for the calculation of the modified MSM or the GER'"
	

	# Get the list of ages (in unit years before present)
	present_year = present_year
	
	if ((is.null(taxon_names) || is.null(taxon_years)) && is.null(fn))
		{
		txt = "STOP ERROR in write_TNT_MSM_input_dates_data(): Either 'fn' or 'taxon_years' and 'taxon_names' must be specified"
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		} # END if (is.null(taxon_years) && is.null(fn))
	
	
	if (is.null(taxon_years) == TRUE)
		{
		#fn = "dq_quote.txt"
		z = read_table_good(fn)

		# number of columns (1 = 1 character for age) and rows (# of taxa)
		s[4] = paste("1 ", nrow(z), "", sep="")

		taxon_names = as.character(z$rownames)
		#for (i in 1:length(taxon_names))
		#	{
		#	taxon_names[i] = substr(taxon_names[i], 2, nchar(taxon_names[i], type="chars"))
		#	}
		taxon_ages = rep(NA, nrow(z))

		taxon_years = rep(NA, nrow(z))
		for (i in 1:nrow(z))
			{
			taxon_years[i] = as.numeric(substr(taxon_names[i], 2, 5))
			} # END for (i in 1:nrow(z))
		} else {
		# number of columns (1 = 1 character for age) and rows (# of taxa)
		s[4] = paste("1 ", length(taxon_names), "", sep="")
		}
		
		
			
	taxon_ages = present_year - taxon_years
	
	# print the uniq taxon ages
	list_of_ages = unique(taxon_ages)
	num_of_ages = length(list_of_ages)
	cat("List of ", num_of_ages, " ages: (note: 32 distinct ages, max!!) \n")
	cat(list_of_ages, "\n")
	if (num_of_ages > 32)
		{
		cat("ERROR!!!!  You have ", num_of_ages, " distinct ages, but TNT and the MSM script only allow 32 character states in a character!\n", sep="")
		}
	
	# Assign character states to the ages
	# 32 characters in TNT = 0 1 2 3 4 5 6 7 8 9 A B C D E F G H I J K L M N O P Q R S T U V
	age_charstates = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V")
	
	
	for (i in 1:length(taxon_ages))
		{
		# add the header length to the output line in the output tnt ages file
		outfn_linenum = i + 4
		
		# Get the tmpage
		tmpage = taxon_ages[i]
		
		# Get the matching code by matching list_of_ages to age_charstates
		tmp_age_charstate = NA
		for (j in 1:length(list_of_ages))
			{
			#cat(tmpage, list_of_ages[j],  age_charstates[j], "\n", sep="	")
			if (as.character(tmpage) == as.character(list_of_ages[j]))
				{
				cat(tmpage, list_of_ages[j],  age_charstates[j], "\n", sep="	")
				tmp_age_charstate = age_charstates[j]
				}
			}
		
		s[outfn_linenum] = paste(taxon_names[i], "	", tmp_age_charstate, "", sep="")
		}
	
		
	s[outfn_linenum + 1] = ";"
	s[outfn_linenum + 2] = ""
	s[outfn_linenum + 3] = "p/;"
	s[outfn_linenum + 4] = ""
	s[outfn_linenum + 5] = ""
	s[outfn_linenum + 6] = "label AGES;"
	s[outfn_linenum + 7] = ""
	s[outfn_linenum + 8] = "/* DEFINE THE AGES OF EACH CHARACTER STATE IN THE FOLLOWING BLOCK"
	s[outfn_linenum + 9] = "	The age should be in mya or any other time unit you want to use."
	s[outfn_linenum + 10] = "*/"
	s[outfn_linenum + 11] = ""
	s[outfn_linenum + 12] = ""
	outfn_linenum = outfn_linenum + 12
	
	for (i in 1:length(list_of_ages))
		{
		outfn_linenum = outfn_linenum + 1
		s[outfn_linenum] = paste("set age[", i-1, "] ", list_of_ages[i], ";", sep="")
	
		}
	
	s[outfn_linenum + 1] = ""
	s[outfn_linenum + 2] = ""
	s[outfn_linenum + 3] = ""
	s[outfn_linenum + 4] = ""
	s[outfn_linenum + 5] = "proc/;"
	s[outfn_linenum + 6] = ""
	s[outfn_linenum + 7] = ""
	
	dim(s) = c(length(s), 1)
	
	# Fix any NAs
	#s[is.na(s)] = ""
	
	write.table(s, file=agefn, quote=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	return(agefn)
	}







# Find the most parsimonious trees, calculate summary statistics
write_TNT_find_MP_trees_cmds <- function(input_chars, weights_cmds, weights_suffixes, rownames_fn, tnt_scripts_fns, Bremer_sup=FALSE, bootstraps=FALSE)
	{
	
	# Cycle through all of the weights
	for (j in 1:length(weights_cmds))
		{
		# Set up weights for this run
		weights_cmd = weights_cmds[j]
		suffix = weights_suffixes[j]
		
		# The base fn for all scripts and results
		fnbase = gsub(".tnt", suffix, input_chars)
		
		# setup output filenames (fns)
		tnt_script_fn = paste(fnbase, "_script.run", sep="")
		tnt_log_fn = paste(fnbase, ".txt", sep="")
		best_trees_tnt_fn = paste(fnbase, "_best_trees.tnttrees", sep="")
		best_trees_nex_fn = paste(fnbase, "_best_trees.nex", sep="")
		best_tree_Oth_brlen_tnt_fn = paste(fnbase, "_best_tree_Oth_brlen.tnttree", sep="")
		best_tree_Oth_brlen_nex_fn = paste(fnbase, "_best_tree_Oth_brlen.nex", sep="")
		
		nelsen_tree_tnt_fn = paste(fnbase, "_nelsen.tnttree", sep="")
		nelsen_tree_nex_fn = paste(fnbase, "_nelsen.nex", sep="")
		
		majority_tree_tnt_fn = paste(fnbase, "_majority.tnttree", sep="")
		majority_tree_nex_fn = paste(fnbase, "_majority.nex", sep="")
		
		comcomp_tree_tnt_fn = paste(fnbase, "_comcomp.tnttree", sep="")
		comcomp_tree_nex_fn = paste(fnbase, "_comcomp.nex", sep="")
		comcomp_tree_w_Bremers_tnt_fn = paste(fnbase, "_comcomp_w_Bremers.tnttree", sep="")
		comcomp_tree_w_Bremers_nex_fn = paste(fnbase, "_comcomp_w_Bremers.nex", sep="")
		
		comcomp_tree_w_bootstraps_tnt_fn = paste(fnbase, "_comcomp_w_bootstraps.tnttree", sep="")
		comcomp_tree_w_bootstraps_nex_fn = paste(fnbase, "_comcomp_w_bootstraps.nex", sep="")
		
		msm_run_fn = paste(fnbase, "_MSM.run", sep="")
		msm_in_dates_fn = paste(fnbase, "_dataTNT.tnt", sep="")
		msm_out_tree_fn = paste(fnbase, "_calibrated_tree.tre", sep="")
		msm_out_fn = paste(fnbase, "_MSM.out", sep="")
		
		
		# comcomp_tree_tnt_fn
		
		MP_treenum_to_chose = 0
		
		# list of lines to print to tnt_script_fn
		l = c()
		# count up the line numbers, starting with 0
		h = 0
		
		# script header
		l[(h=h+1)] = 
		"macro=;		/* Turns on macros for scripting. */
		/*###############################									*/
		/*# Search for most parsimonious trees								*/
		/*###############################									*/
		/*																	*/"
		
		l[(h=h+1)] = paste("/*	cd ", analysisdir, "					*/", sep="")
		
		l[(h=h+1)] = "/*	./tnt.command													*/"
		l[(h=h+1)] = paste("/*	proc ", tnt_script_fn, ";												*/", sep="")
		l[(h=h+1)] = "/*	ipy tnt2newick.py												*/"
		l[(h=h+1)] = "/* Start the log file */"
		l[(h=h+1)] = paste("log ", tnt_log_fn, ";", sep="")
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		# setup to process the data
		l[(h=h+1)] = "/* Turn off ('-') issuing of warnings */"
		l[(h=h+1)] = "warn -;"
		l[(h=h+1)] = "/* Turn on ('=') user-break with 'Esc' or '.', or pause with 'p' */"
		l[(h=h+1)] = "break =;"
		
		
		l[(h=h+1)] = "/* Read in the data */"
		if (weights_cmd != "piwe -;")
			{
			l[(h=h+1)] = "piwe =;"
			}
		l[(h=h+1)] = paste("proc ", input_chars, ";", sep="")
		l[(h=h+1)] = ""
		
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote ==========================================================================;"
		l[(h=h+1)] = paste("quote SCRIPT STARTING ANALYSIS OF: ", input_chars, ";", sep="")
		l[(h=h+1)] = paste("quote DESCRIPTION: ", dataset_description, ".;", sep="")
		l[(h=h+1)] = "quote ==========================================================================;"
		l[(h=h+1)] = ""
		
		
		
		
		l[(h=h+1)] = "/* make tables a simpler parsing format (=); to turn off: - */"
		l[(h=h+1)] = "tables =;"
		l[(h=h+1)] = "/* set number of significant digits to 4 */"
		l[(h=h+1)] = "tables /4;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		l[(h=h+1)] = "quote BEGIN TAXON LABELS;"
		l[(h=h+1)] = "taxlabels;"
		l[(h=h+1)] = "quote END TAXON LABELS;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		l[(h=h+1)] = "quote BEGIN CHARACTER LABELS;"
		l[(h=h+1)] = "cnames;"
		l[(h=h+1)] = "quote END CHARACTER LABELS;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		l[(h=h+1)] = "quote BEGIN INFORMATIVE CHARACTERS;"
		l[(h=h+1)] = "info +;"
		l[(h=h+1)] = "quote END INFORMATIVE CHARACTERS;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote BEGIN UNINFORMATIVE CHARACTERS;"
		l[(h=h+1)] = "info -;"
		l[(h=h+1)] = "quote END UNINFORMATIVE CHARACTERS;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Inactivate characters that are non-parsimony-informative */"
		l[(h=h+1)] = "xinact;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		# Print the outgroup
		l[(h=h+1)] = "/* print outgroup -- should be 1st taxon (taxon 0) */"
		l[(h=h+1)] = "quote BEGIN OUTGROUP;"
		l[(h=h+1)] = "outgroup;"
		l[(h=h+1)] = "quote END OUTGROUP;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		
		l[(h=h+1)] = "/* weights; default is k=3 (decay of relative weighting with */"
		l[(h=h+1)] = "/* increasing number of homoplasy steps, weight<0.01 below 46 steps */"
		l[(h=h+1)] = "/* If 'piwe -;' -- Here, no weights (equal weighting assumption) */"
		
		# only display the relative weights if they are on;
		#  otherwise, display 
		if (weights_cmd == "piwe -;")
			{
			l[(h=h+1)] = weights_cmd
			l[(h=h+1)] = "quote BEGIN DISPLAY OF RELATIVE WEIGHTS BY AMOUNT OF HOMOPLASY;"
			l[(h=h+1)] = "piwe ;"
			l[(h=h+1)] = "quote END WEIGHTS DISPLAY;"
			l[(h=h+1)] = "quote  ;"
			l[(h=h+1)] = "quote  ;"	
			} else {
			l[(h=h+1)] = "piwe =;"
			l[(h=h+1)] = weights_cmd
			l[(h=h+1)] = "quote BEGIN DISPLAY OF RELATIVE WEIGHTS BY AMOUNT OF HOMOPLASY;"
			l[(h=h+1)] = "piwe ;"
			l[(h=h+1)] = "piwe [;"
			l[(h=h+1)] = "quote END WEIGHTS DISPLAY;"
			l[(h=h+1)] = "quote  ;"
			l[(h=h+1)] = "quote  ;"
			}
		
		
		# Do the search
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Run the parsimony search */"
		l[(h=h+1)] = "/* (following online manual recommendations for low numbers of taxa) */"
		l[(h=h+1)] = "/* [NOTE : Stats.run seems to crash if you have more than about 400 trees */"
		l[(h=h+1)] = "/* in memory. i.e. after 'hold 400; mult 100 =tbr drift ;' works, but after */"
		l[(h=h+1)] = "/* hold 500;' it doesn't.]  */"
		l[(h=h+1)] = "hold 400 ;"
		l[(h=h+1)] = "mult 100 =tbr drift ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		# Summary stats
		l[(h=h+1)] = "/* Get CI and RI (retention index) */"
		l[(h=h+1)] = "quote BEGIN SUMMARY STATS FOR MP TREES;"
		l[(h=h+1)] = "stats;"
		l[(h=h+1)] = "quote END SUMMARY STATS FOR MP TREES;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		
		# Tree length
		l[(h=h+1)] = "/* Get tree length */"
		l[(h=h+1)] = "quote BEGIN LENGTH FOR MP TREES;"
		l[(h=h+1)] = "length;"
		l[(h=h+1)] = "quote END LENGTH FOR MP TREES;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		
		# Homoplasy per-character
		l[(h=h+1)] = "quote BEGIN MINMAX FOR MP TREES: reporting minimum & maximum number of steps per-character;"
		l[(h=h+1)] = "minmax;"
		l[(h=h+1)] = "quote END MINMAX FOR MP TREES;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote BEGIN CHOMO FOR MP TREES: reporting homoplasy per-character;"
		l[(h=h+1)] = "chomo;"
		l[(h=h+1)] = "quote END CHOMO FOR MP TREES: reporting homoplasy per-character;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		
		
		# Save the best trees (MP trees = most parsimonious trees)
		l[(h=h+1)] = "/* keep the names in memory: taxname= ; */"
		l[(h=h+1)] = "taxname= ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* These trees all have the same length */"
		l[(h=h+1)] = "/* save trees in TNT format */ "
		l[(h=h+1)] = "/* show tags doesn't work here: ttags ; */"
		l[(h=h+1)] = paste("tsave *", best_trees_tnt_fn, ";", sep="")
		l[(h=h+1)] = paste("export - ", best_trees_nex_fn, ";", sep="")
		# save 1st (0th) tree
		#l[(h=h+1)] = "save 0;"
		l[(h=h+1)] = "save ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Close file receiving saved trees. */"
		l[(h=h+1)] = "tsave/;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* see the trees */"
		l[(h=h+1)] = "tplot ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		# plot one of the trees with branch lengths (# of steps)
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* plot the first of MP trees as a phylogram, WITH branchlengths */"
		l[(h=h+1)] = paste("tchoose ", MP_treenum_to_chose, ";", sep="")
		l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
		l[(h=h+1)] = "ttag- ;"
		l[(h=h+1)] = "/* store the tree tags */"
		l[(h=h+1)] = "ttag= ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote BEGIN PRINTING OF BRANCH LENGTHS OF 0th TREE BY NODE;"
		l[(h=h+1)] = "blength 0;"
		l[(h=h+1)] = "quote END PRINTING OF BRANCH LENGTHS OF 0th TREE BY NODE;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* produce the branch lengths of the 0th tree; hopefully stored as tags */"
		l[(h=h+1)] = "blength *0;"
		l[(h=h+1)] = "/* Save tags in a readable form. */"
		l[(h=h+1)] = "ttags /;"
		l[(h=h+1)] = "/* Show the tree tags (prints to screen I think). */"
		l[(h=h+1)] = "quote BEGIN PRINTING OF BRANCH LENGTHS OF 0th TREE;"
		l[(h=h+1)] = "ttags ;"
		l[(h=h+1)] = "quote END PRINTING OF BRANCH LENGTHS OF 0th TREE;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Save the consensus tree, with tags holding Bremer supports */"
		l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
		l[(h=h+1)] = paste("tsave *", best_tree_Oth_brlen_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* Save the tree(s) with tags ('*') */"
		l[(h=h+1)] = "save *;"
		l[(h=h+1)] = "/* Close file receiving saved trees. */"
		l[(h=h+1)] = "tsave/;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote BEGIN PLOTTING PHYLOGRAM OF 0th TREE WITH BRANCH LENGTHS;"
		l[(h=h+1)] = paste("export -", best_tree_Oth_brlen_nex_fn, ";", sep="")
		l[(h=h+1)] = "quote END PLOTTING PHYLOGRAM OF 0th TREE WITH BRANCH LENGTHS;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		l[(h=h+1)] = "tplot ;"
		l[(h=h+1)] = "quote END PLOTTING PHYLOGRAM WITH BRANCH LENGTHS;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		
		l[(h=h+1)] = "quote BEGIN CALCULATION OF CONSENSUS TREES;"
		
		# Calculate strict (Nelsen) consensus
		l[(h=h+1)] = "/* strict (Nelsen) consensus:								*/"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		l[(h=h+1)] = paste("proc ", best_trees_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
		l[(h=h+1)] = "ttag- ;"
		l[(h=h+1)] = "/* store the tree tags */"
		l[(h=h+1)] = "ttag= ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Calculate strict (Nelsen) consensus tree */"
		l[(h=h+1)] = "nelsen *;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* This saves the last tree (the consensus tree)						*/"
		l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
		l[(h=h+1)] = paste("tsave *", nelsen_tree_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* Save the tree(s) (last tree, '/') with tags ('*') */"
		l[(h=h+1)] = "save /;"
		l[(h=h+1)] = "/* Close file receiving saved trees. */"
		l[(h=h+1)] = "tsave/;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		
		
		# Calculate majority rule consensus
		l[(h=h+1)] = "/* majority rule consensus:								*/"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		l[(h=h+1)] = paste("proc ", best_trees_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
		l[(h=h+1)] = "ttag- ;"
		l[(h=h+1)] = "/* store the tree tags */"
		l[(h=h+1)] = "ttag= ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Calculate majority-rule consensus tree */"
		l[(h=h+1)] = "majority *;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* This saves the last tree (the consensus tree)						*/"
		l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
		l[(h=h+1)] = paste("tsave *", majority_tree_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* Save the tree(s) (last tree, '/') with tags ('*') */"
		l[(h=h+1)] = "save /;"
		l[(h=h+1)] = "/* Close file receiving saved trees. */"
		l[(h=h+1)] = "tsave/;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		# Calculate Bremer consensus (combinable components)
		l[(h=h+1)] = "/* Combinable components (Bremer consensus):								*/
		/* COMCOMP 				 												*/
		/* Calculate combinable component (=Bremer ) consensus tree				*/
		/*    N/L    display consensus for tree(s) N, excluding taxon (taxa) L 	*/
		/*    *N/L   same, but keep consensus as last tree in memory			*/
		/* 	
		*/"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		l[(h=h+1)] = paste("proc ", best_trees_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
		l[(h=h+1)] = "ttag- ;"
		l[(h=h+1)] = "/* store the tree tags */"
		l[(h=h+1)] = "ttag= ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Calculate combinable consensus tree */"
		l[(h=h+1)] = "comcomp *;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* This saves the last tree (the consensus tree)						*/"
		l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
		l[(h=h+1)] = paste("tsave *", comcomp_tree_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* Save the tree(s) (last tree, '/') with tags ('*') */"
		l[(h=h+1)] = "save /;"
		l[(h=h+1)] = "/* Close file receiving saved trees. */"
		l[(h=h+1)] = "tsave/;"
		l[(h=h+1)] = ""
		
		l[(h=h+1)] = "quote END CALCULATION OF CONSENSUS TREES;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = ""
		
		
		
		#Bremer_sup = FALSE
		if (Bremer_sup)
			{
			# Bremer support
			l[(h=h+1)] = "/* Bremer support 													*/"
			l[(h=h+1)] = "/* choose the last tree in memory; presumably the consensus tree 	*/"
			l[(h=h+1)] = "tchoose /;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
			l[(h=h+1)] = "ttag- ;"
			l[(h=h+1)] = "/* store the tree tags */"
			l[(h=h+1)] = "ttags =;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Run Kip Will's modified Bremer support script */"
			l[(h=h+1)] = "/* Note: this only works well if the consensus tree is close to the MP tree; */"
			l[(h=h+1)] = "/* which it usually isn't, except in perfect situations. */"
			l[(h=h+1)] = "KWBremer ;"
			
			# Turn off warnings again, since KWBremer turns them on
			l[(h=h+1)] = "/* Turn off ('-') issuing of warnings */"
			l[(h=h+1)] = "warn -;"
			
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Save tags in a readable form. */"
			l[(h=h+1)] = "ttags /;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Show the tree tags (prints to screen I think). */"
			l[(h=h+1)] = "ttags ;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Save the consensus tree, with tags holding Bremer supports */"
			l[(h=h+1)] = paste("tsave * ", comcomp_tree_w_Bremers_tnt_fn, ";", sep="")
			l[(h=h+1)] = "save *;"
			l[(h=h+1)] = "/* Close file receiving saved trees. */"
			l[(h=h+1)] = "tsave/;"
			l[(h=h+1)] = paste("export - ", comcomp_tree_w_Bremers_nex_fn, ";", sep="")
			l[(h=h+1)] = ""
			l[(h=h+1)] = ""
			l[(h=h+1)] = ""
			
			# Summary stats on just the consensus tree
			l[(h=h+1)] = "/* Select the last tree (the consensus tree), and discard the rest */"
			l[(h=h+1)] = "tchoose /;"
			l[(h=h+1)] = "quote BEGIN SUMMARY STATS ON CONSENSUS TREE;"
			l[(h=h+1)] = "/* Get CI and RI (retention index) */"
			l[(h=h+1)] = "stats;"
			l[(h=h+1)] = "quote END SUMMARY STATS ON CONSENSUS TREE;"
			l[(h=h+1)] = ""
			
			# Tree length
			l[(h=h+1)] = "/* Get tree length */"
			l[(h=h+1)] = "quote BEGIN LENGTH FOR CONSENSUS TREE;"
			l[(h=h+1)] = "length;"
			l[(h=h+1)] = "quote END LENGTH FOR CONSENSUS TREE;"
			l[(h=h+1)] = ""
			
			l[(h=h+1)] = "quote BEGIN MINMAX FOR CONSENSUS TREE: reporting minimum & maximum number of steps per-character;"
			l[(h=h+1)] = "minmax;"
			l[(h=h+1)] = "quote END MINMAX FOR CONSENSUS TREE;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "quote BEGIN CHOMO: reporting homoplasy per-character;"
			l[(h=h+1)] = "chomo;"
			l[(h=h+1)] = "quote END CHOMO: reporting homoplasy per-character;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = ""
			}
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		############################################
		# MSM stuff
		############################################
		
		sourcedir = '/Dropbox/_njm/'
		source4 = 'tnt_R_utils_v1.R'
		source(paste(sourcedir, source4, sep=""))
		
		
		# Write the years to dates...
		###################################################################
		# Write the dates of the taxa to a file (max 32 taxa)
		# for analysis by MSM.run
		###################################################################
		dir = "/Users/nick/Documents/NCSE/_issues/_issues_history/Darwin_quote/__v3/_3_tntall/"
		#dir = "/Users/nick/Documents/NCSE/_issues/_issues_history/Darwin_quote/"
		setwd(dir)
		
		# =========================================================
		# Read in the table with the row names (which have to contain the dates, in this script)
		# =========================================================# version 24 removed all "#" signs from links and elsewhere, these apparently scrambled read.table
		#fn = "2010-08-04_Darwin_quote_db_v34.txt"
		#fn = "input_chars.txt"
		fn = rownames_fn
		fnxls = "2010-08-04_Darwin_quote_db_v34.xls"
		
		dirfn = paste(dir, fn, sep="")
		dirfnxls = paste(dir, fnxls, sep="")
		
		write_TNT_MSM_input_dates_data(fn, msm_in_dates_fn)
		
		
		# Write a new version of the the MSM.run script, changing the filenames to this specific run...
		write_TNT_MSM_runfile(fn_to_copy="MSM.run", msm_in_dates_fn, msm_run_fn, msm_out_tree_fn, msm_out_fn)
		
		
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote BEGIN RUNNING MSM.RUN;"
		l[(h=h+1)] = "/* Commands to run the TNT script to calculate MSM, GER on */"
		l[(h=h+1)] = "/* the various most-parsimonious trees, and pick the best. */"
		l[(h=h+1)] = paste("/* cd ", dir, "      */", sep="")
		l[(h=h+1)] = "/* ./tnt.command        */"
		l[(h=h+1)] = "/* (not used inside bigscript...) proc runMSM_v1.tnt   */"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* */"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		
		
		l[(h=h+1)] = "/* Load the data */"
		
		#l[(h=h+1)] = "proc dataMSM.tnt;"
		l[(h=h+1)] = paste("proc ", msm_in_dates_fn, ";", sep="")
		
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Load the trees */"
		l[(h=h+1)] = paste("proc ", best_trees_tnt_fn, ";", sep="")
		
		#l[(h=h+1)] = "run MSM.run;"
		l[(h=h+1)] = paste("run ", msm_run_fn, ";", sep="")
		l[(h=h+1)] = "quote END RUNNING MSM.RUN;"
		
		
		#bootstraps = FALSE
		if(bootstraps)
			{
			##############################################
			# Calculate bootstraps (kinda slow)
			##############################################
			l[(h=h+1)] = ""
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Reload the raw data */"
			l[(h=h+1)] = paste("proc ", input_chars, ";", sep="")
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Reload the MP trees */"
			l[(h=h+1)] = paste("proc ", comcomp_tree_tnt_fn, ";", sep="")
			l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
			l[(h=h+1)] = "ttag- ;"
			l[(h=h+1)] = "/* store the tree tags */"
			l[(h=h+1)] = "ttag= ;"
			l[(h=h+1)] = ""
			
			l[(h=h+1)] = "/* choose the last tree in memory; presumably the consensus tree */"
			l[(h=h+1)] = "tchoose /;"
			
			l[(h=h+1)] = "ttags =;"
			
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Turn on reporting for time-consuming operations*/"
			l[(h=h+1)] = "report =;"
			
			l[(h=h+1)] = "resample boot replications 100 from 0;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Save tags in a readable form. */"
			l[(h=h+1)] = "ttags /;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Show the tree tags (prints to screen I think). */"
			l[(h=h+1)] = "ttags ;"
			l[(h=h+1)] = ""
			l[(h=h+1)] = "/* Save the consensus tree, with tags holding bootstrap supports */"
			l[(h=h+1)] = paste("tsave * ", comcomp_tree_w_bootstraps_tnt_fn, ";", sep="")
			l[(h=h+1)] = "save *;"
			l[(h=h+1)] = "/* Close file receiving saved trees. */"
			l[(h=h+1)] = "tsave/;"
			l[(h=h+1)] = paste("export - ", comcomp_tree_w_bootstraps_nex_fn, ";", sep="")
			l[(h=h+1)] = ""
			l[(h=h+1)] = ""
			l[(h=h+1)] = ""
			}
		
		
		
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote ==========================================================================;"
		l[(h=h+1)] = paste("quote SCRIPT ENDING ANALYSIS OF: ", input_chars, ";", sep="")
		l[(h=h+1)] = paste("quote DESCRIPTION: ", dataset_description, ".;", sep="")
		l[(h=h+1)] = "quote ==========================================================================;"
		l[(h=h+1)] = ""
		
		
		
		
		
		l[(h=h+1)] = "/* End the log */"
		l[(h=h+1)] = "log/;"
		l[(h=h+1)] = ""
		
		l[(h=h+1)] = "/* Exit... */"
		#l[(h=h+1)] = "procedure/;"
		l[(h=h+1)] = "zzz;"
		l[(h=h+1)] = ""
		
		l[(h=h+1)] = ""
		
		
		# Write the list of commands to a file
		write.table(l, file=paste(analysisdir, tnt_script_fn, sep=""), quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
		
		# Add to the list of script fns
		tnt_scripts_fns = c(tnt_scripts_fns, tnt_script_fn)

		}	

	return(tnt_scripts_fns)
	}





# Basic TNT run:
basic_TNT_run <- function(dataTNT_fn, analysisdir, ordered_chars=NULL, force_cmd=NULL, dataset_description="description here", Bremer_sup=TRUE, bootstraps=FALSE, runstats=FALSE, weights_cmd=NULL, mxram=50, nstates=6)
	{
	suffix = ""
	
	# The base fn for all scripts and results
	fnbase = gsub(".tnt", suffix, dataTNT_fn)
	
	# setup output filenames (fns)
	tnt_script_fn = paste(fnbase, "_script.run", sep="")
	tnt_log_fn = paste(fnbase, "_tntlog.txt", sep="")
	best_trees_tnt_fn = paste(fnbase, "_best_trees.tnttrees", sep="")
	best_trees_nex_fn = paste(fnbase, "_best_trees.nex", sep="")
	best_tree_Oth_brlen_tnt_fn = paste(fnbase, "_best_tree_Oth_brlen.tnttree", sep="")
	best_tree_Oth_brlen_nex_fn = paste(fnbase, "_best_tree_Oth_brlen.nex", sep="")
	
	each_tree_brlen_tnt_fn = paste(fnbase, "_each_tree_brlen", sep="")
	
	
	nelsen_tree_tnt_fn = paste(fnbase, "_nelsen.tnttree", sep="")
	nelsen_tree_nex_fn = paste(fnbase, "_nelsen.nex", sep="")
	
	majority_tree_tnt_fn = paste(fnbase, "_majority.tnttree", sep="")
	majority_tree_nex_fn = paste(fnbase, "_majority.nex", sep="")
	
	comcomp_tree_tnt_fn = paste(fnbase, "_comcomp.tnttree", sep="")
	comcomp_tree_nex_fn = paste(fnbase, "_comcomp.nex", sep="")
	comcomp_tree_w_Bremers_tnt_fn = paste(fnbase, "_comcomp_w_Bremers.tnttree", sep="")
	comcomp_tree_w_Bremers_nex_fn = paste(fnbase, "_comcomp_w_Bremers.nex", sep="")
	
	comcomp_tree_w_bootstraps_tnt_fn = paste(fnbase, "_comcomp_w_bootstraps.tnttree", sep="")
	comcomp_tree_w_bootstraps_nex_fn = paste(fnbase, "_comcomp_w_bootstraps.nex", sep="")
	
	
	MP_treenum_to_chose = 0
	
	# list of lines to print to tnt_script_fn
	l = c()
	# count up the line numbers, starting with 0
	h = 0
	
	# script header
	l[(h=h+1)] = 
	"macro=;		/* Turns on macros for scripting. */
	/*###############################									*/
	/*# Search for most parsimonious trees								*/
	/*###############################									*/
	/*																	*/"
	
	l[(h=h+1)] = paste("/*	cd ", analysisdir, "					*/", sep="")
	
	l[(h=h+1)] = "/*	./tnt.command													*/"
	l[(h=h+1)] = paste("/*	proc ", tnt_script_fn, ";												*/", sep="")
	l[(h=h+1)] = "/*	ipy tnt2newick.py												*/"
	l[(h=h+1)] = "/* Start the log file */"
	l[(h=h+1)] = paste("log ", tnt_log_fn, ";", sep="")
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	# setup to process the data
	l[(h=h+1)] = "/* Turn off ('-') issuing of warnings */"
	l[(h=h+1)] = "warn -;"
	l[(h=h+1)] = "/* Turn on ('=') user-break with 'Esc' or '.', or pause with 'p' */"
	l[(h=h+1)] = "break =;"

	l[(h=h+1)] = "/*   */"
	l[(h=h+1)] = "/* Set the amount of RAM */"
	l[(h=h+1)] = paste("mxram ", mxram, "; nstates ", nstates, ";", sep="")
	l[(h=h+1)] = "/*   */"
	
	
	l[(h=h+1)] = "/* Read in the data; default is equal weighting (piwe -;) */"
	if (is.null(weights_cmd))
		{
		# 
		l[(h=h+1)] = "piwe -;"
		weights_cmd = "piwe -;"
		}
	l[(h=h+1)] = paste("proc ", dataTNT_fn, ";", sep="")
	l[(h=h+1)] = ""
	
	l[(h=h+1)] = ""
	l[(h=h+1)] = "quote ==========================================================================;"
	l[(h=h+1)] = paste("quote SCRIPT STARTING ANALYSIS OF: ", dataTNT_fn, ";", sep="")
	l[(h=h+1)] = paste("quote DESCRIPTION: ", dataset_description, ".;", sep="")
	l[(h=h+1)] = "quote ==========================================================================;"
	l[(h=h+1)] = ""
	

	# ORDERED characters
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Character ordering options; user-input if desired */"
	if (is.null(ordered_chars) == FALSE)
		{
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Input the ordered characters */"
		
		# input 20 characters at a time to avoid line length issues
		starts = seq(1, length(ordered_chars), 20)
		ends = c(-1+starts[2:length(starts)], length(ordered_chars))
		starts
		ends
		
		for (i in 1:length(starts))
			{
			chars_to_input_on_this_line = paste(ordered_chars[starts[i]:ends[i]], collapse=" ")
			l[(h=h+1)] = paste("ccode +", chars_to_input_on_this_line, ";", sep="")
			}
		}
	l[(h=h+1)] = "quote ;"
	l[(h=h+1)] = "quote BEGIN DISPLAY OF CHARACTER ORDERING;"
	l[(h=h+1)] = "ccode ;"
	l[(h=h+1)] = "quote END DISPLAY OF CHARACTER ORDERING;"
	l[(h=h+1)] = "quote ;"
	l[(h=h+1)] = "quote ;"
	
		
	
	l[(h=h+1)] = "/* make tables a simpler parsing format (=); to turn off: - */"
	l[(h=h+1)] = "tables =;"
	l[(h=h+1)] = "/* set number of significant digits to 4 */"
	l[(h=h+1)] = "tables /4;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	
	l[(h=h+1)] = "quote BEGIN TAXON LABELS;"
	l[(h=h+1)] = "taxlabels;"
	l[(h=h+1)] = "quote END TAXON LABELS;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	l[(h=h+1)] = "quote BEGIN CHARACTER LABELS;"
	l[(h=h+1)] = "cnames;"
	l[(h=h+1)] = "quote END CHARACTER LABELS;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	
	l[(h=h+1)] = "quote BEGIN INFORMATIVE CHARACTERS;"
	l[(h=h+1)] = "info +;"
	l[(h=h+1)] = "quote END INFORMATIVE CHARACTERS;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "quote BEGIN UNINFORMATIVE CHARACTERS;"
	l[(h=h+1)] = "info -;"
	l[(h=h+1)] = "quote END UNINFORMATIVE CHARACTERS;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Inactivate characters that are non-parsimony-informative */"
	l[(h=h+1)] = "xinact;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	
	# Set the ingroup (and thus the outgroup by default)
	if (is.null(force_constraints) == FALSE)
		{
		
		
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Specify ingroup */"
		l[(h=h+1)] = "/* Force 1 monophyletic group, the ingroup, and set the others to outgroup */"
		l[(h=h+1)] = force_cmd
		l[(h=h+1)] = "/* Set the outgroup (default is 0 anyway)*/"
		l[(h=h+1)] = "outgroup 0;"
		l[(h=h+1)] = "/* Make sure constraints are on */"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote CHECK CONSTRAINTS"
		l[(h=h+1)] = "constrain =;"
		l[(h=h+1)] = "/* Constraints is ON */"
		l[(h=h+1)] = "/* Double-check the constraints & that they are on. */"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote FORCE CONSTRAINTS ON;"
		l[(h=h+1)] = "force ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "force *;"
		l[(h=h+1)] = ""
		}
	
	# Print the outgroup
	l[(h=h+1)] = "/* print outgroup -- should be 1st taxon (taxon 0) */"
	l[(h=h+1)] = "quote BEGIN OUTGROUP;"
	l[(h=h+1)] = "outgroup;"
	l[(h=h+1)] = "quote END OUTGROUP;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	
	
	l[(h=h+1)] = "/* weights; default is k=3 (decay of relative weighting with */"
	l[(h=h+1)] = "/* increasing number of homoplasy steps, weight<0.01 below 46 steps */"
	l[(h=h+1)] = "/* If 'piwe -;' -- Here, no weights (equal weighting assumption) */"
	
	# only display the relative weights if they are on;
	#  otherwise, display 
	if (weights_cmd == "piwe -;")
		{
		l[(h=h+1)] = weights_cmd
		l[(h=h+1)] = "quote BEGIN DISPLAY OF RELATIVE WEIGHTS BY AMOUNT OF HOMOPLASY;"
		l[(h=h+1)] = "piwe ;"
		l[(h=h+1)] = "quote END WEIGHTS DISPLAY;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"	
		} else {
		l[(h=h+1)] = "piwe =;"
		l[(h=h+1)] = weights_cmd
		l[(h=h+1)] = "quote BEGIN DISPLAY OF RELATIVE WEIGHTS BY AMOUNT OF HOMOPLASY;"
		l[(h=h+1)] = "piwe ;"
		l[(h=h+1)] = "piwe [;"
		l[(h=h+1)] = "quote END WEIGHTS DISPLAY;"
		l[(h=h+1)] = "quote  ;"
		l[(h=h+1)] = "quote  ;"
		}
	
	

	# Do the search
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Run the parsimony search */"
	l[(h=h+1)] = "/* (following online manual recommendations for low numbers of taxa) */"
	l[(h=h+1)] = "/* [NOTE : Stats.run seems to crash if you have more than about 400 trees */"
	l[(h=h+1)] = "/* in memory. i.e. after 'hold 400; mult 100 =tbr drift ;' works, but after */"
	l[(h=h+1)] = "/* hold 500;' it doesn't.]  */"
	l[(h=h+1)] = "hold 400 ;"
	l[(h=h+1)] = "mult 100 =tbr drift ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	# Summary stats
	l[(h=h+1)] = "/* Get CI and RI (retention index) */"
	l[(h=h+1)] = "quote BEGIN SUMMARY STATS FOR MP TREES;"
	l[(h=h+1)] = "stats;"
	l[(h=h+1)] = "quote END SUMMARY STATS FOR MP TREES;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	
	# Tree length
	l[(h=h+1)] = "/* Get tree length */"
	l[(h=h+1)] = "quote BEGIN LENGTH FOR MP TREES;"
	l[(h=h+1)] = "length;"
	l[(h=h+1)] = "quote END LENGTH FOR MP TREES;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	
	# Homoplasy per-character
	l[(h=h+1)] = "quote BEGIN MINMAX FOR MP TREES: reporting minimum & maximum number of steps per-character;"
	l[(h=h+1)] = "minmax;"
	l[(h=h+1)] = "quote END MINMAX FOR MP TREES;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "quote BEGIN CHOMO FOR MP TREES: reporting homoplasy per-character;"
	l[(h=h+1)] = "chomo;"
	l[(h=h+1)] = "quote END CHOMO FOR MP TREES: reporting homoplasy per-character;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	
	
	# Save the best trees (MP trees = most parsimonious trees)
	l[(h=h+1)] = "/* keep the names in memory: taxname= ; */"
	l[(h=h+1)] = "taxname= ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* These trees all have the same length */"
	l[(h=h+1)] = "/* save trees in TNT format */ "
	l[(h=h+1)] = "/* show tags doesn't work here: ttags ; */"
	l[(h=h+1)] = paste("tsave *", best_trees_tnt_fn, ";", sep="")
	l[(h=h+1)] = paste("export - ", best_trees_nex_fn, ";", sep="")
	# save 1st (0th) tree
	#l[(h=h+1)] = "save 0;"
	l[(h=h+1)] = "save ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Close file receiving saved trees. */"
	l[(h=h+1)] = "tsave/;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* see the trees */"
	l[(h=h+1)] = "tplot ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""





	# Save each tree individually with branch lengths
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote ======================================================== ;"
	l[(h=h+1)] = "quote BEGIN PLOT AND SAVE ALL TREES WITH BRANCH LENGTHS;"
	l[(h=h+1)] = "quote ======================================================== ;"
	
	l[(h=h+1)] = ""
	l[(h=h+1)] = "var ="
	l[(h=h+1)] = "   + fns [ ( ntrees + 1 ) ]"
	l[(h=h+1)] = "   + bl"
	l[(h=h+1)] = ";"
	l[(h=h+1)] = "loop 0 ntrees"
	l[(h=h+1)] = "	/* progress command requires */"
	l[(h=h+1)] = "	/* report- ;                 */"
	l[(h=h+1)] = "	/* progress #1 (ntrees+1) Calculating branch lengths... ;*/"
	l[(h=h+1)] = "	/* set bl blength[#1] ;*/"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "	/* blank out the old tree tags to avoid extra commas */"
	l[(h=h+1)] = "	ttag- ;"
	l[(h=h+1)] = "	/* store the tree tags */"
	l[(h=h+1)] = "	ttag= ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "	quote BEGIN PRINTING OF BRANCH LENGTHS OF [#1]th TREE BY NODE;"
	l[(h=h+1)] = "	blength #1;"
	l[(h=h+1)] = "	quote END PRINTING OF BRANCH LENGTHS OF [#1]th TREE BY NODE;"
	l[(h=h+1)] = "	"
	l[(h=h+1)] = "	/* print branch lengths to table */"
	l[(h=h+1)] = "	/* produce the branch lengths of the 0th tree; hopefully stored as tags */"
	l[(h=h+1)] = "	blength *#1;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "	/* Save tags in a readable form. */"
	l[(h=h+1)] = "	ttags /;"
	l[(h=h+1)] = "	/* Show the tree tags (prints to screen I think). */"
	l[(h=h+1)] = "	quote BEGIN PRINTING OF BRANCH LENGTHS OF [#1]th TREE;"
	l[(h=h+1)] = "	ttags ;"
	l[(h=h+1)] = "	quote END PRINTING OF BRANCH LENGTHS OF [#1]th TREE;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = paste("	tsave *", 	each_tree_brlen_tnt_fn, "#1", "tr.tnttree;", sep="")
	l[(h=h+1)] = "	/* Save the tree(s) with tags ('*') */"
	l[(h=h+1)] = "	save *;"
	l[(h=h+1)] = "	/* Close file receiving saved trees. */"
	l[(h=h+1)] = "	tsave/;"
	l[(h=h+1)] = "stop"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote ======================================================== ;"
	l[(h=h+1)] = "quote END PLOT AND SAVE ALL TREES WITH BRANCH LENGTHS;"
	l[(h=h+1)] = "quote ======================================================== ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""

	
	
	



	
	# plot one of the trees with branch lengths (# of steps)
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* plot the first of MP trees as a phylogram, WITH branchlengths */"
	l[(h=h+1)] = paste("tchoose ", MP_treenum_to_chose, ";", sep="")
	l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
	l[(h=h+1)] = "ttag- ;"
	l[(h=h+1)] = "/* store the tree tags */"
	l[(h=h+1)] = "ttag= ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "quote BEGIN PRINTING OF BRANCH LENGTHS OF 0th TREE BY NODE;"
	l[(h=h+1)] = "blength 0;"
	l[(h=h+1)] = "quote END PRINTING OF BRANCH LENGTHS OF 0th TREE BY NODE;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* produce the branch lengths of the 0th tree; hopefully stored as tags */"
	l[(h=h+1)] = "blength *0;"
	l[(h=h+1)] = "/* Save tags in a readable form. */"
	l[(h=h+1)] = "ttags /;"
	l[(h=h+1)] = "/* Show the tree tags (prints to screen I think). */"
	l[(h=h+1)] = "quote BEGIN PRINTING OF BRANCH LENGTHS OF 0th TREE;"
	l[(h=h+1)] = "ttags ;"
	l[(h=h+1)] = "quote END PRINTING OF BRANCH LENGTHS OF 0th TREE;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Save the consensus tree, with tags holding Bremer supports */"
	l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
	l[(h=h+1)] = paste("tsave *", best_tree_Oth_brlen_tnt_fn, ";", sep="")
	l[(h=h+1)] = "/* Save the tree(s) with tags ('*') */"
	l[(h=h+1)] = "save *;"
	l[(h=h+1)] = "/* Close file receiving saved trees. */"
	l[(h=h+1)] = "tsave/;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "quote BEGIN PLOTTING PHYLOGRAM OF 0th TREE WITH BRANCH LENGTHS;"
	l[(h=h+1)] = paste("export -", best_tree_Oth_brlen_nex_fn, ";", sep="")
	l[(h=h+1)] = "quote END PLOTTING PHYLOGRAM OF 0th TREE WITH BRANCH LENGTHS;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	
	l[(h=h+1)] = "tplot ;"
	l[(h=h+1)] = "quote END PLOTTING PHYLOGRAM WITH BRANCH LENGTHS;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	
	
	

	
	
	
	
	
	
	l[(h=h+1)] = "quote BEGIN CALCULATION OF CONSENSUS TREES;"
	
	# Calculate strict (Nelsen) consensus
	l[(h=h+1)] = "/* strict (Nelsen) consensus:								*/"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	l[(h=h+1)] = paste("proc ", best_trees_tnt_fn, ";", sep="")
	l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
	l[(h=h+1)] = "ttag- ;"
	l[(h=h+1)] = "/* store the tree tags */"
	l[(h=h+1)] = "ttag= ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Calculate strict (Nelsen) consensus tree */"
	l[(h=h+1)] = "nelsen *;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* This saves the last tree (the consensus tree)						*/"
	l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
	l[(h=h+1)] = paste("tsave *", nelsen_tree_tnt_fn, ";", sep="")
	l[(h=h+1)] = "/* Save the tree(s) (last tree, '/') with tags ('*') */"
	l[(h=h+1)] = "save /;"
	l[(h=h+1)] = "/* Close file receiving saved trees. */"
	l[(h=h+1)] = "tsave/;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	
	
	
	# Calculate majority rule consensus
	l[(h=h+1)] = "/* majority rule consensus:								*/"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	l[(h=h+1)] = paste("proc ", best_trees_tnt_fn, ";", sep="")
	l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
	l[(h=h+1)] = "ttag- ;"
	l[(h=h+1)] = "/* store the tree tags */"
	l[(h=h+1)] = "ttag= ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Calculate majority-rule consensus tree */"
	l[(h=h+1)] = "majority *;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* This saves the last tree (the consensus tree)						*/"
	l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
	l[(h=h+1)] = paste("tsave *", majority_tree_tnt_fn, ";", sep="")
	l[(h=h+1)] = "/* Save the tree(s) (last tree, '/') with tags ('*') */"
	l[(h=h+1)] = "save /;"
	l[(h=h+1)] = "/* Close file receiving saved trees. */"
	l[(h=h+1)] = "tsave/;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	
	# Calculate Bremer consensus (combinable components)
	l[(h=h+1)] = "/* Combinable components (Bremer consensus):								*/
	/* COMCOMP 				 												*/
	/* Calculate combinable component (=Bremer ) consensus tree				*/
	/*    N/L    display consensus for tree(s) N, excluding taxon (taxa) L 	*/
	/*    *N/L   same, but keep consensus as last tree in memory			*/
	/* 	
	*/"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	l[(h=h+1)] = paste("proc ", best_trees_tnt_fn, ";", sep="")
	l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
	l[(h=h+1)] = "ttag- ;"
	l[(h=h+1)] = "/* store the tree tags */"
	l[(h=h+1)] = "ttag= ;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* Calculate combinable consensus tree */"
	l[(h=h+1)] = "comcomp *;"
	l[(h=h+1)] = ""
	l[(h=h+1)] = "/* This saves the last tree (the consensus tree)						*/"
	l[(h=h+1)] = "/* Open a tree file in parenthetical format ('*') */"
	l[(h=h+1)] = paste("tsave *", comcomp_tree_tnt_fn, ";", sep="")
	l[(h=h+1)] = "/* Save the tree(s) (last tree, '/') with tags ('*') */"
	l[(h=h+1)] = "save /;"
	l[(h=h+1)] = "/* Close file receiving saved trees. */"
	l[(h=h+1)] = "tsave/;"
	l[(h=h+1)] = ""
	
	l[(h=h+1)] = "quote END CALCULATION OF CONSENSUS TREES;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = "quote  ;"
	l[(h=h+1)] = ""
	
	
	
	#Bremer_sup = FALSE
	if (Bremer_sup)
		{
		# Bremer support
		l[(h=h+1)] = "/* Bremer support 													*/"
		l[(h=h+1)] = "/* choose the last tree in memory; presumably the consensus tree 	*/"
		l[(h=h+1)] = "tchoose /;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
		l[(h=h+1)] = "ttag- ;"
		l[(h=h+1)] = "/* store the tree tags */"
		l[(h=h+1)] = "ttags =;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Run Kip Will's modified Bremer support script */"
		l[(h=h+1)] = "/* Note: this only works well if the consensus tree is close to the MP tree; */"
		l[(h=h+1)] = "/* which it usually isn't, except in perfect situations. */"
		l[(h=h+1)] = "KWBremer ;"
		
		# Turn off warnings again, since KWBremer turns them on
		l[(h=h+1)] = "/* Turn off ('-') issuing of warnings */"
		l[(h=h+1)] = "warn -;"
		
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Save tags in a readable form. */"
		l[(h=h+1)] = "ttags /;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Show the tree tags (prints to screen I think). */"
		l[(h=h+1)] = "ttags ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Save the consensus tree, with tags holding Bremer supports */"
		l[(h=h+1)] = paste("tsave * ", comcomp_tree_w_Bremers_tnt_fn, ";", sep="")
		l[(h=h+1)] = "save *;"
		l[(h=h+1)] = "/* Close file receiving saved trees. */"
		l[(h=h+1)] = "tsave/;"
		l[(h=h+1)] = paste("export - ", comcomp_tree_w_Bremers_nex_fn, ";", sep="")
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		
		# Summary stats on just the consensus tree
		l[(h=h+1)] = "/* Select the last tree (the consensus tree), and discard the rest */"
		l[(h=h+1)] = "tchoose /;"
		l[(h=h+1)] = "quote BEGIN SUMMARY STATS ON CONSENSUS TREE;"
		l[(h=h+1)] = "/* Get CI and RI (retention index) */"
		if (runstats)
			{
			l[(h=h+1)] = "stats;"
			} else {
			l[(h=h+1)] = "/* Note: User chose to not run stats.run; perhaps a memory problem, or a poor consensus tree? */"
			}
		l[(h=h+1)] = "quote END SUMMARY STATS ON CONSENSUS TREE;"
		l[(h=h+1)] = ""
		
		# Tree length
		l[(h=h+1)] = "/* Get tree length */"
		l[(h=h+1)] = "quote BEGIN LENGTH FOR CONSENSUS TREE;"
		l[(h=h+1)] = "length;"
		l[(h=h+1)] = "quote END LENGTH FOR CONSENSUS TREE;"
		l[(h=h+1)] = ""
		
		l[(h=h+1)] = "quote BEGIN MINMAX FOR CONSENSUS TREE: reporting minimum & maximum number of steps per-character;"
		l[(h=h+1)] = "minmax;"
		l[(h=h+1)] = "quote END MINMAX FOR CONSENSUS TREE;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "quote BEGIN CHOMO: reporting homoplasy per-character;"
		l[(h=h+1)] = "chomo;"
		l[(h=h+1)] = "quote END CHOMO: reporting homoplasy per-character;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		}
	l[(h=h+1)] = "/* */"
	l[(h=h+1)] = "/* */"
	l[(h=h+1)] = "/* */"
	l[(h=h+1)] = "/* */"
	l[(h=h+1)] = "/* */"
	l[(h=h+1)] = "/* */"
	l[(h=h+1)] = "/* */"
	l[(h=h+1)] = ""
	l[(h=h+1)] = ""
	
	

	#bootstraps = FALSE
	if(bootstraps)
		{
		##############################################
		# Calculate bootstraps (kinda slow)
		##############################################
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Reload the raw data */"
		l[(h=h+1)] = paste("proc ", dataTNT_fn, ";", sep="")
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Reload the MP trees */"
		l[(h=h+1)] = paste("proc ", comcomp_tree_tnt_fn, ";", sep="")
		l[(h=h+1)] = "/* blank out the old tree tags to avoid extra commas */"
		l[(h=h+1)] = "ttag- ;"
		l[(h=h+1)] = "/* store the tree tags */"
		l[(h=h+1)] = "ttag= ;"
		l[(h=h+1)] = ""
		
		l[(h=h+1)] = "/* choose the last tree in memory; presumably the consensus tree */"
		l[(h=h+1)] = "tchoose /;"
		
		l[(h=h+1)] = "ttags =;"
		
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Turn on reporting for time-consuming operations*/"
		l[(h=h+1)] = "report =;"
		
		l[(h=h+1)] = "resample boot replications 100 from 0;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Save tags in a readable form. */"
		l[(h=h+1)] = "ttags /;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Show the tree tags (prints to screen I think). */"
		l[(h=h+1)] = "ttags ;"
		l[(h=h+1)] = ""
		l[(h=h+1)] = "/* Save the consensus tree, with tags holding bootstrap supports */"
		l[(h=h+1)] = paste("tsave * ", comcomp_tree_w_bootstraps_tnt_fn, ";", sep="")
		l[(h=h+1)] = "save *;"
		l[(h=h+1)] = "/* Close file receiving saved trees. */"
		l[(h=h+1)] = "tsave/;"
		l[(h=h+1)] = paste("export - ", comcomp_tree_w_bootstraps_nex_fn, ";", sep="")
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		l[(h=h+1)] = ""
		}
	
	
	
	l[(h=h+1)] = ""
	l[(h=h+1)] = "quote ==========================================================================;"
	l[(h=h+1)] = paste("quote SCRIPT ENDING ANALYSIS OF: ", dataTNT_fn, ";", sep="")
	l[(h=h+1)] = paste("quote DESCRIPTION: ", dataset_description, ".;", sep="")
	l[(h=h+1)] = "quote ==========================================================================;"
	l[(h=h+1)] = ""
	
	
	
	
	
	l[(h=h+1)] = "/* End the log */"
	l[(h=h+1)] = "log/;"
	l[(h=h+1)] = ""
	
	l[(h=h+1)] = "/* Exit... */"
	#l[(h=h+1)] = "procedure/;"
	l[(h=h+1)] = "zzz;"
	l[(h=h+1)] = ""
	
	l[(h=h+1)] = ""
	
	
	# Write the list of commands to a file
	write.table(l, file=paste(analysisdir, tnt_script_fn, sep=""), quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	return(tnt_script_fn)
	}


