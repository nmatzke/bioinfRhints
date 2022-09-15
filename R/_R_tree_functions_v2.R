# =============================================
# _R_tree_functions_v1.R: many useful utility functions
#   for dealing with phylogenetic trees in R
#
# by Nick Matzke
# Copyright 2011-infinity
# matzkeATberkeley.edu
# January 2011
#
# Please link/cite if you use this, email me if you have 
#   thoughts/improvements/corrections.
#
##############################################################
#
# Free to use/redistribute under:
# Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0) 
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the above license, linked here:
# 
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
# Summary:
#
# You are free:
#
#   * to Share -- to copy, distribute and transmit the work
#   * to Remix -- to adapt the work
#
# Under the following conditions:
#
#   * Attribution -- You must attribute the work in the manner 
#     specified by the author or licensor (but not in any way that 
#     suggests that they endorse you or your use of the work).
#   * Noncommercial -- You may not use this work for commercial purposes. 
#
#   * Share Alike -- If you alter, transform, or build upon this work,
#     you may distribute the resulting work only under the same or
#     similar license to this one. 
#
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
###################################################################

# R functions for dealing with trees
# sourcedir = '/Dropbox/_njm/'
# source3 = '_R_tree_functions_v1.R'
# source(paste(sourcedir, source3, sep=""))


# Check the unique states in each character
checkchars <- function(tmpdf2)
	{
	unique_chars_list = rep(NA, nrow(tempdf2))
	for (i in 1:nrow(tempdf2))
		{
		chars = sort(unique(c(unlist(tempdf2[i,]))))
		
		l = length(chars)
		chars = paste(chars, sep="", collapse=",")
		
		unique_chars_list[i] = chars
		
		charstxt = paste(i, ") length=", l, ", chars=", chars, sep="")
		cat(charstxt, "\n", sep="")
		}
	return(unique_chars_list)
	}

postorder_traversal_down <- function(tr2, curnode, list_visited_internal_nodes)
	{
	anc_branch = tr2$edge[tr2$edge[,2] == curnode,]
	anc_node = anc_branch[1]
	
	list_visited_internal_nodes = c(list_visited_internal_nodes, anc_node)
	
	postorder_traveral_up(tr2, curnode, list_visited_internal_nodes)
	
	}

postorder_traveral_up <- function(tr2, curnode, list_visited_internal_nodes)
	{
	print("BLAH!")
	}


# Add a root branch to a standard APE phylo object/tree
add_root_to_tree <- function(tr, root_brlen=100)
	{
	trstr = write.tree(tr, file="")
	new_trstr = sub(pattern=";", replacement=paste("ingroup:", root_brlen, ";", sep=""), x=trstr)
	newtr = read.tree(file="", text=new_trstr, root.edge=TRUE)
	return(newtr)
	}


# Add a root branch to a standard BEAST phylo object/tree
add_root_to_nexus_trees <- function(trfn, out_trfn, root_brlen=100)
	{
	# Find the tree strings
	X <- scan(file = trfn, what = "", sep = "\n", quiet = TRUE)
  
    # Get line numbers for END; and ENDBLOCK;
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    
    # Get line numbers with semicolons
    semico <- grep(";", X)
        
    # Line # for beginning of TREES block
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# WITHOUT TRANSLATE block
	# Start & end lines of the trees block
	start <- i1+1
	end <- endblock[endblock > i1][1] - 1
	treelines <- X[start:end]
	
	for (i in 1:length(treelines))
		{
		trstr = treelines[i]
		
		# Add the root branchlength to each tree
		#new_trstr = sub(pattern=";", replacement=paste("ingroup:", root_brlen, ";", sep=""), x=trstr)
		new_trstr = sub(pattern=";", replacement=paste(":", root_brlen, ";", sep=""), x=trstr)
		
		# Insert back into lines
		X[(start-1) + i] = new_trstr
		}
	
	# Write to outfn
	write.table(X, file=out_trfn, append=FALSE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)
	
	return(out_trfn)
	}


# Try to get axisPhylo to plot all the way to 0 mya, when tips are
# non-contemporaneous
axisPhylo2 <- function(side = 1, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (lastPP$type %in% c("phylogram", "cladogram")) {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
			# Original behavior
			#x <- pretty(lastPP$xx)
			# Better behavior (all the way to 0)
            x <- pretty(c(lastPP$xx, 0))
            
            if (lastPP$direction == "rightwards") maxi <- max(lastPP$xx)
            else {
	            maxi <- min(lastPP$xx)
                x <- -x
            }
        } else {
			# Original behavior
			#x <- pretty(lastPP$yy)
			# Better behavior (all the way to 0)
            x <- pretty(c(lastPP$yy, 0))

            if (lastPP$direction == "upwards") maxi <- max(lastPP$yy)
            else {
                maxi <- min(lastPP$yy)
                x <- -x
            }
        }
    }
    axis(side = side, at = c(maxi - x), labels = abs(x), ...)
}



fixtree <- function(tree)
	{
	# Write a tree to Newick, then read back in
	# (hopefully this fixes whatever formatting issue was a problem)
	tmpfn = "tmp_junktree.tree"
	#write.tree(tree, tmpfn)
	#newtree2 = read.tree(tmpfn)
	
	newtree2 = write.tree(tree)
	
	return(newtree2)
	}

get_newick_str <- function(tree)
	{
	write.tree(tree, "temp.newick")
	str1 = read.table("temp.newick", skip=2)
	str2 = str1[1,1]
	return(str2)
	}


chardtf_w_taxa_as_cols_to_chardtf_w_taxa_as_rows <- function(char_dtf)
	{
	new_dtf = t(char_dtf)
	return(new_dtf)
	}


# Convert from a dataframe, with the taxa as columns, 
# back to a list (i.e., what you get from read.nexus.data
# of taxa, with the list of characters for each taxon
dtf_to_nexus_list <- function(char_dtf_w_taxa_as_cols)
	{
	tmp_nexus_list = NULL
	
	# This assumes that you've put the names into the column headers
	taxon_names = names(char_dtf_w_taxa_as_cols)
	
	for (colnum in 1:ncol(char_dtf_w_taxa_as_cols))
		{
		tmp_nexus_list[[taxon_names[colnum]]] = char_dtf_w_taxa_as_cols[, colnum]
		}
	
	return(tmp_nexus_list)
	}



# This assumes the taxa are rows, the characters are columns (?)
# (as.data.frame(read.nexus.data(filename)) produces the opposite, BTW)
dtf_to_phyDat <- function(char_dtf, datatype="USER", sites_not_site_patterns=TRUE)
	{
	require(phangorn)	# for phyDat
	
	# convert columns to character
	d = apply(char_dtf, 2, as.character)
	d2 = as.matrix(d)
	tmplevels = as.character(sort(unique(unlist(char_dtf))))
	phyd_tmp_site_patterns = phangorn::phyDat(d2, type="USER", levels=tmplevels)
	
	if (sites_not_site_patterns == FALSE)
		{
		phyd = phyd_tmp_site_patterns
		} else {
		# Convert to sites, not site patterns 
		phyd_characters = as.character(phyd_tmp_site_patterns)
	
		phyd = list()
		for (i in 1:nrow(phyd_characters))
			{
			phyd[[i]] = phyd_characters[i,]
			}
		attr(phyd, "class") = "phyDat"
		} # END if (sites_not_site_patterns == FALSE)
		
	attr(phyd, "names") = row.names(char_dtf)
	attr(phyd, "type") = datatype
	
	return(phyd)
	}

phyDat_to_str_dtf <- function(phyDat_seqs)
	{
	#seqs2 = as.character(phyDat_seqs)
	seqs2 = sapply(phyDat_seqs, as.character)
	
	(seqs3 = apply(seqs2, 1, paste, collapse=""))
	
	seqs4 = adf(seqs3)
	names(seqs4) = c("data")
	row.names(seqs4) = names(phyDat_seqs)
	
	return(seqs4)
	}

phyDat_to_seqlist <- function(phyDat_seqs, rows_or_cols="rows", subtract_1=TRUE)
	{
	# Get the length (number of taxa)
	#numseqs = length(phyDat_seqs)
	numseqs = nrow(phyDat_seqs)
	if (is.null(numseqs) == TRUE)
		{
		numseqs = length(phyDat_seqs)
		}
	
	# Number of characters
	numchars = nchar(phyDat_seqs[1])
	
	# convert to character
	# seqs2 = as.character(phyDat_seqs)
	as_numeric_then_char <- function(x)
		{
		as.character(as.numeric(x))
		}
	seqs2 = sapply(phyDat_seqs, as_numeric_then_char)

	tmp_rownames = row.names(phyDat_seqs)
	
	# for each sequence, make it a character list
	seqs3 = vector("list", length=numseqs)
	
	if (subtract_1 == TRUE)
		{
		if (rows_or_cols == "rows")
			{
			for (i in 1:numseqs)
				{
				seqs3[[i]] = as.numeric(seqs2[i, ]) - 1
				}
			} # END if (rows_or_cols == "cols")
		if (rows_or_cols == "cols")
			{
			for (i in 1:numseqs)
				{
				seqs3[[i]] = as.numeric(seqs2[, i]) - 1
				}
			} # END if (rows_or_cols == "cols")
		} # END if (subtract_1 == TRUE)
	if (subtract_1 == FALSE)
		{
		if (rows_or_cols == "rows")
			{
			for (i in 1:numseqs)
				{
				seqs3[[i]] = as.numeric(seqs2[i, ])
				}
			} # END if (rows_or_cols == "cols")
		if (rows_or_cols == "cols")
			{
			for (i in 1:numseqs)
				{
				seqs3[[i]] = as.numeric(seqs2[, i])
				}
			} # END if (rows_or_cols == "cols")
		} # END if (subtract_1 == FALSE)
	
	# Add the names
	attr(seqs3, "names") = tmp_rownames
	
	return(seqs3)
	}


phyDat_to_nexus_good <- function(phyDat_seqs, outfn="tmpnex.nex", printflag=FALSE, rows_or_cols="rows", subtract_1=TRUE,  format="dna", datablock=FALSE)
	{
	defaults='
	phyDat_seqs=phyDat_010
	outfn="phyDat_010.nex"
	printflag=FALSE
	rows_or_cols="cols"
	subtract_1=TRUE
	format="restriction"
	datablock=TRUE
	'
	
	
	#numseqs = length(phyDat_seqs)
	numseqs = nrow(phyDat_seqs)
	
	# convert sequences to sequence list, with names as taxa names
	seqs3 = phyDat_to_seqlist(phyDat_seqs, rows_or_cols=rows_or_cols, subtract_1=subtract_1)
	names(seqs3) = names(phyDat_seqs)
	
	# write to NEXUS using Nick's function, not default	
	#write.nexus.data(seqs3, file=outfn, interleaved=FALSE, gap = "-", missing = "?")
	write_nexus_data2(seqs3, file=outfn, format=format, datablock=datablock, interleaved=FALSE, gap = "-", missing = "?")
	
	if (printflag == TRUE)
		{
		moref(outfn)
		}
	
	printstr = paste("phyDat_to_nexus_good(): ", numseqs, " sequences written to ", outfn, "\n", sep="")
	cat(printstr)

	return(outfn)
	}


dtf_to_phylip <- function(char_dtf, outfn="tmpdata.phylip", txtTF=FALSE, max_namelength="max")
	{
	default='
	char_dtf = t(morph_df)
	outfn="tmpdata.phylip"
	txtTF=TRUE
	max_namelength="max"
	multistate_chars_as="paren"
	'
	
	# write a data.frame to phylip
	# you may want to put columns with different numbers of states into different columns!!
	

	
	# convert columns to character
	d = apply(char_dtf, 2, paste)
	d2 = as.data.frame(d)
	tmplevels = as.character(sort(unique(unlist(char_dtf))))
	numcols = ncol(d2)
	
	
	
	str1 = paste(nrow(d2), numcols, sep=" ")
	list_of_strings = NULL
	list_of_strings = c(list_of_strings, str1)
	
	# Maximum numer of characters
	if (max_namelength == "max")
		{
		max_namelength = max(nchar(row.names(char_dtf))) + 1
		}
	
	for (i in 1:nrow(d2))
		{
		tmpname = row.names(char_dtf)[i]
		numchars = nchar(tmpname)
		if(numchars > max_namelength)
			{
			tmpname = strtrim(tmpname, max_namelength)
			numspaces_to_add = 0
			} else {
			numspaces_to_add = max_namelength-numchars
			}
		
		spaces_to_add = paste(rep(" ", numspaces_to_add), collapse="")
		
		charstring = paste(d2[i, ], collapse="")
		outrowstr = paste(tmpname, spaces_to_add, " ", charstring, sep="")
		
		list_of_strings = c(list_of_strings, outrowstr)
		}
	
	
	if (txtTF == FALSE)
		{
		# Start writing to file
		write(list_of_strings[[1]], file=outfn, ncolumns=1, append=FALSE, sep="")
		
		# Continue writing
		for (i in 2:length(list_of_strings))
			{
			outrowstr = list_of_strings[[i]]
			write(outrowstr, file=outfn, ncolumns=1, append=TRUE, sep="")
			} # END for (i in 1:nrow(d2))
		return(outfn)
		} # END if (txtTF == FALSE)
	
	# Otherwise, paste everything together into a string
	seqs_txt = paste(unlist(list_of_strings), collapse="\n", sep="")
	
	return(seqs_txt)
	}


resample_dtf_matrix <- function(char_dtf)
	{
	nrows = nrow(char_dtf)
	new_dtf = apply(char_dtf, 2, sample, size=nrows, replace=FALSE)
	row.names(new_dtf) = row.names(char_dtf)
	return(new_dtf)
	}







steps_dtf <- function(tree, char_dtf, by_char=FALSE)
	{
	# Calculate CI, given a tree and a data.frame
	# Note: row.names should be correct, order doesn't matter
	# For this to work, for some reason you MUST use non-default (Sankoff) parsimony,
	# but on default (equal weights)
	# method="fitch" seems to mis-calculate on the star trees
	# (calculating g)
	# 	# From IB200a notes:
	# ci = m/s  -----where m = minimum number of steps in a character (number of states -1) 
	#       s = steps actually realized on a given tree 

	
	# For CI by-character
	if (by_char == TRUE)
		{
		ncols = ncol(char_dtf)
		steps_list = rep(NA, ncols)
		
		for (i in 1:ncols)
			{
			# Get one column
			char_col = adf(char_dtf[,i])
			row.names(char_col) = row.names(char_dtf)
			
			# actual steps for that column
			phyd = dtf_to_phyDat(char_col)
			steps = parsimony(tree, phyd, method="sankoff")
			
			steps_list[i] = steps
			}
		
		return(steps_list)
		}
	
	# actual steps (Sankoff parsimony)
	phyd = dtf_to_phyDat(char_dtf)
	steps = parsimony(tree, phyd, method="sankoff")

	return(steps)	
	}


g_dtf <- function(tree, char_dtf, by_char=FALSE)
	{
	# Note: row.names should be correct, order doesn't matter
	# For this to work, for some reason you MUST use non-default (Sankoff) parsimony,
	# but on default (equal weights)
	# method="fitch" seems to mis-calculate on the star trees
	# (calculating g)

	# polytomize tree to a star
	biggest_brlen = max(tree$edge.length)
	star_tree = di2multi(tree, tol=(biggest_brlen+1) )

	# calculate m
	minsteps_list = apply(char_dtf, 2, function(x) length(unique(x)))

	
	# For CI by-character
	if (by_char == TRUE)
		{
		ncols = ncol(char_dtf)
		g_list = rep(NA, ncols)
		
		for (i in 1:ncols)
			{
			# Get one column
			char_col = adf(char_dtf[,i])
			row.names(char_col) = row.names(char_dtf)
			
			# actual steps for that column
			phyd = dtf_to_phyDat(char_col)
			
			# g for that column
			g = parsimony(star_tree, phyd, method="sankoff")
			g_list[i] = g
			}
		
		return(g_list)
		}
	
	# actual steps (Sankoff parsimony)
	phyd = dtf_to_phyDat(char_dtf)

	# g for that column
	g = parsimony(star_tree, phyd, method="sankoff")

	return(g)
	}


minsteps_dtf <- function(tree, char_dtf, by_char=FALSE)
	{
	# calculate m
	minsteps_list = apply(char_dtf, 2, function(x) length(unique(x)))
	
	if (by_char == TRUE)
		{
		return(minsteps_list)
		}
	
	return(sum(minsteps_list))
	}


CI_dtf <- function(tree, char_dtf, by_char=FALSE)
	{
	# Calculate CI, given a tree and a data.frame
	# Note: row.names should be correct, order doesn't matter
	# For this to work, for some reason you MUST use non-default (Sankoff) parsimony,
	# but on default (equal weights)
	# method="fitch" seems to mis-calculate on the star trees
	# (calculating g)
	# 	# From IB200a notes:
	# ci = m/s  -----where m = minimum number of steps in a character (number of states -1) 
	#       s = steps actually realized on a given tree 

	
	# calculate m
	minsteps_list = apply(char_dtf, 2, function(x) length(unique(x)))
	
	# For CI by-character
	if (by_char == TRUE)
		{
		ncols = ncol(char_dtf)
		steps_list = rep(NA, ncols)
		
		for (i in 1:ncols)
			{
			# Get one column
			char_col = adf(char_dtf[,i])
			row.names(char_col) = row.names(char_dtf)
			
			# actual steps for that column
			phyd = dtf_to_phyDat(char_col)
			steps = parsimony(tree, phyd, method="sankoff")
			
			steps_list[i] = steps
			}
		
		CIvals = minsteps_list / steps_list
		return(CIvals)
		}
	
	# For CI totals
	minsteps = sum(minsteps_list)
	
	# actual steps (Sankoff parsimony)
	phyd = dtf_to_phyDat(char_dtf)
	steps = parsimony(tree, phyd, method="sankoff")

	# calculate CI
	CIval = minsteps/steps
	return(CIval)	
	}



RI_dtf <- function(tree, char_dtf, by_char=FALSE)
	{
	# Note: row.names should be correct, order doesn't matter
	# For this to work, for some reason you MUST use non-default (Sankoff) parsimony,
	# but on default (equal weights)
	# method="fitch" seems to mis-calculate on the star trees
	# (calculating g)

	# From IB200a
	# ri = (g - s)/(g - m)  ------where g= minimum steps on the worst tree (=bush) 
	
	# polytomize tree to a star
	biggest_brlen = max(tree$edge.length)
	star_tree = di2multi(tree, tol=(biggest_brlen+1) )

	# calculate m
	minsteps_list = apply(char_dtf, 2, function(x) length(unique(x)))

	
	# For CI by-character
	if (by_char == TRUE)
		{
		ncols = ncol(char_dtf)
		steps_list = rep(NA, ncols)
		g_list = rep(NA, ncols)
		
		for (i in 1:ncols)
			{
			# Get one column
			char_col = adf(char_dtf[,i])
			row.names(char_col) = row.names(char_dtf)
			
			# actual steps for that column
			phyd = dtf_to_phyDat(char_col)
			steps = parsimony(tree, phyd, method="sankoff")
			steps_list[i] = steps
			
			# g for that column
			g = parsimony(star_tree, phyd, method="sankoff")
			g_list[i] = g
			}
		
		
		RIvals = (g_list - steps_list) / (g_list - minsteps_list)
		return(RIvals)
		}
	
	# For RCI totals
	minsteps = sum(minsteps_list)
	
	# actual steps (Sankoff parsimony)
	phyd = dtf_to_phyDat(char_dtf)
	steps = parsimony(tree, phyd, method="sankoff")

	# g for that column
	g = parsimony(star_tree, phyd, method="sankoff")

	# calculate CI
	RIval = (g - steps) / (g - minsteps)
	return(RIval)
	}

parsim_stats_fast <- function(tree, char_dtf, by_char=TRUE)
	{
	# Note: row.names should be correct, order doesn't matter
	# For this to work, for some reason you MUST use non-default (Sankoff) parsimony,
	# but on default (equal weights)
	# method="fitch" seems to mis-calculate on the star trees
	# (calculating g)

	# From IB200a
	# ri = (g - s)/(g - m)  ------where g= minimum steps on the worst tree (=bush) 
	
	# polytomize tree to a star
	biggest_brlen = max(tree$edge.length)
	star_tree = di2multi(tree, tol=(biggest_brlen+1) )

	# calculate m
	minsteps_list = apply(char_dtf, 2, function(x) length(unique(x))-1)

	
	# For CI by-character
	if (by_char == TRUE)
		{
		ncols = ncol(char_dtf)
		steps_list = rep(NA, ncols)
		g_list = rep(NA, ncols)
		
		for (i in 1:ncols)
			{
			# Get one column
			char_col = adf(char_dtf[,i])
			row.names(char_col) = row.names(char_dtf)
			
			# actual steps for that column
			phyd = dtf_to_phyDat(char_col)
			steps = parsimony(tree, phyd, method="sankoff")
			steps_list[i] = steps
			
			# g for that column
			g = parsimony(star_tree, phyd, method="sankoff")
			g_list[i] = g
			}
		
		
		#RIvals = (g_list - steps_list) / (g_list - minsteps_list)
		
		parsim_counts = adf(cbind(minsteps_list, g_list, steps_list))
		names(parsim_counts) = c("minstep", "g", "nstep")
		
		parsim_counts_ttls = colSums(parsim_counts)
		parsim_counts = rbind(parsim_counts, parsim_counts_ttls)
		
		CI = round(parsim_counts$minstep / parsim_counts$nstep, 3)
		RI =  round((parsim_counts$g - parsim_counts$nstep) / (parsim_counts$g - parsim_counts$minstep), 3)
		RCI =  round( RI * CI, 3)
		
		parsim_counts = cbind(parsim_counts, CI, RI, RCI)
		
		row.names(parsim_counts)[nrow(parsim_counts)] = "total"
		
		return(parsim_counts)
		}
	
	# For RCI totals
	minsteps = sum(minsteps_list)
	
	# actual steps (Sankoff parsimony)
	phyd = dtf_to_phyDat(char_dtf)
	steps = parsimony(tree, phyd, method="sankoff")

	# g for that column
	g = parsimony(star_tree, phyd, method="sankoff")

	# calculate CI
	CI = round((minsteps / steps), 3)
	RI = round((g - steps) / (g - minsteps), 3)
	RCI = round(CI * RI, 3)
	
	parsim_ttls = adf(c(minsteps, g, steps, CI, RI, RCI))
	names(parsim_ttls) = c("minstep", "g", "nstep", "CI", "RI", "RCI")
	row.names(parsim_ttls) = "total"
	return(parsim_ttls)
	}



RCIcalc <- function(CIvals, RIvals)
	{
	# Note: row.names should be correct, order doesn't matter
	# For this to work, for some reason you MUST use non-default (Sankoff) parsimony,
	# but on default (equal weights)
	# method="fitch" seems to mis-calculate on the star trees
	# (calculating g)

	# These problems for CI noted above may be overcome by excluding autapomorphies OR calculating a Rescaled 
	# Consistency Index.   
	# 
	# RC = RI*CI   
	# 
	# This removes the impact of any characters that do not
	# contribute to the "fit" of the data to the tree (e.g., 
	# autapomorphies ci=1.0 and ri=0.0) 
	# 
	
	RCIvals = CIvals * RIvals
	return(RCIvals)
	}

CI_phyd <- function(tree, phyd)
	{
	# Note: row.names should be correct, order doesn't matter
	# For this to work, for some reason you MUST use non-default (Sankoff) parsimony,
	# but on default (equal weights)
	# method="fitch" seems to mis-calculate on the star trees
	# (calculating g)

	# Calculate total CI, given a phyDat dataset and a tree
	# To get individual CIs, use CI_dtf
	#
	# From IB200a notes:
	# ci = m/s  -----where m = minimum number of steps in a character (number of states -1) 
	#       s = steps actually realized on a given tree 

	# convert phyDat object to data.frame:
	x = as.character(phyd)
	minsteps_list = apply(x, 2, function(x) length(unique(x)))
	minsteps = sum(minsteps_list)
	
	# actual steps
	steps = parsimony(tree, phyd, method="sankoff")

	# calculate CI
	CIval = minsteps/steps
	return(CIval)	
	}


parsim_stats <- function(tree, char_dtf)
	{
	# Calculate the standard parsimony statistics
	CIval = CI_dtf(tr2, char_dtf)
	CIvals = CI_dtf(tr2, char_dtf, by_char=TRUE)
	
	RIval = RI_dtf(tr2, char_dtf)
	RIvals = RI_dtf(tr2, char_dtf, by_char=TRUE)

	RCIval = RCIcalc(CIval, RIval)
	RCIvals = RCIcalc(CIvals, RIvals)

	steps_list = steps_dtf(tr2, char_dtf, by_char=TRUE)
	steps = sum(steps_list)
	
	minsteps_list = minsteps_dtf(tr2, char_dtf, by_char=TRUE)
	minsteps = sum(minsteps)

	g_list = g_dtf(tr2, char_dtf, by_char=TRUE)
	g_ttl = sum(g_list)

	
	# Assemble the output
	minstep = c(minsteps_list, minsteps)
	g = c(g_list, g_ttl)
	nstep = c(steps_list, steps)
	CI = round(c(CIvals, CIval), 3)
	RI = round(c(RIvals, RIval), 3)
	RCI = round(c(RCIvals, RCIval), 3)
	
	parstats = cbind(minstep, g, nstep, CI, RI, RCI)
	parstats = adf(parstats)
	row.names(parstats) = c(names(char_dtf), "total")
	
	return(parstats)
	
	}



#######################################################
# Calculate the null distribution and p-values for your parsimony stats
#######################################################
parsim_stats_null <- function(tr2, dtf, parsim_stats, numsims=10, return_null_sims=TRUE)
	{
		
	numtips = length(tr2$tip.label)
	null_array = array(data=NA, dim=c(dim(parsim_stats), numsims))
	dim(null_array)
	nrows = dim(null_array)[1]
	ncols = dim(null_array)[2]
	
	for (i in 1:numsims)
	#for (i in 1:1)
		{
		cat("Calculating null #", i, "\n", sep="")
		# Sample tip names WITHOUT replacement
		random_order = sample(x=1:numtips, size=numtips, replace=FALSE)
		dtf_random = dtf
		row.names(dtf_random) = row.names(dtf_random)[random_order]
		
		# Calculate states on the null data
		null_parsim_stats = parsim_stats_fast(tree=tr2, char_dtf=dtf_random, by_char=TRUE)
		# print(null_parsim_stats)
		null_array[1:nrows,1:ncols,i] = matrix(unlist_df2(null_parsim_stats))
		}
	
	# Get the mean, sd, etc
	mean_null = apply(X=null_array, MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
	sd_null = apply(X=null_array, MARGIN=c(1,2), FUN=sd, na.rm=TRUE)
	min_null = apply(X=null_array, MARGIN=c(1,2), FUN=min, na.rm=TRUE)
	max_null = apply(X=null_array, MARGIN=c(1,2), FUN=max, na.rm=TRUE)
	lower05 = apply(X=null_array, MARGIN=c(1,2), FUN=quantile, probs=0.05, na.rm=TRUE)
	upper95 = apply(X=null_array, MARGIN=c(1,2), FUN=quantile, probs=0.95, na.rm=TRUE)
	lower025 = apply(X=null_array, MARGIN=c(1,2), FUN=quantile, probs=0.05, na.rm=TRUE)
	upper975 = apply(X=null_array, MARGIN=c(1,2), FUN=quantile, probs=0.95, na.rm=TRUE)
	
	# format to data.frame with the same column & row names as parsim_stats
	mean_null = adf2(mean_null)
	sd_null = adf2(sd_null)
	min_null = adf2(min_null)
	max_null = adf2(max_null)
	lower05_null = adf2(lower05)
	upper95_null = adf2(upper95)
	lower025_null = adf2(lower025)
	upper975_null = adf2(upper975)
	
	names(mean_null) = names(parsim_stats)
	names(sd_null) = names(parsim_stats)
	names(min_null) = names(parsim_stats)
	names(max_null) = names(parsim_stats)
	names(lower05_null) = names(parsim_stats)
	names(upper95_null) = names(parsim_stats)
	names(lower025_null) = names(parsim_stats)
	names(upper975_null) = names(parsim_stats)
	
	rownames(mean_null) = rownames(parsim_stats)
	rownames(sd_null) = rownames(parsim_stats)
	rownames(min_null) = rownames(parsim_stats)
	rownames(max_null) = rownames(parsim_stats)
	rownames(lower05_null) = rownames(parsim_stats)
	rownames(upper95_null) = rownames(parsim_stats)
	rownames(lower025_null) = rownames(parsim_stats)
	rownames(upper975_null) = rownames(parsim_stats)
	
	# calculate the p-value on the hypothesis that the observed value is located within the 
	# confidence intervals of the null distribution
	pnorm(q=1, mean=0.1, sd=0.1, lower.tail=FALSE)

	# If the sd is 0, due to invariance in the simulations (e.g. autapomorphic characters),
	# replace with a small positive value, so that if observed = mean, 
	# this is different than observed > mean...
	zero_SD_TF = sd_null[,3:6] == 0
	num_zero_SDvals = sum(zero_SD_TF)
	if (num_zero_SDvals > 0)
		{
		cat("Warning: your SD values on the null distribution had ", num_zero_SDvals, " cells where \nstandard deviation (sd) == 0. \nThese are being replaced with 1e-20 so that if observed = mean,\nthis is different than observed > mean.\n\n", sep="")
		
		# Fix the SD
		sd_null[,3:6][zero_SD_TF] = 1e-20
		}
	
	
	
	pvals = mapply(FUN=pnorm, q=parsim_stats, mean=mean_null, sd=sd_null, lower.tail=FALSE, SIMPLIFY=TRUE, USE.NAMES=TRUE)
	pvals = adf2(pvals)
	rownames(pvals) = rownames(parsim_stats)
	# invariant columns don't have meaningful p-values
	pvals$minstep = NA
	pvals$g = NA
	pvals$g = NA
	# for number of steps, the H1 hypothesis is that the observed value is LOWER, not HIGHER, than the null
	pvals$nstep = 1 - pvals$nstep
	
	# Merge everything
	parsim_stats_results = list()
	parsim_stats_results$parsim_stats = parsim_stats
	parsim_stats_results$mean_null = mean_null
	parsim_stats_results$sd_null = sd_null
	parsim_stats_results$min_null = min_null
	parsim_stats_results$max_null = max_null
	parsim_stats_results$lower05_null = lower05_null
	parsim_stats_results$upper95_null = upper95_null
	parsim_stats_results$lower025_null = lower025_null
	parsim_stats_results$upper975_null = upper975_null
	parsim_stats_results$pvals = pvals
	
	
	
	if (return_null_sims == TRUE)
		{
		parsim_stats_results$null_array = null_array
		}
		
	return(parsim_stats_results)
	}













# Get the log likelihood given data, phylogeny, and Q rate matrix
# (diagonals are negative so rows add to 0)
calc_loglike <- function(x, phy, Qmat, uncertain_tips=FALSE, tipstate_probs=NULL)
	{
	nl = nrow(Qmat)
	
	nb.node = phy$Nnode
	nb.tip = length(phy$tip.label)
	
	# likelihoods computed at all of these nodes
	comp <- numeric(nb.tip + nb.node)
	
	liks <- matrix(0, nb.tip + nb.node, nl)
	TIPS <- 1:nb.tip
	
	# initial likelihoods of all the observed states = 1
	# or, you could put site probabilities here
	if (uncertain_tips == TRUE)
		{
		# tipstate probabilities IN THE RIGHT ORDER!!!
		# a numeric matrix with
		# rows = tips (species)
		# cols = state probabilities
		liks[TIPS, ] <- tipstate_probs
		}
	else
		{
		liks[cbind(TIPS, x)] <- 1
		}
	
	# This is CRUCIAL!!
	phy <- reorder(phy, "pruningwise")

	for (i in seq(from = 1, by = 2, length.out = nb.node))
		{
		j <- i + 1
		anc <- phy$edge[i, 1]
		des1 <- phy$edge[i, 2]
		des2 <- phy$edge[j, 2]
		
		#print("Checking Q matrix")
		#cat("p=", p, ", rate=", rate, "\n", sep=" ")
		#print(Q)
		#print(phy$edge.length[i])
		v.l <- matexpo(Qmat * phy$edge.length[i]) %*% liks[des1, 
		  ]
		v.r <- matexpo(Qmat * phy$edge.length[j]) %*% liks[des2, 
		  ]


		#v.l <- expm(Qmat * phy$edge.length[i], method="Ward77") %*% liks[des1, 
		#  ]
		#v.r <- expm(Qmat * phy$edge.length[j], method="Ward77") %*% liks[des2, 
		#  ]


		#v.l <- exp(Qmat * phy$edge.length[i]) %*% liks[des1, 
		#  ]
		#v.r <- exp(Qmat * phy$edge.length[j]) %*% liks[des2, 
		#  ]


		v <- v.l * v.r
		comp[anc] <- sum(v)
		liks[anc, ] <- v/comp[anc]
		}
	
	# why times 2??
	#output_loglike = -2 * sum(log(comp[-TIPS]))
	output_loglike = 2 * sum(log(comp[-TIPS]))
	
	# matexpo sometimes produces log-likelihoods greater than 0 (???)
	if (output_loglike >= 0)
		{
		output_loglike = NaN
		}
	
	return(output_loglike)
	}


calc_loglike_lotsa_sims <- function(y, phy, Qmat)
	{
	nl = nrow(Qmat)
	
	nb.node = phy$Nnode
	nb.tip = length(phy$tip.label)
	
	# likelihoods computed at all of these nodes
	comp <- numeric(nb.tip + nb.node)
	
	liks <- matrix(0, nb.tip + nb.node, nl)
	TIPS <- 1:nb.tip
	
	# This is CRUCIAL!!
	phy <- reorder(phy, "pruningwise")
	
	numsims = dim(y)[3]
	output_loglike_list = array(NA, numsims)
	for (simn in 1:numsims)
		{
		x = y[, , simn]
		
		# initial likelihoods of all the observed states = 1
		liks[cbind(TIPS, x)] <- 1
		

		for (i in seq(from = 1, by = 2, length.out = nb.node))
			{
			j <- i + 1
			anc <- phy$edge[i, 1]
			des1 <- phy$edge[i, 2]
			des2 <- phy$edge[j, 2]
			
			#print("Checking Q matrix")
			#cat("p=", p, ", rate=", rate, "\n", sep=" ")
			#print(Q)
			#print(phy$edge.length[i])
			v.l <- matexpo(Qmat * phy$edge.length[i]) %*% liks[des1, 
			  ]
			v.r <- matexpo(Qmat * phy$edge.length[j]) %*% liks[des2, 
			  ]
			v <- v.l * v.r
			comp[anc] <- sum(v)
			liks[anc, ] <- v/comp[anc]
			}
		#output_loglike = -2 * sum(log(comp[-TIPS]))
		tmp_output_loglike = 2 * sum(log(comp[-TIPS]))
		output_loglike_list[j] = tmp_output_loglike
		}
	
	if (sum(is.na(output_loglike_list)) > 0)
		{
		print(paste("ERROR: ", sum(is.na(output_loglike_list)), " sims produced NAs for log like", sep=""))
		output_loglike = NaN
		}
	else
		{
		output_loglike = sum(output_loglike_list)
		}
	return(output_loglike)
	}






# Fill a rate matrix from CRP (Chinese Restaurant Process)
make_rate_matrix = function(nrows, alpha=1, rel_rate_prior=c(0,10))
	{
	num_rates = nrows * nrows - nrows
	temprates = rchinese(num_rates, alpha)
	
	uniqs = sort(unique(temprates))
	num_uniq = length(uniqs)
	lambdas = rep(NA, num_uniq)
	lambdas = runif(num_uniq, rel_rate_prior[1], rel_rate_prior[2])
	
	outmat = matrix(NA, nrow=nrows, ncol=nrows)
	diag(outmat) = 0
	outmat[is.na(outmat)] = temprates
	
	for (i in 1:num_uniq)
		{
		outmat[outmat==uniqs[i]] = lambdas[i]
		}
	
	diag(outmat) = -rowSums(outmat)
	
	return(outmat)
	}




# Get plot coordinates from the last plotted phylogeny
get_phylo3_plotcoords <- function()
	{
	zz = get(x="last_plot.phylo", envir=.PlotPhyloEnv)
	
	x = zz$xx
	y = zz$yy

	xy = as.data.frame(cbind(x,y))
	return(xy)
	}






## Moved to BioGeoBEARS_add_fossils_randomly_v1.R
# trace_parents_up <- function(nodenum, t, depthtime)
# 	{
# 	# Trace from a node up to its parents etc., a specified distance
# 	parent_node = get_parent_for_trace_parents_up(nodenum, t)
# 	
# 	# print nodenum
# 	#print(nodenum)
# 	#print(parent_node)
# 	#print(depthtime)
# 	
# 	length_to_parent = t$edge.length[t$edge[,2] == nodenum]
# 	#cat("length_to_parent: ", length_to_parent, ", length=", length(length_to_parent), "\n", sep="")
# 	#cat("depthtime: ", depthtime, "\n", sep="")
# 	if (length(length_to_parent) == 0)
# 		{
# 		print("ERROR: trace_parents_up() -- no length_to_parent returned, probably overshot bottom of tree")
# 		return(NA)
# 		}
# 	if (length_to_parent == depthtime)
# 		{
# 		print("ERROR: trace_parents_up() doesn't want to find an EXACT match between depthtime and a node; this will lead to problems in hook addition!")
# 		}
# 	if (length_to_parent > depthtime)
# 		{
# 		# you're done!
# 		return(nodenum)
# 		}
# 	else
# 		{
# 		# burrow up to parents
# 		depthtime = depthtime - length_to_parent
# 		parent_node = trace_parents_up(parent_node, t, depthtime)
# 		return(parent_node)
# 		}
# 	return(nodenum)
# 	}

## Moved to BioGeoBEARS_add_fossils_randomly_v1.R

# get_parent_for_trace_parents_up <- function(nodenum, t, printflag=FALSE)
# 	{
# 	matching_edges = findall(nodenum, t$edge[,2])
# 	parent_nodenum = t$edge[,1][matching_edges][1]
# 	
# 	if (printflag)
# 		{
# 		print(paste("nodenum=", nodenum, " parent_nodenum=", parent_nodenum, sep=""))
# 		}
# 	if (is.na(parent_nodenum))
# 		{
# 		if (printflag)
# 			{
# 			print(paste("get_parent(): node ", nodenum, " has no parent, it's probably the root!\nAnd you missed whatever parent you were actually trying to find!", sep=""))
# 			}
# 		}
# 	return(parent_nodenum)
# 	}
# 	

## Moved to BioGeoBEARS_add_fossils_randomly_v1.R
# 
# dist_between_direct_ancestors <- function(ancestor_node, descendant_node, t, totaldist=0, printflag=FALSE)
# 	{
# 	# Recursive algorithm to get distance between descendent and ancestor
# 	
# 	# Error trap should operate before this
# 	
# 	if (ancestor_node == descendant_node)
# 		{
# 		if (printflag)
# 			{
# 			print("dist_between_direct_ancestors(): ancestor_node == descendant_node")
# 			}
# 		return(totaldist)
# 		}
# 	
# 	parent_node = get_parent(descendant_node, t)
# 	dist_to_parent = t$edge.length[t$edge[,2] == descendant_node]
# 	
# 	totaldist = dist_to_parent + totaldist
# 	
# 	#print(paste(parent_node, ancestor_node, sep=""))
# 	if (parent_node == ancestor_node)
# 		{
# 		return(totaldist)
# 		}
# 	else
# 		{
# 		totaldist = dist_between_direct_ancestors(ancestor_node, parent_node, t, totaldist)
# 		return(totaldist)
# 		}
# 	}


## Moved to BioGeoBEARS_add_fossils_randomly_v1.R
# get_daughter_nodes <- function(nodenum, tr, nodes=NULL)
# 	{
# 	if(is.null(nodes))
# 		{
# 		nodes = vector()
# 		}
# 	daughter_nodes = tr$edge[which(tr$edge[,1]==nodenum),2]
# 	
# 	# Error check, in case the starting nodenum is a tip
# 	if ((length(daughter_nodes) == 0) && (length(nodes)==0))
# 		{
# 		nodes = c(nodes, nodenum)
# 		} else {
# 		nodes = c(nodes, daughter_nodes)
# 		}
# 		
# 	internal_nodes_indices = which(daughter_nodes > length(tr$tip.label))
# 	if(length(internal_nodes_indices) > 0)
# 		{
# 		for (i in 1:length(internal_nodes_indices))
# 			{
# 			nodes = get_daughter_nodes(nodenum=daughter_nodes[internal_nodes_indices[i]], tr=tr, nodes=nodes)
# 			}
# 		}
# 	return(nodes)
# 	}

## Moved to BioGeoBEARS_add_fossils_randomly_v1.R
# get_daughter_tipnums <- function(nodenum, tr)
# 	{
# 	nodes = get_daughter_nodes(nodenum, tr, nodes=NULL)
# 	tips_TF = nodes <= length(tr$tip.label)
# 	tipnums = nodes[tips_TF]
# 	return(tipnums)
# 	}


# get_daughter_tipnames <- function(nodenum, tr)
# 	{
# 	tipnums = get_daughter_tipnums(nodenum, tr)
# 	tipnames = tr$tip.label[tipnums]
# 	return(tipnames)
# 	}

## moved to BioGeoBEARS_stratified_v01.R
# get_daughters <- function(nodenum, t)
# 	{
# 	daughter_edgenums = findall(nodenum, t$edge[,1])
# 	daughter_nodenums = t$edge[,2][daughter_edgenums]
# 	return(daughter_nodenums)
# 	}

# get_parent <- function(nodenum, t)
# 	{
# 	matching_edges = findall(nodenum, t$edge[,2])
# 	parent_nodenum = t$edge[,1][matching_edges][1]
# 	return(parent_nodenum)
# 	}
# 	
# get_level <- function(nodenum, t, tmplevel=0)
# 	{
# 	parent_nodenum = get_parent(nodenum, t)
# 	if (is.na(parent_nodenum))
# 		{
# 		#tmplevel = 0
# 		return(tmplevel)
# 		}
# 	else
# 		{
# 		#print(paste("parent_nodenum: ", parent_nodenum, " level: ", tmplevel, sep=""))
# 		tmplevel = tmplevel + 1
# 		tmplevel = get_level(parent_nodenum, t, tmplevel)
# 		return(tmplevel)
# 		}
# 	# If an error occurs
# 	return(NA)
# 	}

get_node_info <- function(nodenum, t, printthis=FALSE)
	{
	# find all of the edges descending from this node
	edgenums = findall(nodenum, t$edge[,1])
	# find the daughters of these nodes
	daughters = t$edge[edgenums,2]
	# find the length of the branches attaching to these nodes
	brlens = t$edge.length[edgenums]
	
	# find all of the edges mother to this node (hopefully, just one!)
	mother_edgenums = findall(nodenum, t$edge[,2])
	if (length(mother_edgenums) > 1)
		{
		print(paste("ERROR: get_node_info found more than one mother for node ", nodenum, sep=""))
		}
	if (length(mother_edgenums) == 0)
		{
		print(paste("get_node_info found root at node ", nodenum, sep=""))
		mother = ""
		}
	mother = t$edge[mother_edgenums, 1]
	mother_brlens = t$edge.length[mother_edgenums]
	
	if (printthis == TRUE)
		{
		print(paste("edgenums: ", edgenums, sep=""))
		print(paste("daughters: ", daughters, sep=""))
		print(paste("brlens: ", brlens, sep=""))
		print(paste("mother_edgenums: ", mother_edgenums, sep=""))
		print(paste("mother: ", mother, sep=""))
		print(paste("mother_brlens: ", mother_brlens, sep=""))
		}
	
	# Store data in a structure
	nodeinfo = c()
	nodeinfo$numdaughters = length(daughters)
	nodeinfo$edgenums = edgenums
	nodeinfo$daughters = daughters
	nodeinfo$brlens = brlens
	nodeinfo$mother_edgenums = mother_edgenums
	nodeinfo$mother = mother
	nodeinfo$mother_brlens = mother_brlens
	return(nodeinfo)	
	}



# get_edge_times_before_present <- function(t)
# 	{
# 	#height above root
# 	hts_at_end_of_branches_aka_at_nodes = t$edge.length
# 	hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
# 	h = hts_at_end_of_branches_aka_at_nodes
# 
# 	# times before present, below (ultrametric!) tips
# 	# numbers are positive, i.e. in millions of years before present
# 	#                       i.e. mybp, Ma
# 	times_before_present = get_max_height_tree(t) - h
# 
# 	
# 	# fill in the ages of each node for the edges
# 	edge_ages = t$edge
# 	edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
# 	edge_ages[,2] = h[t$edge[,2]]	# top of branch
# 
# 	# fill in the times before present of each node for the edges
# 	edge_times_bp = t$edge
# 	edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
# 	edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
# 	
# 	return(edge_times_bp)
# 	}



get_edge_ages_above_root <- function(t)
	{
	#height above root
	hts_at_end_of_branches_aka_at_nodes = t$edge.length
	hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
	h = hts_at_end_of_branches_aka_at_nodes

	# fill in the ages of each node for the edges
	edge_ages = t$edge
	edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
	edge_ages[,2] = h[t$edge[,2]]	# top of branch
	
	return(edge_ages)
	}



# Trace the height of the tree from the startnode to the tip
find_sisters <- function(t, nodenum)
	{
	# Find the mother node
	mother_edgenums = findall(nodenum, t$edge[,2])
	return(mother_edgenums)
	}

# # this returns the NUMBERS identifying each node
# get_nodenums <- function(t)
# 	{
# 	# get just the unique node numbers from the edge list (left column: start node; right column: end node):
# 	nodenames = unique(c(t$edge))
# 	ordered_nodenames = nodenames[order(nodenames)]
# 	return(ordered_nodenames)
# 	}
# 
# get_nodenum_structural_root <- function(t, print_nodenum=FALSE)
# 	{
# 	#numnodes = length(t$tip.label) + length(t$node.label)
# 	#ordered_nodes = 1:length(numnodes)
# 	
# 	ordered_nodes = get_nodenums(t)
# 
# 	root_nodenums_list = c()
# 	for (n in 1:length(ordered_nodes))
# 		{
# 		tmpnode = ordered_nodes[n]
# 		if (tmpnode %in% t$edge[,2])
# 			{
# 			blah = TRUE
# 			}
# 		else
# 			{
# 			if (print_nodenum == TRUE)
# 				{
# 				cat("get_nodenum_structural_root(): Root nodenum = ", tmpnode, sep="")
# 				}
# 			root_nodenums_list = c(root_nodenums_list, tmpnode)
# 			}
# 		}
# 	return(root_nodenums_list)
# 	}

# This assumes that the root node is the first node in the numbered list 1:Nnode,
#  *after* the tip nodes, which are labeled 1:ntips
get_rootnodenum_implied <- function(t)
	{
	rootnodenum = length(t$tip.label) + 1
	return(rootnodenum)
	}

get_root_name <- function(obj)
	{
	bts = branching.times(obj)
	rootname = names(bts[1])
	cat("get_root_name(obj: Root name = ", rootname, "\n")
	return(rootname)
	}

get_root_as_num <- function(obj)
	{
	rootname = get_nodenum_structural_root(obj)
	root_node_num = as.numeric(rootname)
	return(root_node_num)
	}

# get_TF_tips <- function(obj)
# 	{
# 	# Get TF for nodes being tips
# 	
# 	# BIG CHANGE?
# 	#TF_tips = match_list1_in_list2(1:length(dists_from_root), obj$tip.label)
# 	TF_tips = match_list1_in_list2(1:length(obj$edge), 1:length(obj$tip.label))
# 	#TF_tips = obj$tip.label[TF_tips_indices]
# 	return(TF_tips)
# 	}
# 
# get_node_ages_of_tips <- function(obj)
# 	{
# 	TF_tips = get_TF_tips(obj)
# 	root_node_num = get_nodenum_structural_root(obj)
# 	dists_from_root = dist.nodes(obj)[root_node_num, ]
# 	node_ages_of_tips = dists_from_root[TF_tips]
# 	return(node_ages_of_tips)
# 	}
# 
# get_all_node_ages <- function(obj)
# 	{
# 	node_ages = dist.nodes(obj)[get_nodenum_structural_root(obj), ]
# 	return(node_ages)
# 	}
# 
# get_max_height_tree <- function(obj)
# 	{
# 	max_height = max(get_node_ages_of_tips(obj))
# 	return(max_height)
# 	}

get_node_names_of_tips <- function(obj)
	{
	TF_tips = get_TF_tips(obj)
	root_node_num = get_root_as_num(obj)
	dists_from_root = dist.nodes(obj)[root_node_num, ]
	node_names_of_tips = names(dists_from_root[TF_tips])
	return(node_names_of_tips)
	}

get_TF_tips_alive <- function(obj)
	{
	node_ages_of_tips = get_node_ages_of_tips(obj)
	TF_tips_alive = (round(node_ages_of_tips, digits=5) == round(max(node_ages_of_tips), digits=5))
	return(TF_tips_alive)
	}


# get_indices_of_tip_nodes <- function(obj)
# 	{
# 	tip_indices = 1:length(obj$tip.label)
# 	return(tip_indices)
# 	}
	
# get_indices_of_branches_under_tips <- function(obj)
# 	{
# 	tip_indices = get_indices_of_tip_nodes(obj)
# 	branchnums_under_tips = get_indices_where_list1_occurs_in_list2_noNA(tip_indices, obj$edge[, 2])
# 	return(branchnums_under_tips)
# 	}


midpoint <- function(tree)
	{
	# Root a tree at the midpoint
	# Source: https://stat.ethz.ch/pipermail/r-sig-phylo/2010-September/000750.html
	require(phangorn)
	
	dm = cophenetic(tree)
	tree = unroot(tree)
	rn = max(tree$edge)+1
	maxdm = max(dm)
	ind =  which(dm==maxdm,arr=TRUE)[1,]
	
	# "Ancestors" function requires phangorn
	#e <- simpleError("Error: 'Ancestors' function inside of 'midpoint' requires: 'library(phangorn)' in your script.")
	#tryCatch(Ancestors(tree, ind[1], "parent"), error = function(e) e, finally=print("tryCatch error message"))
	
	# If it works:
	library(phangorn)
	tmproot = Ancestors(tree, ind[1], "parent")
	
	tree = phangorn:::reroot(tree, tmproot)
	edge = tree$edge
	el = tree$edge.length
	children = tree$edge[,2]
	left = match(ind[1], children)
	tmp = Ancestors(tree, ind[2], "all")
	tmp= c(ind[2], tmp[-length(tmp)])
	right = match(tmp, children)
	if(el[left]>= (maxdm/2))
		{
		edge = rbind(edge, c(rn, ind[1]))
		edge[left,2] = rn
		el[left] = el[left] - (maxdm/2)
		el = c(el, maxdm/2)
		}
	else
		{
	    sel = cumsum(el[right])
	    i = which(sel>(maxdm/2))[1]
	    edge = rbind(edge, c(rn, tmp[i]))
	    edge[right[i],2] = rn
	    eltmp =  sel[i] - (maxdm/2)
		#el = c(el, sel[i] - (maxdm/2))
	    el = c(el, el[right[i]] - eltmp)
	    el[right[i]] = eltmp
		}
	tree$edge.length = el
	tree$edge=edge
	tree$Nnode  = tree$Nnode+1
	phangorn:::reorderPruning(phangorn:::reroot(tree, rn))
	}


midpoint2 <- function(tree, counter)
	{
	# Root a tree at the midpoint
	# Source: https://stat.ethz.ch/pipermail/r-sig-phylo/2010-September/000750.html
	#require(phangorn)
	cat("Midpoint rooting tree #", counter, "\n", sep="")
	
	dm = cophenetic(tree)
	tree = unroot(tree)
	rn = max(tree$edge)+1
	maxdm = max(dm)
	ind =  which(dm==maxdm,arr=TRUE)[1,]
	
	# "Ancestors" function requires phangorn
	#e <- simpleError("Error: 'Ancestors' function inside of 'midpoint' requires: 'library(phangorn)' in your script.")
	#tryCatch(Ancestors(tree, ind[1], "parent"), error = function(e) e, finally=print("tryCatch error message"))
	
	# If it works:
	#library(phangorn)
	tmproot = Ancestors(tree, ind[1], "parent")
	
	tree = phangorn:::reroot(tree, tmproot)
	edge = tree$edge
	el = tree$edge.length
	children = tree$edge[,2]
	left = match(ind[1], children)
	tmp = Ancestors(tree, ind[2], "all")
	tmp= c(ind[2], tmp[-length(tmp)])
	right = match(tmp, children)
	if(el[left]>= (maxdm/2))
		{
		edge = rbind(edge, c(rn, ind[1]))
		edge[left,2] = rn
		el[left] = el[left] - (maxdm/2)
		el = c(el, maxdm/2)
		}
	else
		{
	    sel = cumsum(el[right])
	    i = which(sel>(maxdm/2))[1]
	    edge = rbind(edge, c(rn, tmp[i]))
	    edge[right[i],2] = rn
	    eltmp =  sel[i] - (maxdm/2)
		#el = c(el, sel[i] - (maxdm/2))
	    el = c(el, el[right[i]] - eltmp)
	    el[right[i]] = eltmp
		}
	tree$edge.length = el
	tree$edge=edge
	tree$Nnode  = tree$Nnode+1
	phangorn:::reorderPruning(phangorn:::reroot(tree, rn))
	}

# 
# 
# # axisPhylo() with more flexibility in labeling
# axisPhylo2 <- function (side = 1, roundlabels=FALSE, minage=0, ...) 
# 	{
#     lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
#     if (lastPP$type %in% c("phylogram", "cladogram")) {
#         if (lastPP$direction %in% c("rightwards", "leftwards")) {
#             x <- pretty(lastPP$xx)
#             if (lastPP$direction == "rightwards") 
#                 maxi <- max(lastPP$xx)
#             else {
#                 maxi <- min(lastPP$xx)
#                 x <- -x
#             }
#         }
#         else {
#             x <- pretty(lastPP$yy)
#             if (lastPP$direction == "upwards") 
#                 maxi <- max(lastPP$yy)
#             else {
#                 maxi <- min(lastPP$yy)
#                 x <- -x
#             }
#         }
#     }
#     if (roundlabels)
# 		{
# 		axis(side = side, at = c(maxi - x), labels = abs(minage+x), ...)
# 		}
# 	else
# 		{
# 		axis(side = side, at = c(maxi - x), labels = round(abs(minage+x), digits=roundlabels), ...)
# 		}
# 	}




# http://www.mail-archive.com/r-sig-phylo@r-project.org/msg00731.html
axisPhyloFact <- function (side = 1, fact = 1, ...) 
	{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (lastPP$type %in% c("phylogram", "cladogram"))
    	{
        if (lastPP$direction %in% c("rightwards", "leftwards"))
        	{
            x <- pretty(lastPP$xx)
            if (lastPP$direction == "rightwards") 
                maxi <- max(lastPP$xx)
            else
            	{
                maxi <- min(lastPP$xx)
                x <- -x
            	}
        	}
        else
        	{
            x <- pretty(lastPP$yy)
            if (lastPP$direction == "upwards") 
                maxi <- max(lastPP$yy)
            else
            	{
                maxi <- min(lastPP$yy)
                x <- -x
            	}
        	}
    	}
    axis(side = side, at = c(maxi - x), labels = abs(x * fact), ...)
	}





extract_TMRCA_from_newick_files <- function(newick_fns, specifiers, TMRCAs_Rdata_fn="list_of_TMRCAs.Rdata", clades_Rdata_fn="list_of_clades.Rdata", runthis=TRUE)
	{
	defaults='
	newick_fns=slashslash(list.files(path=treedir, pattern=".newick", full.names=TRUE))
	specifiers=c("H.sapiens", "P.troglodytes")
	TMRCAs_Rdata_fn="list_of_TMRCAs.Rdata"
	clades_Rdata_fn="list_of_clades.Rdata"
	TMRCAs_Rdata_fn="list_of_morphology_TMRCAs.Rdata"
	clades_Rdata_fn="list_of_morphology_clades.Rdata"
	runthis = TRUE
	'
	# Get the number of input trees
	numtrees = length(newick_fns)
	
	if (runthis)
		{
		# Extract the TMRCA of H.sapiens and chimps
		list_of_TMRCAs = rep(NA, numtrees)
		list_of_clades = rep(NA, numtrees)
		
		for (i in 1:length(newick_fns))
			{
			cat("\nGetting TMRCA(", specifiers[1], ",", specifiers[2], ") for tree #", i, ": ", sep="")
			
			morphtree = read.tree(newick_fns[i])
			morphtree$tip.label
		
			dtf = prt(morphtree, printflag=FALSE, get_tipnames=TRUE)
			
			# Find lists of tipnames with commas, with "H_sapiens" and "P_troglodytes"
			tipnames = dtf$tipnames
			keepTF1 = grepl(pattern=",", x=tipnames)
			keepTF2 = grepl(pattern=specifiers[1], x=tipnames)
			keepTF3 = grepl(pattern=specifiers[2], x=tipnames)
			keepTF = (keepTF1 + keepTF2 + keepTF3) == 3
			tipnames_that_match = tipnames[keepTF]
			
			# Get the number of taxa in each matching clade; i.e. number of commas + 1
			number_of_taxa_in_each_clade = 1 + count_chars(X=tipnames_that_match, char=",")
			
			# The MRCA is the smallest one
			MRCA_clade = tipnames_that_match[number_of_taxa_in_each_clade == min(number_of_taxa_in_each_clade)]
			list_of_clades[i] = MRCA_clade
			#print(MRCA_clade)
			
			# The dtf row is 
			rownums = (1:nrow(dtf))
			row_for_MRCA_node = rownums[dtf$tipnames == MRCA_clade]
			
			time_bp = dtf[row_for_MRCA_node,]$time_bp
			list_of_TMRCAs[i] = time_bp
			cat(time_bp)
			}
		
		# Save results
		save(list_of_TMRCAs, file=TMRCAs_Rdata_fn)
		save(list_of_clades, file=clades_Rdata_fn)
		} else {
		load(file=TMRCAs_Rdata_fn)
		load(file=clades_Rdata_fn)
		}
	results = NULL
	results$list_of_TMRCAs = list_of_TMRCAs
	results$list_of_clades = list_of_clades
	
	return(results)
	}


DNAbin_to_list_of_strings <- function(tmpdata)
	{
	# Get the names of the taxa
	tmpdata_rownames = attr(tmpdata, "dimnames")[[1]]
	
	# Convert DNAbin to character matrix
	tmpdata2 = as.character(tmpdata)

	# Convert each row to a list of characters
	tmpdata3 = lapply(seq_len(nrow(tmpdata2)), function(i) tmpdata2[i,])
	
	# Add the names back in:
	names(tmpdata3) = tmpdata_rownames
	
	return(tmpdata3)
	}


NEXUS_fn_to_full_DNAbin <- function(nexus_fn)
	{
	tmpdata = read.nexus.data(file=nexus_fn)
	numrows = length(tmpdata)
	numcols = length(tmpdata[[1]])

	tmpdata1a = matrix( unlist(tmpdata, use.names=FALSE), nrow=numrows, ncol=numcols, byrow=TRUE)
	rownames(tmpdata1a) = names(tmpdata)
	tmpdata2 = as.DNAbin(tmpdata1a)
	attributes(tmpdata2)
	
	return(tmpdata2)
	}

#######################################################
# Test molecular clock
#######################################################
#
# Based on:
# Lemey, P. and D. Posada (2009). "Molecular clock analysis." The phylogenetic handbook: a practical approach to phylogenetic analysis and hypothesis testing. 
# M. Salemi, A.-M. Vandamme and P. Lemey. Cambridge, UK ; New York, Cambridge University Press: xxvi, 723 p. 
# Link: http://www.kuleuven.be/aidslab/phylogenybook/home.html
# 
test_molclock <- function(data_fn="", datatype="nexus", jModelTest_dir="/Applications/jmodeltest-2.1.3/", paup_exe="/Applications/paup", jModelTest_AIC_or_BIC="AIC", printlevel=1, returnwhat="all")
	{
	defaults='
	fasta_fn="/Users/nickm/Desktop/__projects/_hominin_phylo_dating/molecular_clock_test/mtDNA_homs_mod2_noGaps.fasta"
	data_fn=fasta_fn; datatype="fasta"; jModelTest_dir="/Applications/jmodeltest-2.1.3/"; paup_exe="/Applications/paup"; jModelTest_AIC_or_BIC="AIC"; printlevel=1; returnwhat="all"
	
	nexus_fn="/Users/nickm/Desktop/__projects/_hominin_phylo_dating/molecular_clock_test/mtDNA_homs_mod2_noGaps_tmp.nexus"
	data_fn=nexus_fn; datatype="nexus"; jModelTest_dir="/Applications/jmodeltest-2.1.3/"; paup_exe="/Applications/paup"; jModelTest_AIC_or_BIC="AIC"; printlevel=1; returnwhat="all"
	'
	
	require(ape)
	
	# Check the file types
	if (datatype == "nexus")
		{
		# Get filename and make FASTA filename
		nexus_data_fn = data_fn
		fasta_data_fn = paste(get_fn_prefix(nexus_data_fn), "_tmp.fasta", sep="")
		
		# Read the NEXUS data file (slower)
		tmpdata2 = NEXUS_fn_to_full_DNAbin(nexus_fn=nexus_data_fn)
		ntaxa = attr(tmpdata2, "dim")[1]
		
		# Make the FASTA file
		write.dna(x=tmpdata2, file=fasta_data_fn, format="fasta", append=FALSE, nbcol=-1, colsep="")
		}
	
	if (datatype == "fasta")
		{
		# Get filename and make FASTA filename
		fasta_data_fn = data_fn
		nexus_data_fn = paste(get_fn_prefix(fasta_data_fn), "_tmp.nexus", sep="")

		
		tmpdata = read.dna(file=fasta_data_fn, format="fasta")
		ntaxa = attr(tmpdata, "dim")[1]
		
		tmpdata3 = DNAbin_to_list_of_strings(tmpdata)
		write.nexus.data(x=tmpdata3, file=nexus_data_fn, format="dna", datablock=FALSE, interleaved=FALSE)
		}
	
	
	# Now, run jModelTest
	cmd = jModelTest(fasta_fn=fasta_data_fn, logfn="", cmdpart3=" -g 4 -i -f -AIC -BIC -a -w", jModelTest_dir="/Applications/jmodeltest-2.1.3/", returnwhat="cmd")
	logfn = jModelTest(fasta_fn=fasta_data_fn, logfn="", cmdpart3=" -g 4 -i -f -AIC -BIC -a -w", jModelTest_dir="/Applications/jmodeltest-2.1.3/")
	
	
	# Would you like to select the model using AIC or BIC?
	# (This assumes you ran both, in cmdpart3)
	
	if (jModelTest_AIC_or_BIC == "BIC")
		{		
		# Extract the PAUP block (AIC selection):
		tmplines = readLines(logfn)
		lines_to_keep = extract_lines_startstr_to_endstr(lines=tmplines, string_to_start_at="\\[!", string_to_end_at="END;", printflag=TRUE, include_endstring=TRUE, instance_to_find=2)
		cat(lines_to_keep, sep="\n")
		}
	
	if (jModelTest_AIC_or_BIC == "AIC")
		{	
		# Extract the PAUP block (BIC selection):
		tmplines = readLines(logfn)
		lines_to_keep = extract_lines_startstr_to_endstr(lines=tmplines, string_to_start_at="\\[!", string_to_end_at="END;", printflag=TRUE, include_endstring=TRUE, instance_to_find=1)
		cat(lines_to_keep, sep="\n")
		}

	
	# Make the PAUP commands file for running the molecular clock test!
	if (printlevel >= 1)
		{
		PAUPnexus_cmds_fn = PAUP_test_clock_cmds(nexus_data_fn=nexus_data_fn, PAUPblock_lines=lines_to_keep, PAUP_logfn="", PAUPnexus_cmds_fn="", returnwhat="PAUPnexus", paup_exe="/Applications/paup")
		moref(PAUPnexus_cmds_fn)
		}
	
		
	# Run the molecular clock test!
	PAUP_logfn = PAUP_test_clock_cmds(nexus_data_fn=nexus_data_fn, PAUPblock_lines=lines_to_keep, PAUP_logfn="", PAUPnexus_cmds_fn="", returnwhat="log", paup_exe="/Applications/paup")
	
	if (printlevel >= 1)
		{
		moref(PAUP_logfn)
		}
	
	# Extract the log-likelihood values from the PAUP output
	# Get the lines of the PAUP log file that have "-ln L" in them
	PAUP_lines = readLines(PAUP_logfn)
	matchlines_TF = grepl(pattern="-ln L", x=PAUP_lines)
	matchlines_TF
	
	LnL_lines = PAUP_lines[matchlines_TF]
	LnL_lines
	
	library(gdata)	# for trim
	LnL_lines = LnL_lines[c(2,4)]
	LnL_lines = gsub(pattern="-ln L", replacement="", LnL_lines)
	LnL_lines = trim(LnL_lines)
	LnL_vals = as.numeric(LnL_lines)
	if (printlevel >= 1)
		{
		cat("LnL_vals:\n", sep="")
		print(LnL_vals)
		}
	
	# Get the actual LnL values of the models
	# null = clock (fewer parameters)
	# alternative = non-clock (more parameters)
	alt_LnL = -1 *  LnL_vals[1]
	null_LnL = -1 * LnL_vals[2]
	
	# Number of degrees of freedom in a test of the strict clock
	# In a non-clock tree, 2n-3 branches (of an unrooted tree) have to be estimated
	# In a clock tree, you have a rooted tree, and only n-1 branch lengths to infer
	# Thus, the difference in the number of parameters is ntaxa-2, i.e. n-2
	ntaxa = ntaxa
	
	num_free_branches_nonclock = (2*ntaxa) - 3
	num_free_branches_clock = ntaxa - 1
	
	degrees_of_freedom = num_free_branches_nonclock - num_free_branches_clock

	# Difference in number of parameters
	# (this is the number of degrees of freedom)
	# delta_numparameters = numparams1 - numparams2
	
	if (returnwhat == "pval")
		{
		LRT_pval = lrttest(LnL_1=alt_LnL, LnL_2=null_LnL, numparams1=num_free_branches_nonclock, numparams2=num_free_branches_clock, returnwhat="pval")			
		# p = 0.1328429
		# not significant!
		if (printlevel >= 1)
			{
			cat("LnL_vals:\n", sep="")
			print(LnL_vals)
			}
		return(LRT_pval)
		}

	if (returnwhat == "all")
		{
		# Run Likelihood Ratio Test (LRT)
		LRT_pval = lrttest(LnL_1=alt_LnL, LnL_2=null_LnL, numparams1=num_free_branches_nonclock, numparams2=num_free_branches_clock, returnwhat="pval")			

		# Run Likelihood Ratio Test (LRT)

		LRT_AIC_results = AICstats_2models(LnL_1=alt_LnL, LnL_2=null_LnL, numparams1=num_free_branches_nonclock, numparams2=num_free_branches_clock)
		LRT_AIC_results$alt = "nonclock"
		LRT_AIC_results$null = "clock"
		
		# p = 0.1328429
		# not significant!
		if (printlevel >= 1)
			{
			cat("LRT_AIC_results:\n", sep="")
			print(LRT_AIC_results)
			}
		return(LRT_AIC_results)
		}
	return("ERROR")		
	}





# Getting stats from a BEAST consensus tree
# From phyloch
'
file = trfn
file = confn
'
extractBEASTstats_orig <- function (file, digits=4, printflag=FALSE) 
	{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    tab <- unlist(strsplit(X, "\\["))[-1]
    tab <- gsub("&|;|\\]", "", tab)
    tab <- gsub(":.+$", "", tab)


    # This extracts ~19 objects delimited by [brackets
    # In the consensus tree output
    foo <- function(x)
    	{
        x <- unlist(strsplit(x, ","))
        x
    	}
    tab <- lapply(tab, foo)
    

	#for (i in seq(along = tab))
	# same as
	for (i in 1:length(tab))
    	{
    	# in each element in a tab(le) item, get the
    	# indices of anything containing "{"
        ind <- grep("[{]", tab[[i]])
        
        # Change the name text of those elements (remove e.g. {)
        names <- gsub("=.+$", "", tab[[i]][ind])
        tab[[i]][ind] <- gsub("[{]", "", tab[[i]][ind])
        
        # Replace = with _MIN=
        tab[[i]][ind] <- gsub("=", "_MIN=", tab[[i]][ind])
        
        # Remove the } bracket also
        tab[[i]][ind + 1] <- gsub("[}]", "", tab[[i]][ind + 1])
        
        # And call this item MAX=
        tab[[i]][ind + 1] <- paste(paste(names, "MAX=", sep = "_"), tab[[i]][ind + 1])
    	}


	# Empty data frame
    ttab <- data.frame()
    
    # Take everything before the = signs, and uniqify it
    stats <- unique(gsub("=.+$", "", unlist(tab)))
    
    # Remove "!rotate" (A FigTree code)
    stats = stats[stats != "!rotate"]
    stats
    
    # Go through all the tabs/nodes
    # Put them in a big table for each node
    for (i in seq(along = tab))
    	{
    	
    	# Go through all the summary statistics by name
        for (j in seq(along = stats))
        	{
        	# Get the index of the summary statistic you are looking for
            ind <- grep(paste("^", stats[j], "=", sep = ""), tab[[i]])
            
            # If the summary statistic is found at this node,
            if (length(ind) > 0)
            	{
            	# Remove the name information, and make the number numeric
                v <- as.numeric(gsub(paste(stats[j], "=", sep = ""), "", tab[[i]][ind]))
                
                if (printflag == TRUE)
                	{
	                cat(i, j, ind, v, "\n", sep="	")
	                }
                
                ttab[i, j] <- v
				}
			}
		}
		
	# Name the columns of the table
    colnames(ttab) <- stats
    
    # Figure out which are the tips (doesn't matter, really)
    tip <- which(is.na(ttab$posterior))
    ttab
    
    
# 	phy <- read.beast(file, digits = digits)
# 	int <- phy$Nnode
#     tips <- length(phy$tip.label)
#     node <- (tips + 1):(tips + int)
#     M <- cbind(node, ttab)
#  	   
#     return(M)
	}

# Getting stats from a BEAST consensus tree
# From phyloch
read.beast.table_original <- function (file, digits = 2) 
	{
    phy <- read.beast_original(file, digits = digits)
    int <- phy$Nnode
    stats <- phy[-(1:4)]
    M <- matrix(unlist(stats), nrow = int, byrow = FALSE)
    colnames(M) <- names(stats)
    tips <- length(phy$tip.label)
    node <- (tips + 1):(tips + int)
    M <- cbind(node, M)
    M
	}





'
file=confn
'

extractBEASTstats3 <- function (file) 
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    phy <- read.nexus(file)
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("^.*tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    X <- gsub("[!]rotate=true,*|[!]rotate=false,*", "", X)
    X <- gsub(";$", "", X)
    vals <- unlist(strsplit(X, "\\][[:alnum:]:.)]*\\[*&*"))
    foo <- function(x) {
        x <- gsub(",([[:lower:]])", "xxx\\1", x)
        x <- unlist(strsplit(x, "xxx"))
        names(x) <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[1]
        })
        x <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[2]
        })
        x <- gsub("[{}]", "", x)
        x <- strsplit(x, ",")
        lapply(x, as.numeric)
    }
    vals <- lapply(vals, foo)
    edges <- gsub("\\[(&[[:alnum:]_=%!.,{}-]+)\\]", "", X)
    edges <- unlist(strsplit(gsub("[()]*", "", edges), ":"))
    tips <- gsub("^[0-9]+.[0-9]+,*", "", edges)


	re1='(,)'								# Any Single Character 1
	#re2='(?:[a-z][a-z]*[0-9]+[a-z0-9]*)'	# Alphanum 1
	re2='(?:[a-zA-Z0-9]*)'					# Alphanum 1
    regstr = paste(re1, re2, sep="")
    
    #gsub(regstr, "", "asd,234g")


    # remove everything after comma
    # original
    # edges <- as.numeric(gsub("\\[,$|[[:alpha:]].+$\\]", "", edges))
    # modified
    edges <- as.numeric(gsub(pattern=regstr, "", edges))
    
    tab <- cbind(edges, nodes = c("root", head(tips, -1)), id = rep(NA, 
        length(edges)))
    tab[tab[, 2] == "", 2] <- "internal"
    internal <- grep("internal|root", tab[, 2])
    tab[-internal, 3] <- seq(along = phy$tip.label)
    intnodes <- match(tab[internal, 1], phy$edge.length[phy$edge[, 
        2] > 50], nomatch = 0) + 51
    tab[internal, 3] <- intnodes
    names(vals)[internal] <- paste("NODE", tab[internal, 3], 
        sep = ":")
    names(vals)[-internal] <- paste("TIP", tab[-internal, 3], 
        sep = ":")
    tips <- vals[-internal]
    nodes <- vals[internal]
}




extractMrBayesStats4 <- function (file) 
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    phy <- read.nexus(file)
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("^.*tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    X <- gsub("[!]rotate=true,*|[!]rotate=false,*", "", X)
    X <- gsub(";$", "", X)
    vals <- unlist(strsplit(X, "\\][[:alnum:]:.)]*\\[*&*"))
    foo <- function(x) {
        x <- gsub(",([[:lower:]])", "xxx\\1", x)
        x <- unlist(strsplit(x, "xxx"))
        names(x) <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[1]
        })
        x <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[2]
        })
        x <- gsub("[{}]", "", x)
        x <- strsplit(x, ",")
        lapply(x, as.numeric)
    }
    vals <- lapply(vals, foo)
    edges <- gsub("\\[(&[[:alnum:]_=%!.,{}-]+)\\]", "", X)
    edges <- unlist(strsplit(gsub("[()]*", "", edges), ":"))
    tips <- gsub("^[0-9]+.[0-9]+,*", "", edges)


	re1='(,)'								# Any Single Character 1
	#re2='(?:[a-z][a-z]*[0-9]+[a-z0-9]*)'	# Alphanum 1
	re2='(?:[a-zA-Z0-9]*)'					# Alphanum 1
    regstr = paste(re1, re2, sep="")
    
    #gsub(regstr, "", "asd,234g")


    # remove everything after comma
    # original
    # edges <- as.numeric(gsub("\\[,$|[[:alpha:]].+$\\]", "", edges))
    # modified
    edges <- as.numeric(gsub(pattern=regstr, "", edges))
    
    tab <- cbind(edges, nodes = c("root", head(tips, -1)), id = rep(NA, 
        length(edges)))
    tab[tab[, 2] == "", 2] <- "internal"
    internal <- grep("internal|root", tab[, 2])
    tab[-internal, 3] <- seq(along = phy$tip.label)
    intnodes <- match(tab[internal, 1], phy$edge.length[phy$edge[, 
        2] > 50], nomatch = 0) + 51
    tab[internal, 3] <- intnodes
    names(vals)[internal] <- paste("NODE", tab[internal, 3], 
        sep = ":")
    names(vals)[-internal] <- paste("TIP", tab[-internal, 3], 
        sep = ":")
    tips <- vals[-internal]
    nodes <- vals[internal]
}





# Reading a BEAST tree, with the stats
# from phyloch
setup='
file=confn
digits = NULL
'
read.beast_original <- function (file, digits = NULL, printflag=FALSE) 
	{
	# Scan the files in
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats_orig(file, printflag=FALSE)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)
        
    # Nodes without posterior values are tip nodes (always 100%!)
    interior <- which(!is.na(tab$posterior))

    # Right is the lines containing ]
    # -- basically the tree lines
    RIGHT <- grep("\\]", X)    
		
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT))
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    
    # Get line numbers with semicolons
    semico <- grep(";", X)
    
    # Line # for beginning of TREES block
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Find and use the TRANSLATE block, if it exists
    # Line # for beginning of TRANSLATE block
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    
    # No TRANSLATE block
    if (length(i2) != 0)
    	{    
		# Last line with TRANSLATE items
		end <- semico[semico > i2][1]
		
		# Lines with TRANSLATE items
		x <- X[(i2 + 1):end]
		x <- unlist(strsplit(x, "[,; \t]"))
		x <- x[nzchar(x)]
		
		# Get the translation information
		TRANS <- matrix(x, ncol = 2, byrow = TRUE)
		TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
		n <- dim(TRANS)[1]

		# Start & end lines of the trees block
		start <- semico[semico > i2][1] + 1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]
		} else {
		# WITHOUT TRANSLATE block
		# Start & end lines of the trees block
		start <- i1+1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]

		}
    
    
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS files...)
    tree <- gsub("^.*= *", "", tree)
    
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    
    # Make an empty vector, 1 cell for each stat
    nodestats <- vector(mode = "list", length = dim(tab)[2])

	# Go through each nodestat
    for (i in seq(along = nodestats))
    #for (i in 1:1)
    	{
        newtree <- tree
        
        # Paste together the nodestat & branchlength
        val <- tab[, i]
        ggg <- paste(val, brl, sep = "")
        ggg[length(ggg)] <- paste(tail(val, 1), ";", sep = "")
        
        # Go through the interior nodes, and replace
        # the branchlength with the node statistic in question
        for (j in interior)
        	{
        	newtree <- gsub(brl[j], ggg[j], newtree)
        	}
        dt <- read.tree(text = newtree)
        
        # Then get these as node labels
        z <- dt$node.label
        z[z == "NA"] <- 9999
        z <- as.numeric(z)
        z[z == 9999] <- NA
        
        # Tabulate the nodestats
        nodestats[[i]] <- z
        names(nodestats)[i] <- colnames(tab)[i]
    	}
    
    
    # Read the input file as plain NEXUS
    tr <- read.nexus(file)
    
    # Add the stats to the standard tree object
    tr <- c(tr[1:length(tr)], nodestats[1:length(nodestats)])
    class(tr) <- ("phylo")
    
    # Add the origin attribute
    attr(tr, "origin") <- file
    tr
	}



# Reading a BEAST tree, with the stats
# from phyloch
# but then putting into a table with prt
# adding printflags, etc.

# read_beast_prt() produces these columns in the prt table:
# names(beast_nodestats_dtf) = c("internal_nodenums_from1", "internal_nodenums", "rate_range_MIN", "rate_range_MAX", "height_95%_HPD_MIN", "height_95%_HPD_MAX", "height_median", "rate", "height", "rate_median", "height_range_MIN", "height_range_MAX", "rate_95%_HPD_MIN", "rate_95%_HPD_MAX", "posterior", "height_HPD_width", "height_range", "rate_HPD_width", "height_CV", "rate_CV")


setup='
file=confn
digits = NULL
get_tipnames=TRUE # If get_tipnames==TRUE, add a row containing the tipnames in the clade, in alphabetical order
'
read_beast_prt <- function (file, digits = NULL, get_tipnames=TRUE, printflag=FALSE) 
	{
	# Scan the files in
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats_orig(file, printflag=FALSE)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)
        
    # Nodes without posterior values are tip nodes (always 100%!)
    interior <- which(!is.na(tab$posterior))

    # Right is the lines containing ]
    # -- basically the tree lines
    RIGHT <- grep("\\]", X)    
		
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT))
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    
    # Get line numbers with semicolons
    semico <- grep(";", X)
    
    # Line # for beginning of TREES block
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)


 
 	# Find and use the TRANSLATE block, if it exists
    # Line # for beginning of TRANSLATE block
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    
    # No TRANSLATE block
    if (length(i2) != 0)
    	{    
		# Last line with TRANSLATE items
		end <- semico[semico > i2][1]
		
		# Lines with TRANSLATE items
		x <- X[(i2 + 1):end]
		x <- unlist(strsplit(x, "[,; \t]"))
		x <- x[nzchar(x)]
		
		# Get the translation information
		TRANS <- matrix(x, ncol = 2, byrow = TRUE)
		TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
		n <- dim(TRANS)[1]

		# Start & end lines of the trees block
		start <- semico[semico > i2][1] + 1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]
		} else {
		# WITHOUT TRANSLATE block
		# Start & end lines of the trees block
		start <- i1+1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]

		}
    

     
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS files...)
    tree <- gsub("^.*= *", "", tree)
    
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    
    # Make an empty vector, 1 cell for each stat
    nodestats <- vector(mode = "list", length = dim(tab)[2])

	# Go through each nodestat
    for (i in seq(along = nodestats))
    #for (i in 1:1)
    	{
        newtree <- tree
        
        # Paste together the nodestat & branchlength
        val <- tab[, i]
        ggg <- paste(val, brl, sep = "")
        ggg[length(ggg)] <- paste(tail(val, 1), ";", sep = "")
        
        # Go through the interior nodes, and replace
        # the branchlength with the node statistic in question
        for (j in interior)
        	{
        	newtree <- gsub(brl[j], ggg[j], newtree)
        	}
        dt <- read.tree(text = newtree)
        
        # Then get these as node labels
        z <- dt$node.label
        z[z == "NA"] <- 9999
        z <- as.numeric(z)
        z[z == 9999] <- NA
        
        # Tabulate the nodestats
        nodestats[[i]] <- z
        names(nodestats)[i] <- colnames(tab)[i]
    	}
    
    
    # Read the input file as plain NEXUS
    tr <- read.nexus(file)
    
    # Add the stats to the standard tree object
    tr <- c(tr[1:length(tr)], nodestats[1:length(nodestats)])
    class(tr) <- ("phylo")

    # Add the origin attribute
    attr(tr, "origin") <- file
    tr
    
    
    #######################################################
	# Get the prt table    
	#######################################################
    tr_table = prt(tr, printflag=FALSE, relabel_nodes = FALSE, time_bp_digits=7, get_tipnames=get_tipnames)
    
    # make the internal node numbers
    internal_nodenums = seq(from=length(tr$tip.label)+1, to=length(tr$tip.label)+tr$Nnode, by=1)
	internal_nodenums_from1 = seq(from=1, to=tr$Nnode, by=1)

    # Make a table of internal node stats from the BEAST stats
    height_HPD_width = tr$"height_95%_HPD_MAX" - tr$"height_95%_HPD_MIN"
    height_range = tr$"height_range_MAX" - tr$"height_range_MIN"
    rate_HPD_width = tr$"rate_95%_HPD_MAX" - tr$"rate_95%_HPD_MIN"
    
    # Calculate the CVs of height and rate
    height_CV = (height_HPD_width/(1.96*2)) / tr$height
    rate_CV = (rate_HPD_width/(1.96*2)) / tr$rate
    
    
    
    beast_nodestats_table_temp = cbind(internal_nodenums_from1, internal_nodenums, tr$rate_range_MIN, tr$rate_range_MAX, tr$"height_95%_HPD_MIN", tr$"height_95%_HPD_MAX", tr$height_median, tr$rate, tr$height, tr$rate_median, tr$height_range_MIN, tr$height_range_MAX, tr$"rate_95%_HPD_MIN", tr$"rate_95%_HPD_MAX", tr$posterior, height_HPD_width, height_range, rate_HPD_width, height_CV, rate_CV)
    
    # Add tip rows as NAs
    tiprows = matrix(data=NA, nrow=length(tr$tip.label), ncol=ncol(beast_nodestats_table_temp))
    
    beast_nodestats_table = rbind(tiprows, beast_nodestats_table_temp)
    beast_nodestats_dtf = adf2(beast_nodestats_table)
    names(beast_nodestats_dtf) = c("internal_nodenums_from1", "internal_nodenums", "rate_range_MIN", "rate_range_MAX", "height_95%_HPD_MIN", "height_95%_HPD_MAX", "height_median", "rate", "height", "rate_median", "height_range_MIN", "height_range_MAX", "rate_95%_HPD_MIN", "rate_95%_HPD_MAX", "posterior", "height_HPD_width", "height_range", "rate_HPD_width", "height_CV", "rate_CV")
    
	# head(beast_nodestats_dtf)
	# tail(beast_nodestats_dtf)
    
    
    prt_beast_nodestats = cbind(tr_table, beast_nodestats_dtf)
    
	
	# Return a list with the tree and the megatable
	beastcon = NULL
	beastcon$tr = tr
	beastcon$prt_beast_nodestats = prt_beast_nodestats
	
	return(beastcon)
	}









# Reading a BEAST tree, with the stats
# from phyloch
setup='
trfn="temp.nexus"
digits = NULL
'
read_beasttree_rates <- function (trfn, digits = NULL) 
	{
	# Scan the trfns in
    X <- scan(file = trfn, what = "", sep = "\n", quiet = TRUE)
    X_orig = X
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    lines_with_brackets = grep("\\[", X)
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats2(trfn)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)

	#interior <- which(!is.na(tab$posterior))
    
    # Right is the lines containing ]
    # -- basically the tree lines
	RIGHT <- grep("\\]", X)    
	
	
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT) > 0)
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
                # also remove from X_org
                X_orig <- X_orig[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
   	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]

	
	# Otherwise do this...
	#  If there isn't an endblock, then the last tree is the last line of interest
	#line_of_last_tree = 
	
	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]
	
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS trfns...)
    tree <- gsub("^.*= *", "", tree)
    

	#====================================
	# Get branch lengths in text trfn order    
	#====================================
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    ########################################
    
    ########################################
    # Make an empty vector, 1 cell for each stat
    nodestats <- vector(mode = "list", length = dim(tab)[2])
    ########################################


	# Go through ALL the nodes, and replace
	# the branchlength with the rate statistic
	#
	# This produces a rateogram

	rateogram <- tree
	
	for (j in 1:nrow(tab))
		{		
		new_branch_length = paste(":", as.character(tab[j,]), sep="")
		
		# Replace JUST THE FIRST match
		matches = gregexpr(brl[j], rateogram)[[1]]
		matches_end = matches-1+attr(matches,"match.length")
		
		newtree_chars = strsplit(rateogram, "")[[1]]
		firstchars = newtree_chars[1:(matches[1]-1)]
		lastchars = newtree_chars[(matches_end[1]+1):length(newtree_chars)]
		
		middle_chars = new_branch_length
		
		
		# Merge back together
		rateogram = paste(c(firstchars, middle_chars, lastchars), sep="", collapse="")
		
		#newtree <- gsub(brl[j], new_branch_length, newtree)
		}
	rateogram_tr <- read.tree(text = rateogram)
	orig_tr <- read.nexus(trfn)
	
	# Put the edge lengths into the original tree
	orig_tr$edge.length = rateogram_tr$edge.length
	rateogram_tr = orig_tr
	
	
	return(rateogram_tr)
	}










read_beasttree_states <- function (trfn, digits = NULL) 
	{
	# Scan the trfns in
    X <- scan(file = trfn, what = "", sep = "\n", quiet = TRUE)
    X_orig = X
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    lines_with_brackets = grep("\\[", X)
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats2(trfn)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)

	#interior <- which(!is.na(tab$posterior))
    
    # Right is the lines containing ]
    # -- basically the tree lines
	RIGHT <- grep("\\]", X)    
	
	
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT) > 0)
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
                # also remove from X_org
                X_orig <- X_orig[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
   	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]

	
	# Otherwise do this...
	#  If there isn't an endblock, then the last tree is the last line of interest
	#line_of_last_tree = 
	
	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]
	
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS trfns...)
    tree <- gsub("^.*= *", "", tree)
    

	#====================================
	# Get branch lengths in text trfn order    
	#====================================
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    ########################################
    
    ########################################
    # Make an empty vector, 1 cell for each stat
    nodestats <- vector(mode = "list", length = dim(tab)[2])
    ########################################

	return(nodestats)
	}









setup='
fn = trfn
fn = confn
regexp = "(tree)( )(STATE)(_)(\\d+)"
'
extractBEASTstats2 <- function (fn, regexp = "(tree)( )(STATE)(_)(\\d+)") 
	{
	
	# Scan in the NEXUS file
    X1 <- scan(file = fn, what = "", sep = "\n", quiet = TRUE)
    
    # phyloch original search string (assumes consensus tree)
    #Y <- X[grep("tree [[:space:]]+=", X)]
	
	# use your desired regular expression for the beginning of each tree line
	# The default gets e.g. "tree STATE_2323092"
	#regexp = "(tree)( )(STATE)(_)(\\d+)"

    X2 <- X1[grep(regexp, X1)]
        
    #Y <- X[grep("tree", X)]
    
    
    tmpre = NULL
	ri = 0
	tmpre[[(ri=ri+1)]] = '(tree)'	# Word 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 1
	tmpre[[(ri=ri+1)]] = '(STATE)'	# Word 2
	tmpre[[(ri=ri+1)]] = '(_)'	# Any Single Character 1
	tmpre[[(ri=ri+1)]] = '(\\d+)'	# Integer Number 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 2
	tmpre[[(ri=ri+1)]] = '(\\[.*?\\])'	# Square Braces 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 3
	tmpre[[(ri=ri+1)]] = '(=)'	# Any Single Character 2
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 4
	tmpre[[(ri=ri+1)]] = '(\\[&R\\])'	# Square Braces 2
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 5
	
	tmpre_string = list2str(tmpre, spacer="")
	(tmpre_string)

    # Remove the header   
    X3 <- gsub(tmpre_string, "", X2)
    
	# Extract the tree header
	treeheader_txt = extract_regexp(tmpstr=X2, tmpre_string)
    
    # Parse the tree header
    words = strsplit(treeheader_txt, split=" ")[[1]]
    
    # get the STATE
    tmpstate = gsub(pattern="STATE_", replacement="", words[2])
    
    # Extract everything within square braces
	regexp_braces = "(\\[.*?\\])"
    brackets_strings = extract_regexp(tmpstr=X3, tmpre_string=regexp_braces)
    names(brackets_strings) = NULL
    brackets_strings = unlist(brackets_strings)
	smm(brackets_strings)
	
	# Get the strings after each "["    
    #tab <- unlist(strsplit(X3, "\\["))[-1]
    brackets_strings1 = gsub("\\[&rate=", "", brackets_strings)
    brackets_strings2 = gsub("\\]", "", brackets_strings1)
    
    
    tab = as.numeric(brackets_strings2)
    
    
    
    # PARSE BRACKETS STRINGS & THIS WILL WORK
    
    
    
    tab = matrix(tab, ncol=1)
    tab = adf(tab)
    names(tab) = c("rates")
    
    #tab2 <- gsub("&|;|\\]", "", tab)
    #tab <- gsub(":.+$", "", tab)
    #foo <- function(x)
    #	{
    #    x <- unlist(strsplit(x, ","))
    #    x
   	# 	}
    #tab <- lapply(tab, foo)
    #for (i in seq(along = tab))
    #	{
    #    ind <- grep("[{]", tab[[i]])
    #    names <- gsub("=.+$", "", tab[[i]][ind])
    #    tab[[i]][ind] <- gsub("[{]", "", tab[[i]][ind])
    #    tab[[i]][ind] <- gsub("=", "_MIN=", tab[[i]][ind])
    #    tab[[i]][ind + 1] <- gsub("[}]", "", tab[[i]][ind + 1])
    #    tab[[i]][ind + 1] <- paste(paste(names, "MAX=", sep = "_"), tab[[i]][ind + 1])
	#    }
    #ttab <- data.frame()
    #stats <- unique(gsub("=.+$", "", unlist(tab)))
    #for (i in seq(along = tab)) 
    #	{
     #   for (j in seq(along = stats))
      #  	{
       #     ind <- grep(paste("^", stats[j], "=", sep = ""), 
        #        tab[[i]])
         #   if (length(ind) > 0)
          #  	{
           #     v <- as.numeric(gsub(paste(stats[j], "=", sep = ""), 
            #      "", tab[[i]][ind]))
             #   ttab[i, j] <- v
            	#}
   #     	}
    #	}
    #colnames(ttab) <- stats
    #tip <- which(is.na(ttab$posterior))
    return(tab)
	}


# The getAllSubTrees function below is a necessary subfunction that atomizes a
# tree into each individual subclade and was provided compliments of Luke Harmon.
#
# http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
getAllSubtrees <- function(tr, minSize=2)
	{
	require(geiger)
	require(ape)

	# Make empty list
	res = list()
	
	# Start the count
	count = 1
	
	# Get the number of tips
	ntip<-length(phy$tip.label)
	
	# For each internal node
	for (i in 1:phy$Nnode)
		{
		# Get leaves
		l = node.leaves(phy, ntip+i)
		bt = match(phy$tip.label, l)
		
		if (sum(is.na(bt)) == 0)
			{
			st = phy
			}
		else
			{
			st = drop.tip(phy, phy$tip.label[is.na(bt)])
			}
			
		if(length(st$tip.label) >= minSize)
			{
			res[[count]] = st
			count = count+1
			}
		}
	return(res)
	}











#######################################################
# Some functions for analyzing BEAST file outputs
#######################################################

# read_beast_prt() produces these columns in the prt table:
# names(beast_nodestats_dtf) = c("internal_nodenums_from1", "internal_nodenums", "rate_range_MIN", "rate_range_MAX", "height_95%_HPD_MIN", "height_95%_HPD_MAX", "height_median", "rate", "height", "rate_median", "height_range_MIN", "height_range_MAX", "rate_95%_HPD_MIN", "rate_95%_HPD_MAX", "posterior", "height_HPD_width", "height_range", "rate_HPD_width", "height_CV", "rate_CV")

# Plot chart of corresponding nodes, giving number of dataset 1, dataset 2, and column

plot_tr_beastcon_pairs <- function(tr_beastcon1=NULL, tr_beastcon2=NULL, dataset1=1, dataset2=3, coltxt, fns, analysis_names, analysis_colors, runspchs, test_one_to_one_slope=FALSE, linktxt=" - ", use_coltxt_label=TRUE, labelline=0, axissize=0.8, legend_x=0.05, legend_y=0.95, printflag=FALSE)
	{
	defaults='
	dataset1=1
	dataset2=3
	coltxt="posterior"
	'
	
	# Load the consensus trees and extract the data, if you haven't already
	#tr1 = read.beast_original(file=confn1)
	#tr_table1 = extractBEASTstats_orig(file=confn1)
	if (is.null(tr_beastcon1))
		{
		confn1 = fns[dataset1]
		tr_beastcon1 = read_beast_prt(file=confn1, get_tipnames=TRUE, printflag=printflag)
		}

	#tr2 = read.beast_original(file=confn2)
	#tr_table2 = extractBEASTstats_orig(file=confn2)
	if (is.null(tr_beastcon2))
		{
		confn2 = fns[dataset2]
		tr_beastcon2 = read_beast_prt(file=confn2, get_tipnames=TRUE, printflag=printflag)
		}
	
	# Get matching tipnames
	# Names from 1 found in 2
	names1 = tr_beastcon1$prt_beast_nodestats$tipnames
	names2 = tr_beastcon2$prt_beast_nodestats$tipnames
	
	corresp_rownums	= get_corresponding_rownums(names1, names2)
	rownums1 = corresp_rownums[,1]
	rownums2 = corresp_rownums[,2]
	
	# Node numbers only roughly correspond
	#plot(corresp_rownums[,1], corresp_rownums[,2])
	
	# Get the numbers of interest
	names(tr_beastcon1$prt_beast_nodestats)
	
	cmdstr1 = paste("xvals = tr_beastcon1$prt_beast_nodestats$'", coltxt, "'[rownums1]", sep="")
	eval(parse(text=cmdstr1))
	cmdstr2 = paste("yvals = tr_beastcon2$prt_beast_nodestats$'", coltxt, "'[rownums2]", sep="")
	eval(parse(text=cmdstr2))

	xmin = 0
	xmax = max(xvals, na.rm=TRUE)
	ymin = 0
	ymax = max(yvals, na.rm=TRUE)
	maxboth = max(xmax, ymax)


	if (use_coltxt_label == TRUE)
		{
		# xlabel = paste(coltxt, ", ", analysis_names[dataset1], sep="")
		# ylabel = paste(coltxt, ", ", analysis_names[dataset2], sep="")
		xlabel = paste(analysis_names[dataset1], linktxt, coltxt, sep="")
		ylabel = paste(analysis_names[dataset2], linktxt, coltxt, sep="")
		} else {
		xlabel = paste(analysis_names[dataset1], sep="")
		ylabel = paste(analysis_names[dataset2], sep="")
		}
	plot(xvals, yvals, xlim=c(xmin,maxboth), ylim=c(ymin,maxboth), xlab="", ylab="")
	# Indent xlabel and ylabel
	# xlabel
	mtext(text=xlabel, side=1, line=labelline, cex=axissize)
	# ylabel
	mtext(text=ylabel, side=2, line=labelline, cex=axissize)
	
	
	#title()
	# Plot a 1:1 line
	lines(x=c(xmin,maxboth), y=c(ymin,maxboth), lty="dotted")
	
	if (test_one_to_one_slope == TRUE)
		{
		# Take out the 1:1 trend, we only care about deviations from the 1:1 trend.
		yvals = yvals - xvals
		
		# Make a linear model
		xydata = adf(cbind(xvals, yvals))
		lm_result = lm(formula=yvals~xvals, data=xydata)
		
		
		# 
		R2 = summary(lm_result)$r.squared
		slope = summary(lm_result)$coefficients[2,1]
		intercept = summary(lm_result)$coefficients[1,1]

		slopeSE = summary(lm_result)$coefficients[2,2]
		slope95 = 1.96*slopeSE
		slope_pval = summary(lm_result)$coefficients[2,4]
		intercept_pval = summary(lm_result)$coefficients[1,4]
		
		
		R2txt = paste("R2=", format(R2, digits=3), sep="")
		slopetxt = paste("m=", format(slope, digits=3), "; int=", format(intercept, digits=3), sep="")#, " +/- ", format(slope95, digits=3), sep="")
		slp_pvaltxt = paste("slpP=", format(slope_pval, digits=3), sep="")
		int_pvaltxt = paste("intP=", format(intercept_pval, digits=3), sep="")

		txt_to_plot = paste("Dev. from 1:1:", slopetxt, slp_pvaltxt, int_pvaltxt, sep="\n")
		
		# http://stackoverflow.com/questions/3761410/how-can-i-plot-my-r-squared-value-on-my-scatterplot-using-r
		# bty suppresses box
		legend_x = legend_x * xmax
		legend_y = legend_y * ymax
		
		# print(legend_x)
		# print(legend_y)
		
		legend(x=legend_x, y=legend_y, bty="n", legend=txt_to_plot, cex=axissize)

		
		return(lm_result)
		}
	
	return(cbind(xvals, yvals))
	}


check_for_difference_from_one_to_one_slope <- function(tr_beastcon1=NULL, tr_beastcon2=NULL, dataset1, dataset2, coltxt, fns, analysis_names, analysis_colors, runspchs)
	{

	# Load the consensus trees and extract the data, if you haven't already
	#tr1 = read.beast_original(file=confn1)
	#tr_table1 = extractBEASTstats_orig(file=confn1)
	if (is.null(tr_beastcon1))
		{
		confn1 = fns[dataset1]
		tr_beastcon1 = read_beast_prt(file=confn1, get_tipnames=TRUE)
		}

	#tr2 = read.beast_original(file=confn2)
	#tr_table2 = extractBEASTstats_orig(file=confn2)
	if (is.null(tr_beastcon2))
		{
		confn2 = fns[dataset2]
		tr_beastcon2 = read_beast_prt(file=confn2, get_tipnames=TRUE)
		}
	
	# Get matching tipnames
	# Names from 1 found in 2
	names1 = tr_beastcon1$prt_beast_nodestats$tipnames
	names2 = tr_beastcon2$prt_beast_nodestats$tipnames
	
	corresp_rownums	= get_corresponding_rownums(names1, names2)
	rownums1 = corresp_rownums[,1]
	rownums2 = corresp_rownums[,2]
	
	# Node numbers only roughly correspond
	#plot(corresp_rownums[,1], corresp_rownums[,2])
	
	# Get the numbers of interest
	names(tr_beastcon1$prt_beast_nodestats)
	
	cmdstr1 = paste("xvals = tr_beastcon1$prt_beast_nodestats$'", coltxt, "'[rownums1]", sep="")
	eval(parse(text=cmdstr1))
	cmdstr2 = paste("yvals = tr_beastcon2$prt_beast_nodestats$'", coltxt, "'[rownums2]", sep="")
	eval(parse(text=cmdstr2))

	xmin = 0
	xmax = max(xvals, na.rm=TRUE)
	ymin = 0
	ymax = max(yvals, na.rm=TRUE)
	maxboth = max(xmax, ymax)

	# xlabel = paste(coltxt, ", ", analysis_names[dataset1], sep="")
	# ylabel = paste(coltxt, ", ", analysis_names[dataset2], sep="")
	# xlabel = paste(analysis_names[dataset1], " - ", coltxt, sep="")
	# ylabel = paste(analysis_names[dataset2], " - ", coltxt, sep="")
	# plot(xvals, yvals, xlim=c(xmin,maxboth), ylim=c(ymin,maxboth), xlab=xlabel, ylab=ylabel)
	
	#title()
	# Plot a 1:1 line
	# lines(x=c(xmin,maxboth), y=c(ymin,maxboth), lty="dotted")

	#######################################
	# Test for one-to-one slope	
	#######################################
	# Take out the 1:1 trend, we only care about deviations from the 1:1 trend.
	yvals = yvals - xvals
	
	# Make a linear model
	xydata = adf(cbind(xvals, yvals))
	lm_result = lm(formula=yvals~xvals, data=xydata)
	
	
	# 
	R2 = summary(lm_result)$r.squared
	slope = summary(lm_result)$coefficients[2,1]
	slopeSE = summary(lm_result)$coefficients[2,2]
	slope95 = 1.96*slopeSE
	pval = summary(lm_result)$coefficients[2,4]
	
	
	R2txt = paste("R2=", format(R2, digits=3), sep="")
	slopetxt = paste("m=", format(slope, digits=3), " +/- ", format(slope95, digits=3), sep="")
	pvaltxt = paste("p=", format(pval, digits=3), sep="")

	txt_to_plot = paste(R2txt, slopetxt, pvaltxt, sep="\n")
		
	return(lm_result)
	}
























# Tree-to-tree distances
# (returns branch-score differences, i.e. treedist[2]; for apply to a list of trees)
# BUT: JUST USE *mapply* INSTEAD...
t2t_dist <- function(k, subset_of_trees1, subset_of_trees2)
	{
	tmpdist = treedist( subset_of_trees1[[k]], subset_of_trees2[[k]])[2]
	return(tmpdist)
	}

# defaults='UPtr = read.tree(((((((Syn7803:0.007398,Syn7805:0.009548):0.02046,Syn9311:0.044607):0.005999,(Syn9917:0.032992,Syn9916:0.025461):0.009812):0.013396,(((Pro9303:0.004518,Pro9313:0.004096):0.045163,((Pro9211:0.055912,Pro1375:0.061717):0.023249,((ProNAT2:0.002222,ProNAT1:0.002492):0.073484,((Pro9515:0.01412,ProMED4:0.018279):0.024564,(((Pro9301:0.006455,Pro9601:0.004998):0.004017,Pro9215:0.011527):0.004599,Pro9312:0.012784):0.024095):0.117301):0.029427):0.037262):0.016632,(Syn5701:0.065547,(((Syn7002:0.15763,(SynJ23B:0.022676,SynJ33A:0.027358):0.211298):0.059181,(Syn6301:0.001198,Syn7942:0.000226):0.107629):0.163424,SynR307:0.064874):0.02067):0.032984):0.016042):0.020577,(Syn8102:0.031365,Syn9605:0.024734):0.005393):0.030896,SynB107:0.004734,Syn9902:0.008388);)'


gammadist_equiv_trees <- function(symdist, MASTdist, nsplits, ntaxa)
	{
	# This calculates Ge's gamma for trees with the same number of splits and taxa
	#
	# Calculate gamma(tree1, tree2), from Ge, F., L.-S. Wang and J. Kim (2005). 
	# "The Cobweb of Life Revealed by Genome-Scale Estimates of 
	# Horizontal Gene Transfer." PLoS Biol 3(10): e316.
	# http://dx.doi.org/10.1371/journal.pbio.0030316 
	#
	# m = # splits in tree #1 (here, 28)
	# n = # splits in tree #2 (here, 28)
	# x = # of taxa (here, 27)
	#
	# gammastat = ((symdist - abs(m-n)) / (2* min(c(m, n)))) - ( MASTdist / (x-3) )
	
	x = ntaxa
	gammastat = ((symdist) / (2* nsplits)) - ( MASTdist / (x-3) )
	return(gammastat)
	}



gammadist_between_trees <- function(symdist, MASTdist, nsplits_tr1, nsplits_tr2, ntaxa)
	{
	# This calculates Ge's gamma for trees with a different number of splits (same # of taxa)
	#
	# Calculate gamma(tree1, tree2), from Ge, F., L.-S. Wang and J. Kim (2005). 
	# "The Cobweb of Life Revealed by Genome-Scale Estimates of 
	# Horizontal Gene Transfer." PLoS Biol 3(10): e316.
	# http://dx.doi.org/10.1371/journal.pbio.0030316 
	#
	# m = # splits in tree #1 (here, 28)
	# n = # splits in tree #2 (here, 28)
	# x = # of taxa (here, 27)
	#
	# gammastat = ((symdist - abs(m-n)) / (2* min(c(m, n)))) - ( MASTdist / (x-3) )
	m = nsplits_tr1
	n = nsplits_tr2
	x = ntaxa
	gammastat = ((symdist - abs(m-n)) / (2* min(c(m, n)))) - ( MASTdist / (x-3) )
	return(gammastat)
	}


defaults = '
gammastat_between = gammastats_between_ordered[190]
'
gammadist_pval <- function(gammastat_between, gammastats_within_ordered)
	{
	# Calculate average p-value
	sum_rank = sum(gammastat_between >= gammastats_within_ordered)
	gamma_pval = sum_rank / length(gammastats_within_ordered)
	return(gamma_pval)
	}


min_phyd <- function(phyd)
	{
	###############################################################################
	# min_phyd: find the minimums of a phylogenetic distance matrix, for each tip
	###############################################################################
	# replace diagonal 0s with NAs
	# also replace the upper Triangle
	tmpphyd = phyd
	diag(tmpphyd) = NA
	upperTriangle(tmpphyd) = NA
	
	maxrow = nrow(phyd)
	
	# leave out the last column of tmpphyd, which is just NA
	mindists = colmin(tmpphyd[, 1:(ncol(tmpphyd)-1) ])
	
	list_of_min_pairs = c()
	for (i in 1:ncol(tmpphyd))
		{
		# avoid repeating pairs...with phyd[(i+1):maxrow, i]
		tmp_closest_tips = which(tmpphyd[, i] == mindists[i])
		
		tmplist = cbind(rep(i, length(tmp_closest_tips)), tmp_closest_tips)
		
		list_of_min_pairs = rbind(list_of_min_pairs, tmplist)
		}
	return(list_of_min_pairs)
	}


#######################################################
# fix Psychotria tree so tips all come to 0
#######################################################

#######################################################
# make_all_tips_come_to_0
#######################################################
#' Take a tree, ensure all tips end at 0
#' 
#' Makes tree precisely ultrametric (unlike e.g. the default <i>Psychotria</i> tree).
#'
#' This function ADDS the time_before_present
#' 
#' @param trfn The Newick tree filename
#' @param outfn The filename for the resulting file
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be treated as fossils and not averaged into the tips.
#' This is currently set to 0, on the assumption that you want to make ALL tips come to zero, including fossils.
#' @return \code{phy} The corrected phylogeny
#' @export
#' @seealso \code{\link[ape]{read.tree}}, \code{\link{prt}}, \code{\link{average_tr_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
make_all_tips_come_to_0 <- function(trfn, outfn, fossils_older_than=0.0)
	{
	phy = read.tree(trfn)


	print("Input tree:")
	print("is.ultrametric(phy)")
	print(is.ultrametric(phy))
	
	prt_table = prt(phy, fossils_older_than=fossils_older_than)
	
	
	prt_table$time_bp
	
	for (i in 1:nrow(prt_table))
		{
		if ((prt_table$time_bp[i] > 0) && (prt_table$node.type[i] == "tip") && (prt_table$fossils[i]==FALSE))
			{
			edgeTF = phy$edge[,2] == i
			phy$edge.length[edgeTF] = phy$edge.length[edgeTF] + prt_table$time_bp[i]
			}
		}
	
	print("Output tree:")
	print("is.ultrametric(phy)")
	print(is.ultrametric(phy))
	
	prt(phy)
	
	#write.tree(phy, file="Psychotria_5.2_ultrametric.newick")
	write.tree(phy, file=outfn)
	return(phy)
	}










extend_tips_to_ultrametricize <- function(obj, age_of_root, tips_end_at_this_date=NA)
	{
	
	#t3 = 
	
	print("node ages of tips:")
	tip_ages = age_of_root + get_node_ages_of_tips(obj)
	#print(tip_ages)
	
	
	if (is.na(tips_end_at_this_date))
		{
		tips_end_at_this_date = max(tip_ages)
		}
	
	nums_to_add_to_tip_to_ultrametricize = tips_end_at_this_date - tip_ages
	
	indices_of_branches_under_tips = get_indices_of_branches_under_tips(obj)

	obj$edge.length[indices_of_branches_under_tips] = obj$edge.length[indices_of_branches_under_tips] + nums_to_add_to_tip_to_ultrametricize
	
	return(obj)
	}



edges_existing_at_correct_time_bp_TF <- function(time_slice, edge_times_bp, roundto=5)
	{
	# find the edges that exist in the right time
	
	# (note: round to a default of 6 digits (a single year) with roundto; this is 
	#  important for whether or not lineages exist at time=0 before present)

	# timepoint is younger or equal to the oldest end of the branch
	edges_that_start_below_time = round(edge_times_bp[, 1], digits=roundto) > time_slice
	
	# timepoint is older than the youngest end of the branch	
	edges_that_end_after_time = round(edge_times_bp[, 2], digits=roundto) <= time_slice
	edges_that_exist_in_the_right_time = edges_that_start_below_time + edges_that_end_after_time == 2
	return(edges_that_exist_in_the_right_time)
	}



################################################
# Tree drawing
################################################
drawtree_branches_heatmap <- function(tr, nums_for_edges, titletxt="", rescaled_colors="", tmp_br_wid=8)
	{
	# heat map
	# scale to 0-1, times 130, round down, +1, = scale from 20 to 120
	tmp_colors = rev(heat.colors(130, alpha=1))
	
	if (rescaled_colors == "")
		{
		rescaled_colors = 15+floor( 99*((nums_for_edges-min(nums_for_edges)) / (max(nums_for_edges)-min(nums_for_edges))) )
		}
	br.col = tmp_colors[rescaled_colors]
	plot(tr, edge.color=br.col, edge.width=tmp_br_wid, label.offset=0.05)
	
	# Legend
	value_range = pretty(nums_for_edges)
	value_range[1] = round(min(nums_for_edges), digits=2)
	value_range[length(value_range)] = round(max(nums_for_edges), digits=2)

	value_range_rescaled = 15+floor( 99*((value_range-min(nums_for_edges)) / (max(nums_for_edges)-min(nums_for_edges))) )
	
	color_range = tmp_colors[value_range_rescaled]
	legend(x="bottomleft", legend=value_range, bty="0", fill=color_range, y.intersp=0.75, pt.cex=2, cex=1, xjust=0, yjust=0, title=titletxt)

	}
	
drawtree_branches_bluered <- function(tr, nums_for_edges, titletxt="", rescaled_colors="", tmp_br_wid=8)
	{
	# blue-red
	# scale to 0-1, times 100, round down, +1, = scale from 1 to 100
	tmp_colors = rev(rainbow(110, start=0, end=4.5/6, alpha=1))
	if (rescaled_colors == "")
		{
		rescaled_colors = 5+floor( 99*((nums_for_edges-min(nums_for_edges)) / (max(nums_for_edges)-min(nums_for_edges))) )
		}
	br.col = tmp_colors[rescaled_colors]
	plot(tr, edge.color=br.col, edge.width=tmp_br_wid, label.offset=0.05)

	# Legend
	value_range = pretty(nums_for_edges)
	value_range[1] = round(min(nums_for_edges), digits=2)
	value_range[length(value_range)] = round(max(nums_for_edges), digits=2)

	value_range_rescaled = 5+floor( 99*((value_range-min(nums_for_edges)) / (max(nums_for_edges)-min(nums_for_edges))) )
	
	color_range = tmp_colors[value_range_rescaled]
	legend(x="bottomleft", legend=value_range, bty="o", fill=color_range, y.intersp=0.75, pt.cex=2, cex=1, xjust=0, yjust=0, title=titletxt)

	}




#############################################################################
# LTT plots (lineages-thru-time plots)
#############################################################################
LTT_fossil <- function(times, tr)
	{
	# get lineages-thru-time for a specified series of time points before the 
	# present, by looking at the actual number of lineages at the time units,
	# not by counting bifurcations (as the standard LTT function does)

	# get a list of the edge start/stops in the phylogeny	
	edge_times_bp = get_edge_times_before_present(tr)
	
	# set up the output array
	lineages_count = rep(NA, length(times))

	for (i in 1:length(times))
		{
		time_slice = times[i]
		branches_existing = edges_existing_at_correct_time_bp_TF(time_slice, edge_times_bp)

		# count the number of lineages crossing that time point
		# according to this particular phylogeny
		lineages_count[i] = sum(branches_existing)
		}
	return(lineages_count)
	}


LTTplot_fossil <- function(times, counts)
	{
	plot(times, counts, pch="", xlab="Time (mya)", ylab="# of lineages crossing time point", xlim=rev(range(times)))
	title("Lineage-through-time plot including fossils")

	add_LTT_line_fossil(times, counts)
	
	}

add_LTT_line_fossil <- function(times, counts, tmpcol="blue")
	{
	# Lines of average counts & 95% confidence intervals
	lines(times, counts, type="l", lty=1, lwd=2, col=tmpcol, xpd=FALSE)
	
	}


################################################################################
# TREE MODIFICATION FUNCTIONS (e.g. adding hooks, choosing certain branches)
################################################################################
add_hooks <- function(tr, list_of_times_before_present, brlen_of_side_branch=0.0000001, plottree=FALSE)
	{
	# Take a list of ages, add hooks to any branch existing at that age

	# OK, RE-FUCKING DO
	# Gather a list of tip labels, and then record the distance below those tips that
	# you would go down to attach a hook
	hooktree = tr
	list_of_daughter_tipnames_to_add_hooks_below = c()
	list_of_ages_below_daughter = c()   # assumes ultrametric
	ntips = length(hooktree$tip.label)
	for (i in 1:length(list_of_times_before_present))
		{
		# Get the edges that exist at the time_slice in question
		time_slice = as.numeric(list_of_times_before_present[i])
		edge_times_bp = get_edge_times_before_present(hooktree)
		edges_that_exist_in_the_right_time = edges_existing_at_correct_time_bp_TF(time_slice, edge_times_bp)
		
		# get the nodes daughter to the branches that match
		nodenums_to_add_hooks_to = hooktree$edge[,2][edges_that_exist_in_the_right_time]

		# calculate the times parent to these daughters at which to insert the hooks		
		#times_before_daughter_nodes = time_slice - edge_times_bp[edges_that_exist_in_the_right_time, 2]
		
		# trace these nodes to their tips in the (UNALTERED ORIGINAL) tree
		for (j in 1:length(nodenums_to_add_hooks_to))
			{
			nodenum = nodenums_to_add_hooks_to[j]
			
			# If the node is a tip
			if (nodenum <= ntips)
				{
				daughter_tipname_to_add_hooks_below = hooktree$tip.label[nodenum]
				
				edgenums = 1:nrow(edge_times_bp)
				edgenum = edgenums[tr$edge[,2] == nodenum]
				node_time_of_daughter = edge_times_bp[edgenum, 2]
				
				time_before_daughter_nodes = time_slice - node_time_of_daughter
				
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, daughter_tipname_to_add_hooks_below)
				list_of_ages_below_daughter = c(list_of_ages_below_daughter, time_before_daughter_nodes)			
				} else {
				# Get *a* tip to compare time slice to
				temp_tipnames = get_all_daughter_tips_of_a_node(nodenum, hooktree)
				temp_tipname = temp_tipnames[1]
				namenums = 1:length(hooktree$tip.label)
				temp_tipnum = namenums[hooktree$tip.label == temp_tipname]
				
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, hooktree$tip.label[temp_tipnum])

				#
				edgenums = 1:nrow(edge_times_bp)
				edgenum = edgenums[tr$edge[,2] == temp_tipnum]
				node_time_of_daughter = edge_times_bp[edgenum, 2]
				
				time_before_daughter_nodes = time_slice - node_time_of_daughter

				list_of_ages_below_daughter = c(list_of_ages_below_daughter, time_before_daughter_nodes)
				}
			print(list_of_ages_below_daughter)
			}
		}
		
	
	# Now, attach the hooks
	for (i in 1:length(list_of_daughter_tipnames_to_add_hooks_below))
		{
		#print(paste("i=", i, sep=""))
		tipname = list_of_daughter_tipnames_to_add_hooks_below[i]
		depthtime = as.numeric(list_of_ages_below_daughter[i])
		
		print(depthtime)
		print(list_of_ages_below_daughter)
		
		hooktree = add_hook(hooktree, tipname, depthtime, plottree=plottree)
		}
	
	return(hooktree)
	}


add_hooks_old <- function(tr, list_of_times_before_present, brlen_of_side_branch=0.0000001, plottree=FALSE)
	{
	# Take a list of ages, add hooks to any branch existing at that age

	# OK, RE-FUCKING DO
	# Gather a list of tip labels, and then record the distance below those tips that
	# you would go down to attach a hook
	hooktree = tr
	list_of_daughter_tipnames_to_add_hooks_below = c()
	list_of_ages_below_daughter = c()   # assumes ultrametric
	ntips = length(hooktree$tip.label)
	for (i in 1:length(list_of_times_before_present))
		{
		# Get the edges that exist at the time_slice in question
		time_slice = as.numeric(list_of_times_before_present[i])
		edge_times_bp = get_edge_times_before_present(hooktree)
		edges_that_exist_in_the_right_time = edges_existing_at_correct_time_bp_TF(time_slice, edge_times_bp)
		
		# get the nodes daughter to the branches that match
		nodenums_to_add_hooks_to = hooktree$edge[,2][edges_that_exist_in_the_right_time]

		# calculate the times parent to these daughters at which to insert the hooks		
		times_before_daughter_nodes = time_slice - edge_times_bp[edges_that_exist_in_the_right_time, 2]
		
		# trace these nodes to their tips in the (UNALTERED ORIGINAL) tree
		for (j in 1:length(nodenums_to_add_hooks_to))
			{
			nodenum = nodenums_to_add_hooks_to[j]
			if (nodenum <= ntips)
				{
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, hooktree$tip.label[nodenum])
				list_of_ages_below_daughter = c(list_of_ages_below_daughter, times_before_daughter_nodes[j])			
				}
			else
				{
				temp_tips = get_all_daughter_tips_of_a_node(nodenum, hooktree)
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, temp_tips[1])
				list_of_ages_below_daughter = c(list_of_ages_below_daughter, times_before_daughter_nodes[j])
				}
			}
		}
		
	
	# Now, attach the hooks
	for (i in 1:length(list_of_daughter_tipnames_to_add_hooks_below))
		{
		#print(paste("i=", i, sep=""))
		tipname = list_of_daughter_tipnames_to_add_hooks_below[i]
		depthtime = as.numeric(list_of_ages_below_daughter[i])
		
		hooktree = add_hook(hooktree, tipname, depthtime, plottree=plottree)
		}
	
	return(hooktree)
	}


add_hook2 <- function(t, list_of_times_before_present, tipname_to_grepl, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=FALSE)
	{
	# Take a list of ages, add hooks to any branch existing at that age

	# OK, RE-FUCKING DO
	# Gather a list of tip labels, and then record the distance below those tips that
	# you would go down to attach a hook
	hooktree = t
	list_of_daughter_tipnames_to_add_hooks_below = c()
	list_of_ages_below_daughter = c()   # assumes ultrametric
	ntips = length(hooktree$tip.label)
	for (i in 1:length(list_of_times_before_present))
		{
		# Get the edges that exist at the time_slice in question
		time_slice = as.numeric(list_of_times_before_present[i])
		edge_times_bp = get_edge_times_before_present(hooktree)
		edges_that_exist_in_the_right_time = edges_existing_at_correct_time_bp_TF(time_slice, edge_times_bp)
		
		# get the nodes daughter to the branches that match
		nodenums_to_add_hooks_to = hooktree$edge[,2][edges_that_exist_in_the_right_time]

		# calculate the times parent to these daughters at which to insert the hooks		
		times_before_daughter_nodes = time_slice - edge_times_bp[edges_that_exist_in_the_right_time, 2]
		
		# trace these nodes to their tips in the (UNALTERED ORIGINAL) tree
		for (j in 1:length(nodenums_to_add_hooks_to))
			{
			nodenum = nodenums_to_add_hooks_to[j]
			if (nodenum <= ntips)
				{
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, hooktree$tip.label[nodenum])
				list_of_ages_below_daughter = c(list_of_ages_below_daughter, times_before_daughter_nodes[j])			
				} else
				{
				temp_tips = get_all_daughter_tips_of_a_node(nodenum, hooktree)
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, temp_tips[1])
				list_of_ages_below_daughter = c(list_of_ages_below_daughter, times_before_daughter_nodes[j])
				}
			}
		}
		
	# pick one daughter node
	tmptable = cbind(list_of_daughter_tipnames_to_add_hooks_below, list_of_ages_below_daughter)
	#tmptable = tmptable[order(list_of_ages_below_daughter), ]
	index_of_closest_to_present_branch = match(TRUE, grepl(tipname_to_grepl, tmptable[,1]))
	
	# pull out just this one branch; if no suitable branch found, ignore it and continue
	if (is.na(index_of_closest_to_present_branch))
		{
		return(hooktree)
		} else {
		list_of_daughter_tipnames_to_add_hooks_below = tmptable[index_of_closest_to_present_branch, 1]
		list_of_ages_below_daughter = tmptable[index_of_closest_to_present_branch, 2]
		}
	
	# Now, attach the hooks
	for (i in 1:length(list_of_daughter_tipnames_to_add_hooks_below))
		{
		if (printflag)
			{
			print(paste("i=", i, sep=""))
			}
		tipname = list_of_daughter_tipnames_to_add_hooks_below[i]
		depthtime = as.numeric(list_of_ages_below_daughter[i])
		
		hooktree = add_hook(hooktree, tipname, depthtime, plottree=plottree)
		}
	
	return(hooktree)
	}


## Moved to BioGeoBEARS_add_fossils_randomly_v1.R
# 
# add_hook <- function(t, tipname, depthtime, brlen_of_side_branch=0.0000001, plottree = FALSE, printflag=FALSE, newtipname="default")
# 	{
# 	# Add a hook (a small side tip) to a phylogeny
# 	#
# 	# e.g.:
# 	# Do spatial variogram by doing points from many different species
# 	# add tips to tree
# 	#cat("owls(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);", file = "ex.tre", sep = "\n")
# 	#t <- read.tree("ex.tre")
# 	#prt(t)
# 	#newtree = add_hook(t, brlen_of_side_branch=0.0000001, plottree = TRUE)
# 
# 	if (printflag)
# 		{
# 		cat("add_hook(): Adding below tipname ",  tipname, "\n", sep="")
# 		}
# 
# 	newtree = t
# 	#daughter_nodenum_to_add_hook_below = 4
# 	#height_below_daughter_at_which_to_add_it = 1.0
# 	#brlen_of_side_branch=0.0000001
# 	
# 	# find the node, below this you will add the 
# 	tip_nodenum = which(t$tip.label == tipname)
# 	
# 	if (printflag)
# 		{
# 		print(paste("addhook(): tip_nodenum = ", tip_nodenum, sep=""))
# 		}
# 
# 	daughter_nodenum_to_add_hook_below = trace_parents_up(tip_nodenum, t, depthtime)
# 	
# 	# Error trap in case of e.g. overrun of bottom of tree
# 	if (is.na(daughter_nodenum_to_add_hook_below))
# 		{
# 		cat("add_hook() error: daughter_nodenum_to_add_hook_below is NA, probably overshot bottom of tree (?); not adding hook\n")
# 		return(t)
# 		}
# 	
# 	if (printflag)
# 		{
# 		print(paste("addhook(): internal nodenum to add below = ", daughter_nodenum_to_add_hook_below, sep=""))
# 		}
# 	
# 	tip_to_ancestor_node_dist = dist_between_direct_ancestors(daughter_nodenum_to_add_hook_below, tip_nodenum, t, totaldist=0)
# 	
# 	
# 	
# 	height_below_daughter_at_which_to_add_it = depthtime - tip_to_ancestor_node_dist
# 	
# 	
# 	# add a new tip to the list of tips (this is the hook)
# 	new_tip_nodenum = get_nodenum_structural_root(t)
# 	
# 	# bump up all of the nodenums above the node tip by 1
# 	newtree$edge[t$edge >= new_tip_nodenum] = t$edge[t$edge >= new_tip_nodenum] + 1
# 	
# 	# add a new internal node at the end
# 	new_inNode = max(newtree$edge) + 1
# 	
# 	# add two new edges, and replace the old edge
# 	#print(t$edge[,2])
# 	#print(daughter_nodenum_to_add_hook_below)
# 	#print(t$edge[,2] == daughter_nodenum_to_add_hook_below)
# 	old_edge_num = which(t$edge[,2] == daughter_nodenum_to_add_hook_below)
# 	
# 	# extract the edgenums before and after this insertion (exceptions for in case
# 	# if the modified row is the first or last row)
# 	if (old_edge_num == 1)
# 		{
# 		first_old_edges_rownums = NULL
# 		} else {
# 		first_old_edges_rownums = 1:(old_edge_num-1) #newtree$edge[1:(old_edge_num-1), ]
# 		}
# 	if (old_edge_num == nrow(t$edge))
# 		{
# 		second_old_edges_rownums = NULL
# 		} else {
# 		second_old_edges_rownums = (old_edge_num+1):nrow(t$edge) # newtree$edge[, ]
# 		}
# 	
# 	
# 	
# 	# replace the edge, keeping the old parent (which may have increased by 1! use newtree!!), put the new internal node as the daughter)
# 	replacement_edge_row = newtree$edge[old_edge_num, ]
# 	replacement_edge_row[2] = new_inNode
# 	
# 	# subtract the distance below the daughter, from the top
# 	replacement_edge_length = t$edge.length[old_edge_num] - height_below_daughter_at_which_to_add_it
# 	
# 	
# 	# make the new edge, which goes below the old daughter node
# 	# you have to bump the daughter_nodenum_to_add_hook_below if it is
# 	# >= to the new_tip_nodenum
# 	if (daughter_nodenum_to_add_hook_below >= new_tip_nodenum)
# 		{
# 		daughter_nodenum_to_add_hook_below = daughter_nodenum_to_add_hook_below + 1
# 		}
# 	new_edge_below_old_daughter_node = c(new_inNode, daughter_nodenum_to_add_hook_below)
# 	new_edge_below_old_daughter_node_edge_length = height_below_daughter_at_which_to_add_it
# 	
# 	# make the new edge, which goes below the new tip: c(parent, daughter)
# 	new_edge_below_new_tip = c(new_inNode, new_tip_nodenum)
# 	new_edge_below_new_tip_edge_length = brlen_of_side_branch
# 	
# 	
# 	# add the edge rows before the one that is replaced, then the replaced edge, then the other old edges, then the 2 new edges
# 	new_edge_table = rbind(newtree$edge[first_old_edges_rownums, ], replacement_edge_row, newtree$edge[second_old_edges_rownums, ], new_edge_below_old_daughter_node, new_edge_below_new_tip)
# 	
# 	new_edgelength_list = c(t$edge.length[first_old_edges_rownums], replacement_edge_length, t$edge.length[second_old_edges_rownums], new_edge_below_old_daughter_node_edge_length, new_edge_below_new_tip_edge_length)
# 	
# 	# it MAY be important that the node numbers be INTEGER, not NUMERIC
# 	newtree$edge = matrix(as.integer(new_edge_table), ncol=2, byrow=FALSE)
# 	
# 	#row.names(newtree$edge) = NULL
# 	#newtree$edge[,1] = as.integer(newtree$edge[,1])
# 	#newtree$edge[,2] = as.integer(newtree$edge[,2])
# 	
# 	newtree$edge.length = new_edgelength_list
# 	#row.names(newtree$edge.length) = NULL
# 	
# 	# update number of internal nodes
# 	newtree$Nnode = t$Nnode + 1
# 	
# 	# add the new tip to the end of the list of tips
# 	if (newtipname == "default")
# 		{
# 		newtipname = paste("hook", new_tip_nodenum, sep="")
# 		} else {
# 		newtipname = newtipname
# 		}
# 	newtree$tip.label = c(t$tip.label, newtipname)
# 	
# 	cat("Adding ",  newtipname, "\n", sep="")
# 	
# 	
# 	# some crap to fix the tree formatting somehow
# 	# I mean, really, the tree was fucking logically correct, but
# 	# hanging plot and dist.nodes and probably anything
# 	# using reorder(phy, "pruningwise"), but I couldn't figure out why
# 	# I guess the order of the tips is important for some reason?
# 	# like maybe leftmost tip is leftmost in branching?
# 	# wtf kind of data architecture is this?
# 	# anyway, FUCK IT, writing to Newick and reading back fixes it.
# 	newtree = reorder(newtree)
# 	tmpfn = "tmp_junktree.tree"
# 	write.tree(newtree, tmpfn)
# 	newtree2 = read.tree(tmpfn)
# 
# 	
# 	# plot, if desired:
# 	if (plottree == TRUE)
# 		{
# 		cat("add_hook(): plotting/printing the resulting tree...\n", sep="")
# 		prt(newtree)
# 		plot(newtree)
# 		}
# 	
# 	return(newtree2)
# 	}



# Add short extinct tips to phylogeny
junk = '
brlen_of_side_branch=0.0000001
# calc_loglike default min is 0.000001
'
add_EX_hooks <- function(tr, brlen_of_side_branch=0.0000001)
	{
	# Add short extinct tips to phylogeny;
	# These will be short branches added to any lineage that ends before the present (0 mya)
	# Note: add these BEFORE adding hooks representing fossils; those hooks don't go
	# extinct, they just represent the lineage at that timepoint.
	
	tr_table = prt(tr, printflag=FALSE)
	
	# hist(tr_table$time_bp)
	
	c1 = tr_table$time_bp > 0
	c2 = tr_table$node.type == "tip"
	tip_labels_to_make_EX_TF = ( (c1 + c2) == 2 )
	tip_labels_to_make_EX = tr_table$label[tip_labels_to_make_EX_TF]
	
	
	# But, it's not an extinction if the taxon is described as an "-anc" -- because it continued into a "-desc"!!
	tip_labels_to_make_EX = tip_labels_to_make_EX[ grepl("-anc", tip_labels_to_make_EX) == FALSE ]
	extinct_tip_ages = tr_table$time_bp[match(tip_labels_to_make_EX, tr_table$label)]
	
	# make sure they're a little older than the actual tip
	extinct_tip_ages = extinct_tip_ages + brlen_of_side_branch
	
	hooktree2 = tr
	
	# Add the hooks
	for (i in 1:length(tip_labels_to_make_EX))
		{
		tipname_to_grepl = tip_labels_to_make_EX[i]
		list_of_times_before_present = extinct_tip_ages[i]

		old_hooktree2 = hooktree2
		
		# print(paste(z, ": ", tipname_to_grepl, sep=""))
		hooktree2 = add_hook2(hooktree2, list_of_times_before_present, tipname_to_grepl)
		
		# identify new tip, if any:
		newtip = hooktree2$tip.label[hooktree2$tip.label %in% old_hooktree2$tip.label == FALSE]
		#hooktree$tip.label
		
		# print(newtip)
		# print(len(newtip))
		
		if (len(newtip) > 0 )
			{
			cat(tipname_to_grepl, ": yes, extinct tip added\n", sep="")
			#print("yes")
			
			# Update the tipname
			new_tipname = paste(tipname_to_grepl, "_EX", sep="")
			hooktree2$tip.label[hooktree2$tip.label == newtip] = new_tipname
			
			} else {
			# no new tips added; nothing to do
			cat(tipname_to_grepl, ": no new extinct tip added\n", sep="")		
			}
		}
	
	return(hooktree2)
	
	}



chainsaw <- function(tr, timepoint=10)
	{
	# Take a tree and saw it off evenly across a certain timepoint.
	# This removes any tips above the timepoint, and replaces them 
	# with a single tip representing the lineage crossing
	# the timepoint (with a new tip name).

	# Get the tree in a table
	tr_table = prt(tr, printflag=FALSE)
	
	# Find the tips that are less than 10 my old and drop them
	TF_exists_more_recently_than_10mya = tr_table$time_bp < timepoint
	
	# Get the corresponding labels
	labels_for_tips_existing_more_recently_than_10mya = tr_table$label[ TF_exists_more_recently_than_10mya == TRUE ]
	
	###########################################
	# Draft chainsaw function
	###########################################
	# loop through the branches that cross 10 mya
	
	# get a list of the edge start/stops in the phylogeny's edges
	edge_times_bp = get_edge_times_before_present(tr)
	
	# which of these branches cross 10 mya?
	edges_start_earlier_than_10mya = edge_times_bp[, 1] > timepoint
	edges_end_later_than_10mya = edge_times_bp[, 2] <= timepoint
	edges_to_chainsaw = edges_start_earlier_than_10mya + edges_end_later_than_10mya == 2
	
	# then, for each of these edges, figure out how many tips exist descending from it
	nodes_to_chainsaw = tr$edge[, 2][edges_to_chainsaw]
	# Take only internal nodes (? why ?)
	nodes_to_chainsaw = nodes_to_chainsaw[nodes_to_chainsaw > length(tr$tip.label)]
	
	# create a copy of the tree to chainsaw
	tree_to_chainsaw = tr
	
	for (i in 1:length(nodes_to_chainsaw))
		{
		tmp_subtree = extract.clade(tr, nodes_to_chainsaw[i])
		#print(tmp_subtree$tip.label)
		
		tmp_number_of_tips = length(tmp_subtree$tip.label)
		#print(tmp_number_of_tips)
		
		# number of tips to drop = (numtips -1)
		numtips_to_drop = tmp_number_of_tips - 1 
		
		# tips_to_drop
		tmp_labels = tmp_subtree$tip.label
		
		labels_to_drop = tmp_labels[1:numtips_to_drop]
		
		# new label
		label_kept_num = length(tmp_labels)
		label_kept = tmp_labels[label_kept_num]
		new_label = paste("CA_", label_kept, "+", numtips_to_drop, "_tips", sep="")
		tree_to_chainsaw$tip.label[tree_to_chainsaw$tip.label == label_kept] = new_label
		
		# chop off e.g. 2 of the 3 tips
		tree_to_chainsaw = drop.tip(tree_to_chainsaw, labels_to_drop)
		
		}
	#plot(tree_to_chainsaw)
	#axisPhylo()
	
	tree_to_chainsaw_table = prt(tree_to_chainsaw, printflag=FALSE)
	
	tree_to_chainsaw_table_tips_TF_time_bp_LT_10my = tree_to_chainsaw_table$time_bp < timepoint
	
	
	tmp_edge_lengths =  tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
	
	times_bp_for_edges_to_chainsaw = tree_to_chainsaw_table$time_bp[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
	
	adjustment = times_bp_for_edges_to_chainsaw - timepoint
	
	revised_tmp_edge_lengths = tmp_edge_lengths + adjustment
	
	tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my] = revised_tmp_edge_lengths
	
	# revised
	ordered_nodenames = get_nodenums(tree_to_chainsaw)
	parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, tree_to_chainsaw$edge[,2])
	
	NA_false = is.not.na(tree_to_chainsaw_table$edge.length)
	
	tree_to_chainsaw$edge.length[parent_branches[NA_false]] = tree_to_chainsaw_table$edge.length[NA_false]
	
	return(tree_to_chainsaw)
	}







add_species_to_skeleton_tree <- function(t1, diversity_of_each_tip)
	{
	# For input into diversitree, you might want to put real diversities/population numbers
	# into each tip of a skeleton tree.  This does that, in the correct format.
	t1$names = c(t1$names, "clades")
	ntaxa = length(t1$tip.label)
	
	

	# blank clades list:
	#t1$clades = c()
	
	clades = c()

	#max_numsp = max(diversity_of_each_tip)
	max_numsp = sum(diversity_of_each_tip)
	#print(max_numsp)
	numdigits = nchar(as.character(ceiling(max_numsp)))
	#sprintf_format_string = "%04.0f"
	sprintf_format_string = paste("%0", numdigits, ".0f", sep="")
	
	old_numsp = double(1)
	for (i in 1:ntaxa)
		{
		tipname = t1$tip.label[i]
		
		numsp = diversity_of_each_tip[i]
		if (numsp < 2)
			{
			numsp = 1
			old_numsp = old_numsp + 1
			new_numsp = old_numsp + numsp
			}
		else
			{
			old_numsp = old_numsp + 1
			new_numsp = old_numsp + numsp			
			}
		#print(old_numsp)
		
		new_numsp = old_numsp + numsp
		
		spnums = sprintf(sprintf_format_string, old_numsp:new_numsp)
		tipname_splist = paste("sp", spnums, sep="")
		old_numsp = new_numsp

		#
		cmdstr = paste("clades$", tipname, " = tipname_splist", sep="")
		eval(parse(text = cmdstr))
		}
	# I guess we need this (?)
	#t1@order = "cladewise"
	#t1$order = "cladewise"
	attr(t1, "order") = "cladewise"
	
	t1$clades = clades
	
	class(t1) <- c("clade.tree", "phylo")
	
	# return to main
	return(t1)		
	}

















###########################################
# PHYLO-BETA DIVERSITY STUFF
###########################################

get_nonphylo_data <- function(obj, input_states, diversity_of_each_tip)
	{
	outlist = c()
	tmpnames = names(obj$clades)
	ntaxa = length(obj$tip.label)
	
	clades = c()

	old_numsp = 0
	for (i in 1:ntaxa)
		{
		tipname = t1$tip.label[i]
		
		numsp = diversity_of_each_tip[i]
		if (numsp < 2)
			{
			numsp = 1
			old_numsp = old_numsp + 1
			new_numsp = old_numsp + numsp			
			}
		else
			{
			old_numsp = old_numsp + 1
			new_numsp = old_numsp + numsp			
			}
		# To make a repetitive list of the grouping taxa (e.g. genera)
		# taxa_to_add = rep(tmpnames[i], length(old_numsp : new_numsp)
		# outlist = c(outlist, taxa_to_add)

		# To make a repetitive list of the character states in the grouping taxa (e.g. genera)
		states_to_add = rep(input_states[i], length(old_numsp : new_numsp))
		outlist = c(outlist, states_to_add)
		}
	
	return(outlist)
	}



# Read in the consensus tree from MrBayes .con file
get_treestring_from_mb_confile <- function(mb_con_fn, findthis="tree con_50_majrule = ", treenum_to_return=1)
	{
	# Read in the text file:
	lines = scan(mb_con_fn, what="character", sep="\n")

	# Initialize the treecount
	treecount = 0
	
	# Go through the lines
	for (i in 1:length(lines))
		{
		# Trim whitespace (trim comes from the library(gregmisc))
		line = trim(lines[i])
		#cat(line, "\n", sep="")
		
		# Take the first line which starts with "tree con_50_majrule = "
		#findthis = "tree con_50_majrule = "
		result = find_instring(findthis, string=line)

		if (result == TRUE)
			{
			treecount = treecount + 1
			
			if (treecount == treenum_to_return)
				{
				return_string = strsplit(line, findthis)[[1]][2]
				cat("\n")
				cat("get_treestring_from_mb_confile() found tree:", "\n", sep="")
				cat(return_string, "\n")
				return(return_string)
				}
			else
				{
				next()
				}
			}
		}
	cat("get_treestring_from_mb_confile(): ERROR: No tree found!!\n")
	return
	}

# Read in consensus tree from MrBayes to APE tree structure (phylo structure)
get_contree_from_mb_confile <- function(mb_con_fn, findthis="tree con_50_majrule = ", treenum_to_return=1)
	{
	treestr = get_treestring_from_mb_confile(mb_con_fn)
	
	# Read in the tree
	temptree = read.tree(text=treestr)
	
	return(temptree)
	}


# Plot NMMDS (Non-metric, multidimensional scaling, from NMDS package) plots of tree distances
plot_NMMDS <- function(nmds_results_conf)
	{
	ticklength = 0
	plot(nmds_results_conf, col="grey", xlab="", ylab="", tcl=ticklength, xaxt="n", yaxt="n")
	points(nmds_results_conf[1,1], nmds_results_conf[1,2], pch="*", cex=3, col="black")
	return()
	}



# plot all conf in a NMDS_results_list
defaults='
UPindices=riboprots_indices
i=1
j=1
'
plot_NMMDS_list <- function(NMDS_results_list, points_to_color, UPindices=NULL)
	{
	for (i in 1:length(NMDS_results_list))
		{
		nmds_results = NMDS_results_list[[i]]
		
		# Go through each of the confs
		for (j in 1:length(nmds_results$conf))
			{
			plot_NMMDS(nmds_results$conf[j][[1]])
			plot_NMMDS_points(nmds_results$conf[j][[1]], points_to_color)
			plot_NMMDS_species_tree(nmds_results$conf[j][[1]])
			
			#print(UPindices)
			if (is.null(UPindices))
				{
				next()
				}
			else
				{
				# Plot the UP proteins as small black dots
				points(nmds_results$conf[j][[1]][UPindices,1], nmds_results$conf[j][[1]][UPindices,2], col="black", pch=19, cex=0.5)
				}
			}
		}
	return()
	}



plot_NMMDS_list2 <- function(NMDS_results_list, points_to_color, UPindices=NULL)
	{
	for (i in 1:length(NMDS_results_list))
		{
		nmds_results = NMDS_results_list[[i]]
		
		# Go through each of the confs
#		for (j in 1:length(nmds_results$conf))
		for (j in 3:3)
			{
			plot_NMMDS(nmds_results$conf[j][[1]])
			plot_NMMDS_points(nmds_results$conf[j][[1]], points_to_color)
			plot_NMMDS_species_tree(nmds_results$conf[j][[1]])
			
			#print(UPindices)
			if (is.null(UPindices))
				{
				next()
				}
			else
				{
				# Plot the UP proteins as small black dots
				points(nmds_results$conf[j][[1]][UPindices,1], nmds_results$conf[j][[1]][UPindices,2], col="black", pch=19, cex=0.5)
				}
			}
		}
	return()
	}




plot_NMMDS_species_tree <- function(nmds_results_conf)
	{
	points(nmds_results_conf[1,1], nmds_results_conf[1,2], pch="*", cex=3, col="black")
	return()
	}


defaults='
i=1
nmds_results_conf=nmds_results$conf[j][[1]]
tmptable = nmds_results_conf
'
plot_NMMDS_points <- function(nmds_results_conf, points_to_color)
	{
	rows_to_keep = points_to_color$index
	# plot each point individually
	for (i in 1:length(rows_to_keep))
		{
		rows_to_plot = gather_rows_with_indices(nmds_results_conf, rows_to_keep[i])
		points(rows_to_plot[1], rows_to_plot[2], col=points_to_color$color[i], pch=points_to_color$cog[i])
		#cat(points_to_color$color[i], points_to_color$cog[i], sep="	")
		
		}
	return()
	}


# do an NMMDS plot of a distance matrix, colored (heat map) to match the rank order of a list
# of numbers that is the length of nrow(distance matrix)
plot_NMMDS_list_heatcolors <- function(NMDS_results_list, value_to_colorize, title1, title2, UPindices=NULL)
	{
	# set up a blank 3-column matrix to hold the color data
	points_to_color = data.frame(matrix(data=NA, nrow=length(value_to_colorize), ncol=3))
	names(points_to_color) = c("index", "cog", "color")
	points_to_color$index = seq(1, nrow(points_to_color), 1)
	points_to_color$cog = as.numeric(1)
	
	# set up the colors list (you could replace heat.colors for some other colors)
	colslist = heat.colors(nrow(points_to_color))
	
	# Just to check, here is a plot of the colors in order
	# (red should equal low)
	plot(1:length(colslist), rep(1, length(colslist)), col=colslist)
	
	
	# set the values to use to rank the color map
	# value_to_colorize = align_stats$avgdiffsbetw / align_stats$avgdiffswithin
	value_to_colorize = value_to_colorize
	
	# add the values to color to a sequence of indices, 1:nrow(points_to_color)
	tmptable = data.frame(cbind(value_to_colorize, seq(1, nrow(points_to_color), 1)))
	names(tmptable) = c("value_to_colorize", "index")
	
	# order by the value
	tmptable = tmptable[order(tmptable$value_to_colorize), ]

	# now that they are ordered in ascending order, 
	colortable = cbind(tmptable, colslist)
	names(colortable) = c(names(tmptable), "color")
	colortable = colortable[order(colortable$index), ]
	points_to_color$color = colortable$color
	

	# Subplots	
	par(mfrow=c(4,3))
	# Subplot margins: c(bottom, left, top, right) 
	par(mar=c(0, 0, 0, 0))
	# outer margins
	par(oma=c(4, 4, 4, 2))	

	# mtext params:
	# outer = TRUE means outer margin
	# adj = parallel outer margin adjust, 0-1
	# padj = perpendic outer marging adjust, 0-1 within outer margin

	plot_NMMDS_list(NMDS_results_list[1:4], points_to_color, UPindices=UPindices)
	title(title1, outer=TRUE, line=1, cex.main=1.5)
	mtext("run1", outer=TRUE, side=1, adj=0.15, padj=1)
	mtext("run2", outer=TRUE, side=1, adj=0.5, padj=1)
	mtext("run3", outer=TRUE, side=1, adj=0.85, padj=1)
	mtext("RF", outer=TRUE, side=2, adj=0.88, padj=-0.5, las=0)
	mtext("Sym", outer=TRUE, side=2, adj=0.62, padj=-0.5, las=0)
	mtext("FPN", outer=TRUE, side=2, adj=0.36, padj=-0.5, las=0)
	mtext("Euc", outer=TRUE, side=2, adj=0.11, padj=-0.5, las=0)
		
	plot_NMMDS_list(NMDS_results_list[5:8], points_to_color, UPindices=UPindices)
	title(title2, outer=TRUE, line=1, cex.main=1.5)
	mtext("run1", outer=TRUE, side=1, adj=0.15, padj=1)
	mtext("run2", outer=TRUE, side=1, adj=0.5, padj=1)
	mtext("run3", outer=TRUE, side=1, adj=0.85, padj=1)
	mtext("RF", outer=TRUE, side=2, adj=0.88, padj=-0.5, las=0)
	mtext("Sym", outer=TRUE, side=2, adj=0.62, padj=-0.5, las=0)
	mtext("FPN", outer=TRUE, side=2, adj=0.36, padj=-0.5, las=0)
	mtext("Euc", outer=TRUE, side=2, adj=0.11, padj=-0.5, las=0)

	
	return()
	}






# do an NMMDS plot of a distance matrix, colored (heat map) to match the rank order of a list
# of numbers that is the length of nrow(distance matrix)
plot_NMMDS_list_heatcolors2 <- function(NMDS_results_list, value_to_colorize, title1, title2, UPindices=NULL)
	{
	# set up a blank 3-column matrix to hold the color data
	tmp_points_to_color = data.frame(matrix(data=NA, nrow=length(value_to_colorize), ncol=3))
	names(tmp_points_to_color) = c("index", "cog", "color")
	tmp_points_to_color$index = seq(1, nrow(tmp_points_to_color), 1)
	tmp_points_to_color$cog = as.numeric(1)
	
	# set up the colors list (you could replace heat.colors for some other colors)
	colslist = heat.colors(nrow(tmp_points_to_color))
	
	# Just to check, here is a plot of the colors in order
	# (red should equal low)
	plot( 1:length(colslist), rep(1, length(colslist)), col=colslist, xaxt="n")

	# Make x-axis
	x_axis_indices_to_tick = pretty(1:length(colslist))
	# Put the literal 1st and last value in
	x_axis_indices_to_tick[1] = 1
	x_axis_indices_to_tick[length(x_axis_indices_to_tick)] = max(length(x_axis_indices_to_tick))
	x_axis_labels_for_ticks = formatC(sort(value_to_colorize)[x_axis_indices_to_tick], format="g", digits=4)
	axis(side=1, at=x_axis_indices_to_tick, labels=x_axis_labels_for_ticks )
	
	# set the values to use to rank the color map
	# value_to_colorize = align_stats$avgdiffsbetw / align_stats$avgdiffswithin
	value_to_colorize = value_to_colorize
	
	# add the values to color to a sequence of indices, 1:nrow(tmp_points_to_color)
	tmptable = data.frame(cbind(value_to_colorize, seq(1, nrow(tmp_points_to_color), 1)))
	names(tmptable) = c("value_to_colorize", "index")
	
	# order by the value
	tmptable = tmptable[order(tmptable$value_to_colorize), ]

	# now that they are ordered in ascending order, 
	colortable = cbind(tmptable, colslist)
	names(colortable) = c(names(tmptable), "color")
	colortable = colortable[order(colortable$index), ]
	tmp_points_to_color$color = colortable$color
	

	# Subplots	
	par(mfrow=c(4,3))
	# Subplot margins: c(bottom, left, top, right) 
	par(mar=c(0, 0, 0, 0))
	# outer margins
	par(oma=c(4, 4, 4, 2))	

	# mtext params:
	# outer = TRUE means outer margin
	# adj = parallel outer margin adjust, 0-1
	# padj = perpendic outer marging adjust, 0-1 within outer margin

	plot_NMMDS_list(NMDS_results_list[1:4], tmp_points_to_color, UPindices=UPindices)
	title(title1, outer=TRUE, line=1, cex.main=1.5)
	mtext("run1", outer=TRUE, side=1, adj=0.15, padj=1)
	mtext("run2", outer=TRUE, side=1, adj=0.5, padj=1)
	mtext("run3", outer=TRUE, side=1, adj=0.85, padj=1)
	mtext("RFdist", outer=TRUE, side=2, adj=0.88, padj=-0.5, las=0)
	mtext("RFdist TL=1", outer=TRUE, side=2, adj=0.62, padj=-0.5, las=0)
	mtext("RFdist maxht=1", outer=TRUE, side=2, adj=0.36, padj=-0.5, las=0)
	mtext("RFdist r8s ht=1", outer=TRUE, side=2, adj=0.11, padj=-0.5, las=0)
	
	return()
	}







########################################################
# TRACING THE ANNOYING R/APE PHYLO OBJECT FORMAT
########################################################

edges_of_droptip <- function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0) 
{	
	#nums = 1:dim(phy$edge)[1]
	##new_edges = phy$edge
	indlist = list()
	if (class(phy) != "phylo") 
		stop("object \"phy\" is not of class \"phylo\"")
	phy <- new2old.phylo(phy)
	phy$edge = cbind(phy$edge, 1:dim(phy$edge)[1])

	if (subtree) {
		trim.internal <- TRUE
		edge.bak <- phy$edge
	}
	tmp <- as.numeric(phy$edge)
	nb.tip <- max(tmp)
	nodes <- setdiff(tmp, 1:nb.tip)
	nobr <- is.null(phy$edge.length)
	if (is.numeric(tip)) 
		tip <- phy$tip.label[tip]
	del <- phy$tip.label %in% tip
	ind <- which(phy$edge[, 2] %in% as.character(which(del)))
	print("Saving ind1 to indlist")
	length_indlist = length(indlist) + 1
	indlist[[length_indlist]] = ind
  	phy$edge <- phy$edge[-ind, ]
  	#new_edges <- #new_edges[-ind, ]`
  	#nums <- #nums[-ind]
	if (!nobr) 
		phy$edge.length <- phy$edge.length[-ind]
	phy$tip.label <- phy$tip.label[!del]
	if (trim.internal) {
		if (root.edge) {
		    seq.nod <- list()
		    for (i in phy$edge[, 2][as.numeric(phy$edge[, 2]) > 
		        0]) {
		        vec <- i
		        j <- i
		        while (j != "-1") {
		          ind <- which(phy$edge[, 2] == j)
					print("Saving ind2 to indlist")
					length_indlist = length(indlist) + 1
					indlist[[length_indlist]] = ind

		          j <- phy$edge[ind, 1]
		          vec <- c(vec, j)
		        }
		        seq.nod[[i]] <- vec
		    }
		    sn <- lapply(seq.nod, rev)
		    i <- 1
		    x <- unlist(lapply(sn, function(x) x[i]))
		    while (length(unique(x)) == 1) {
		        x <- unlist(lapply(sn, function(x) x[i]))
		        i <- i + 1
		    }
		    MRCA <- sn[[1]][i - 2]
		    newrootedge <- if (is.null(phy$root.edge)) 
		        0
		    else phy$root.edge
		    for (i in 1:root.edge) {
		        ind <- which(phy$edge[, 2] == MRCA)
				print("Saving ind3 to indlist")
				length_indlist = length(indlist) + 1
				indlist[[length_indlist]] = ind
		        newrootedge <- newrootedge + phy$edge.length[ind]
		        MRCA <- phy$edge[ind, 1]
		        if (MRCA == "-1" && i < root.edge) {
		          newrootedge <- newrootedge
		          break
		        }
		    }
		    phy$root.edge <- newrootedge
		}
		else {
		    if (!is.null(phy$root.edge)) 
		        phy$root.edge <- NULL
		}
		while (!all(phy$edge[, 2][as.numeric(phy$edge[, 2]) < 
		    0] %in% phy$edge[, 1])) {
		    temp <- phy$edge[, 2][as.numeric(phy$edge[, 2]) < 
		        0]
		    k <- temp %in% phy$edge[, 1]
		    ind <- phy$edge[, 2] %in% temp[!k]
		    phy$edge <- phy$edge[!ind, ]
		    #nums <- #nums[-ind]

				print("Saving ind4 to indlist")
				length_indlist = length(indlist) + 1
				indlist[[length_indlist]] = ind
				#new_edges <- #new_edges[!ind, ]

		    if (!nobr) 
		        phy$edge.length <- phy$edge.length[!ind]
		}
	}
	else {
		k <- nodes %in% phy$edge[, 1]
		ind <- phy$edge[, 2] %in% nodes[!k]
		phy$edge[which(ind), 2] <- as.character(nb.tip + (1:sum(ind)))
		if (is.null(phy$node.label)) 
		    new.tip.label <- rep("NA", sum(ind))
		else new.tip.label <- phy$node.label[!k]
		phy$tip.label <- c(phy$tip.label, new.tip.label)
	}
	useless.nodes <- names(which(table(phy$edge[, 1]) == 1))
	if (subtree) {
		if (!nobr) 
		    mnbr <- mean(phy$edge.length)
		if (length(useless.nodes) == 1) 
		    n <- length(tip)
		else {
		    seq.nod <- list()
		    wh <- numeric(0)
		    for (i in as.character(which(del))) {
		        vec <- i
		        j <- i
		        while (!(j %in% useless.nodes)) {
		          ind <- which(edge.bak[, 2] == j)
		          wh <- c(wh, ind)
		          j <- edge.bak[ind, 1]
		          vec <- c(vec, j)
		        }
		        seq.nod[[i]] <- vec
		    }
		    n <- table(unlist(lapply(seq.nod, function(x) rev(x)[1])))
		}
		new.lab <- paste("[", n, "_tips]", sep = "")
		for (i in 1:length(useless.nodes)) {
		    wh <- which(phy$edge[, 1] == useless.nodes[i])
		    phy$tip.label <- c(phy$tip.label, new.lab[i])
		    if (wh == dim(phy$edge)[1]) {
		        phy$edge <- rbind(phy$edge, c(useless.nodes[i], 
		          as.character(nb.tip + i)))
		        if (!nobr) 
		          phy$edge.length <- c(phy$edge.length, mnbr)
		    }
		    else {
		        phy$edge <- rbind(phy$edge[1:wh, ], c(useless.nodes[i], 
		          as.character(nb.tip + i)), phy$edge[(wh + 1):dim(phy$edge)[1], 
		          ])
		        if (!nobr) 
		          phy$edge.length <- c(phy$edge.length[1:wh], 
		            mnbr, phy$edge.length[(wh + 1):(dim(phy$edge)[1] - 
		              1)])
		    }
		}
	}
	else {
		for (i in useless.nodes) {
		    ind1 <- which(phy$edge[, 1] == i)
		    ind2 <- which(phy$edge[, 2] == i)
		    phy$edge[ind2, 2] <- phy$edge[ind1, 2]
		    phy$edge <- phy$edge[-ind1, ]
		    #nums = #nums[-ind1]
		    if (!nobr) {
		        phy$edge.length[ind2] <- phy$edge.length[ind2] + 
		          phy$edge.length[ind1]
		        phy$edge.length <- phy$edge.length[-ind1]
		    }
		}
	}
	tmp <- as.numeric(phy$edge)
	if (!is.null(phy$node.label)) {
		x <- unique(tmp)
		x <- x[x < 0]
		phy$node.label <- phy$node.label[-x]
	}
	n <- length(tmp)
	nodes <- tmp < 0
	ind.nodes <- (1:n)[nodes]
	ind.tips <- (1:n)[!nodes]
	new.nodes <- -as.numeric(factor(-tmp[nodes]))
	new.tips <- as.numeric(factor(tmp[!nodes]))
	tmp[ind.nodes] <- new.nodes
	tmp[ind.tips] <- new.tips
	dim(tmp) <- c(n/2, 2)
	mode(tmp) <- "character"
	phy$edge <- tmp
	phy <- old2new.phylo(phy)
	if (!trim.internal || subtree) {
		S <- write.tree(phy)
		phy <- if (nobr) 
		    clado.build(S)
		else tree.build(S)
	}
	#return(list(nodes, tmp, new.nodes, new.tips, save_ind1, save_ind3, save_ind4))
	#return(indlist)
	##new_edges = cbind(#new_edges, 1:dim(#new_edges)[1])
	return(phy)
}



trace_edges <- function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0) 
{
	phy2 <- phy
	
	if (class(phy) != "phylo") 
		stop("object \"phy\" is not of class \"phylo\"")
	phy <- new2old.phylo(phy)
	phy2 <- new2old.phylo(phy2)
	phy2$edge[,2] = 1:dim(phy2$edge)[1]
	
	if (subtree) {
		trim.internal <- TRUE
		edge.bak <- phy$edge
	}
	tmp <- as.numeric(phy$edge)
	nb.tip <- max(tmp)
	nodes <- setdiff(tmp, 1:nb.tip)
	nobr <- is.null(phy$edge.length)
	if (is.numeric(tip)) 
		tip <- phy$tip.label[tip]
	del <- phy$tip.label %in% tip
	ind <- which(phy$edge[, 2] %in% as.character(which(del)))
	phy$edge <- phy$edge[-ind, ]
	phy2$edge <- phy2$edge[-ind, ]

	if (!nobr) 
		phy$edge.length <- phy$edge.length[-ind]
	phy$tip.label <- phy$tip.label[!del]
	if (trim.internal) {
		if (root.edge) {
		    seq.nod <- list()
		    for (i in phy$edge[, 2][as.numeric(phy$edge[, 2]) > 
		        0]) {
		        vec <- i
		        j <- i
		        while (j != "-1") {
		          ind <- which(phy$edge[, 2] == j)
		          j <- phy$edge[ind, 1]
		          vec <- c(vec, j)
		        }
		        seq.nod[[i]] <- vec
		    }
		    sn <- lapply(seq.nod, rev)
		    i <- 1
		    x <- unlist(lapply(sn, function(x) x[i]))
		    while (length(unique(x)) == 1) {
		        x <- unlist(lapply(sn, function(x) x[i]))
		        i <- i + 1
		    }
		    MRCA <- sn[[1]][i - 2]
		    newrootedge <- if (is.null(phy$root.edge)) 
		        0
		    else phy$root.edge
		    for (i in 1:root.edge) {
		        ind <- which(phy$edge[, 2] == MRCA)
		        newrootedge <- newrootedge + phy$edge.length[ind]
		        MRCA <- phy$edge[ind, 1]
		        if (MRCA == "-1" && i < root.edge) {
		          newrootedge <- newrootedge
		          break
		        }
		    }
		    phy$root.edge <- newrootedge
		    phy2$root.edge <- newrootedge

		}
		else {
		    if (!is.null(phy$root.edge)) 
		        phy$root.edge <- NULL
		        phy2$root.edge <- NULL

		}
		while (!all(phy$edge[, 2][as.numeric(phy$edge[, 2]) < 
		    0] %in% phy$edge[, 1])) {
		    temp <- phy$edge[, 2][as.numeric(phy$edge[, 2]) < 
		        0]
		    k <- temp %in% phy$edge[, 1]
		    ind <- phy$edge[, 2] %in% temp[!k]
		    phy$edge <- phy$edge[!ind, ]
		    phy2$edge <- phy2$edge[!ind, ]
		    if (!nobr) 
		        phy$edge.length <- phy$edge.length[!ind]
		}

	}
	else {
		k <- nodes %in% phy$edge[, 1]
		ind <- phy$edge[, 2] %in% nodes[!k]
		phy$edge[which(ind), 2] <- as.character(nb.tip + (1:sum(ind)))
		phy2$edge[which(ind), 2] <- as.character(nb.tip + (1:sum(ind)))

		if (is.null(phy$node.label)) 
		    new.tip.label <- rep("NA", sum(ind))
		else new.tip.label <- phy$node.label[!k]
		phy$tip.label <- c(phy$tip.label, new.tip.label)
	}
	useless.nodes <- names(which(table(phy$edge[, 1]) == 1))
	if (subtree) {
		if (!nobr) 
		    mnbr <- mean(phy$edge.length)
		if (length(useless.nodes) == 1) 
		    n <- length(tip)
		else {
		    seq.nod <- list()
		    wh <- numeric(0)
		    for (i in as.character(which(del))) {
		        vec <- i
		        j <- i
		        while (!(j %in% useless.nodes)) {
		          ind <- which(edge.bak[, 2] == j)
		          wh <- c(wh, ind)
		          j <- edge.bak[ind, 1]
		          vec <- c(vec, j)
		        }
		        seq.nod[[i]] <- vec
		    }
		    n <- table(unlist(lapply(seq.nod, function(x) rev(x)[1])))
		}
		new.lab <- paste("[", n, "_tips]", sep = "")
		for (i in 1:length(useless.nodes)) {
		    wh <- which(phy$edge[, 1] == useless.nodes[i])
		    phy$tip.label <- c(phy$tip.label, new.lab[i])
		    if (wh == dim(phy$edge)[1]) {
		        phy$edge <- rbind(phy$edge, c(useless.nodes[i], 
		          as.character(nb.tip + i)))
		        phy2$edge <- rbind(phy2$edge, c(useless.nodes[i], 
		          as.character(nb.tip + i)))

		        if (!nobr) 
		          phy$edge.length <- c(phy$edge.length, mnbr)
		    }
		    else {
		        phy$edge <- rbind(phy$edge[1:wh, ], c(useless.nodes[i], 
		          as.character(nb.tip + i)), phy$edge[(wh + 1):dim(phy$edge)[1], 
		          ])
		        phy2$edge <- rbind(phy2$edge[1:wh, ], c(useless.nodes[i], 
		          as.character(nb.tip + i)), phy2$edge[(wh + 1):dim(phy$edge)[1], 
		          ])
		        if (!nobr) 
		          phy$edge.length <- c(phy$edge.length[1:wh], 
		            mnbr, phy$edge.length[(wh + 1):(dim(phy$edge)[1] - 
		              1)])
		    }
		}
	}
	else {
		for (i in useless.nodes) {
		    ind1 <- which(phy$edge[, 1] == i)
		    ind2 <- which(phy$edge[, 2] == i)
		    phy$edge[ind2, 2] <- phy$edge[ind1, 2]
		    phy2$edge[ind2, 2] <- phy2$edge[ind1, 2]

		    phy$edge <- phy$edge[-ind1, ]
		    phy2$edge <- phy2$edge[-ind1, ]
		    if (!nobr) {
		        phy$edge.length[ind2] <- phy$edge.length[ind2] + 
		          phy$edge.length[ind1]
		        phy$edge.length <- phy$edge.length[-ind1]
		    }
		}
	}
	tmp <- as.numeric(phy$edge)
	tmp2 <- as.numeric(phy2$edge)

	if (!is.null(phy$node.label)) {
		x <- unique(tmp)
		x <- x[x < 0]
		phy$node.label <- phy$node.label[-x]
	}
	n <- length(tmp)
	nodes <- tmp < 0
	ind.nodes <- (1:n)[nodes]
	ind.tips <- (1:n)[!nodes]
	new.nodes <- -as.numeric(factor(-tmp[nodes]))
	new.tips <- as.numeric(factor(tmp[!nodes]))
	tmp[ind.nodes] <- new.nodes
	tmp[ind.tips] <- new.tips
	dim(tmp) <- c(n/2, 2)
	mode(tmp) <- "character"
	phy$edge <- tmp
	phy <- old2new.phylo(phy)

	

	if (!trim.internal || subtree) {
		S <- write.tree(phy)
		phy <- if (nobr) 
		    clado.build(S)
		else tree.build(S)
	}
	old_edge_nums = as.numeric(phy2$edge[,2])
	return(old_edge_nums)
}



#######################################################
# Test strick molecular clock with PAUP
#######################################################



PAUP_test_clock_cmds <- function(nexus_data_fn, PAUPblock_lines, PAUP_logfn="", PAUPnexus_cmds_fn="", returnwhat="log", paup_exe="/Applications/paup")
	{
	defaults='
	PAUPblock_lines = lines_to_keep
	'# end defaults
	
	if (PAUP_logfn == "")
		{
		prefix = get_fn_prefix(fn=nexus_data_fn)
		PAUP_logfn = slashslash(paste(prefix, "_PAUPlog.txt", sep=""))
		}

	if (PAUPnexus_cmds_fn == "")
		{
		prefix = get_fn_prefix(fn=nexus_data_fn)
		PAUPnexus_cmds_fn = slashslash(paste(prefix, "_PAUPcmds.txt", sep=""))
		}
	
	# Does the model utilize pinvar (proportion of invariant sites?)
	lines_match_TF = grepl(pattern="Lset ", x=PAUPblock_lines)
	line_to_check = PAUPblock_lines[lines_match_TF]
	if (grepl(pattern="pinvar=0;", x=line_to_check) == TRUE)
		{
		pinvar=FALSE
		pinvar_txt = "pinv=0"
		} else {
		pinvar=TRUE
		pinvar_txt = "pinv=est"
		}

	cmds = NULL
	inum = 0
	
	# Write the PAUP block to make an NJ tree under the jModeltest model,
	# Then optimize params & branch lengths under no clock
	# And again, midpoint-rooted with clock
	cmds[[(inum=inum+1)]] = "BEGIN PAUP;"
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[ start logfile ]"
	cmds[[(inum=inum+1)]] = paste("log start replace=yes file=", PAUP_logfn, " ;", sep="")
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[ Load data in NEXUS format ]"
	cmds[[(inum=inum+1)]] = paste("execute ", nexus_data_fn, " ;", sep="")
	cmds[[(inum=inum+1)]] = ""
	
	
	# Load in the jModelTest lines; leave out BEGIN/END
	cmds[[(inum=inum+1)]] = "[ PAUP block from jModelTest to determine model ]"
	for (i in 1:length(PAUPblock_lines))
		{
		# Exclude certain lines:
		if ( (PAUPblock_lines[i] == "BEGIN PAUP;") || (PAUPblock_lines[i] == "END;"))
			{
			next()
			}
		
		# Otherwise, add them
		cmds[[(inum=inum+1)]] = PAUPblock_lines[i]
		}
	cmds[[(inum=inum+1)]] = ""

	cmds[[(inum=inum+1)]] = "[  ]"
	cmds[[(inum=inum+1)]] = ""
	
	# Estimate a distance tree, using the ML parameters for distance
	cmds[[(inum=inum+1)]] = "[ Estimate a distance tree, using the ML parameters for distance ]"
	cmds[[(inum=inum+1)]] = "set criterion=distance;"
	cmds[[(inum=inum+1)]] = "[  ]"
	cmds[[(inum=inum+1)]] = "[ Set the distance criterion to the ML parameters ]"
	cmds[[(inum=inum+1)]] = "dset distance=ML;"
	cmds[[(inum=inum+1)]] = "showdist;"
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[ Give the tree an implied root at the midpoint ]"
	cmds[[(inum=inum+1)]] = "[ Then calculate the NJ tree, to get the topology ]"
	cmds[[(inum=inum+1)]] = "set root=midpoint;"
	cmds[[(inum=inum+1)]] = "nj;"
	cmds[[(inum=inum+1)]] = "showtrees;"
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[ Calculate the likelihood (assuming GTR+G, I is optional)]"
	cmds[[(inum=inum+1)]] = "set criterion=likelihood;"
	cmds[[(inum=inum+1)]] = paste("lset nst=6 rmat=est basefreq=est rates=gamma shape=est ", pinvar_txt, ";", sep="")
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[ This is the -lnL of the alternative hypothesis (all branches vary independently) ]"
	cmds[[(inum=inum+1)]] = "lscores;"
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[Now, we will change the likelihood settings to enforce a molecular clock:]"
	cmds[[(inum=inum+1)]] = "set root=midpoint;"
	cmds[[(inum=inum+1)]] = "roottrees;"
	cmds[[(inum=inum+1)]] = "showtrees;"
	cmds[[(inum=inum+1)]] = paste("lset nst=6 rmat=est basefreq=est rates = gamma shape = est ", pinvar_txt, " clock=yes;", sep="")
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[ This is the -lnL of the null hypothesis (strict clock) ]"
	cmds[[(inum=inum+1)]] = "lscores;"
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "[ Stop the log ]"
	cmds[[(inum=inum+1)]] = "log stop;"
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = ""
	cmds[[(inum=inum+1)]] = "END;"
	
	
	# Now, write the commands to text...
	write.table(x=cmds, file=PAUPnexus_cmds_fn, append=FALSE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)
	
	cmdstr = paste(paup_exe, " -n ", PAUPnexus_cmds_fn, sep="")
	
	# Run, or just return the commands file?
	if (returnwhat == "log")
		{
		system(cmdstr)
		return(PAUP_logfn)
		} 
	if (returnwhat == "PAUPnexus")
		{
		return(PAUPnexus_cmds_fn)
		} 
	
	
	return(PAUPnexus_cmds_fn)
	}

jModelTest <- function(fasta_fn="", logfn="", cmdpart3=" -g 4 -i -f -AIC -BIC -a -w", jModelTest_dir="/Applications/jmodeltest-2.1.3/", returnwhat="logfn")
	{
	defaults='
	fasta_fn=""
	logfn=""
	cmdpart3=" -g 4 -i -f -AIC -BIC -a -w"
	jModelTest_dir="/Applications/jmodeltest-2.1.3/"
	'
	
	# Example data comes from here:
	# http://code.google.com/p/jmodeltest2/wiki/QuickStartConsole
	#
	# java -jar jModelTest.jar -d example-data/aP6.fas -g 4 -i -f -AIC -BIC -a
	# 
	
	# Check inputs
	if (fasta_fn == "")
		{
		fasta_fn = slashslash(paste(jModelTest_dir, "/example-data/aP6.fas", sep=""))
		cat("\nWARNING: Running default jModelTest example dataset: ", fasta_fn, "\n", sep="")
		}
	
	if (logfn == "")
		{
		prefix = get_fn_prefix(fn=fasta_fn)
		logfn = slashslash(paste(prefix, "_log.txt", sep=""))
		}
	
	
	cmdpart1 = "java -jar "
	cmdpart2 = "jModelTest.jar -d "
	# cmdpart3=" -g 4 -i -f -AIC -BIC -a -w"
	
	jModelTest_cmd = paste(cmdpart1, jModelTest_dir, cmdpart2, fasta_fn, cmdpart3, " > ", logfn, sep="")
	
	cat("\nResults in:\n", logfn, sep="")
	
	
	if (returnwhat == "logfn")
		{
		cat("\nRunning:\n", jModelTest_cmd, "\n", sep="")
		system(jModelTest_cmd)
		cat("\n", sep="")
		cat("\n", sep="")
		cat("PAUP run done.\n", sep="")
		cat("\n", sep="")
		cat("\nResults in:\n", logfn, sep="")
		return(logfn)
		}
	
	if (returnwhat == "cmd")
		{
		cat("\n", sep="")
		cat("\n", sep="")
		cat("\nReturning jModelTest_cmd:\n", jModelTest_cmd, "\n", sep="")
		cat("\n", sep="")
		return(jModelTest_cmd)
		}
	
	}



#########################
# Tree stats
#########################
SPR_distance_with_PAUP <- function(tree1, tree2)
	{
	blah='
	cd /Users/nick/Desktop/_LGT/_cyano/_04_treedists/ 
	execute /Users/nick/Desktop/_LGT/_cyano/_02_mb_runs_v1_noParres_plus_sp/A_sptree_prc_simp.nexus.sub1.t;

	[! filter trees so you get only the last 50]
	filter NUMGE=51;
	
	filter NUMGE=195;
	treedist metric=symdiff;
	treedist metric=agd1;
	
	[MAST (Maximum agreement subtree) after Goddard et al. 1994.]
	[FD = no means suppress Frequency Distribution]
	treedist metric=agreement;
	
	
	[load more trees]
	gettrees FILE=/Users/nick/Desktop/_LGT/_cyano/_02_mb_runs_v1_noParres_plus_sp/ProMED4_637448992.nexus.sub1.t MODE=7;
	'
	}


continuous_tree_NJ <- function(params, colnames_to_use, normalize=TRUE, printflag=FALSE)
	{
	
	# subset parameters
	cols_tokeep = names(params) %in% colnames_to_use == TRUE
	tmp_params = params[, cols_tokeep]
	
	prflag(tmp_params, printflag)
	
	if (normalize == TRUE)
		{
		tmp_params = normalize_by_colmax(tmp_params)
		}

	prflag(tmp_params, printflag)

	
	distmat = dist(tmp_params)
	params_tr = nj(distmat)
	
	return(params_tr)
	}


shuffle_tips_dists_to_tree <- function(reftree, params_tr, N=100, topodist=TRUE, bhvdist=TRUE, rfdist=FALSE, treedists=TRUE)
	{
	# N=1000; topodist=TRUE; bhvdist=TRUE; rfdist=FALSE; treedists=TRUE
	
	reftree = unroot(reftree)
	tmp_params_tr = params_tr
	
	# make a null hypothesis
	tiplabs = params_tr$tip.label
	len_tiplabs = length(tiplabs)
	topodist_null = c()
	BHVdist_null = c()
	rfdist_null = c()
	treedists_null = c()
	for (i in 1:N)
		{
		# randomly shuffle tips, with replacement
		tmp_params_tr = params_tr
		
		newtips = sample(tiplabs, size=len_tiplabs, replace=FALSE)
		tmp_params_tr$tip.label = newtips
		
		if (topodist == TRUE)
			{
			tmp_topodist = dist.topo(tmp_params_tr, reftree, method="PH85")
			topodist_null = c(topodist_null, tmp_topodist)
			}
		if (bhvdist == TRUE)
			{
			tmp_BHVdist = dist.topo(tmp_params_tr, reftree, method="BHV01")
			BHVdist_null = c(BHVdist_null, tmp_BHVdist)
			}
		if (rfdist == TRUE)
			{
			tmp_rfdist = RF.dist(tmp_params_tr, reftree)
			rfdist_null = c(rfdist_null, tmp_rfdist)
			}
		if (treedists == TRUE)
			{
			tmp_treedists = treedist(tmp_params_tr, reftree)
			treedists_null = rbind(treedists_null, tmp_treedists)
			}
		}
	
	# make a big table
	dists_table = cbind(topodist_null, BHVdist_null, rfdist_null, treedists_null)
	
	# blank out the row.names (which are identical)
	row.names(dists_table)=c()
	
	dists_table = adf(dists_table)

	# devise the data.frame names
	tmpnames = c()
	if (topodist == TRUE)
		{
		tmpnames = c(tmpnames, "topodist")
		}
	if (bhvdist == TRUE)
		{
		tmpnames = c(tmpnames, "bhvdist")
		}
	if (rfdist == TRUE)
		{
		tmpnames = c(tmpnames, "rfdist")
		}
	if (treedists == TRUE)
		{
		tmpnames = c(tmpnames, names(tmp_treedists))
		}

	names(dists_table) = tmpnames
	
	
	# observed distances added as the last row
	dists_observed = c()
	if (topodist == TRUE)
		{
		tmpdist = dist.topo(params_tr, reftree, method="PH85")
		dists_observed = c(dists_observed, tmpdist)
		}
	if (bhvdist == TRUE)
		{
		tmpdist = dist.topo(params_tr, reftree, method="BHV01")
		dists_observed = c(dists_observed, tmpdist)
		}
	if (rfdist == TRUE)
		{
		tmpdist = RF.dist(params_tr, unroot(reftree))
		dists_observed = c(dists_observed, tmpdist)
		}
	if (treedists == TRUE)
		{
		tmpdist = treedist(params_tr, reftree)
		dists_observed = c(dists_observed, tmpdist)
		}
	
	dists_table = rbind(dists_table, dists_observed)
	
	return(dists_table)
	}


hist_obs_null_diffs <- function(reftree, params_tr, dists_infn=NULL, dists_outfn=NULL, N=100, reftreetxt="", paramstr_txt="", datatxt="")
	{
	# Calculate the distances, if necesary

	numdists = N
	if (is.null(dists_infn))
		{
		dists_table_final = shuffle_tips_dists_to_tree(reftree, params_tr, N=numdists)
		
		# write out distances
		write_table_good(dists_table_final, dists_outfn)
		} else {
		# read in observed distances
		dists_table_final = read_table_good(dists_infn)
		}

	# the observed distances are the last row of the null differences
	dists_observed = unlist(dists_table_final[nrow(dists_table_final), ])

	# do subplots for the two trees
	par(mfrow=c(1,2))
	
	maintxt = reftreetxt
	plot(unroot(reftree), type="unrooted", lab4ut="axial")
	title(maintxt, cex=0.8)
	
	maintxt = paste(paramstr_txt, "\n", datatxt, sep="")
	plot(params_tr, type="unrooted", lab4ut="axial")
	title(maintxt, cex=0.8)
	
	# do histograms of the six different distances
	# set outer margins & number of subplots
	par(oma=c(1,1,5,1), mfrow=c(2,3))
	# histogram all the distances (rfdist = symmetric distance, so exclude)
	for (i in 1:ncol(dists_table_final))
		{
		x = dists_table_final[, i]
		observed_val = dists_observed[i]
		#pval = empirical_pval_slow(x, observed_val)
		pval = round(empirical_pval_fast(x, observed_val), 4)
		tmp_title = paste(names(dists_table_final)[i], "\nnon-parametric p-value = ", pval, sep="")
		
		h = hist(x, main=tmp_title, breaks=50, cex.main=1.2, xlab="distance")
		
		# add the arrow
		arrowtop = 0.3*max(h$counts, na.rm=TRUE)
		arrowbot = 0.03*max(h$counts, na.rm=TRUE)
		arrows(observed_val, arrowtop, observed_val, arrowbot, col="blue", lwd=4)	
		}
	maintxt = paste("Observed distance between ", reftreetxt, " and ", paramstr_txt, "\ncompared to null distribution of distances with tips randomized\n", datatxt, ", N=", N, sep="")
	mtext(maintxt, outer=TRUE)

	}




# Extract a newick string from a NEXUS newick string
extract_newickstr_from_nexusstr <- function(nexusstr, delimiter = " = ")
	{
	words = strsplit(nexusstr, split=delimiter)[[1]]
	newickstr = words[2]
	return(newickstr)
	}




scalebar_loc = function(plotinfo, loc="bottomright", spacer=0.1)
	{
	xmin = plotinfo$x.lim[1]
	xmax = plotinfo$x.lim[2]
	ymin = plotinfo$y.lim[1]
	ymax = plotinfo$y.lim[2]
	
	
	if(loc == "bottomright")
		{
		xloc = (1-spacer) * (xmax-xmin) + xmin
		yloc = (0+spacer) * (ymax-ymin) + ymin
		}

	if(loc == "topright")
		{
		xloc = (1-spacer) * (xmax-xmin) + xmin
		yloc = (1-spacer) * (ymax-ymin) + ymin
		}

	if(loc == "bottomleft")
		{
		xloc = (0+spacer) * (xmax-xmin) + xmin
		yloc = (0+spacer) * (ymax-ymin) + ymin
		}

	if(loc == "topleft")
		{
		xloc = (0+spacer) * (xmax-xmin) + xmin
		yloc = (1-spacer) * (ymax-ymin) + ymin
		}
	
	coords = c(xloc, yloc)
	
	return(coords)
	}

add_scalebar_good = function(plotinfo, loc="bottomright", spacer=0.1, ...)
	{
	coords = scalebar_loc(plotinfo, loc, spacer, ...)
	add.scale.bar(x=coords[1], y=coords[2])
	}


get_treelength <- function(tr)
	{
	treelength = sum(tr$edge.length)
	return(treelength)
	}

rescale_tree <- function(tree, treelength=1)
	{
	tr2 = tree
	tr2$edge.length = (tr2$edge.length / sum(tr2$edge.length) * treelength)
	return(tr2)
	}

setup = '
treelength = 1
'
rescale_newickstr_tree <- function(tmpstr, treelength=1)
	{
	# Take a standard Newick string (tipnames and branchlengths, no other numbers)
	# and rescale the numbers so they add to 1
	
	# Extract just the numbers, with (numbers) (dot) (numbers) pattern
	match_indicies = gregexpr("(?:(\\d+)+\\.+(\\d+))+", tmpstr)[[1]]
	match_indicies_end = match_indicies-1+attr(match_indicies,"match.length")
	
	# Use mapply (multiple apply) to extract those numbers
	x = as.character(mapply(substr, tmpstr, match_indicies, match_indicies_end))
	
	#len(match_indicies)
	#[1] 53
	#len(match_indicies_end)
	#[1] 53
	#len(x)
	#[1] 53
	
	# Convert from strings to numeric
	brlens = as.character(x)
	sumtl = sum(as.numeric(brlens))
	
	brlens_numeric = as.numeric(brlens)
	
	# number of characters in a brlen number
	brlens_numchar = nchar(brlens)
	#unlist(lapply(brlens, nchar))
	


	# Normalize to sum to 1
	brlens2 = brlens_numeric / sumtl * treelength
	# sum(brlens2)

	brlens_numchar_before_decimal = nchar(floor(brlens2))
	brlens_numchar_sprintf_str = mapply(paste, "%0.", brlens_numchar-1-brlens_numchar_before_decimal, "f", sep="", collapse="")

	# create new string of numbers (appropriate length for each number; this assumes 
	# e.g. 0.003 or whatever.
	brlens3 = mapply(sprintf, brlens_numchar_sprintf_str, brlens2)

	# format again to deal with changes from e.g. 10.23 to 0.23
	brlens_numchar_before_decimal2 = nchar(floor(as.numeric(brlens3)))
	brlens_numchar_sprintf_str2 = mapply(paste, "%0.", brlens_numchar-1-brlens_numchar_before_decimal2, "f", sep="", collapse="")


	brlens4 = mapply(sprintf, brlens_numchar_sprintf_str2, as.numeric(brlens3))

	
	#formatstr_before_decimal1 = mapply(rep, "0", brlens_numchar_before_decimal-1)
	#formatstr_before_decimal2 = mapply(paste, formatstr_before_decimal1, brlens_numchar_before_decimal, sep="", collapse="")
	
	#formatstr_before_decimal3 = mapply(paste, "%0", brlens_numchar_before_decimal, ".", brlens_numchar-1-brlens_numchar_before_decimal, "f", sep="", collapse="")
	
	
	#brlens_numchar_sprintf_str = paste("%", formatstr_before_decimal2, ".", brlens_numchar-(brlens_numchar_before_decimal+1), "f", sep="")
	
	#brlens_numchar_sprintf_str = formatstr_before_decimal3
	
	#brlens3_numchars = nchar(brlens3)
	

	brlens4_str = paste(brlens4, sep="", collapse="")
	brlens4_chars = strsplit(brlens4_str, split="")[[1]]
	#length(brlens3_chars)
	
	
	#cbind(unlist(lapply(brlens, nchar)) , unlist(lapply(brlens4, nchar)) )
	
	# Put them back into the string as array
	tmpchars = strsplit(tmpstr, split="")[[1]]
	tmpchars2 = tmpchars
	
	indicies = unlist(mapply(seq, match_indicies, match_indicies_end, 1))
	
	tmpchars2[indicies]
	#length(tmpchars2[indicies])
	
	tmpchars2[indicies] = brlens4_chars
	
	tmpstr2 = paste(tmpchars2, sep="", collapse="")
	
	return(tmpstr2)
	}








# Source: https://svn.mpl.ird.fr/ape/dev/ape/R/drop.tip.R
# Fixed bug identified here:
# https://stat.ethz.ch/pipermail/r-sig-phylo/2010-November/000850.html
# 
## drop.tip.R (2010-11-24)

##   Remove Tips in a Phylogenetic Tree

## Copyright 2003-2010 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

extract.clade <- function(phy, node, root.edge = 0, interactive = FALSE)
{
    Ntip <- length(phy$tip.label)
    ROOT <- Ntip + 1
    Nedge <- dim(phy$edge)[1]
    wbl <- !is.null(phy$edge.length)
    if (interactive) node <- identify(phy)$nodes else {
        if (length(node) > 1) {
            node <- node[1]
            warning("only the first value of 'node' has been considered")
        }
        if (is.character(node)) {
            if (is.null(phy$node.label))
                stop("the tree has no node labels")
            node <- which(phy$node.label %in% node) + Ntip
        }
        if (node <= Ntip)
            stop("node number must be greater than the number of tips")
    }
    if (node == ROOT) return(phy)
    phy <- reorder(phy) # insure it is in cladewise order
    root.node <- which(phy$edge[, 2] == node)
    start <- root.node + 1 # start of the clade looked for
    anc <- phy$edge[root.node, 1] # the ancestor of 'node'
    next.anc <- which(phy$edge[-(1:start), 1] <= anc) # find the next occurence of 'anc' or an 'older' node

    keep <- if (length(next.anc)) start + 0:(next.anc[1] - 1) else start:Nedge

    if (root.edge) {
        NewRootEdge <- phy$edge.length[root.node]
        root.edge <- root.edge - 1
        while (root.edge) {
            if (anc == ROOT) break
            i <- which(phy$edge[, 2] ==  anc)
            NewRootEdge <- NewRootEdge + phy$edge.length[i]
            root.edge <- root.edge - 1
            anc <- phy$edge[i, 1]
        }
        if (root.edge && !is.null(phy$root.edge))
            NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
    }

    phy$edge <- phy$edge[keep, ]
    if (wbl) phy$edge.length <- phy$edge.length[keep]
    TIPS <- phy$edge[, 2] <= Ntip
    tip <- phy$edge[TIPS, 2]
    phy$tip.label <- phy$tip.label[sort(tip)] # <- added sort to avoid shuffling of tip labels (2010-07-21)
    ## keep the ordering so no need to reorder tip.label:
    phy$edge[TIPS, 2] <- order(tip)
    if (!is.null(phy$node.label))
        phy$node.label <- phy$node.label[sort(unique(phy$edge[, 1])) - Ntip]
    Ntip <- length(phy$tip.label)
    phy$Nnode <- dim(phy$edge)[1] - Ntip + 1L
    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format -- same as in root()
    newNb <- integer(Ntip + phy$Nnode)
    newNb[node] <- Ntip + 1L
    sndcol <- phy$edge[, 2] > Ntip
    ## executed from right to left, so newNb is modified before phy$edge:
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <-
        (Ntip + 2):(Ntip + phy$Nnode)
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    phy
}

drop.tip2 <-
    function(phy, tip, trim.internal = TRUE, subtree = FALSE,
             root.edge = 0, rooted = is.rooted(phy), interactive = FALSE)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')

    Ntip <- length(phy$tip.label)
    ## find the tips to drop:
    if (interactive) {
        cat("Left-click close to the tips you want to drop; right-click when finished...\n")
        xy <- locator()
        nToDrop <- length(xy$x)
        tip <- integer(nToDrop)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        for (i in 1:nToDrop) {
            d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
            tip[i] <- which.min(d)
        }
    } else {
        if (is.character(tip))
            tip <- which(phy$tip.label %in% tip)
    }
    if (any(tip > Ntip))
        warning("some tip numbers were higher than the number of tips")

    if (!rooted && subtree) {
        phy <- root(phy, (1:Ntip)[-tip][1])
        root.edge <- 0
    }

    phy <- reorder(phy)
    NEWROOT <- ROOT <- Ntip + 1
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    if (subtree) {
        trim.internal <- TRUE
        tr <- reorder(phy, "pruningwise")
        N <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
                as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]),
                as.integer(Nedge), double(Ntip + Nnode),
                DUP = FALSE, PACKAGE = "ape")[[6]]
    }
    wbl <- !is.null(phy$edge.length)
    edge1 <- phy$edge[, 1] # local copies
    edge2 <- phy$edge[, 2] #
    keep <- !logical(Nedge)

    ## delete the terminal edges given by `tip':
    keep[match(tip, edge2)] <- FALSE

    if (trim.internal) {
        ints <- edge2 > Ntip
        ## delete the internal edges that do not have anymore
        ## descendants (ie, they are in the 2nd col of `edge' but
        ## not in the 1st one)
        repeat {
            sel <- !(edge2 %in% edge1[keep]) & ints & keep
            if (!sum(sel)) break
            keep[sel] <- FALSE
        }
        if (subtree) {
            ## keep the subtending edge(s):
            subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
            keep[subt] <- TRUE
        }
        if (root.edge && wbl) {
            degree <- tabulate(edge1[keep])
            if (degree[ROOT] == 1) {
                j <- integer(0) # will store the indices of the edges below the new root
                repeat {
                    i <- which(edge1 == NEWROOT & keep)
                    j <- c(i, j)
                    NEWROOT <- edge2[i]
                    degree <- tabulate(edge1[keep])
                    if (degree[NEWROOT] > 1) break
                }
                keep[j] <- FALSE
                if (length(j) > root.edge) j <- 1:root.edge
                NewRootEdge <- sum(phy$edge.length[j])
                if (length(j) < root.edge && !is.null(phy$root.edge))
                    NewRootEdge <- NewRootEdge + phy$root.edge
                phy$root.edge <- NewRootEdge
            }
        }
    }

    if (!root.edge) phy$root.edge <- NULL

    ## drop the edges
    phy$edge <- phy$edge[keep, ]
    if (wbl) phy$edge.length <- phy$edge.length[keep]

    ## find the new terminal edges (works whatever 'subtree' and 'trim.internal'):
    TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])

    ## get the old No. of the nodes and tips that become tips:
    oldNo.ofNewTips <- phy$edge[TERMS, 2]

    ## in case some tips are dropped but kept because of 'subtree = TRUE':
    if (subtree) {
        i <- which(tip %in% oldNo.ofNewTips)
        if (length(i)) {
            phy$tip.label[tip[i]] <- "[1_tip]"
            tip <- tip[-i]
        }
    }

    n <- length(oldNo.ofNewTips) # the new number of tips in the tree

    ## the tips may not be sorted in increasing order in the
    ## 2nd col of edge, so no need to reorder $tip.label
    phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
    phy$tip.label <- phy$tip.label[-tip]

    ## make new tip labels if necessary:
    if (subtree || !trim.internal) {
        ## get the numbers of the nodes that become tips:
        node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
        new.tip.label <- if (subtree) {
            paste("[", N[node2tip], "_tips]", sep = "")
        } else {
            if (is.null(phy$node.label)) rep("NA", length(node2tip))
            else phy$node.label[node2tip - Ntip]
        }
        if (!is.null(phy$node.label))
            phy$node.label <- phy$node.label[-(node2tip - Ntip)]
        phy$tip.label <- c(phy$tip.label, new.tip.label)
    }

    ## update node.label if needed:
    if (!is.null(phy$node.label))
        phy$node.label <- phy$node.label[sort(unique(phy$edge[, 1])) - Ntip]

    phy$Nnode <- dim(phy$edge)[1] - n + 1L # update phy$Nnode

    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format -- same as in root()
    newNb <- integer(n + phy$Nnode)
    newNb[NEWROOT] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    ## executed from right to left, so newNb is modified before phy$edge:
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <-
        (n + 2):(n + phy$Nnode)
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    storage.mode(phy$edge) <- "integer"
    collapse.singles(phy)
}







# Simplify NEXUS names to match universal standards
# 
# Names file:
# column 1: original names
# column 2: new names
# 
simplify_NEXUS_names <- function(nameslist=NULL, fn=NULL, maxchar=31)
	{
	if (is.null(nameslist) && is.null(fn))
		{
		errortxt = "\n\nERROR: you must specify either nameslist or fn (filename)\n\n"
		stop(errortxt)
		} # END if (is.null(nameslist) && is.null(fn))
		
	if (!is.null(nameslist) && !is.null(fn))
		{
		warningtxt = "\n\nYou specified both nameslist and fn (filename); taking nameslist.\n\n"
		
		cat(warningtxt)
		
		# Set filename to NULL
		fn = NULL
		} # END if (is.null(nameslist) && is.null(fn))
	
	if ( is.null(nameslist) )
		{
		# Tab-delimited read (but you should only have 1 column)
		df = read.table(file=fn, header=FALSE, sep="\n", strip.white=TRUE, stringsAsFactors=FALSE, quote=NULL, as.is=TRUE)
		nameslist = df[,1]
		} # END if ( is.null(nameslist) )
	
	
	# Make a new nameslist
	new_nameslist = rep("", length=length(nameslist))
	
	# Remove a bunch of ANNOYING characters
	for (i in 1:length(nameslist))
		{
		tmpname = nameslist[i]
		tmpname = gsub("'", "", x=tmpname)
		
		tmpname = gsub("\\*", "", x=tmpname)
		tmpname = gsub("\\(", "", x=tmpname)
		tmpname = gsub("\\)", "", x=tmpname)
		tmpname = gsub("\\?", "", x=tmpname)
		tmpname = gsub("\\+", "_", x=tmpname)
		tmpname = gsub("\\-", "_", x=tmpname)
		tmpname = gsub("\\.", "_", x=tmpname)
		tmpname = gsub("cf_", "", x=tmpname)
		tmpname = gsub("__", "_", x=tmpname)
		tmpname = gsub("__", "_", x=tmpname)
		tmpname = gsub("__", "_", x=tmpname)

		new_nameslist[i] = tmpname
		}

	# Split the names up
	new_nameslist2 = new_nameslist
	for (i in 1:length(new_nameslist))
		{
		tmpname = new_nameslist[i]
		
		words = strsplit(tmpname, split="_")[[1]]
		
		if (length(words) == 1)
			{
			genus_sp = words
			}
		
		if (length(words) >= 2)
			{
			genus_sp = paste(words[1], "_", words[2], sep="")
			}
		
		if (length(words) > 2)
			{
			suffix = ""
			for (j in 3:length(words))
				{
				suffix = paste(suffix, "_", words[j], sep="")
				}
			} else {
			suffix = ""
			}
		
		tmpname2 = paste(genus_sp, suffix, sep="")
		
		# Check the length, subset if too long
		if (nchar(tmpname2) > maxchar)
			{
			tmpname2 = substr(x=tmpname2, start=1, stop=maxchar)
			}
		new_nameslist2[i] = tmpname2
		}
	
	# Check for duplicates
	uniq_names = unique(new_nameslist2)
	
	if (length(uniq_names) < length(new_nameslist2))
		{
		for (u in 1:length(uniq_names))
			{
			uniq_name = uniq_names[u]
			TF = new_nameslist2 == uniq_name
			if (sum(TF) > 1)
				{
				nums = 1:sum(TF)
				max_nchar = max(nchar(as.character(nums)))
				sprintf_fmt = paste("%0", max_nchar, ".f", sep="")
				for (v in 1:length(new_nameslist[TF]))
					{
					numtxt = sprintf(fmt=sprintf_fmt, v)
					length_numtxt = 2+nchar(numtxt)
					
					# Remove "_" if at the end
					current_string = remove_ending_character(new_nameslist2[TF][v], char="_")
					
					length_modstr = nchar(current_string) + length_numtxt
					amount_over = length_modstr - maxchar
					if (amount_over > 0)
						{
						stopval = nchar(current_string) - amount_over
						modname = remove_ending_character(substr(x=current_string, start=1, stop=stopval), char="_")
						} else {
						modname = remove_ending_character(current_string, char="_")
						}
					modname2 = paste(modname, "_n", numtxt, sep="")
					new_nameslist2[TF][v] = modname2
					}
				}
			}
		}	


	# Remove a bunch of "_" characters
	for (i in 1:length(new_nameslist2))
		{
		tmpname = new_nameslist2[i]
		tmpname = gsub("__", "_", x=tmpname)
		tmpname = gsub("__", "_", x=tmpname)
		tmpname = gsub("__", "_", x=tmpname)
		tmpname = remove_ending_character(tmpname, char="_")
		new_nameslist2[i] = tmpname
		}
	
	# Make a table to return
	orig2new_names = cbind(nameslist, new_nameslist2)
	tmp_rownames = 1:nrow(orig2new_names)
	orig2new_names = as.data.frame(orig2new_names, row.names=tmp_rownames, stringsAsFactors=FALSE)

	names(orig2new_names) = c("orig", "new")
	
	return(orig2new_names)
	}


# This function converts ambiguous DNA codes to their single-character IUPAC codes.
nexus_file_DNA_to_IUPAC <- function(file)
	{
	defaults='
	file = "/Dropbox/_njm/__packages/BEASTmasteR_setup/inst/extdata/Venerids_morphDNA_v1/z_old/v130527_DNA_simp3.nex"
	'
	#######################################################
	# If desired, use GREP to convert all ambiguous characters to their IUPAC codes, 
	# according to:
	# http://www.boekhoff.info/?pid=data&dat=fasta-codes
	#######################################################
	strs = NULL
	sn = 0
	
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\)/M/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\)/R/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ T\\)/W/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\)/S/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ T\\)/Y/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(G\\ T\\)/K/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\)/V/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ T\\)/H/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\ T\\)/D/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\ T\\)/B/g;" -pi $(find ', file, ' -type f)', sep="")
	strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\ T\\)/N/g;" -pi $(find ', file, ' -type f)', sep="")
	
	cat("\n\n")
	for (i in 1:sn)
		{
		cat(paste("Running: ", strs[[i]], "...\n", sep=""))
		system(strs[[i]])
		}
	
	}



read_nexus_data_w_Errors <- function (file, convert_ambiguous_to_IUPAC=FALSE) 
	{
	defaults='
	file = "/Dropbox/_njm/__packages/BEASTmasteR_setup/inst/extdata/Venerids_morphDNA_v1/z_old/v130527_DNA_simp3.nex"
	'
	#######################################################
	# If desired, use GREP to convert all ambiguous characters to their IUPAC codes, 
	# according to:
	# http://www.boekhoff.info/?pid=data&dat=fasta-codes
	#######################################################
	if (convert_ambiguous_to_IUPAC == TRUE)
		{
		strs = NULL
		sn = 0
		
		# Copy the file to a new filename
		
		
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\)/M/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\)/R/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ T\\)/W/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\)/S/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ T\\)/Y/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(G\\ T\\)/K/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\)/V/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ T\\)/H/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\ T\\)/D/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\ T\\)/B/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\ T\\)/N/g;" -pi $(find ', file, ' -type f)', sep="")
		
		cat("\n\n")
		for (i in 1:sn)
			{
			cat(paste("Running: ", strs[[i]], "...\n", sep=""))
			system(strs[[i]])
			}

		}
	
	
	#######################################################
	# "find.ntax"
	#######################################################
	"find.ntax" <- function(x)
		{ 
		for (i in 1:NROW(x))
			{
			if (any(f <- grep("\\bntax", x[i], ignore.case = TRUE)))
			{
				ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
				  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
				break
				}
			}
		ntax
		}

	#######################################################
	# "find.nchar"
	#######################################################
	"find.nchar" <- function(x)
		{
		for (i in 1:NROW(x))
			{
			if (any(f <- grep("\\bnchar", x[i], ignore.case = TRUE)))
			{
				nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)", 
				  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
				break
				}
			}
		nchar
		}


	#######################################################
	# "find.matrix.line"
	#######################################################
	"find.matrix.line" <- function(x)
		{
		for (i in 1:NROW(x))
			{
			if (any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE)))
			{
				matrix.line <- as.numeric(i)
				break
				}
			}
		matrix.line
		}


	#######################################################
	# "trim.whitespace"
	#######################################################
	"trim.whitespace" <- function(x)
		{
		gsub("\\s+", "", x)
		}

	#######################################################
	# "trim.semicolon"
	#######################################################
	"trim.semicolon" <- function(x)
		{
		gsub(";", "", x)
		}

	X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE, comment.char = "[", strip.white = TRUE)
	ntax <- find.ntax(X)
	nchar <- find.nchar(X)
	matrix.line <- find.matrix.line(X)
	start.reading <- matrix.line + 1
	Obj <- list()
	length(Obj) <- ntax
	i <- 1
	pos <- 0
	tot.nchar <- 0
	tot.ntax <- 0
	
	for (j in start.reading:NROW(X))
		{
		Xj <- trim.semicolon(X[j])
		if (Xj == "")
			{
			break
			}
		if (any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, 
			ignore.case = TRUE)))
			{
			break
			}
		ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
		if (length(ts) > 2)
			{
			print("printing ts with length > 2:")
			print(ts)
			stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
			}
		if (length(ts) != 2)
			{
			print("printing ts with length != 2:")
			print(ts)
			stop("nexus parser failed to read the sequences (ts!=2)")
			}
		Seq <- trim.whitespace(ts[2])
		Name <- trim.whitespace(ts[1])
		nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
		if (any(l <- grep(nAME, names(Obj))))
			{
			tsp <- strsplit(Seq, NULL)[[1]]
			for (k in 1:length(tsp))
			{
				p <- k + pos
				Obj[[l]][p] <- tsp[k]
				chars.done <- k
				}
			}
		else
			{
			names(Obj)[i] <- Name
			tsp <- strsplit(Seq, NULL)[[1]]
			for (k in 1:length(tsp))
			{
				p <- k + pos
				Obj[[i]][p] <- tsp[k]
				chars.done <- k
				}
			}
		tot.ntax <- tot.ntax + 1
		if (tot.ntax == ntax)
			{
			i <- 1
			tot.ntax <- 0
			tot.nchar <- tot.nchar + chars.done
			if (tot.nchar == nchar * ntax)
			{
				print("ntot was more than nchar*ntax")
				break
				}
			pos <- tot.nchar
			}
		else
			{
			i <- i + 1
			}
		}

	if (tot.ntax != 0)
		{
		cat("ntax:", ntax, "differ from actual number of taxa in file?\n")
		stop("nexus parser did not read names correctly (tot.ntax!=0)")
		}
	for (i in 1:length(Obj))
		{
		if (length(Obj[[i]]) != nchar)
			{
			cat(names(Obj[i]), "has", length(Obj[[i]]), "characters\n")
			stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
			}
		}
	Obj <- lapply(Obj, tolower)
	Obj
	}








# List the unique characters at every column
list_unique_chars <- function(tmpdf2, startchar=1, endchar=ncol(tmpdf2), printall=TRUE)
	{
	numcols = endchar - startchar + 1
	uniq_chars_list = as.list(rep(NA, numcols))
	
	for (i in 1:nrow(tempdf2))
		{
		chars = sort(unique(c(unlist(tempdf2[i,startchar:endchar]))))
	
		# Print to screen if desired
		if (printall == TRUE)
			{
			l = length(chars)
			c = paste(chars, sep="", collapse=",")
			charstxt = paste(i, ") #unique=", l, ", chars: ", c, sep="")
			cat(charstxt, "\n", sep="")
			} # END if (printall == TRUE)
			
		# Store the unique characters
		uniq_chars_list[[i]] = chars
		} # END for (i in 1:nrow(tempdf2))
	
	return(uniq_chars_list)
	} # END list_unique_chars





# For Debugging NEXUS input -- see TNT R stuff 
# for better functions.
defaults='
file=nexus_fn
check_ambig_chars=TRUE

file=infn
check_ambig_chars=TRUE
convert_ambiguous_to=NULL
printall="short"
convert_ambiguous_to_IUPAC=FALSE
'
read_nexus_data2 <- function(file, check_ambig_chars=TRUE, convert_ambiguous_to=NULL, printall="short", convert_ambiguous_to_IUPAC=FALSE) 
	{
	
	
	if (printall != "none")
		{
		cat("\n\nRunning read_nexus_data2(). This function modifies APE's read.nexus.data() in order to\nsuccessfully parse ambiguous morphological characters and other issues.  It might not \nwork perfectly, as NEXUS is a very complex format. Also, just because it reads the file\ndoes not mean that functions in other packages will magically be able to handle\nyour weird morphology data.\n")

		cat("\nFor issues, please email the maintainer, Nick Matzke, at:  matzke@nimbios.org\n")
		
		cat("\nTIPS:")
		cat("\n- You should use Mesquite to export your data to 'simplified\n  NEXUS' format before using this function.")
		cat("\n- Taxon/OTU names should include *NO* spaces or ' characters!")		
		cat("\n- Filenames should also have *NO* spaces or ' characters!")		
		cat("\n- If convert_ambiguous_to==TRUE (default), the script will check for\n  spaces in the character rows in the data matrix, e.g. due to '(0 1)'.", sep="")
		cat("\n- But your names still need to have no spaces, no apostrophes (')...\n\n", sep="")
		cat("\n\n====================================================================\n", sep="")
		cat("Processing ", file, "\n", sep="")
		cat("(notes may follow; set printall='none' to turn off)\n", sep="")
		cat("====================================================================\n\n", sep="")
		} # END if (printall != "none")


	#######################################################
	# If desired, use GREP to convert all ambiguous characters to their IUPAC codes, 
	# according to:
	# http://www.boekhoff.info/?pid=data&dat=fasta-codes
	#######################################################
	if (convert_ambiguous_to_IUPAC == TRUE)
		{
		strs = NULL
		sn = 0
		
		prefix = get_fn_prefix(file) 
		suffix = get_fn_suffix(file) 
		new_fn = paste0(prefix, "_noAmbig", ".", suffix)
		file = new_fn
		
		# After running this, you (SHOULDN'T!) need to check ambiguous
		# characters in DNA
		check_ambig_chars = FALSE
		
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\)/M/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\)/R/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ T\\)/W/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\)/S/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ T\\)/Y/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(G\\ T\\)/K/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\)/V/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ T\\)/H/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\ T\\)/D/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\ T\\)/B/g;" -pi $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\ T\\)/N/g;" -pi $(find ', file, ' -type f)', sep="")
		
		if (printall != "none")
			{
			cat("\n\n")
			cat("read_nexus_data2() is running perl to convert DNA from non-IUPAC to IUPAC...\ne.g. '(A C)' to 'M'...\n\n")
			} # END if (printall != "none")
		
		for (i in 1:sn)
			{
			if (printall != "none")
				{
				cat(paste("Running: ", strs[[i]], "...\n", sep=""))
				}
			system(strs[[i]])
			}

		} # END if (convert_ambiguous_to_IUPAC == TRUE)
	
	if (printall != 'none')
		{
		cat("\n\nReading in data...\n\n")
		}
	
	# Find the number of taxa
    "find.ntax" <- function(x)
    	{
        for (i in 1:NROW(x))
        	{
            if (any(f <- grep("\\bntax", x[i], ignore.case = TRUE)))
            	{
                ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
                  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            	} # END if (any(f <- grep("\\bntax", x[i], ignore.case = TRUE)))
        	} # END for (i in 1:NROW(x))
        return(ntax)
    	} # END "find.ntax" <- function(x)

	# Find the number of characters
    "find.nchar" <- function(x)
    	{
        for (i in 1:NROW(x))
        	{
            if (any(f <- grep("\\bnchar", x[i], ignore.case = TRUE)))
            	{
                nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)", 
                  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
				} # END if (any(f <- grep("\\bnchar", x[i], ignore.case = TRUE)))
			} # END  for (i in 1:NROW(x))
		return(nchar)
		}

	# Find the data matrix startline
    "find.matrix.line" <- function(x)
    	{
        for (i in 1:NROW(x))
        	{
            if (any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE)))
            	{
                matrix.line <- as.numeric(i)
                break
				} # END if (any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE)))
			} # END for (i in 1:NROW(x))
        matrix.line
	    } # END "find.matrix.line" <- function(x)
    
    
    "trim.whitespace" <- function(x)
    	{
        gsub("\\s+", "", x)
    	}
    
    "trim.semicolon" <- function(x)
    	{
        gsub(";", "", x)
    	}
    
    # Check if you can find the file
    if (file.access(file, mode = 4))
    	{
        stop("file could not be found")
    	}
    
    # Scan in the datafile
    X <- scan(file=file, what=character(), sep = "\n", quiet = TRUE, 
        comment.char = "[", strip.white = TRUE)
    ntax <- find.ntax(X)
    nchar <- find.nchar(X)
    matrix.line <- find.matrix.line(X)
    start.reading <- matrix.line + 1
    Obj <- list()
    length(Obj) <- ntax
    i <- 1
    pos <- 0
    tot.nchar <- 0
    tot.ntax <- 0
    





	# Some lines that have caused problems in the past.
	bugcheck = '
	Xj = "y2008Wola 2(0 1)002???91011101001002100511?200?1??01?11?0?30030100000100?1??1???10201?10?1?10?1030101201000?10?1???10?1?1?1?????0???1??00&00201?"
	# What if you have ambiguous characters with more than one state?
	Xj = "Gafrarium_tumidum_FMNH307858+FMNH312294+IM200741731+IM200741732+IM200741733 0000(1 2)0000(0 1)001000022110(1 3 4)0(0 2 3)0(0 2)13000-01(0 2)11(0 2)1120-10010101011210000100??????(0 1 2)"
	Xj = "H._ergaster         (0 1)(0 1)(0 1)(0 1)???????00110111111011110?101101??????????????????????????1????1???1??????11??1???1?1?????????????0??????????0????0?(0 1)????11?????????????????????????????????????????????????0??????????????????????????????????????1?0??10(0 1)10?1?00??0??01?????00100010???0001(0 1)110000001011(0 1)(0 1)01001?0(0 1)110?11100????????????0100(0 1)(0 1)0101100001000001001011011000011(0 1)00(0 1)(0 1)01(0 1)000(0 1)(0 1)(0 1)10(0 1)1(0 1)10(0 1)(0 1)(0 1)0112?(1 2)???0211(0 1)0100(0 1 2)(0 1)0(1 2)(1 2)(0 2)2?00013200?22101?????2???????????????????????????????1?21002112???2122????1?1221102?021222(1 2)2121(0 1)1(0 1)21122?(0 1 2)(0 1)(1 2)0(0 1)0(0 1 2)02(0 1)0(0 1 2)(1 2)(1 2)12331021(0 1)0(0 2)13(1 2)11"
	
	# Squiggely brackets:
	Xj = "Byronosaurus_jaffei                 ?????101???101?1100110101011?00??????20120????????????1?100??????0000001??11????00021??01?0?0???????010110????????0??02??????????????????????????????????????????????????????????????1????0?0??0????????1??????????0????0????100???0?{23}0?100?101100??0100???10?00110?????????????????????????????????????????????????????????????11????????????????????0??11?????0??????????00?00??01?1"

	ts <- strsplit(Xj, "[ \t]+")[[1]]
	Name = ts[1]
	
	tmpSeq = list2str_fast_nosplit(ts[-1], spacer=" ")
	'
	
	
	# Go through the characters, check for ambiguities if desired
    if ((check_ambig_chars == TRUE) && (printall != "none"))
    	{
    	cat("\n\nread_nexus_data2() is checking for ambiguous characters in line...\n", sep="")
    	}
    
    for (j in start.reading:NROW(X))
    	{
    	if ( (check_ambig_chars == TRUE) && (printall != "none") )
    		{
		    cat(j, ", ", sep="")
		    }
        Xj <- trim.semicolon(X[j])
        if (Xj == "")
        	{
            break
        	}
        # Stop once you hit "END"
        if (any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, ignore.case = TRUE)))
            {
            break
	        }
        
        # Convert character string line with e.g. (0 1) to Z (or whatever)
        # count the number of spaces
        # number_of_spaces = count_chars(char_to_count=" ", Xj)
        
	
		
		need_to_correct_ambiguous = FALSE
        if (check_ambig_chars == TRUE)
        	{
        	# split on whitespace
        	ts <- strsplit(Xj, "[ \t]+")[[1]]
        	Name = ts[1]
        	
        	# Take everything that wasn't the name and merge back into old string
        	tmpSeq = list2str_fast_nosplit(ts[-1], spacer=" ")
        	tmpSeq_orig = tmpSeq
        	
        	
        	# Convert any "{" or "}" to "(" and ")"
        	if (grepl(pattern="\\{", x=tmpSeq) == TRUE)
        		{
				tmpSeq = gsub(pattern="\\{", replacement="(", x=tmpSeq)
				tmpSeq = gsub(pattern="\\}", replacement=")", x=tmpSeq)
				}
        	
        	# convert the ambiguous characters to "~" for later processing
        	# locate e.g. (0 1), (0 1 2), etc...
        	ambig_chars_paren_locations = gregexpr('(\\(.*?\\))', tmpSeq, perl=TRUE)[[1]]
			
			# If there ARE ambiguous characters, convert them to ~
			# otherwise, don't
			ambig_chars_TF = FALSE
			if ( length(ambig_chars_paren_locations) > 1)
				{
				ambig_chars_TF = TRUE
				} else {
				if (ambig_chars_paren_locations != -1)
					{
					ambig_chars_TF = TRUE
					} # END if (ambig_chars_paren_locations != -1)
				} # END if ( length(ambig_chars_paren_locations) > 1)
			
			if (ambig_chars_TF == TRUE)
				{
				need_to_correct_ambiguous = TRUE
				ambig_chars_paren_list = NULL
				
				# Count down backwards, to avoid the changing of the sequence length
				for (k in length(ambig_chars_paren_locations):1)
					{
					
					startnum = ambig_chars_paren_locations[k]
					ambig_length = attr(ambig_chars_paren_locations, "match.length")[k]
					endnum = startnum + ambig_length - 1
					tmpSeq_split = strsplit(tmpSeq, split="")[[1]]

					# print(k)
					# print(startnum)
					# print(ambig_length)
					# print(endnum)
					# print(tmpSeq_split)

					
					# Do you want to put the full ambiguous coding into the NEXUS
					# object, or replace with e.g. a question mark ("?") ?
					if (is.null(convert_ambiguous_to))
						{
						# Put in original ambiguous coding into NEXUS object
						ambig_charstring = list2str_fast(tmpSeq_split[startnum:endnum])
						
						# Check for ambig_charstring of length 4,
						# e.g. "(23)"
						# Things SHOULD have spaces, but maybe
						# {23} doesn't need them
						if (nchar(ambig_charstring) == 4)
							{
							ambig_charstring_chars = strsplit(ambig_charstring, split="")[[1]]
							ambig_charstring_chars
							ambig_charstring_chars2 = c(ambig_charstring_chars[1:2], " ", ambig_charstring_chars[3:4])
							ambig_charstring = list2str_fast(ambig_charstring_chars2, spacer="")
							}
						
						
						tmprow = c(startnum, endnum, ambig_length, ambig_charstring)
						} else {
						# Put ? into NEXUS object
						tmprow = c(startnum, endnum, ambig_length, "?")						
						}
			
					# ambiguous character = ~
					#print(tmpSeq)
					#print(ambig_chars_paren_locations)
					#cat("\n", startnum, endnum, startnum-1, "\n")
					
					
					# If you are NOT at the end of the string -- catches where end of string is e.g. (0 1)
					if (endnum < length(tmpSeq_split))
						{
						# This catches the case where the first character is e.g. (0 1)
						if (startnum-1 > 0)
							{
							newSeq = list2str_fast(c(tmpSeq_split[1:(startnum-1)], "~", tmpSeq_split[(endnum+1):length(tmpSeq_split)] ))
							} else {
							newSeq = list2str_fast(c("~", tmpSeq_split[(endnum+1):length(tmpSeq_split)] ))
							}
						} else {
						# If you ARE at the end of the string: -- catches where end of string is e.g. (0 1)
						# This catches the case where the first character is e.g. (0 1)
						if (startnum-1 > 0)
							{
							newSeq = list2str_fast(c(tmpSeq_split[1:(startnum-1)], "~") )
							} else {
							newSeq = list2str_fast(c("", "~") )
							}
						}
	
					# add to the list of ambiguous characters
					ambig_chars_paren_list = rbind(ambig_chars_paren_list, tmprow)
					
					tmpSeq = newSeq
					}
				
				# We will have to store this in a permanent place somewhere
				ambig_chars_paren_list = adf(ambig_chars_paren_list)
				names(ambig_chars_paren_list) = c("startnum", "endnum", "ambig_length", "ambig_charstring")
        		}
        	# locate e.g. [0 1], [0 1 2], etc...
        	#ambig_locations2 = gregexpr('(\\[(.*?)\\])', tmpSeq, perl=TRUE)[[1]]

        	# locate e.g. 0&1, etc...
        	#ambig_locations3 = gregexpr('((.&.))', tmpSeq, perl=TRUE)[[1]]
        	
        	# store the revised sequence in Seq
        	Seq = tmpSeq
        	}
        else
        	{
			ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
			if (length(ts) > 2)
				{
				cat("\n\nERROR: length of ts (data string, name-space-characters) is >2.  printing ts:\n", sep="")
				print(ts)
				stop("\nOriginal error message: nexus parser does not handle spaces in sequences or taxon names (ts>2)")
				}
			if (length(ts) != 2)
				{
				cat("\n\nERROR: length of ts (data string, name-space-characters) is !=2.  printing ts:\n", sep="")
				print(ts)
				stop("\nOriginal error message: nexus parser failed to read the sequences (ts!=2)")
				}
			Seq <- trim.whitespace(ts[2])
			Name <- trim.whitespace(ts[1])
			}
        nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
		
		# Here, this if() didnt run on a standard simplified MrBayes NEXUS file
		if (any(l <- grep(nAME, names(Obj))))
        	{
            tsp <- strsplit(Seq, NULL)[[1]]
            
			# Correct for ambiguous sequences
            if (need_to_correct_ambiguous == TRUE)
            	{
				tsp[tsp=="~"] = ambig_chars_paren_list$ambig_charstring

				# Tell the user what happened if ambiguous characters are detected
				if (is.null(convert_ambiguous_to))
					{
					if (printall != "none")
						{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", tmpSeq_orig, "\n", sep="")
							} # END if (printall == "short")
						} # END if (printall != "none")
					} else {
					if (printall != "none")
						{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "'\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "':", "\n", tmpSeq_orig, " --> \n", sep="")
							cat(list2str_fast_nosplit(tsp), "\n", sep="")				
							} # END if (printall == "short")
						} # END if (printall != "none")
					} # END if (is.null(convert_ambiguous_to))
				}
            
            # input the characters into the object
            for (k in 1:length(tsp))
            	{
                p <- k + pos
                Obj[[l]][p] <- tsp[k]
                chars.done <- k
	            }
    	    }
        else {
			# Here, this else() DID run on a standard simplified MrBayes NEXUS file
            names(Obj)[i] <- Name
            tsp <- strsplit(Seq, NULL)[[1]]

			# Correct for ambiguous sequences
            if (need_to_correct_ambiguous == TRUE)
            	{
				tsp[tsp=="~"] = ambig_chars_paren_list$ambig_charstring

				# Tell the user what happened if ambiguous characters are detected
				if (is.null(convert_ambiguous_to))
					{
				    if (printall != "none")
    					{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", tmpSeq_orig, "\n", sep="")
							}
						} # END if (printall != "none")
					} else {
					if (printall == "short")
						{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "'\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "':", "\n", tmpSeq_orig, " --> \n", sep="")
							cat(list2str_fast_nosplit(tsp), "\n", sep="")				
							}
						} # END if (printall != "none")
					}
				}
            
            for (k in 1:length(tsp))
            	{
                p <- k + pos
                Obj[[i]][p] <- tsp[k]
                chars.done <- k
            	}
	        }
        tot.ntax <- tot.ntax + 1
        if (tot.ntax == ntax) {
            i <- 1
            tot.ntax <- 0
            tot.nchar <- tot.nchar + chars.done
            
            # Make sure the total number of characters = total # of taxa * # of characters per
            if (tot.nchar == nchar * ntax) {
                print("Warning: ntot was more than nchar*ntax")
                break
            }
            pos <- tot.nchar
        }
        else {
            i <- i + 1
        }
    }
    if (tot.ntax != 0) {
        cat("ntax:", ntax, "differ from actual number of taxa in file?\n")
        stop("nexus parser did not read names correctly (tot.ntax!=0)")
    }
    for (i in 1:length(Obj)) {
        if (length(Obj[[i]]) != nchar) {
            cat(names(Obj[i]), "has", length(Obj[[i]]), "characters\n")
            stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
        }
    }
    Obj <- lapply(Obj, tolower)
    
    if (printall != "none")
    	{
		cat("\n...read_nexus_data2() done.\n\n")
		} # END if (printall != "none")
    return(Obj)
	} # END read_nexus_data2




# Check for all ? columns
check_for_all_Qs <- function(nexd, printcounts=FALSE)
	{
	
	rownums = 1:nrow(nexd)
	lengths = apply(nexd, 2, length)
	numQs = apply(nexd, 2, countQs)
	allQs_TF = lengths == numQs
	
	countQs_df = cbind(rownums, lengths, numQs, allQs_TF)
	countQs_df = adf2(countQs_df)
	
	return(countQs_df)
	}

# Count question marks
countQs <- function(tmprow)
	{
	numQs = sum(tmprow == "?", na.rm=TRUE)
	return(numQs)
	}



# Pull out and code observed ambiguous states for each collection of
# observed character states
# for morphology only
# Hannah's data: 69 morphological characters
# return_missing_chars = "list", "correct", or "none"
# (for listing the sites with missing character states, correcting them, or neither)
get_numstates_per_char <- function(nexd, ambig_to_remove=c("\\(", "\\)", " ", ","), return_missing_chars="list", printall="short")
	{
	# Make a list of all ? (allQs)
	ntaxa = ncol(nexd)
	allQs_list = rep(FALSE, times=nrow(nexd))
	
	# Convert data to data.frame
	nexdf = adf(nexd)
	

	#column_strings_morph = NULL
	numstates_morph_list = NULL

	if (printall != "none")
		{
		cat("\n\nget_numstates_per_char() is checking for missing charstates:\n", sep="")
		} # END if (printall != "none")
	missing_charstates_list = NULL
	tmp_charstate_codes = charstate_codes()
	missing_charstates = FALSE
	
	
	# Go through each column, count the number of states
	# (but remove "(", ")" )
	for (rownum in 1:nrow(nexdf))
		{
		# Check for all "?" characters
		if (countQs(nexdf[rownum,]) == ntaxa)
			{
			if (printall != "none")
				{
				cat("\n\nWARNING: Data column #", rownum, " is ALL ? CHARACTERS, i.e. all missing data. PLEASE CORRECT.\n", sep="")
				cat("  (...adding to list 'allQs_colnums'...)   \n", sep="")
				} # END if (printall != "none")
			# Add to list of allQs rows
			allQs_list[rownum] = TRUE
			}
		
		
		# is the column a "standard" character (contains some digits, i.e. "0s" (or 1s!) == "\\d" covers all digits)
		standard_TF = sum(grepl("\\d", nexdf[rownum,])) > 0
		
		#cat("\n")
		if (standard_TF != TRUE)
			{
			cat("\nERROR: rownum that is not a character of type 'standard': ", rownum, "\n", sep="")
			nexdf[rownum,]
			numstates = 0
			numstates_morph_list = c(numstates_morph_list, numstates)
			
			} else {
			# are there any ambiguous characters?
			ambiguous_TF = grepl("\\(", nexdf[rownum,])
			
			# If there are ambiguous characters, convert them 
			# from e.g. (0 1) to 01
			if (sum(ambiguous_TF) > 0)
				{
				ambig_chars = nexdf[rownum,][ambiguous_TF]
				
				for (i in 1:length(ambig_to_remove))
					{
					# replace "(", ")", " "
					#gsub("\\(", "", c("(0 1)", "(0 1)", "(0 1)"))
					ambig_chars = gsub(ambig_to_remove[i], "", ambig_chars)
					ambig_chars = gsub(ambig_to_remove[i], "", ambig_chars)
					ambig_chars = gsub(ambig_to_remove[i], "", ambig_chars)
					}
				
				nexdf[rownum,][ambiguous_TF] = ambig_chars
				}
			# add these to the list
			unique_characters = unique(strsplit(list2str_fast(nexdf[rownum,]), split="")[[1]])
			
			# drop "-" and "?"
			unique_characters = unique_characters[unique_characters != "-"]
			unique_characters = unique_characters[unique_characters != "?"]
			numstates = length(unique_characters)
			
			# Check for characters missing states (the highest character
			# should equal the one you get from counting up states)
			unique_characters = unique_characters[order(unique_characters)]

			if ( (numstates > 1) && (unique_characters[length(unique_characters)] != tmp_charstate_codes[length(unique_characters)]) )
				{
				if (printall != "none")
					{
					cat("\n\nWarning: in rownum #", rownum, " the maximum charstate is\nnot the same as the maximum charstate derived by counting up charstates...\nmax charstate: ", unique_characters[length(unique_characters)], "\nmax code: ", tmp_charstate_codes[length(unique_characters)], "\n", sep="")
					print(unique_characters)
					#print(tmp_charstate_codes)
					cat("Old row #", rownum, "\n", sep="")
					print(as.character(nexdf[rownum,]))
					} # END if (printall != "none")

				
				missing_charstates_list = c(missing_charstates_list, rownum)
				missing_charstates = TRUE
				if (printall != "none")
					{
					cat(rownum, ", ", sep="")
					} # END if (printall != "none")

				
				# If you want to correct these, do it here
				if (return_missing_chars == "correct")
					{
					#cat("correcting, ", sep="")
					TFs_to_change_list = NULL
					TFs_to_change_num = 0
					unique_chars_to_change_list = NULL
					new_character_values_list = NULL
					
					for (j in 1:numstates)
						{
						unique_char = unique_characters[j]
						if (unique_char != tmp_charstate_codes[j])
							{
							unique_chars_to_change_list = c(unique_chars_to_change_list, unique_char)
							new_character_values_list = c(new_character_values_list, tmp_charstate_codes[j])
							TFs_to_change_list[[(TFs_to_change_num=TFs_to_change_num+1)]] = grepl(unique_char, nexdf[rownum,])
							} # END if (unique_char != tmp_charstate_codes[j])
						} # END for (j in 1:numstates)
					# Now fix these, ONLY in the appropriate cells
					
					if (printall != "none")
						{
						cat("\n...correcting: ", rownum, "...\n", sep="")
						} # END if (printall != "none")
					for (j in 1:length(unique_chars_to_change_list))
						{
						unique_char = unique_chars_to_change_list[j]
						TF_to_change = TFs_to_change_list[[j]] 
						nexdf[rownum, ][TF_to_change] = gsub(unique_char, new_character_values_list[j], nexdf[rownum,][TF_to_change])
						} # END for (j in 1:length(unique_chars_to_change_list))
					if (printall != "none")
						{
						cat("New row #", rownum, "\n", sep="")
						} # END if (printall != "none")

					if (printall != "none")
						{
						print(as.character(nexdf[rownum,]))
						} # END if (printall != "none")
					} # END if (return_missing_chars == "correct")
				} # END if ( (numstates > 1) && (unique_characters...
			numstates_morph_list = c(numstates_morph_list, numstates)			
			} # END if (standard_TF != TRUE)
		# end forloop
		} #END for (rownum in 1:nrow(nexdf))
		
	# If missing character states were detected, 
	# print "...done. No missing character states detected."
	if (missing_charstates == FALSE)
		{
		if (printall != "none")
			{
			cat("...done. No missing character states detected.\n", sep="")
			cat("\n")
			} # END if (printall != "none")
		} # END if (missing_charstates == FALSE)
		

	# Also return the columns that have missing states, if desired.
	if (return_missing_chars == "list")
		{
		results = NULL
		results$numstates_morph_list = numstates_morph_list
		results$missing_charstates_list = missing_charstates_list
		results$allQs_list = allQs_list
		return(results)
		}

	if (return_missing_chars == "correct")
		{
		results = NULL
		results$nexdf = nexdf
		results$numstates_morph_list = numstates_morph_list
		results$missing_charstates_list = missing_charstates_list
		results$allQs_list = allQs_list
		return(results)
		}
	
	retrieve_cmds='
	numstates_morph_list = res$numstates_morph_list
	missing_charstates_list = res$missing_charstates_list
	allQs_list = res$allQs_list
	morph_df2_corrected = res$nexdf
	'
			
	return(numstates_morph_list)
	}





# Ambiguous single-letter codes for morphology
charstate_codes <- function(max_numstates=32)
	{
	# TNT can do 32 character states
	statenums = 1:max_numstates
	
	alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
	
	alphabet = LETTERS
	
	statecodes = c(0:9, LETTERS[1:22])
	return(statecodes)
	}



# Ambiguous single-letter codes for morphology
charstate_code_combinations_all <- function(max_numstates=5)
	{
	# TNT can do 32 character states
	statenums = 1:max_numstates
	
	alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
	
	alphabet = LETTERS
	
	statecodes = c(0:9, LETTERS[1:22])
	
	# the number of states for which ALL ambiguities are TOTALLY available is only 5
	# (5 states = 
	# length(charstate_code_combinations_all(5))
	# = 31 )
	#
	# (6 states = 
	# length(charstate_code_combinations_all(6))
	# = 63 )
	
	combinations_list = NULL
	current_statecodes = statecodes[1:max_numstates]
	for (i in 1:max_numstates)
		{
		tmp_combinations = combinations(length(current_statecodes), i, current_statecodes)
		
		# For single draws
		if ( ncol(tmp_combinations) == 1 )
			{
			combinations_list = c(combinations_list, c(tmp_combinations))	
			} else {
			# For multiple draws
			tmp_combinations = apply(tmp_combinations, 1, paste, collapse="")
			combinations_list = c(combinations_list, tmp_combinations)
			}
		}
	return(combinations_list)
	}












###############################################
# Tree plotting
###############################################

plot_cont_vals_on_tree <- function(chtr2, tmptipvals, tmpnodevals, texttxt)
	{
	# Plot some continuous values on a tree in pie chart form
	#
	tree_ht = get_max_height_tree(chtr2)
	numtips = length(chtr2$tip.label)
	numnodes = chtr2$Nnode
	label_offset = 0.05*tree_ht
	plot(chtr2, show.node.label=FALSE, cex=0.9, no.margin = TRUE, x.lim=1.6*tree_ht, label.offset=label_offset)
	
	param_vals = c(tmptipvals, tmpnodevals)
	maxval = max(param_vals)
	param_vals_normed = param_vals / maxval
	
	tipvals = param_vals_normed[1:numtips]
	nodevals = param_vals_normed[(numtips+1) : (numtips+numnodes)]
	tmp_colors = c("darkgray", "white")
	tiplabels(pch = ".", pie=tipvals, piecol=tmp_colors, cex = 1) # adj = 0)
	nodelabels(pie=nodevals, cex=1, piecol=tmp_colors)
	
	legend(0, 3, legend=c("fraction of max value"), fill=tmp_colors)
	text(0, 4, labels=texttxt, pos=4)
	
	}


defaults = '
x = elev
phy = tr3
xaxt = "s"
underscore = FALSE
show.names = TRUE
show.xaxis.values = TRUE
method = "pic"
'

traitgram2 <- function (x, phy, xaxt = "s", underscore = FALSE, show.names = TRUE, show.xaxis.values = TRUE, method = c("ace", "pic"), ...) 
    {
    require(ape)
    Ntaxa = length(phy$tip.label)
    Ntot = Ntaxa + phy$Nnode
    phy = node.age(phy)
    ages = phy$ages[match(1:Ntot, phy$edge[, 2])]
    ages[Ntaxa + 1] = 0

    if (class(x) %in% c("matrix", "array"))
    	{
        xx = as.numeric(x)
        names(xx) = row.names(x)
        } else {
    	xx = x
    	}

	# Check to make sure the data have names
	if (NA %in% names(xx) || is.null(names(xx)))
		{
		cat("traitgram2 WARNING: applying phylogeny tipnames to data, they\nbetter be in the same order as the tips!\n")
		names(xx) = phy$tip.label
		}


    if (!is.null(names(xx)))
    	{
        umar = 0.1
        if (!all(names(xx) %in% phy$tip.label))
        	{
            print("trait and phy names do not match")
            return()
	        }
        xx = xx[match(phy$tip.label, names(xx))]
	    } else {
	    umar = 0.1
	    }
    lmar = 0.2
    if (xaxt == "s") 
        if (show.xaxis.values)
        	{
            lmar = 1
            } else
            {
            lmar = 0.5
            }
    if (method[1] == "ace")
    	{
        xanc = ace(xx, phy)$ace
        } else {
        xanc = pic3(xx, phy)[, 3]
        }
    xall = c(xx, xanc)
    a0 = ages[phy$edge[, 1]]
    a1 = ages[phy$edge[, 2]]
    x0 = xall[phy$edge[, 1]]
    x1 = xall[phy$edge[, 2]]
    tg = par(bty = "n", mai = c(lmar, 0.1, umar, 0.1))
    if (show.names)
    	{
        maxNameLength = max(nchar(names(xx)))
        ylim = c(min(ages), max(ages) * (1 + maxNameLength/50))
        if (!underscore)
        	{
            names(xx) = gsub("_", " ", names(xx))
            }
	    } else {
	    ylim = range(ages)
	    }
    plot(range(c(x0, x1)), range(c(a0, a1)), type = "n", xaxt = "n", 
        yaxt = "n", xlab = "", ylab = "", bty = "n", ylim = ylim, 
        cex.axis = 0.8)
    if (xaxt == "s") 
        if (show.xaxis.values) 
            axis(1, labels = TRUE)
        else axis(1, labels = FALSE)
    segments(x0, a0, x1, a1)
    if (show.names)
    	{
        text(sort(xx), max(ages), labels = names(xx)[order(xx)], 
            adj = -0, srt = 90)
	    }
    on.exit(par(tg))
	}


# write.nexus2 == unnecessary; just make sure you use the FILE= argument to distinguish the trees object
# from the filename object








# Just the standard plot.phylo (for phylo3 APE objects)
# In addition to the standard outputs, this function returns xx and yy,
# the plotting coordinates for those nodes.
plot_phylo3_nodecoords <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = "horizontal", tip.color = "black", 
    plot = TRUE, rotate.tree = 0, ...) 
{
    Ntip <- length(x$tip.label)
    if (Ntip == 1) {
        warning("found only one tip in the tree")
        return(NULL)
    }
    if (any(tabulate(x$edge[, 1]) == 1)) 
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C("node_height", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
        DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C("node_depth", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
        edge.length) .C("node_depth_edgelength", as.integer(Ntip), 
        as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
            2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
        "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
        "upwards", "downwards"))
    if (is.null(x$edge.length)) 
        use.edge.length <- FALSE
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
        root.edge <- FALSE
    if (type == "fan" && root.edge) {
        warning("drawing root edge with type = 'fan' is not yet supported")
        root.edge <- FALSE
    }
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge
    if (phyloORclado) {
        phyOrder <- attr(x, "order")
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x)
            if (!identical(x$edge, xe)) {
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
            }
        }
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = "pruningwise")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) 
                node.pos <- 2
        }
        if (node.pos == 1) 
            yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
        else {
            ans <- .C("node_height_clado", as.integer(Ntip), 
                as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 
                  2]), as.integer(Nedge), double(Ntip + Nnode), 
                as.double(yy), DUP = FALSE, PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
                  1
            xx <- max(xx) - xx
        }
        else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                z$edge.length)
        }
    }
    else {
        rotate.tree <- 2 * pi * rotate.tree/360
        switch(type, fan = {
            TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
            xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
            theta <- double(Ntip)
            theta[TIPS] <- xx
            theta <- c(theta, numeric(Nnode))
            theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, 
                theta)
            if (use.edge.length) {
                r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                  Nedge, z$edge.length)
            } else {
                r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
                r <- 1/r
            }
            theta <- theta + rotate.tree
            xx <- r * cos(theta)
            yy <- r * sin(theta)
        }, unrooted = {
            nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
            xx <- XY$M[, 1] - min(XY$M[, 1])
            yy <- XY$M[, 2] - min(XY$M[, 2])
        }, radial = {
            X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            X[X == 1] <- 0
            X <- 1 - X/Ntip
            yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
            Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
            xx <- X * cos(Y + rotate.tree)
            yy <- X * sin(Y + rotate.tree)
        })
    }
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") 
                xx <- xx + x$root.edge
            if (direction == "upwards") 
                yy <- yy + x$root.edge
        }
    }
    if (no.margin) 
        par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par("pin")[1]
                strWi <- strwidth(x$tip.label, "inches")
                xx.tips <- xx[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * xx.tips + 
                  strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(xx.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(xx.tips + strWi/alp)
                  else max(xx.tips)
                }
                x.lim[2] <- tmp
            }
            else x.lim <- c(1, Ntip)
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
            } else x.lim <- c(min(xx), max(xx))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            } else x.lim <- c(0, max(xx))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            } else x.lim <- c(-1, 1)
        })
    }
    else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            x.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.03 * cex)
            else -1
    }
    if (phyloORclado && direction == "leftwards") 
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
                y.lim <- c(1, Ntip)
            else {
                y.lim <- c(0, NA)
                pin2 <- par("pin")[2]
                strWi <- strwidth(x$tip.label, "inches")
                yy.tips <- yy[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * yy.tips + 
                  strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(yy.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(yy.tips + strWi/alp)
                  else max(yy.tips)
                }
                y.lim[2] <- tmp
            }
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            } else y.lim <- c(min(yy), max(yy))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            } else y.lim <- c(0, max(yy))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    }
    else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            y.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
            else -1
    }
    if (phyloORclado && direction == "downwards") 
        yy <- max(yy) - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
            x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
            y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 
        1
    else NA
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, 
        axes = FALSE, asp = asp, ...)
    if (plot) {
        if (is.null(adj)) 
            adj <- if (phyloORclado && direction == "leftwards") 
                1
            else 0
        if (phyloORclado && show.tip.label) {
            MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
            loy <- 0
            if (direction == "rightwards") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == "leftwards") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                  adj)
            }
            if (!horizontal) {
                psr <- par("usr")
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                  psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == "downwards") {
                  loy <- -loy
                  srt <- 180 + srt
                }
            }
        }
        if (type == "phylogram") {
            phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                edge.color, edge.width, edge.lty)
        }
        else {
            if (type == "fan") {
                ereorder <- match(z$edge[, 2], x$edge[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
                circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                  r, edge.color, edge.width, edge.lty)
            }
            else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                edge.lty)
        }
        if (root.edge) 
            switch(direction, rightwards = segments(0, yy[ROOT], 
                x$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
                yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]), 
                upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge), 
                downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                  yy[ROOT] + x$root.edge))
        if (show.tip.label) {
            if (is.expression(x$tip.label)) 
                underscore <- TRUE
            if (!underscore) 
                x$tip.label <- gsub("_", " ", x$tip.label)
            if (phyloORclado) 
                text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
                  adj = adj, font = font, srt = srt, cex = cex, 
                  col = tip.color)
            if (type == "unrooted") {
                if (lab4ut == "horizontal") {
                  y.adj <- x.adj <- numeric(Ntip)
                  sel <- abs(XY$axe) > 0.75 * pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    1.05
                  sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
                    pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    (2 * abs(XY$axe)[sel]/pi - 0.5)
                  sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
                  y.adj[sel] <- strheight(x$tip.label)[sel]/2
                  sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
                  y.adj[sel] <- -strheight(x$tip.label)[sel] * 
                    0.75
                  text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + 
                    y.adj * cex, x$tip.label, adj = c(adj, 0), 
                    font = font, srt = srt, cex = cex, col = tip.color)
                }
                else {
                  adj <- abs(XY$axe) > pi/2
                  srt <- 180 * XY$axe/pi
                  srt[adj] <- srt[adj] - 180
                  adj <- as.numeric(adj)
                  xx.tips <- xx[1:Ntip]
                  yy.tips <- yy[1:Ntip]
                  if (label.offset) {
                    xx.tips <- xx.tips + label.offset * cos(XY$axe)
                    yy.tips <- yy.tips + label.offset * sin(XY$axe)
                  }
                  font <- rep(font, length.out = Ntip)
                  tip.color <- rep(tip.color, length.out = Ntip)
                  cex <- rep(cex, length.out = Ntip)
                  for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                    cex = cex[i], x$tip.label[i], adj = adj[i], 
                    font = font[i], srt = srt[i], col = tip.color[i])
                }
            }
            if (type %in% c("fan", "radial")) {
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                angle <- atan2(yy.tips, xx.tips)
                if (label.offset) {
                  xx.tips <- xx.tips + label.offset * cos(angle)
                  yy.tips <- yy.tips + label.offset * sin(angle)
                }
                s <- xx.tips < 0
                angle <- angle * 180/pi
                angle[s] <- angle[s] + 180
                adj <- as.numeric(s)
                font <- rep(font, length.out = Ntip)
                tip.color <- rep(tip.color, length.out = Ntip)
                cex <- rep(cex, length.out = Ntip)
                for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                  x$tip.label[i], font = font[i], cex = cex[i], 
                  srt = angle[i], adj = adj[i], col = tip.color[i])
            }
        }
        if (show.node.label) 
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
                x$node.label, adj = adj, font = font, srt = srt, 
                cex = cex)
    }
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, edge = xe, xx = xx, yy = yy)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
        envir = .PlotPhyloEnv)
    invisible(L)
}










# Run ClustalW
# https://stat.ethz.ch/pipermail/r-help/2007-June/134427.html
run_clustalw <- function(dna=NULL, infn="temp_unaligned.fasta", outfn="temp_aligned.fasta", clustalw_full_path="clustalw", options_txt="-outorder=input")
	{
	
	# If DNA is not null, write it to FASTA for input into clustalw
	if (dna != NULL)
		{
		write.dna(dna, file=infn, format="fasta")
		}
	
	cmdstr = paste(clustalw_full_path, " -infile=", infn, " -outfile=", outfn, " ", options_txt, sep="")
	system(cmdstr)
	
	return(outfn)
	}
