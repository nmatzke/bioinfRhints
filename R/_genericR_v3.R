# =============================================
# genericR_v1.R: many useful utility functions
#   for R
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
#   * to Share — to copy, distribute and transmit the work
#   * to Remix — to adapt the work
#
# Under the following conditions:
#
#   * Attribution — You must attribute the work in the manner 
#     specified by the author or licensor (but not in any way that 
#     suggests that they endorse you or your use of the work).
#   * Noncommercial — You may not use this work for commercial purposes. 
#
#   * Share Alike — If you alter, transform, or build upon this work,
#     you may distribute the resulting work only under the same or
#     similar license to this one. 
#
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# 
###################################################################

# Generic utility functions for R
# Load with:
# sourcedir = '/drives/Dropbox/_njm/'
# source3 = '_genericR_v1.R'
# source(paste(sourcedir, source3, sep=""))

# NOTE: MANY OF THESE FUNCTIONS ARE BEING MOVED TO:

#sourceall(path="/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

# for e.g. calc_loglike
# sourcedir = '/Dropbox/_njm/'
# source3 = '_R_tree_functions_v1.R'
# source(paste(sourcedir, source3, sep=""))





# source all .R files in a directory, except "compile" and "package" files
sourceall <- function(path=path, pattern="\\.R", ...)
	{
	Rfiles = list.files(path=path, pattern="\\.R", ...)
	
	# Files to remove
	Rfiles_remove_TF1 = grepl("compile", Rfiles)
	Rfiles_remove_TF2 = grepl("package", Rfiles)
	Rfiles_remove_TF = (Rfiles_remove_TF1 + Rfiles_remove_TF2) >= 1
	
	Rfiles = Rfiles[Rfiles_remove_TF == FALSE]

	cat("\nSourcing Rfiles in ", path, "...\n", sep="")

	
	for (Rfile in Rfiles)
		{
		cat("Sourcing Rfile: ", Rfile, "\n", sep="")
		fullfn = slashslash(paste(path, Rfile, sep=""))
		source(fullfn, chdir=TRUE)
		}

	cat("\nDone sourcing Rfiles in ", path, "...\n", sep="")
	return(path)
	}


# Get a list of directories
# http://stackoverflow.com/questions/4749783/how-to-obtain-a-list-of-directories-within-a-directory-like-list-files-but-i
list_dirs <- function(path=".", pattern=NULL, all.dirs=FALSE, full.names=FALSE, ignore.case=FALSE)
	{
	all <- list.files(path, pattern, all.dirs, full.names, recursive=FALSE, ignore.case)
	all[file.info(all)$isdir]
	}









#######################################################
# slashslash:
#######################################################
#' Remove double slash (slash a slash)
#' 
#' Shortcut for: \code{gsub(pattern="//", replacement="/", x=tmpstr)}
#' 
#' This function is useful for removing double slashes that can
#' appear in full pathnames due to inconsistencies in trailing
#' slashes in working directories etc.
#'
#' @param tmpstr a path that you want to remove double slashes from
#' @return outstr a string of the fixed path
#' @export
#' @seealso \code{\link[base]{getwd}}, \code{\link[base]{setwd}}, \code{\link[base]{gsub}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' tmpstr = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/modiscloud/extdata/2002_raw_MYD//MYD03.A2002185.0645.005.2009192031332.hdf"
#'
#' outstr = slashslash(tmpstr)
#' outstr
#'
slashslash <- function(tmpstr)
	{
	outstr = gsub(pattern="//", replacement="/", x=tmpstr)
	return(outstr)
	}
















##################################################################
# abbreviations
##################################################################
# Abbreviation of paste; default separation is "" instead of " "
p <- function(x, sep="", collapse=NULL)
	{
	return(paste(x, sep, collapse))
	}

# Return a row with lots of standard summary stats
summstats <- function(values, db="", descrip="")
	{
	# Calculate a variety of standard estimates
	num_estimates = length(values)
	minval = min(values)
	maxval = max(values)
	rangewidth = maxval - minval
	q025 = quantile(values, 0.025)
	q975 = quantile(values, 0.975)
	HPD95width = q975 - q025
	meanTT = mean(values, na.rm=TRUE)
	sdTT = sd(values, na.rm=TRUE)
	varTT = var(values, na.rm=TRUE)
	sd95_2tailed = 1.96*sdTT
	CIupper = meanTT + sd95_2tailed
	CIlower = meanTT - sd95_2tailed
	CIwidth = CIupper - CIlower
	
	db = db
	descrip = descrip
	
	tmprow = c(db, descrip, num_estimates, meanTT, sdTT, varTT, sd95_2tailed, CIlower, CIupper, CIwidth, q025, q975, HPD95width, minval, maxval, rangewidth)
	
	return(tmprow)
	}

# Return the names of the columns in summstats
summstats_colnames <- function()
	{
	tmpnames = c("db", "timerange", "num_estimates", "meanTT", "sdTT", "varTT", "sd95_2tailed", "CIlower", "CIupper", "CIwidth", "q025", "q975", "HPD95width", "minval", "maxval", "rangewidth")
	
	return(tmpnames)
	}


# Print the first numlines part of the item to screen
hp <- function(item, numlines=15)
	{
	txt = capture.output(print(item))
	
	cat("\n")
	cat("hp() -- 'head print' of first ", numlines, " lines of item:\n\n")
	for (i in 1:numlines)
		{
		cat(txt[i], "\n", sep="")
		}
	return()
	}


#######################################################
# Common manipulations of R objects
#######################################################
unlist_df_factor <- function(df)
	{
	outdf <- data.frame(lapply(df, function(x) factor(unlist(x))))
	}


# length, like in Python
len <- function(x)
	{
	return(length(x))
	}

# the number of distinct values may be an issue,
# let's count those
count_uniq <- function(x)
	{
	length(unique(x))
	}
count_zeros <- function(x)
	{
	sum(x==0)
	}


# Return the number of matching characters
count_chars <- function(char_to_count, tmpstr)
	{
	tmpchars = strsplit(tmpstr, "")
	matches = tmpchars == char_to_count
	number_that_match = sum(matches)
	return(number_that_match)
	}

# Return current date and time in a format suitable for labeling a file
datetime <- function()
	{
	txt = as.character(Sys.time())
	# "2012-03-18 19:51:32 PDT"
	
	n = nchar(txt)
	txt2 = substr(txt, start=1, stop=n-0)
	txt3 = gsub(" ", "_", txt2)
	txt4 = gsub("\\:", "", txt3)
	
	return(txt4)
	}

# Convert # of seconds since the beginning of 1970
# (numeric, or derived from POSIX_ct times)
# to a POSIXlt date object.
#
#
# Note: R assumes seconds are in your computer's local timezone,
# and if your output timezone is GMT, the time will be altered.
#
# I.e., if your computer is in PST, the UTC/GMT time that results
# will be +2 more hours; so -2 is the correct adjustment, even though
# in real life, PST = UTC-8 (weird I know).
datetime_from_secs <- function(POSIX_ct_dates, local_timezone_shift_from_UTC=-2)
	{
	cat("\n\nNote: datetime_from_secs() is assuming that the input seconds are\nfrom a timezone at UTC", as.character(local_timezone_shift_from_UTC), ".\n", sep="")
	
	# Convert, and correct for the addition R makes to convert PST to UTC
	POSIX_lt_dates = as.POSIXlt(POSIX_ct_dates+(local_timezone_shift_from_UTC*60*60), tz="UTC", origin="1970-01-01 00:00.00 UTC")
	
	return(POSIX_lt_dates)
	}


# 2014-01-30_NJM:
# moving to BioGeoBEARS readwrite
# for process_DIVA_output

# Trim surrounding spaces, AND remove any tabs
trim2 <- function(tmpstr)
	{
	require(gdata) # for trim
	
	# Convert any tabs to space; leading/trailing tabs will be removed by trim
	tmpstr = gsub(pattern="\t", replacement=" ", tmpstr)
	tmpstr = trim(tmpstr)
	return(tmpstr)
	}






# Get fns with suffix, from either directory or fns list
# Do NOT use dot in the suffix
get_fns_matching_txt <- function(tmpdir=NULL, fns=NULL, text_to_match = NULL, returnfullnames=TRUE)
	{
	# If tmpdir is NOT null, get those files from list.files.
	if (is.null(tmpdir) == FALSE)
		{
		fns = list.files(tmpdir, full.names=returnfullnames)
		# Remove //, change to /
		fns = slashslash(fns)
		
		fns_without_paths = list.files(tmpdir, full.names=FALSE)
		}
	
	# Return true/false for matched text
	TF = grepl(pattern=text_to_match, x=fns_without_paths)
	
	matching_fns = fns[TF]
	return(matching_fns)
	}



# Get fns with suffix, from either directory or fns list
# Do NOT use dot in the suffix
get_fns_with_suffix <- function(tmpdir=NULL, fns=NULL, suffix_to_match = "newick")
	{
	# If tmpdir is NOT null, get those files from list.files.
	if (is.null(tmpdir) == FALSE)
		{
		fns = list.files(tmpdir)
		}
	suffixes = get_fn_suffixes(fns=fns)
	
	nums = match_item_in_list2(item=suffix_to_match, list2=suffixes, return_indices=TRUE)
	
	matching_fns = fns[nums]
	return(matching_fns)
	}


# Get the suffixes of a list of filenames
get_fn_suffixes <- function(tmpdir=NULL, fns=NULL)
	{
	if (is.null(tmpdir) == FALSE)
		{
		fns = list.files(tmpdir)
		}
	
	# Split the filenames on "."
	splitfns = strsplit(fns, "\\.")
	
	# For each item, return last
	suffixes = mapply(return_last_nosingles, splitfns)
	
	return(suffixes)
	}


# Return last item
return_last <- function(tmplist)
	{
	return(tmplist[length(tmplist)])
	}

# Return last item, NULL if only one item
return_last_nosingles <- function(tmplist)
	{
	if (length(tmplist) == 1)
		{
		return(NULL)
		} else {
		return(tmplist[length(tmplist)])
		}
	}


##################################
# head/tail lines to file
# Instead of grep, use sed
##################################
# http://www.unix.com/shell-programming-scripting/14477-how-extract-range-lines-file.html
# cat combined_75000_trees.trees |  | cat

# eliminate all spaces, tabs, whitespace,
# don't allow "" in result
strsplit_nowhitespace <- function(tmpstr)
	{
	# split on any whitespace
	words = strsplit(tmpstr, "[[:space:]+]")[[1]]
	
	# remove any blanks that remain
	words = words[words != ""]
	
	return(words)
	}


defaults='
numlines_str = "   75091 combined_75000_trees.trees"
# numlines
# numlines()
'
linecount <- function(fn)
	{
	# SLLLLLOOWW as this does a full word count
	# cmdstr = paste("wc -l ", fn, sep="")
	# numlines_str = system(cmdstr, intern=TRUE)
	# words = strsplit_nowhitespace(numlines_str)
	
	#print(words)
	# numlines = as.numeric(words[1])
	
	# Faster Line Count With Grep Rather Than Wc
	# http://www.dzone.com/snippets/faster-line-count-grep-rather
	cmdstr = paste("grep -c '\n' ", fn, sep="")
	numlines_str = system(cmdstr, intern=TRUE)
	numlines = as.numeric(numlines_str)
	
	
	# Starting at 5:20:
	
	# Using these commands:
	# http://www.theunixschool.com/2010/06/5-different-ways-to-count-total-lines.html
	# cat /Users/nickm/Downloads/AT1_MrB_simple3.time.trees.txt | wc -l
	# sed -n '$=' "/Users/nickm/Downloads/AT1_MrB_simple3.time.trees.txt"
	# grep -c '\n' /Users/nickm/Downloads/AT1_MrB_simple3.time.trees.txt
	# awk 'END {print NR}' /Users/nickm/Downloads/AT1_MrB_simple3.time.trees.txt
	# 
	# LOL...all identical -- 5-6 minutes!
	
	return(numlines)
	}




# Extract the string (strings?) matching the regexp
extract_regexp <- function(tmpstr, tmpre_string)
	{
	matches = gregexpr(tmpre_string, tmpstr)[[1]]
	matches_end = matches-1+attr(matches,"match.length")
	x = mapply(substr, tmpstr, matches, matches_end)
	return(x)
	}


# Get the first digit.  If it is a :, it's not a tip...
get_firstchar <- function(tmpstr)
	{
	tmpchars = strsplit(tmpstr, "")[[1]]
	first_char = tmpchars[1]
	return(first_char)
	}


# Normalize values to a large positive number
normalize_vals <- function(char_dtf, maxvals = NULL, baseval=0, normval=100)
	{
	
	# Normalize values to between 0 and 1
	# set these vals to 1
	if (is.null(maxvals))
		{
		maxvals = apply(char_dtf, 2, max)
		}
	char_dtf_fraction_of_one = char_dtf
	for (rownum in 1:nrow(char_dtf))
		{
		char_dtf_fraction_of_one[rownum, ] = char_dtf[rownum, ] / maxvals
		} 
	
	# Multiply by 100
	char_dtf_fraction_of_100 = char_dtf_fraction_of_one * 100
	
	# Add 1000
	char_dtf_1000 = char_dtf_fraction_of_100 + 1000
	
	# to reverse, subtract 1000, divide by 100, multiply by maxval
	}



# Convert 0-1 values to logit
Pvals_to_logit <- function(x, minP=0.0001)
	{
	x1 = x
	
	maxP = 1-minP
	
	vals_LT_min_TF = x < minP
	vals_GT_max_TF = x > maxP
	
	cat("\nRunning Pvals_to_logit():\n")
	cat("Warning: ", sum(vals_LT_min_TF), " of ", length(x), " values converted to ", minP, "\n", sep="")
	cat("Warning: ", sum(vals_GT_max_TF), " of ", length(x), " values converted to ", 1-minP, "\n", sep="")
	
	x1[vals_LT_min_TF] = minP
	x1[vals_GT_max_TF] = 1-minP
	
	logitx = log(x / (1-x))
	
	
	#########################################################
	# Error check (weird floating-point behavior
	# seems to sometimes produce logit values outside the allowed range...
	#########################################################
	min_logit = log(minP / (1-minP))
	max_logit = log(maxP / (1-maxP))

	logitx[logitx < min_logit] = min_logit
	logitx[logitx > max_logit] = max_logit
	
	return(logitx)
	}
	
# Convert logit values to 0-1
logit_to_Pvals <- function(x)
	{
	x1 = x
	
	#vals_LT_min_TF = x < minP
	#vals_GT_max_TF = x > (1-minP)
	
	#cat("\nRunning Pvals_to_logit():\n")
	#cat("Warning: ", sum(vals_LT_min_TF), " of ", length(x), " values converted to ", minP, "\n", sep="")
	#cat("Warning: ", sum(vals_GT_max_TF), " of ", length(x), " values converted to ", 1-minP, "\n", sep="")
	
	#x1[vals_LT_min_TF] = minP
	#x1[vals_GT_max_TF] = 1-minP
	
	invlogitx = exp(x) / (1 + exp(x))
	
	return(invlogitx)
	}
	

# Long summarize object
smml <- function(x, mxl=NA)
	{
	smm(x, mxl=NA)
	}


# Short summarize object
smm <- function(x, mxl=10)
	{
	cat("\n")
	#cat("smm(): PROVIDING OVERALL SUMMARY OF OBJECT...\n")
	cat("CLASS:	", class(x), "\n", sep="")

	dims = dim(x)
	

	if (length(dims) > 1)
		{
		# if maxcols is NA, print all cols
		if (is.na(mxl))
			{
			mxl = dims[2]
			}
		# if maxcols is GT dims2, use dims2
		if (mxl > dims[2])
			{
			mxl = dims[2]
			}
		
		# Cols to print
		colstp = 1:mxl
		}

	
	
	total_length = 1
	for (d in 1:length(dims))
		{
		total_length = total_length * dims[d]
		}
	
	cat("#DIMS:	", length(dims), "\n", sep="")

	if (!is.null(dims))
		{
		if (length(dims) == 1)
			{
			cat("LENGTH:	",  length(x), "\n", sep="")	
			}
		if (length(dims) == 2)
			{
			cat("DIMS:	nrow=", dims[1], " ncol=", dims[2], "\n", sep="")
			}
		if (length(dims) == 3)
			{
			cat("DIMS:	nrow=", dims[1], " ncol=", dims[2], " nz  =	", dims[3], "\n", sep="")
			}
		if (length(dims) > 3)
			{
			cat("DIMS:	", dims, "\n", sep="")
			}
		}

	cat("TTLEN:	", total_length, "\n", sep="")

	cat("NAMES: ")
	if (length(dims > 1))
		{
		cat(names(x)[1:mxl])
		} else {
		cat(names(x))
		}
	cat("\n")
		

	cat("HEAD/TAIL:\n")

	# if dim is null, just do length etc.
	if (is.null(dims))
		{
		cat("LENGTH:	",  length(x), "\n", sep="")
		cat("1ST4: \n")
		print(x[1:4])
	
		cat("LAST4: \n")
		maxval = length(x)
		print(x[(maxval-3) : maxval])
		}
	
	
	if (length(dims) == 3)
		{
		cat("TOP4: \n")
		print(x[1:4, colstp, 1])
		print(x[1:4, colstp, 2])
	
		cat("BOT4: \n")
		print(x[(dims[1]-3):(dims[1]-0), colstp, 1])
		print(x[(dims[1]-3):(dims[1]-0), colstp, 2])
		}
	if (length(dims) == 2)
		{
		cat("TOP4: \n")
		print(x[1:4, colstp])
	
		cat("BOT4: \n")
		print(x[(dims[1]-3):(dims[1]-0), colstp])
		}
	if (length(dims) == 1)
		{
		cat("1ST4: \n")
		print(x[1:4])
	
		cat("LAST4: \n")
		print(x[(dims[1]-3):(dims[1]-0)])
		}
	if (length(dims) > 3)
		{
		cat("1ST4: \n")
		print(x[1:4])
	
		cat("LAST4: \n")
		print(x[(dims[1]-3):(dims[1]-0)])
		}
	}

# Summarize object
summ <- function(x)
	{
	cat("\n")
	cat("SUMM(): PROVIDING OVERALL SUMMARY OF OBJECT...\n")
	cat("\nCLASS OF OBJECT:\n")
	print(class(x))

	cat("\nDIMENSIONS OF OBJECT:\n")
	print(dim(x))

	cat("\nLENGTH OF OBJECT:\n")
	print(length(x))

	cat("\nATTRIBUTES OF OBJECT:\n")
	print(attributes(x))

	cat("\nSUMMARY() OF OBJECT:\n")
	print(summary(x))
	
	cat("\ncls.df(): print the classes of each column (if it's a data.frame):\n")
	cls.df(x, printout=TRUE)
	
	}

print.default2 <- function (x, digits = NULL, quote = TRUE, na.print = NULL, print.gap = NULL, right = FALSE, max = NULL, useSource = TRUE, ...)
	{
	if (is.numeric(x))
		{
		x <- as.numeric(sprintf("%7.3f", x))
		}
	noOpt <- missing(digits) && missing(quote) && missing(na.print) && missing(print.gap) && missing(right) && missing(max) && missing(useSource) && length(list(...)) == 0L
	
	.Internal(print.default(x, digits, quote, na.print, print.gap, right, max, useSource, noOpt))
	}



setup = "
startline=1
endline=100
startstr_to_find='Iteration'
"

get_startline <- function(fn, startline=1, endline=100, startstr_to_find)
	{
	# Read part of a large file, and find the starting string
	skipnum = startline - 1
	tmp_numlines = endline - startline +1
	
	# Read lines
	# readlines
	tmplines = scan(file=fn, what="character", skip=skipnum, nlines=tmp_numlines, sep="\n")

	
	for (i in 1:length(tmplines))
		{
		tmpline = tmplines[i]
		#print(tmpline)
		if (grepl(startstr_to_find, tmpline))
			{
			outline_num = skipnum + i
			return(outline_num)
			}
		}
	# If startstring not found, return NA
	return(NA)
	}


# Run something like this:
# grep -m 2 --line-number "Iteration[[:space:]]" params_ln_continuous_BayesTraits.txt.log.txt
grep_startline <- function(fn, hitnum=1, startstr_to_find="Iteration[[:space:]]")
	{
	# Use grep to return the line number of the FIRST line containing the matching string
	cmdstr = paste("grep -m ", hitnum, " --line-number ", startstr_to_find, " ", fn, sep="")

	# intern = TRUE reports the result to R
	grep_return = system(cmdstr, intern=TRUE) 
	
	# split on whitespace (spaces and tabs)
	grep_list = strsplit(grep_return, ":")[[1]]
	
	linenum = as.numeric(grep_list[1])	
	return(linenum)
	
	# If startstring not found, return NA
	return(NA)
	}
	
	
setup="
linenum = 100
fn = 'params_ln_continuous_BayesTraits.txt.log.txt'
"
sed_a_line_by_number <- function(fn, linenum)
	{
	# retrieve a particular line by line number (fastest)
	cmdstr = paste("sed '", linenum, "q;d' ", fn, sep="")
	sed_result = system(cmdstr, intern=TRUE) 
	
	return(sed_result)
	}


get_numlines <- function(fn)
	{
	# Get the number of lines in a file
	
	# Check if the file exists
	TF = file.exists(fn)
	if (TF == FALSE)
		{
		txt = paste0("\nWARNING in get_numlines(): file '", fn, "' does not exist in directory\n'", getwd(), "'. Returning linecount=0.")
		cat(txt)
		linecount = 0
		return(linecount)
		}
	
	# intern = TRUE reports the result to R
	cmdstr = paste("wc -l ", fn, sep="")
	wc_return = system(cmdstr, intern=TRUE) 
	
	# split on whitespace (spaces and tabs)
	wc_list = strsplit_whitespace(wc_return)
	
	linecount = as.numeric(wc_list[1])
	return(linecount)
	} # END get_numlines <- function(fn)

ppdf <- function(tmp_pdffn=pdffn, a="")
	{
	# open whatever pdffn you have been working on
	#dev.off()
	
	# Turn off all
	if (a == "all")
		{
		graphics.off()
		}
	
	cmdstr = paste("open ", tmp_pdffn, sep="")
	system(cmdstr)
	}
	
# Merge PDFs using Ghostscript
# http://hints.macworld.com/article.php?story=2003083122212228
merge_pdfs <- function(pdffns, merge_pdffn, wd=getwd())
	{
	setwd(wd)
	pdfs_string = paste(pdffns, sep="", collapse=" ")
	
	gs_string = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile="
	
	merge_pdffn = paste(wd, "/", merge_pdffn, sep="")
	merge_pdffn = slashslash(merge_pdffn)
	
	cmdstr = paste(gs_string, merge_pdffn, " ", pdfs_string, sep="")
	print(cmdstr)
	system(cmdstr)
	}


# prflag <- function(x, printflag=TRUE)
# 	{
# 	# A standard function to print (or not) certain variables,
# 	#   based on a master printflag
# 	# This avoids having to comment in/out various code chunks
# 	#   while debugging.
# 	if (printflag == TRUE)
# 		{
# 		# CAT instead of PRINT if it's a string or numeric
# 		if (is.character(x))
# 			{
# 			cat(x, "\n", sep="")
# 			}
# 		if (is.numeric(x))
# 			{
# 			cat(x, "\n", sep="")
# 			} else {
# 			print(x)
# 			}
# 		}
# 	else
# 		{
# 		pass="BLAH"
# 		}
# 	}


checkwd <- function(tmpwd, lastchar_should_be = "/")
	{
	# Error check: tmpwd must end in slash.
	# (for Windows, \\ )
	
	characters = strsplit(tmpwd, "")[[1]]
	if (characters[length(characters)] != "/")
		{
		tmpwd = paste(tmpwd, "/", sep="")
		} else {
		tmpwd
		}
	return(tmpwd)
	}




#######################################################
# extract_fn_from_path:
#######################################################
#' Get the filename from a path
#'
#' The filename is split on slashes, and the last item is taken; this should be just
#' the filename.
#' 
#' @param fn_with_path The filename, with partial or full path
#' @return \code{fn} The extracted filename
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' fn_with_path = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/modiscloud/extdata/2002_raw_MYD/MYD35_L2.A2002185.1910.005.2007206043609.hdf"
#' extract_fn_from_path(fn_with_path)
#'
extract_fn_from_path <- function(fn_with_path)
	{
	words = strsplit(fn_with_path, split="/")[[1]]
	fn = words[length(words)]
	return(fn)
	}

# get path
# getpath_
# get_path_

# get_suffix
# get suffix_

# Backup file, if it exists
fn_bkup <- function(tmpdir, fn)
	{
	# Check if the old log file exists
	files_list = list.files(path=tmpdir)
	if (fn %in% files_list)
		{
		cmdstr = paste("mv ", fn, " ", paste(fn, ".bkup", sep=""), sep="")
		system(cmdstr)
		}
	}

# printall <- function(dtf, chunksize_toprint = 40, printflag=TRUE)
# 	{
# 	# Print everything in a data frame, in chunks of e.g. 50 rows
# 	if (nrow(dtf) <= chunksize_toprint)
# 		{
# 		prflag(dtf, printflag=printflag)
# 		return(dtf)
# 		}
# 	rows_toprint = seq(1, nrow(dtf), chunksize_toprint)
# 	
# 	if (printflag == TRUE)
# 		{
# 		for (i in 1 : (length(rows_toprint)-1) )
# 			{
# 			tmp_rows_toprint_start = rows_toprint[i]
# 			tmp_rows_toprint_end = rows_toprint[i+1]
# 			prflag(dtf[tmp_rows_toprint_start:tmp_rows_toprint_end, ])
# 			}
# 		
# 		# Then print the end
# 		tmp_rows_toprint_start = rows_toprint[length(rows_toprint)]
# 		tmp_rows_toprint_end = nrow(dtf)
# 		prflag(dtf[tmp_rows_toprint_start:tmp_rows_toprint_end, ])
# 		}	
# 	}


printall2 <- function(dtf, maxprint=99999)
	{
	old_maxprint = getOption("max.print")
	
	options(max.print = maxprint)
	print(dtf)
	options(max.print = old_maxprint)
	}



# moved to biogeobears stratified
# # Get indices of all matches to a list
# findall <- function(what, inlist)
# 	{
# 	TFmatches = inlist == what
# 	indices = 1:length(inlist)
# 	matching_indices = indices[TFmatches]
# 	return(matching_indices)
# 	}

# print to screen the header of a file
headf <- function(fn, nlines=5)
	{
	#lines = scan(file=fn, what="character", sep="\n")
	for (i in 1:nlines)
		{
		line = scan(file=fn, what="character", sep="\n", skip=(i-1), nlines=1)
		print(line)
		}
	}

# print to screen the tail of a file
tailf <- function(fn)
	{
	#lines = scan(file=fn, what="character", sep="\n")
	#numlines = length(lines)
	numlines = linecount(fn)
	for (i in (numlines-5):numlines)
		{
		#print(lines[i])
		line = scan(file=fn, what="character", sep="\n", skip=(i-1), nlines=1)
		print(line)
		}
	}

tailfast <- function(fn)
	{
	cmdstr = paste("tail ", fn, sep="")
	system(cmdstr)
	}


# remove punctuation
remove_punct <- function(tempstr)
	{
	# Do find-replace to convert the names to all underscores, no spaces, no ' or "
	# spaces to _
	tempstr = gsub(" ", "_", tempstr)
	# - to _
	tempstr = gsub("-", "_", tempstr)
	# . to _
	tempstr = gsub("\\.", "_", tempstr)
	# "'" to ""
	tempstr = gsub("'", "", tempstr)
	# '"' to ''
	tempstr = gsub('"', '', tempstr)
	
	return(tempstr)
	}

# remove punctuation from a list
remove_punct_list <- function(templist)
	{
	for (k in 1:length(templist))
		{
		tempstr = templist[[k]]
		templist[[k]] = remove_punct(tempstr)
		}
	return(templist)
	}

#######################################
# remove "'" from a file 
#######################################
# remove apostrophes
remove_apostrophes <- function(tempstr)
	{
	# Do find-replace to convert the names to all underscores, no spaces, no ' or "
	# spaces to _
	tempstr = gsub("'", "", tempstr)
	return(tempstr)
	}

# (e.g., MrBayes, read_nexus_data2
# remove punctuation from a list
remove_apostrophes_from_list <- function(templist)
	{
	for (k in 1:length(templist))
		{
		tempstr = templist[[k]]
		templist[[k]] = remove_apostrophes(tempstr)
		}
	return(templist)
	}

remove_apostrophes_from_file <- function(fn, outfn=paste(fn, "_noApostrophes", sep=""))
	{
	# Get the list of strings from the file
	templist = scan(fn, what="character", sep="\n")
	
	# Remove the apostrophes
	templist2 = remove_apostrophes_from_list(templist)
	
	write_lines_good(templist2, outfn)
	return(outfn)
	}



# split out: Taxon (subtaxon)
split_out_subtaxon <- function(templist)
	{
	templist_col1 = templist
	templist_col2 = templist
	for (k in 1:length(templist))
		{
		tempstr = templist[[k]]
		words = strsplit(tempstr, " \\(")[[1]]
		if (length(words) == 1)
			{
			templist_col1[[k]] = words[1]
			templist_col2[[k]] = ""
			}
		if (length(words) == 2)
			{
			templist_col1[[k]] = words[1]
			templist_col2[[k]] = gsub("\\)", "", words[2])
			}
		
		}
	
	newcols = cbind(templist_col1, templist_col2)
	return(newcols)
	}



# http://r.789695.n4.nabble.com/Remove-empty-list-from-list-td806706.html
delete.NULLs <- function(x.list)
	{
	# delele null/empty entries in a list
    x.list[unlist(lapply(x.list, length) != 0)]
	} 

remove_len_zero <- function(tmparr)
	{
	newarr = tmparr[lapply(tmparr, nchar) != 0]
	return(newarr)
	}




# match using 2 criteria
match_items <- function(tmplist, crit1, crit2)
	{
	outlist = c()
	
	crit1_list = c()
	for (i in 1:length(tmplist))
		{
		if (grepl(crit1, tmplist[i]))
			{
			crit1_list = c(crit1_list, tmplist[i])
			}
		}

	crit2_list = c()	
	for (i in 1:length(crit1_list))
		{
		if (grepl(crit2, crit1_list[i]))
			{
			crit2_list = c(crit2_list, crit1_list[i])
			}
		
		}
	return(crit2_list)
	}



# match using crit1, but skip (leave out) if they have crit2
match_item_nomatch <- function(tmplist, crit1, crit2)
	{
	outlist = c()
	
	crit1_list = c()
	for (i in 1:length(tmplist))
		{
		if (grepl(crit1, tmplist[i]))
			{
			crit1_list = c(crit1_list, tmplist[i])
			}
		}

	crit2_list = c()	
	for (i in 1:length(crit1_list))
		{
		if (grepl(crit2, crit1_list[i]))
			{
			blah = "blah"
			}
		else
			{
			crit2_list = c(crit2_list, crit1_list[i])
			}
		
		}
	return(crit2_list)
	}




# Get corresponding rownums for matching cells, in order
get_corresponding_rownums <- function(names1, names2)
	{
	names_from_1_in_2_TF = names1 %in% names2
	names_from_2_in_1_TF = names2 %in% names1
	
	# names_from_1_in_2 = tr_beastcon1$prt_beast_nodestats$tipnames[names_from_1_in_2_TF]
	# names_from_2_in_1 = tr_beastcon2$prt_beast_nodestats$tipnames[names_from_2_in_1_TF]
	names_from_1_in_2 = names1[names_from_1_in_2_TF]
	names_from_2_in_1 = names2[names_from_2_in_1_TF]
	
	rownums1 = (1:length(names_from_1_in_2_TF))[names_from_1_in_2_TF]
	rownums2 = (1:length(names_from_2_in_1_TF))[names_from_2_in_1_TF]
	
	alpha_order_kept_names_1 = order(names_from_1_in_2)
	alpha_order_kept_names_2 = order(names_from_2_in_1)
	
	corresp_node_rownums1 = rownums1[alpha_order_kept_names_1]
	corresp_node_rownums2 = rownums2[alpha_order_kept_names_2]
	
	corresp_rownums = cbind(corresp_node_rownums1, corresp_node_rownums2)
	return(corresp_rownums)
	}






replace_str_file <- function(fn, oldstr, newstr)
	{
	# Read a text file into a list of strings
	lines = scan(fn, what="character", sep="\n")

	# Replace output date-calibrated parsimony trees file
	newlines = replace_str_lines(lines, oldstr, newstr)
	
	# Write the list of lines to a file
	write.table(newlines, file=fn, quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	cat("In: ", fn, ", '", oldstr, "' replaced with '", newstr, "'\n", sep="")
	
	return(fn)
	}


replace_str_lines <- function(lines, oldstr, newstr)
	{

	# Replace output date-calibrated parsimony trees file
	lines = gsub(oldstr, newstr, lines)
	
	return(lines)
	}


# replace each in a text file that matches, with newstr
replace_each_matching_line <- function(fn, str_to_find, newstr)#, blankstr="___")
	{
	# Read a text file into a list of strings
	lines = scan(fn, what="character", blank.lines.skip=FALSE, sep="\n")

	for (i in 1:length(lines))
		{
	#	if (lines[i] == "")
	#		{
	#		lines[i] = blankstr
	#		}
	#	else
	#		{
		if (grepl(str_to_find, lines[i]))
			{
			#print(str_to_find)
			#print(lines[i])
			lines[i] = newstr
			cat("In: ", fn, ", line #", i, " replaced with '", newstr, "'\n", sep="")
			}
	#		}
		}
	# Write the list of lines to a file
	newfn = paste(fn, "_new", sep="")
	write.table(lines, file=newfn, quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	return(newfn)
	}


# Null lines that match
null_lines <- function(list_of_lines, str_to_match)
	{
	
	list_of_lines2 = c()
	
	for (i in 1:length(list_of_lines))
		{
		if (grepl(str_to_match, list_of_lines[i]))
			{
			# This, in effect, deletes the row containing the string
			blah = 0
			}
		else
			{
			list_of_lines2 = c(list_of_lines2, list_of_lines[i])
			}
		}
	return(list_of_lines2)
	}





# Concatenate a list of files into one big file
cat_files <- function(infiles_list, outfn, tmpskip=0)
	{
	new_lines = c()
	list_num_files = c()
	
	for (i in 1:length(infiles_list))
		{
		# input file
		infn = infiles_list[i]
		
		lines = scan(infn, what="character", sep="\n", skip=tmpskip)
		
		list_num_files = c(list_num_files, length(lines))

		new_lines = c(new_lines, lines)
		}
	
	new_lines = write.table(new_lines, file=outfn, quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	return(list_num_files)
	}




print_table_sigs <- function(x, numsig=4,  printout=FALSE)
	{
	new.table = signif_digits_df(dfnums_to_numeric(x), numsig, printout)
	if (printout == TRUE)
		{
		cat("\n")
		print(new.table)
		cat("\n")
		}
	return(new.table)
	}

print_table_row <- function(x, numsig=4,  printout=FALSE)
	{
	require(gdata)	# for "trim" function
	
	if (length(numsig) == 1)
		{
		new.table = signif(as.numeric(x), numsig)
		}
	else
		{
		tmplist = rep(NA, length(numsig))
		for (i in 1:length(numsig))
			{
			cmdstr = paste('tmplist[i] = sprintf("%1.', numsig[i], 'f", x[i])', sep="")
			eval(parse(text = cmdstr))
			}
		}
	
	outstr = ""
	for (item in tmplist)
		{
		outstr = paste(outstr, trim(item), sep="	")
		}
	
	print(paste(tmplist, sep="	"), quote=FALSE)
	#cat("\n")
	return(tmplist)
	}





df_everything_to_char <- function(dtf, max_NAs=0.5)
	{
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	
	for (i in 1:numcols)
		{
			# Get one column, convert to character:
			cmdstr = paste("dtf$", dtf_names[i], " = as.character(dtf$", dtf_names[i], ")", sep="")
			eval(parse(text = cmdstr))			
		}
	
	return(dtf)
	}



df_everything_to_factor <- function(dtf, max_NAs=0.5)
	{
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	
	for (i in 1:numcols)
		{
			# Get one column, convert to character:
			cmdstr = paste("dtf$", dtf_names[i], " = as.factor(dtf$", dtf_names[i], ")", sep="")
			eval(parse(text = cmdstr))			
		}
	
	return(dtf)
	}


df_factors_to_char <- function(dtf, max_NAs=0.5)
	{
	dtf_classes = cls.df(dtf, printout=TRUE)
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	for (i in 1:numcols)
		{
		# Get one column:
		cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], "')", sep="")
		eval(parse(text = cmdstr))
		
		#cat(i, ": ", dtf_names[i], "	=	", cls_col, "\n", sep="")
		cls_col_list[i] = cls_col
		}
	
	for (i in 1:numcols)
		{
		if (cls_col_list[i] == "factor")
			{
			# Get one column, convert to character:
			cmdstr = paste("newcol = as.character(dtf$", dtf_names[i], ")", sep="")
			eval(parse(text = cmdstr))			
			
			cmdstr = paste("dtf$", dtf_names[i], " = newcol", sep="")
			eval(parse(text = cmdstr))				
			}
		}
	tmp_classes = cls.df(dtf)
	dtf_classes$newclasses = tmp_classes[,ncol(tmp_classes)]
	cat("\n")
	cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, ") reports: dataframe 'dtf_classes' has ", nrow(dtf_classes), " rows, ", ncol(dtf_classes), " columns.\n", sep="")
	cat("...names() and classes() of each column below...\n", sep="")
	cat("\n")
	print(dtf_classes)
	
	return(dtf)
	}


# Remove factors from the character-like columns of a data.frame
# (leaves the numbers as numbers)
df_nonum_factors_to_char <- function(dtf, max_NAs=0.5)
	{
	dtf_classes = cls.df(dtf, printout=TRUE)
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	for (i in 1:numcols)
		{
		# Get one column:
		cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], "')", sep="")
		eval(parse(text = cmdstr))
		
		#cat(i, ": ", dtf_names[i], "	=	", cls_col, "\n", sep="")
		cls_col_list[i] = cls_col
		}
	
	for (i in 1:numcols)
		{
		if (cls_col_list[i] == "factor")
			{
			# Get one column, convert to character:
			cmdstr = paste("newcol = as.character(dtf$", dtf_names[i], ")", sep="")
			eval(parse(text = cmdstr))			
			
			cmdstr = paste("dtf$", dtf_names[i], " = newcol", sep="")
			eval(parse(text = cmdstr))				
			}
		}
	tmp_classes = cls.df(dtf)
	dtf_classes$newclasses = tmp_classes[,ncol(tmp_classes)]
	cat("\n")
	cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, ") reports: dataframe 'dtf_classes' has ", nrow(dtf_classes), " rows, ", ncol(dtf_classes), " columns.\n", sep="")
	cat("...names() and classes() of each column below...\n", sep="")
	cat("\n")
	print(dtf_classes)
	
	return(dtf)
	}



# dput on an S4 object, or an S4 object within an S3 object, or something
# Source:
# http://stackoverflow.com/questions/3466599/dputting-an-s4-object
dput2 <- function (x,
                   file = "",
                   control = c("keepNA", "keepInteger", "showAttributes")){
    if (is.character(file))
        if (nzchar(file)) {
            file <- file(file, "wt")
            on.exit(close(file))
        }
        else file <- stdout()
    opts <- .deparseOpts(control)
    if (isS4(x)) {
        cat("new(\"", class(x), "\"\n", file = file, sep = "")
        for (n in slotNames(x)) {
            cat("    ,", n, "= ", file = file)
            dput2(slot(x, n), file = file, control = control)
        }
        cat(")\n", file = file)
        invisible()
    } else if(length(grep('@',capture.output(str(x)))) > 0){
      if(is.list(x)){
        cat("list(\n", file = file, sep = "")
        for (i in 1:length(x)) {
          if(!is.null(names(x))){
            n <- names(x)[i]
            if(n != ''){
              cat("    ,", n, "= ", file = file)
            }
          }
          dput2(x[[i]], file = file, control = control)
        }
        cat(")\n", file = file)
        invisible()
      } else {
        stop('S4 objects are only handled if they are contained within an S4 object or a list object')
      }
    }
    else .Internal(dput(x, file, opts))
}


# dput on an S4 object, or an S4 object within an S3 object, or something
# Source:
# http://stackoverflow.com/questions/3466599/dputting-an-s4-object
# can't get this modified version to insert commas appropriately, so we just
# have reformatted function here
dput2a <- function (x, file = "", control = c("keepNA", "keepInteger", "showAttributes"))
	{
    if (is.character(file))
    	{
        if (nzchar(file))
        	{
            file <- file(file, "wt")
            on.exit(close(file))
	        }
	    else
	    	{
	    	file <- stdout()
	    	}
	    }

    opts <- .deparseOpts(control)
    if (isS4(x)) # if #1
    	{
        cat("new(\"", class(x), "\"\n", file = file, sep = "")
        for (n in slotNames(x))
        	{
            cat("    ,", n, "= ", file = file)
            dput2a(slot(x, n), file = file, control = control)
        	}
        cat(")\n", file = file)
        invisible()
        } else if(length(grep('@',capture.output(str(x)))) > 0) # if #2
        {
        if(is.list(x))
        	{
			cat("list(\n", file = file, sep = "")
			for (i in 1:length(x))
				{
		        if(!is.null(names(x)))
		        	{
					n <- names(x)[i]
					if(n != '')
						{
						cat("    ,", n, "= ", file = file)
						}
			        }
          		dput2a(x[[i]], file = file, control = control)
          		}
			cat(")\n", file = file)
			invisible()
			} else {
			stop('S4 objects are only handled if they are contained within an S4 object or a list object')
			}
	    } else { #if #3
    	.Internal(dput(x, file, opts))
    	}
	}

# multiline sed; the example fixed dput2 output on S4 objects
multiline_sed <- function(fn, patternstr, newstr, outfn="sed_result.txt")
	{
	# R requires \\n here (for system() command)
	# UNIX/Terminal will just want \n
	
	#patternstr = ')\\nnew(\"Polygons\"'
	#newstr = '),new("Polygons"'
	
	# Do find/replace on mulitline (removes \n but who cares)
	
	# This sed works, although it removes the \n:
	#sed -n '1h;1!H;${;g;s|)\nnew(\"Polygons\"|),new("Polygons"|g;p;}' tmppoly2.txt > tmppoly2_edited.txt;

	k = 0
	cmds = NULL
	
	# Starter; from here:
	# http://austinmatzko.com/2008/04/26/sed-multi-line-search-and-replace/
	cmds[[(k=k+1)]] = "sed -n '1h;1!H;${;g;"
	
	# preface to starting string
	cmds[[(k=k+1)]] = "s|"
	
	# pattern string; user specifies escapes etc.
	#cmds[[(k=k+1)]] = ')\\nnew(\"Polygons\"'
	cmds[[(k=k+1)]] = patternstr
	
	# transition to replace string
	cmds[[(k=k+1)]] = "|"
	
	# replacement string; user specifies escapes etc.
	#cmds[[(k=k+1)]] = '),new("Polygons"'
	cmds[[(k=k+1)]] = newstr
	
	# End sed
	cmds[[(k=k+1)]] = "|g;p;}'"
	
	# input filename and output filename
	#outfn = gsub(pattern="\\.", replacement="_edited.", x=fn)
	#cmds[[(k=k+1)]] = paste(" ", fn, " > sed_result.txt;", sep="")
	cmds[[(k=k+1)]] = paste(" ", fn, " > ", outfn, ";", sep="")
	
	cmdstr = paste(cmds, collapse="", sep="")
	return(cmdstr)
	}



moving_average <- function(xvals, yvals, byvals)
	{
	# moving averages
	intervals = seq(from=min(xvals), to=max(xvals), by=byvals)
	xmeans = c()
	ymeans = c()
	for (i in 1:(length(intervals)-1))
		{
		starti = intervals[i]
		endi = intervals[i+1]
		
		indices1 = xvals >= starti
		indices2 = xvals < endi
		indices_in_range = (indices1 + indices2) > 1
		
		vals = yvals[indices_in_range]
		if (length(vals) > 0)
			{
			xmean = mean(c(starti, endi))
			ymean = mean(vals)
			xmeans = c(xmeans, xmean)
			ymeans = c(ymeans, ymean)
			}
		}
	tmp = cbind(xmeans, ymeans)
	moving_avgs = data.frame(tmp)
	return(moving_avgs)
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


# This will only extract 1 group between 1st startflag & 1st doneflag
extract_lines_startstr_to_endstr <- function(lines, string_to_start_at, string_to_end_at, printflag=TRUE, include_endstring=FALSE, instance_to_find=1)
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



# Extract the strings that have a certain word at a certain position
extract_lines_with_word_at_wordnum <- function(lines, word_to_find, wordnum=1, delimiter=" ", printflag=FALSE)
	{
	linenums_to_keep = c()
	
	for (i in 1:length(lines))
		{
		tmpline = lines[i]
		if (tmpline == "")
			{
			next
			}
		words = strsplit(tmpline, split=delimiter)[[1]]
		#print(words)
		if (words[wordnum] == word_to_find)
			{
			#print(words[wordnum])
			linenums_to_keep[length(linenums_to_keep)+1] = i
			if(printflag)
				{			
				print(paste("Storing line #", i, sep=""))
				}

			}
		}
	#print(linenums_to_keep)
	newlines = lines[linenums_to_keep]
	return(newlines)
	}

# This is super-slow
array_to_text_sucky <- function(inarray, spacer)
	{
	tmpstr = inarray[1]
	for (i in 2:length(inarray))
		{
		tmpstr = paste(tmpstr, inarray[i], sep=spacer)
		}
	return(tmpstr)
	}

array_to_text <- function(inarray, spacer)
	{
	tmpstr = paste(as.character(inarray), sep="", collapse=spacer)
	return(tmpstr)
	}


minmax_pretty <- function(x)
	{
	minmax = c( min(pretty(x)), max(pretty(x)) )
	return(minmax)
	}

minmax_literal <- function(x)
	{
	minmax = c( min(x, na.rm=TRUE), max(x, na.rm=TRUE) )
	return(minmax)
	}




convert_nums_to_circular_lat <- function(nums, maxval=90)
	{
	mini = floor(min(nums))
	maxi = ceiling(max(nums))
	nums5 = rep(0, length(nums))

	# go through 0-maxval, maxval-2*maxval, 2*maxval-3*maxval, 3*maxval-4*maxval
	# translates:0-maxval, maxval-0,   0- -maxval,  -maxval-0
	edgeval = maxval * 4

	# get the nums in terms of edgevals
	nums1 = nums/edgeval
	
		
	# Don't change numbers between -maxval and maxval, but change the rest of them...

	##############################################	
	# degrees > 0 (positive)
	##############################################
	posvals = nums > 0
	
	# remove all cycles except the last
	nums2 = nums1[posvals] - floor(nums1[posvals])
	
	# convert back to lat:
	nums3 = nums2 * edgeval

	nums4 = nums3

	vals_to_change = (nums3 > maxval) + (nums3 <= 2*maxval) == 2
	nums4[vals_to_change] = 2*maxval - nums3[vals_to_change]

	vals_to_change = (nums3 > 2*maxval) + (nums3 <= 3*maxval) == 2
	nums4[vals_to_change] = -1*(nums3[vals_to_change] - 2*maxval)

	vals_to_change = (nums3 > 3*maxval) + (nums3 <= 4*maxval) == 2
	nums4[vals_to_change] = 4*maxval - nums3[vals_to_change]

	nums5[posvals] = nums4


	##############################################
	# degrees < 0 (negative)
	##############################################
	negvals = nums < 0
	
	# remove all cycles except the last
	nums2 = nums1[negvals] - ceiling(nums1[negvals])
	
	# convert back to lat:
	nums3 = nums2 * edgeval

	nums4 = nums3


	vals_to_change = (nums3 < -maxval) + (nums3 >= -2*maxval) == 2
	nums4[vals_to_change] = -2*maxval - nums3[vals_to_change]

	vals_to_change = (nums3 < -2*maxval) + (nums3 >= -3*maxval) == 2
	nums4[vals_to_change] = -1*(nums3[vals_to_change] - -1*2*maxval)

	vals_to_change = (nums3 < -3*maxval) + (nums3 >= -4*maxval) == 2
	nums4[vals_to_change] = -4*maxval - nums3[vals_to_change]

	nums5[negvals] = nums4

		
	return(nums5)
	}



convert_nums_to_long <- function(nums, maxval=180)
	{
	# get the nums in terms of 180s
	nums0 = nums
	nums1 = nums
	
	# For nums > 0
	TFnums_gt0 = nums > 0
	nums1[TFnums_gt0] = nums0[TFnums_gt0]/maxval
	
	# remove all cycles except the last
	nums1[TFnums_gt0] = nums1[TFnums_gt0] - floor(nums1[TFnums_gt0])
	
	# convert back to lat:
	nums1[TFnums_gt0] = nums1[TFnums_gt0] * maxval
	

	# For nums < 0
	TFnums_lt0 = nums < 0
	nums1[TFnums_lt0] = nums0[TFnums_lt0]/maxval
	
	# remove all cycles except the last
	nums1[TFnums_lt0] = nums1[TFnums_lt0] - ceiling(nums1[TFnums_lt0])
	
	# convert back to lat:
	nums1[TFnums_lt0] = nums1[TFnums_lt0] * maxval
	
	
	
	return(nums1)
	}








# Basic correlation plot
basic_xy_plot <- function(x, y, titlestart_txt="title\n", xlab=names(x), ylab=names(y), xlim=minmax_pretty(x), ylim=minmax_pretty(y))
	{
	cortest_xy = cor.test(x, y, xlim=c(-1,1), ylim=c(-1,1))
	lm_xy = lm(y ~ x)
	slope = lm_xy$coefficients[2]
	intercept = lm_xy$coefficients[1]
	slopetxt = round(lm_xy$coefficients[2], 4)
	intercepttxt = round(lm_xy$coefficients[1], 4)
	corval = round(cortest_xy$estimate, 2)
	#pval = round(cortest_xy$p.value, 4)
	pval = cortest_xy$p.value
	
	
	titletxt = paste(titlestart_txt, "cor=", corval, "; slope=", slopetxt, "; int=", intercepttxt, "; p=", pval, sep="")
	plot(x, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
	
	# make segment
	minx = min(x, na.rm=TRUE)
	miny = min(y, na.rm=TRUE)
	maxx = max(x, na.rm=TRUE)
	maxy = max(y, na.rm=TRUE)
	
	x0 = minx
	x1 = maxx
	y0 = slope*x0 + intercept
	y1 = slope*x1 + intercept
	
	segments(x0, y0, x1, y1)
	
	title(titletxt)

	}

plot_basic_xy <- function(x, y, titlestart_txt="title\n", xlab=names(x), ylab=names(y), xlim=minmax_pretty(x), ylim=minmax_pretty(y))
	{
	cortest_xy = cor.test(x, y, xlim=c(-1,1), ylim=c(-1,1))
	lm_xy = lm(y ~ x)
	slope = lm_xy$coefficients[2]
	intercept = lm_xy$coefficients[1]
	slopetxt = round(lm_xy$coefficients[2], 4)
	intercepttxt = round(lm_xy$coefficients[1], 4)
	corval = round(cortest_xy$estimate, 2)
	pval = cortest_xy$p.value
	
	titletxt = paste(titlestart_txt, "cor=", corval, "; slope=", slopetxt, "; int=", intercepttxt, "; p=", pval, sep="")
	plot(x, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
	
	# make segment
	minx = min(x, na.rm=TRUE)
	miny = min(y, na.rm=TRUE)
	maxx = max(x, na.rm=TRUE)
	maxy = max(y, na.rm=TRUE)
	
	x0 = minx
	x1 = maxx
	y0 = slope*x0 + intercept
	y1 = slope*x1 + intercept
	
	segments(x0, y0, x1, y1)
	
	title(titletxt)

	}



add_lm_line_to_plot <- function(x, y, ...)
	{
	lm_xy = lm(y ~ x)
	
	# DON'T ROUND FOR PLOTTING!!
	#slope = round(lm_xy$coefficients[2], 2)
	#intercept = round(lm_xy$coefficients[1], 2)
	slope = lm_xy$coefficients[2]
	intercept = lm_xy$coefficients[1]

	# make segment
	minx = min(x, na.rm=TRUE)
	miny = min(y, na.rm=TRUE)
	maxx = max(x, na.rm=TRUE)
	maxy = max(y, na.rm=TRUE)
	
	x0 = minx
	x1 = maxx
	y0 = slope*x0 + intercept
	y1 = slope*x1 + intercept
	
	segments(x0, y0, x1, y1, ...)
	
	return(lm_xy)
	}
	




# Fancy correlation plots
#http://www.personality-project.org/R/r.useful.html
#some fancy graphics   -- adapted from help.cor
panel.cor.scale <- function(x, y, digits=2, prefix="", cex.cor)
	{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r = (cor(x, y,use="pairwise"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex * abs(r))
	}
	
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
	{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r = (cor(x, y,use="pairwise"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex )
	}
	
panel.hist <- function(x, ...)
	{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
	}

panel.points <- function(x, tmppch=20, ...)
	{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    par(pch=tmppch)
    pts <- points(x, pch=".")
	}


pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE, tmppch=NULL) 
	{
	oldpar = par("pch")
	par(pch = tmppch)
	if (smooth )
	 	{
		if (scale)
			{
			pairs(x, diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.smooth)
		    #pairs(x, diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=function(x) panel.points(x) )
		    #pairs(x, diag.panel=panel.hist, upper.panel=panel.cor.scale, lower.panel=function(x) points(x, pch=tmppch))
		    #pairs(x, diag.panel=panel.hist, upper.panel=panel.cor.scale, lower.panel=panel.points)#, lower.panel=panel.points(x, pch=tmppch))
			}
		else
			{
		    #pairs(x, diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=function(x) panel.points(x) )
			pairs(x, diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
			#pairs(x, diag.panel=panel.hist, upper.panel=panel.cor, lower.panel=function(x) points(x, pch=tmppch))
		    #pairs(x, diag.panel=panel.hist, upper.panel=panel.cor, lower.panel=panel.points)#, lower.panel=panel.points(x, pch=tmppch))

			} #else  {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
		}
	else      #smooth is not true
		{
		if (scale)
			{
			pairs(x, diag.panel=panel.hist,upper.panel=panel.cor.scale)
			}
		else
			{
			pairs(x, diag.panel=panel.hist,upper.panel=panel.cor)
			}
		} #end of else (smooth)

	par(pch = oldpar)
	}   #end of function

 
pairs.panels.lm <- function (x,y,tmplm=TRUE,scale=FALSE, tmppch=NULL) 
	{
	oldpar = par("pch")
	par(pch = tmppch)
	if (tmplm)
		{
		if (scale)
			{
			pairs(x, diag.panel=panel.hist, upper.panel=panel.cor.scale, lower.panel=panel.smooth, span=0.5, f=0.5)
		   }
		else
			{
			pairs(x, diag.panel=panel.hist, upper.panel=panel.cor, lower.panel=panel.smooth)
	   	} #else  {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
		}
	else      #smooth is not true
		{
		if (scale)
			{
			pairs(x, diag.panel=panel.hist,upper.panel=panel.cor.scale)
			}
		else
			{
			pairs(x, diag.panel=panel.hist,upper.panel=panel.cor) }
			} #end of else (smooth)
	par(pch = oldpar)
	}   #end of function






#function to replace upper triangle of matrix with unattenuated correlations, and the diagonal with reliabilities
#input is a correlation matrix and a vector of reliabilities
	
 correct.cor <- function(x,y) { n=dim(x)[1]   
	   { x[1,1] <- y[1,1]
	   for (i in 2:n) {
		  x[i,i] <- y[1,i]   #put reliabilities on the diagonal
		  k=i-1
		  for (j in 1:k) {
		     x[j,i] <- x[j,i]/sqrt(y[1,i]*y[1,j])  }   #fix the upper triangular part of the matrix
	   
		    }
		  return(x)  }}
		 
	
 #difference of dependent (paired) correlations (following Steiger,J., 1980, Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87, 245-251)
 paired.r <- function(xy,xz,yz,n) {
	  diff <- xy-xz
	  determin=1-xy*xy - xz*xz - yz*yz + 2*xy*xz*yz
	  av=(xy+xz)/2
	  cube= (1-yz)*(1-yz)*(1-yz)
	  t2 = diff * sqrt((n-1)*(1+yz)/(((2*(n-1)/(n-3))*determin+av*av*cube)))
	  return(t2)
	   }
	
 fisherz <- function(rho)  {0.5*log((1+rho)/(1-rho)) }   #converts r to z     



###################################################
# LaTeX-related functions
###################################################
#
# to run tex2div
library(tools)

# You also need a well-set-up .Rprofile file to make it work (see old
# or current .Rprofile, just one line)

# Might not need much more than:
############################################################################
# See /Dropbox/_njm/_njm/LaTeX_functioning/ for a detailed example of LaTeX / Sweave
# Rnw Snw etc. working...
#############################################################################




## Set working directory
##wd = "/Users/nick/Desktop/_IB286_Marshall_paleobio/"
#wd = "/Users/nick/Desktop/_2010-10-12_work/_Marshall_PBDB/_IB286_Marshall_paleobio/"
#wd = "/Dropbox/_njm/_njm/LaTex_functioning/"
#setwd(wd)

runlatex <- function(fn, wd=getwd())
	{
	# Setup directories for LaTex files and the working directories
	# Old, 2011 version -- still seems to work (2013-02-13, NJM)
	swwds = c("/usr/local/texlive/2011basic/texmf-dist/tex/latex/base/", wd, "/usr/texbin/")
	
	
	# New, 2012 version
	#swwds = c("/usr/local/texlive/2012basic/texmf-dist/tex/latex/base/", wd, "/usr/texbin/")

	# New version, installed from:
	# /Users/nickm/Downloads/install-tl-unx.tar.gz 
	# cd /Users/nickm/Downloads/install-tl-20130128 
	# sudo ./install-tl
	# # This installed to: 
	#  Add /usr/local/texlive/2012/texmf/doc/man to MANPATH, if not dynamically determined.
	# Add /usr/local/texlive/2012/texmf/doc/info to INFOPATH.
	#  Most importantly, add /usr/local/texlive/2012/bin/universal-darwin
	#  to your PATH for current and future sessions.
	#
	# cd ~
	# bbedit .bash_profile
	# bbedit .Rprofile
	
	# These commands build and open the Sweave tex and PDF files...
	# permfn = "example-2.Snw"
	# permfn = "hw2_v2.Snw"
	permfn = fn
	Sweave(permfn)
	
	texfn = sub(".Snw", ".tex", permfn)
	texi2dvi(texfn, pdf="TRUE", quiet=FALSE, texinputs=swwds)
	
	pdffn = sub(".Snw", ".pdf", permfn)
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	
	return(cmdstr)
	}





##############################################################
# TreeRogue 0.1: TreeThief-like tree grabbing using x,y 
# coordinates digitized from an image of a phylogenetic tree.
##############################################################
# GOAL: to process x, y coordinates into a Newick-format tree
##############################################################
# by Nick Matzke
# Copyright 2010
# matzke@berkeley.edu
# 10/27/2010
##############################################################
#
# Free to use/redistribute under GNU General Public License, version 2 
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
# http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
# 
###################################################################
#
# #===================================================
# Run with these commands
# ===================================================
# library(ape)
# 
# put the text files (at bottom) into your working directory
# wd = "/Users/nick/Desktop/__projects/2010-11-01_cone_snails/"
# setwd(wd)
# 
# xy2 = treerogue_read_files()
# xy3 = treerogue_associate_branch_bottoms_with_nodes(xy2)
# tr = build_tree_using_corners(xy3)
# plot(tr)
#
####################################################################
#
# NOTES
# *Heavily* modified from a very limited script posted here by 
# bbolker@gmail.com:
# https://stat.ethz.ch/pipermail/r-sig-phylo/2008-November/000173.html
#
# Background: I worked up this script after getting frustrated at
# (1) the failure of the famous "TreeThief" to work on any modern machine;
# (2) my failure to get the newer "TreeSnatcher" to work on Mac OS X 10.4.11
# (some weirdness about Java and X11, as far as I can tell), and
# (3) no other good options.
#
# Summary: This script takes in x,y coordinates (with the lower left as the origin)
# of nodes, tips, and branch bottoms ("corners"), and builds a tree out of it.
# 
# It assumes:
# (1) Your tree is horizontal left-to-right, with the tips on the right
# (2) Your tree is a "square" tree (i.e. no diagonal/curved branches)
#
# I captured my x,y coordinates using GraphClick 3.0, available for $8 at:
# http://www.arizona-software.ch/graphclick/
#
# (there is a free trial download, but only lets you export 10 coordinates at a
# time, so it is pointless)
#
# REQUIRED INPUTS:
#   (for each, digitize the relevant points in GraphClick, and
#    File-->Export to a text file):
# 
# (Note: all text files should have a header line)
#
# 1. A tab-delimited text file with x and y for each internal node
#
# 2. A tab-delimited text file with x and y for each tip
#
# 2a. A text file with tip names, in order from top-to-bottom
# 
# 3. A tab-delimited text file with x and y for each tip for each "corner"
#    (i.e., the bottom of each branch).
#
# 4. For now, do NOT digitize the bottom of the root of the tree, if the
#    image you are digitizing has one.  You could add the root length later
#    manually, if you like.
#
# 5. The tree must be fully dichotomous (if the image you are digitizing is not,
#    you can "fake it" by resolving polytomies by clicking digitization points
#    to, in effect, manually resolve these polytomies with very short branches.
#    Note that you will have to add a corner for each internal node you add (!).
#
#    The R APE package can collapse short branches to polytomies later, if you like.
#
# Trees do not have to be ultrametric, and digitization does not have to be 
# exact -- the script will attempt to match the most likely branch bottoms
# to the nodes (a graphic is produced by R so that you can check the results
# and tinker if necessary).
#
# Requires the APE library.
# 
#############################################################

library(ape)


# Assumes default filenames
treerogue_read_files <- function(internal_nodes_file="internal.txt", tip_nodes_file="tips.txt", corner_nodes_file="corners.txt", tipnames_file="tipnames.txt", txt_file_type="txt", points_files_have_column_headings=TRUE, tipnames_file_has_column_heading=TRUE)
	{
	defaults='
	# Mac / Nick / Raiff
	internal_nodes_file="internal.txt"
	tip_nodes_file="tips.txt"
	corner_nodes_file="corners.txt"
	tipnames_file="tipnames.txt"
	txt_file_type="txt"
	points_files_have_column_headings=TRUE
	tipnames_file_has_column_heading=TRUE
	
	# Windows/Raji
	internal_nodes_file="Nodes.csv"
	tip_nodes_file="Tips.csv"
	corner_nodes_file="Corners.csv"
	tipnames_file="tipnames.txt"
	txt_file_type="csv"
	points_files_have_column_headings=FALSE
	tipnames_file_has_column_heading=TRUE
	'
	
	
	# The text files types can be tab-delimited text (.txt), or comma-delimited (csv)
	if (txt_file_type == "txt")
		{
		internal = read.table(file=internal_nodes_file, header=points_files_have_column_headings, skip=0, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		tips = read.table(file=tip_nodes_file, header=points_files_have_column_headings, skip=0, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		corners = read.table(file=corner_nodes_file, header=points_files_have_column_headings, skip=0, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		tipnames = read.table(file=tipnames_file, header=tipnames_file_has_column_heading, skip=0, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		}

	# The text files types can be tab-delimited text (.txt), or comma-delimited (csv)
	if (txt_file_type == "csv")
		{
		internal = read.table(file=internal_nodes_file, header=points_files_have_column_headings, skip=0, sep=",", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		tips = read.table(file=tip_nodes_file, header=points_files_have_column_headings, skip=0, sep=",", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		corners = read.table(file=corner_nodes_file, header=points_files_have_column_headings, skip=0, sep=",", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		tipnames = read.table(file=tipnames_file, header=tipnames_file_has_column_heading, skip=0, sep=",", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		}

	if (points_files_have_column_headings == FALSE)
		{
		names(internal) = c("x","y")
		names(tips) = c("x","y")
		names(corners) = c("x","y")
		names(tipnames) = "tipnames"
		}
	if (tipnames_file_has_column_heading == FALSE)
		{
		names(tipnames) = "tipnames"
		}
	
	# combine all the points for plotting purposes
	allpoints = rbind(corners, tips, internal)

	# Setup the plot
	plot(allpoints, pch=".", col="white")
	
	# Plot the points
	points(tips, xlim=c(0,7), ylim=c(0,9), pch="t", col="red")
	title("Plotting your tips (t), internal (i), corners (c). Look especially for\noutlying and missing points.")
	points(internal, pch="i", col="blue")
	points(corners, pch="c", col="green")
	points(corners, xlim=c(0,7), ylim=c(0,9), pch=".", col="black")
	
	# sort the tips from top to bottom in y
	tips = tips[order(tips$y, decreasing = TRUE), ]
	
	# sort the internals from left to right in x
	internal = internal[order(internal$x, decreasing=FALSE), ]
	
	if (nrow(tips) != nrow(tipnames))
		{
		print("ERROR: the number of tips must equal the length of the tipnames!")
		print(paste("Instead, nrow(tipnames) =", nrow(tipnames), "and nrow(tips) =", nrow(tips), sep=" "))
		}
	
	if ((nrow(tips)-1) != nrow(internal))
		{
		print("ERROR: the number of tips-1 must equal the number of the internal nodes!")
		print(paste("Instead, nrow(tips) =", nrow(tips), "and nrow(internal) =", nrow(internal), sep=" "))
		}
	
	nodetypes = c(rep("tip", nrow(tipnames)), rep("internal", nrow(internal)))
	nodenames = unlist(c(tipnames, rep("", nrow(internal))))
	xy = rbind(tips, internal)
	xy2 = cbind(xy, nodetypes, nodenames)
	xy2 = as.data.frame(xy2)
	names(xy2) = c("x", "y", "nodetypes", "nodenames")
	
	xy2 = df_nonum_factors_to_char(xy2, max_NAs=0.5)
	
	xy2$nodetypes[xy2$x == min(xy2$x)] = "root"
	
	if (nrow(corners) != (nrow(xy2)-1))
		{
		print("ERROR: the number of corners must equal the number of nodes-1  !")
		print(paste("Instead, length(nodes) =", nrow(xy2), "and nrow(corners) =", nrow(corners), sep=" "))
		}

	# sort file so that the tips are first
	# tip nodes in order from top to bottom:
	xytips = xy[xy$tipname != "", ]
	
	tips_in_order = xy2$tipname[xy2$tipname != ""]

	return(xy2)

	}


treerogue_associate_branch_bottoms_with_nodes <- function(xy2, corner_nodes_file="corners.txt", points_files_have_column_headings=TRUE, txt_file_type="txt")
	{
	defaults='
	# Mac / Nick / Raiff
	corner_nodes_file="corners.txt"
	points_files_have_column_headings=TRUE
	
	
	# Windows/Raji
	corner_nodes_file="Corners.csv"
	points_files_have_column_headings=FALSE
	'

	# Load the coordinates of the corners

	# The text files types can be tab-delimited text (.txt), or comma-delimited (csv)
	if (txt_file_type == "txt")
		{
		corners = read.table(file=corner_nodes_file, header=points_files_have_column_headings, skip=0, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		}
	if (txt_file_type == "csv")
		{
		corners = read.table(file=corner_nodes_file, header=points_files_have_column_headings, skip=0, sep=",", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
		}

	if (points_files_have_column_headings == FALSE)
		{
		names(corners) = c("x","y")
		}

	# Process the corners and attach them to nodes
	bots = corners
	names(bots) = c("x", "y")

	# Get xy data
	xx = xy2$x
	yy = xy2$y
	
	dtf = associate_branch_bottoms_with_nodes(xx, yy, bots)

	xy3 = cbind(xy2, dtf$chosen_node_bottoms_xx, dtf$chosen_node_bottoms_yy)
	names(xy3) = c(names(xy2), "cx", "cy")
	
	write.table(xy3, "linked_corners.txt", quote=FALSE, sep="	", row.names = FALSE,
		   col.names = TRUE)

	return(xy3)
	}



check_for_duplicate_points <- function(bots, pointtype="default")
	{
	# Get the list of distances
	distbots_matrix = as.matrix(dist(bots))

	# Remove zero distances on diagonal (self-to-self)
	diag(distbots_matrix) = 100
	
	bots_distances = c(distbots_matrix)
	
	
	n1 = 1:nrow(bots)
	n2 = 1:nrow(bots)
	
	closest10 = sort(bots_distances)[1:10]
	
	tmpstr = paste("\nChecking for identical or near-identical ", pointtype, ". Most likely:\n", sep="")

	tmpstr = paste("\nIf you have any distances super-near zero, check for duplicate points.\n", sep="")

	cat(tmpstr)



	closest_bots = NULL
	# For each of the closest 10 distances,
	for (i in 1:length(closest10))
		{
		# Find bots distances that match EXACTLY
		TF = bots_distances == closest10[i]
		TFmat = matrix(TF, nrow=nrow(distbots_matrix), ncol=ncol(distbots_matrix), byrow=FALSE)
		
		sumTF_rows = apply(TFmat, 1, sum)
		sumTF_cols = apply(TFmat, 2, sum)
		
		# This causes a problem, if there are identical bots_distances, and thus 2 instead of 1 in the sum
		#for (j in 1:sum(sumTF_rows))
		for (j in 1:sum(sumTF_rows>0))
			{
			cat(i, j, "\n", sep=", ")
			print(n1)
			print(sumTF_rows)
			print(n1[sumTF_rows > 0])
			
			nums1 = n1[sumTF_rows > 0][j]
			#nums2 = n2[sumTF_cols > 0][j]
			
			bots1 = bots[sumTF_rows > 0, ][j,]
			bots1 = unlist(bots1)
			
			tmprow = c(round(closest10[i], 6), nums1, bots1)
			closest_bots = rbind(closest_bots, tmprow)
			}

		}
		closest_bots
	closest_bots = adf2(closest_bots)
	
	tmpnames = c("dist", "rownum", "x", "y")
	names(closest_bots) = tmpnames
	
	return(closest_bots)
	}
	


# Associate branch bottom coordinates with nodes, and plot the results;
# The user may then edit the output associations if they so desire.
associate_branch_bottoms_with_nodes <- function(xx, yy, bots)
	{
	# There should be one less branch bottom than there are internal nodes
	# (because the root node (should) have no digitized "corner" beneath it)
	nodes = 1:length(xx)
	if (length(nodes) != nrow(bots) +1)
		{
		print("ERROR: the number of corners must equal the number of nodes-1  !")
		print(paste("Instead, length(nodes) =", length(nodes), "and nrow(bots) =", nrow(bots), sep=" "))
		} else {
		
		#######################################
		# Check for duplicate points
		#######################################
		closest_bots = check_for_duplicate_points(bots=bots, pointtype="corners")
		print(closest_bots)
		
		
		closest_tops = check_for_duplicate_points(bots=cbind(xx,yy), pointtype="tops")
		print(closest_tops)
		
		# OK, find bottom of branch to go with each top of branch
		# an array saying which branches have a bottom
		node_with_no_bottom = rep(TRUE, length(nodes))

		# these are the remaining branch bottoms that have not been associated yet
		bots_with_no_top = rep(TRUE, nrow(bots))
		
		bxx = bots$x
		byy = bots$y
		botsnodes = 1:nrow(bots)
		
		# Empty values to hold nodes
		chosen_node_bottoms_xx = rep(NA, length(xx))
		chosen_node_bottoms_yy = rep(NA, length(yy))
		i=0
		while (sum(node_with_no_bottom) > 1)
			{
			i=i+1
			#print(i)
			# look for branch bottom coordinates in a narrow window to the left of the node
			# basically (a) smallest slope and (b) closest distance in x
			
			## find next node to include (the rightmost internal node)
			maxvals = xx * node_with_no_bottom
			nextnode <- which( maxvals == max(maxvals, na.rm=TRUE) )[1]
			
			####################################
			# Find the best matching branch bottom/corner for each node
			####################################
			# This is trial-and-error, you may have to plink to find a function 
			# that works.
			# That, or do manual edits to the tree later...
			####################################
			ydist <- byy - yy[nextnode]
			xdist <- bxx - xx[nextnode]

			# Rank of the y distances
			rank_ydist = rank(abs(ydist))

			# calculate the slops
			xyslopes <- abs(ydist/xdist)
			
			# the best ancestor will have a low slope to the branch bottom, and a short negative distance in x
			xdist_neg = xdist
			xdist_neg[xdist > 0] = 0
			xdist_neg[xdist < 0] = -1 * xdist_neg[xdist < 0]
			# normalize to units of minimum absolute distance
			min_dist = (min(abs(xdist[xdist!=0]), na.rm=TRUE))
			xdist_neg_norm = (xdist_neg / min_dist)

			# short positive distances are less good (half as good) than short negative distances
			xdist_pos = xdist
			xdist_pos[xdist < 0] = 0
			xdist_pos[xdist > 0] = xdist_pos[xdist > 0]
			xdist_pos_norm = (xdist_pos / min_dist) * 100

			rank_xdist = rank_ydist
			rank_xdist[xdist <= 0] = 1
			rank_xdist[xdist > 0] = 10000000
			
			###########################
			# Plink here especially...
			###########################
			rank_slope = (xyslopes^2)
			#final_rank = rank_ydist * abs(ydist) + xyslopes^0.5 * xdist_neg_norm + xdist_pos_norm
			
			# This worked on dinos pt1
			#final_rank = 10000*(rank_ydist * abs(ydist)) + xyslopes^0.5 *  xdist_neg_norm * xdist_pos_norm

			# This worked on dinos pt1
			final_rank = 10000*(rank_ydist * abs(ydist)) + xyslopes^0.5 *  xdist_neg_norm * xdist_pos_norm

			###########################

			branch_bot_fits = final_rank			
			best_fit = which(branch_bot_fits == min(branch_bot_fits, na.rm=TRUE))[1]
			
			#bottom_to_add = botsnodes[bots_with_no_top][best_fit]
			bottom_to_add = botsnodes[best_fit]
			
			chosen_node_bottoms_xx[nextnode] = bxx[bottom_to_add]
			chosen_node_bottoms_yy[nextnode] = byy[bottom_to_add]
			
			bxx = bxx[-bottom_to_add]
			byy = byy[-bottom_to_add]
			
			# remove the node from the list needing branch bottoms
			node_with_no_bottom[nextnode] = FALSE
			bots_with_no_top[bottom_to_add] = FALSE
			
			
			}
		}
	
	tmpdata = cbind(nodes, xx, yy, chosen_node_bottoms_xx, chosen_node_bottoms_yy)
	#print(tmpdata)
	
	dtf = as.data.frame(tmpdata)
	plot(dtf$y, dtf$chosen_node_bottoms_yy)
	title("Y-coord of branch tops vs. branch bottoms.  *Should* be linear.")
	
	plot(c(xx, bots$x), c(yy, bots$y), pch="")
	points(xx, yy, pch="n")
	points(bots$x, bots$y, pch="b")
	title("Use this plot to check if branch bottoms match branch tops (tips/nodes)")
	segments(xx, yy, chosen_node_bottoms_xx, chosen_node_bottoms_yy)
	
	plot(c(xx, bots$x), c(yy, bots$y), pch="")
	text(xx, yy, label=paste("n", 1:length(xx), sep=""))
	text(bots$x, bots$y, labe=paste("b", 1:nrow(bots), sep=""))
	title("This plot has node numbers so you can manually hack it.  However problems are\n almost always due to missing nodes or duplicate nodes while digitizing.")
	segments(xx, yy, chosen_node_bottoms_xx, chosen_node_bottoms_yy)
	

	return(dtf)
	}



build_tree_using_corners <- function(xy3)
	{
	# define the tip.labels
	tip.labels = xy3$nodenames
	tip.labels = tip.labels[tip.labels != ""]
	if (!missing(tip.labels))
		{
		ntips <- length(tip.labels)
		}
	
	
	xx = xy3$x
	yy = xy3$y
	cx = xy3$cx
	cy = xy3$cy
	
	nodes <- 1:length(xx)
	is.tip <- nodes <= ntips

	# keep track of the nodes which are unlinked
	unlinked_nodes = rep(TRUE, length(nodes))

	
	# Checks (kinda) if internal nodes are ordered from left-to-right
	if (which.min(xx) != ntips+1)
		{
		## 
		print("ERROR: reorder nodes the way ape/phylo expects! (tips first, then internals in order from left-to-right.")
		#yy[internal] <- rev(yy[!is.tip])[order(xx[!is.tip])]
		#xx[internal] <- rev(yy[!is.tip])[order(xx[!is.tip])]
		}

	edges <- matrix(nrow=0,ncol=2)
	edge.length <- numeric(0)
	nnode <- length(xx)-ntips

	while (sum(unlinked_nodes) > 1)
		{
		## find next node to include (the rightmost internal node)
		nextnode <- which(!is.tip & xx==max(xx[!is.tip]))[1]

		## find daughters
		
		# get the distance (in y) to all of the other corners
		ydist <- yy-yy[nextnode]
		xdist <- xx-xx[nextnode]

		# Check if it's the root
		if ( is.na(cy[nextnode]) )
			{
			cy[nextnode] = yy[nextnode]
			cx[nextnode] = 0			# leftmost coordinate must be 0!
			}
		
		cydist <- yy-cy[nextnode]
		cxdist <- xx-cx[nextnode]

		
		# find the TWO tips closest to this internal node, which are RIGHT of this node
		# this only finds the CLOSEST tip in Y, we want the TWO closest tips!!
		#daughters <- which(is.tip & dist==min(dist[is.tip]))
		
		# rank the ydistances in the y direction
		rank_of_ydists = order(cydist)
		
		# rank the xdistances in the x direction
		rank_of_xdists = order(cxdist)
		
		# get the node numbers in order; delete from this list as they are eliminated
		nodes_to_keep = nodes
		
		# daughter nodes must be to the right (in x) of the nextnode
		# (and they must be unlinked)
		nodes_to_keep = nodes_to_keep[unlinked_nodes][xdist[unlinked_nodes] > 0]
		
		# daughter nodes should be the two closest corners to nextnode (in y, mostly AND x)
		absolute_dist_from_node = 100*abs(cydist[nodes_to_keep]) + 1*abs(cxdist[nodes_to_keep])
		
		# sort the distances
		absolute_dist_from_node_sorted = sort(absolute_dist_from_node)
		
		# take the 2nd smallest absolute distance
		y_abs_dist_tokeep = absolute_dist_from_node_sorted[2]
		
		nodes_to_keep_final = nodes_to_keep[absolute_dist_from_node <= y_abs_dist_tokeep]
		print(paste("Linking: #", nodes_to_keep_final[1], " ", tip.labels[nodes_to_keep_final[1]], ", #", nodes_to_keep_final[2], " ", tip.labels[nodes_to_keep_final[2]], sep=""))
		
		#daughters <- which(is.tip & dist==min(dist[is.tip]))
		daughters = nodes_to_keep_final

		## be careful with numeric fuzz?
		edges <- rbind(edges,
					   nodes[c(nextnode,daughters[1])],
					   nodes[c(nextnode,daughters[2])])
		edge.length <- c(edge.length,xx[daughters]-xx[nextnode])

		# add nextnode to the list of tips (which are not available nodes for the nextnode)
		is.tip[nextnode] <- TRUE

		# remove the daughters & coordinates from the list of available nodes
		unlinked_nodes[daughters] = FALSE
		print(sum(unlinked_nodes))
		
		#xx <- xx[-daughters]
		#yy <- yy[-daughters]

		
		# remove the daughters from the list of possible nodes to link
		#unlinked_nodes
		
		#is.tip <- is.tip[-daughters]
		#nodes <- nodes[-daughters]
		}
	tr <- list(tip.labels=tip.labels,
			 edge=edges,
			 edge.length=edge.length,
			 Nnode=nnode)

	class(tr) <- "phylo"
	tr <- reorder(tr)
	tr$tip.labels = tip.labels
	return(tr)
	}



# This attempts to build the tree without corners.
# This works OK for young nodes, but gets increasingly bad
# with older nodes, unless you have a perfectly symmetrical tree.
build_tree_without_using_corners <- function(xx, yy, tip.labels, poly=numeric(0), debug=FALSE)
	{
	# define the tips
	if (!missing(tip.labels))
		{
		ntips <- length(tip.labels)
		}
	nodes <- 1:length(xx)
	is.tip <- nodes <= ntips
	
	# keep track of the nodes which are unlinked
	unlinked_nodes = rep(TRUE, length(nodes))
	
	if (which.min(xx) != ntips+1)
		{
		## 
		print("ERROR: reorder nodes the way ape/phylo expects! (tips first, then internals in order from left-to-right.")
		#yy[internal] <- rev(yy[!is.tip])[order(xx[!is.tip])]
		#xx[internal] <- rev(yy[!is.tip])[order(xx[!is.tip])]
		}

	edges <- matrix(nrow=0,ncol=2)
	edge.length <- numeric(0)
	nnode <- length(xx)-ntips

	
	


	while (sum(unlinked_nodes) > 1)
		{
		## find next node to include (the rightmost internal node)
		nextnode <- which(!is.tip & xx==max(xx[!is.tip]))[1]

		## find daughters
		
		# get the distance (in y) to all of the other nodes
		ydist <- yy-yy[nextnode]
		xdist <- xx-xx[nextnode]
		
		# find the TWO tips closest to this internal node, which are RIGHT of this node
		# this only finds the CLOSEST tip in Y, we want the TWO closest tips!!
		#daughters <- which(is.tip & dist==min(dist[is.tip]))
		
		# rank the ydistances in the y direction
		rank_of_ydists = order(ydist)
		
		# rank the xdistances in the x direction
		rank_of_xdists = order(xdist)
		
		# get the node numbers in order; delete from this list as they are eliminated
		nodes_to_keep = nodes
		
		# daughter nodes must be to the right (in x) of the nextnode
		# (and they must be unlinked)
		nodes_to_keep = nodes_to_keep[unlinked_nodes][xdist[unlinked_nodes] > 0]
		
		# daughter nodes should be the two closest to nextnode (in y)
		absolute_dist_from_node = abs(ydist[nodes_to_keep])
		
		# sort the distances
		absolute_dist_from_node_sorted = sort(absolute_dist_from_node)
		
		# take the 2nd smallest absolute distance
		y_abs_dist_tokeep = absolute_dist_from_node_sorted[2]
		
		nodes_to_keep_final = nodes_to_keep[absolute_dist_from_node <= y_abs_dist_tokeep]
		print(paste("Linking: #", nodes_to_keep_final[1], " ", tip.labels[nodes_to_keep_final[1]], ", #", nodes_to_keep_final[2], " ", tip.labels[nodes_to_keep_final[2]], sep=""))
		
		#daughters <- which(is.tip & dist==min(dist[is.tip]))
		daughters = nodes_to_keep_final

		## be careful with numeric fuzz?
		edges <- rbind(edges,
					   nodes[c(nextnode,daughters[1])],
					   nodes[c(nextnode,daughters[2])])
		edge.length <- c(edge.length,xx[daughters]-xx[nextnode])

		# add nextnode to the list of tips (which are not available nodes for the nextnode)
		is.tip[nextnode] <- TRUE

		# remove the daughters & coordinates from the list of available nodes
		unlinked_nodes[daughters] = FALSE
		
		#xx <- xx[-daughters]
		#yy <- yy[-daughters]

		
		# remove the daughters from the list of possible nodes to link
		unlinked_nodes
		
		#is.tip <- is.tip[-daughters]
		#nodes <- nodes[-daughters]
		}
	zz <- list(tip.labels=tip.labels,
			 edge=edges,
			 edge.length=edge.length,
			 Nnode=nnode)

	class(zz) <- "phylo"
	zz <- reorder(zz)
	zz$tip.labels = tip.labels
	return(zz)
	}


treethief <- function(xypts)
	{
	
	## from ?plot.tree:
	cat("(((Strix_aluco:4.2,Asio_otus:4.2):3.1,",
		"Athene_noctua:7.3):6.3,Tyto_alba:13.5);",
		file = "ex.tre", sep = "\n")
	tree.owls <- read.tree("ex.tre")
	#plot(tree.owls)
	unlink("ex.tre") # delete the file "ex.tre"
	
	#plot(tree.owls)
	xy <- get("last_plot.phylo",envir=.PlotPhyloEnv)
	xx <- xy$xx
	yy <- xy$yy
	points(xx,yy,col="white",pch=16,cex=2)
	text(xx,yy,col=2,1:length(xx))
	
	
	newtree <- build.tree(xx,yy,tree.owls$tip.label)
	
	data(bird.orders)
	plot(bird.orders,show.node.label=TRUE)
	xy <- get("last_plot.phylo",envir=.PlotPhyloEnv)
	points(xx,yy,col="white",pch=16,cex=2)
	text(xx,yy,col=2,1:length(xx))
	
	xx <- xy$xx
	yy <- xy$yy
	newtree2 <- build.tree(xx,yy,bird.orders$tip.label)
	}








####################################
# charnames_to_nexus_text
####################################
# Take a list of names, add numbers & quotes, to produce NEXUS character matrix input, e.g.:
#
# 		1 'it_is1', 2 'it_is1_blank', 3 'it_is1a', 4 'it_is1b', 5 'it_is1c', 6 not1, 7 'not1_blank', 8 not1a, 9 not1b, 10 the1, 11 'the1_blank', 12 intstr1, 13 'intstr1_blank', 14 intstr1b, 15 intstr1c, 16 intstr2, 17 'intstr2_blank', 18 intstr2b, 19 intstr2c, 20 of1, 21 the2, 22 'species1_blank', 23 species1a, 24 species1b, 25 that1, 26 'that1_blank', 27 survive1, 28 'survive1_blank', 29 survive1a, 30 survive1b, 31 of1, 32 'it_is2', 33 not2, 34 'not2_blank', 35 the3, 36 species2, 37 species2a, 38 species2b, 39 that2, 40 'that2_blank', 41 that2a, 42 survives2a, 43 'survives2_blank', 44 but1, 45 'but1_blank', 46 but1a, 47 rather, 48 'rather_blank', 49 rathera, 50 'it_is3_blank', 51 'it_is3a', 52 'it_is3b', 53 'it_is3c', 54 the4, 55 'the4_blank', 56 species3, 57 'species3_blank', 58 species3a, 59 that3, 60 'that3_blank', 61 survives3, 62 is1, 63 'is1_blank', 64 the5, 65 'the5_blank', 66 one, 67 that4, 68 'that4_blank', 69 is2, 70 'able_best', 71 'able_best_blank', 72 'able_best_a', 73 'able_best_b', 74 to1, 75 adapt1, 76 'adapt1_blank', 77 adapt1a, 78 adapt1b, 79 to2, 80 'to2_blank', 81 and1, 82 'and1_blank', 83 to3, 84 adjust1, 85 'adjust1_blank', 86 adjust1a, 87 adjust1b, 88 best1, 89 'best1_blank', 90 to4, 91 the6, 92 'the6_blank', 93 changing, 94 'changing_blank', 95 'long_end', 96 'long_end_blank', 97 'long_enda', 98 'long_endb', 99 'long_endc', 100 'long_endd', 101 comma1, 102 comma2, 103 comma3, 104 comma4, 105 'intel_before_strong', 106 'General_Field', 107 paraphrase, 108 'Darwin_attribution', 109 'Darwin_subsource_simp', 110 'if_evo_theory', 111 'if_Darwins_theory', 112 'if_Origin', 113 'cited_source' ; 

charnames_to_nexus_text <- function(inarray, quotemark="")
	{
	tmpstr = c()
	for (i in 1 : (length(inarray)-1) )
		{
		tmptmp = paste(i, " ", quotemark, inarray[i], quotemark, ", ", sep="")
		tmpstr = paste(tmpstr, tmptmp, sep="")
		}
	# add the last one, without comma
	tmpstr = paste(tmpstr, i, " ", quotemark, inarray[i+1], quotemark, " ", sep="")
	return(tmpstr)
	}


chars_to_nexus_fmt <- function(charmatrix, namewidth=12, require_question_marks=FALSE, ambiguity_check=TRUE)
	{
	require(gdata)	# for "trim" function
	
	outstr = ""
	width = ncol(charmatrix)
	
	# Debugging
	#print(namewidth)
	
	for (i in 1:nrow(charmatrix))
		{
		
		# This line caused a bug
		#tmpstr1 = as.character(charmatrix$rownames[i])
		tmpstr1 = as.character(row.names(charmatrix)[i])


		len = nchar(tmpstr1, type="chars")
		numspaces_to_add = namewidth - len
		
		# Debugging
		#print(i)
		#print(tmpstr1)
		#print(len)
		#print(numspaces_to_add)
		spaces_list = rep(" ", numspaces_to_add)
		tmpstr2 = array_to_text(c(tmpstr1, spaces_list), spacer="")
		
		
		chars = charmatrix[i,2:width]
		tmpstr3 = array_to_text(inarray=chars, spacer="")
		
		if (ambiguity_check)
			{
			
			# check for "&" chars: convert 1&2 to (1 2)
			if (grepl("&", tmpstr3))
				{
				for (j in 1:(width-1))
					{
					if (grepl("&", chars[j]))
						{
						if (require_question_marks == FALSE)
							{
							# Code the character with the possible states
							char = strsplit(as.character(chars[j]), "&")
							newchars_str = array_to_text(char[[1]], spacer=" ")
							newchar = paste("(", newchars_str, ")", sep="")
							chars[j] = newchar
							cat("Ambiguous character detected & converted...", charmatrix[i,1], ", char #", j, ": ", newchar, "\n")
							}
						else
							{
							# Code the character as "?"
							chars[j] = "?"
							}
						}
					}
				tmpstr3 = array_to_text(inarray=chars, spacer="")
				}		
			} # End ambiguity check
		tmpline = paste("	", tmpstr2, tmpstr3, "\n", sep="") 
		outstr = paste(outstr, tmpline, sep="")
		}
	outstr = trim(outstr)
	#outstr = paste("	", outstr, sep="")
	return(outstr)
	}






charmatrix_to_nexus_file <- function(charmatrix, nexus_outfn, quotemark="", namewidth=11, require_question_marks=FALSE, ambiguity_check=TRUE)
	{
	# initialize array of strings
	# l = lines
	l = c()
	
	l[1] = "#NEXUS"
	l[2] = "[ written by Nick Matzke R script 'quotescr_v11.R', written 9/1/2010 ]"
	l[3] = ""
	l[4] = "BEGIN TAXA;"
	l[5] = "	TITLE Taxa;"
	l[6] = paste("	DIMENSIONS NTAX=", nrow(charmatrix), ";", sep="")
	l[7] = "TAXLABELS"
	tmpstr = array_to_text(charmatrix$rownames, spacer=" ")
	l[8] = paste("		", tmpstr, sep="")
	l[9] = "	;"
	l[10] = ""
	l[11] = "END;"
	l[12] = ""
	l[13] = ""
	l[14] = "BEGIN CHARACTERS;"
	l[15] = "	TITLE Character_Matrix;"
	l[16] = paste("	DIMENSIONS NCHAR=", (ncol(charmatrix)-1), ";", sep="")
	
	tmpstr = "	FORMAT DATATYPE = STANDARD GAP = - MISSING = ?;"
	
	# maybe we can skip this?
	#SYMBOLS =   0 1 2 3 4 5 6 7 8 9 A B C D";"
	l[17] = tmpstr
	l[18] = "	CHARSTATELABELS "
	tmpstr = charnames_to_nexus_text(inarray = names(charmatrix)[2:ncol(charmatrix)] , quotemark="")
	
	l[19] = paste("		", tmpstr, ";", sep="")
	
	l[20] = "	MATRIX"
	
	# assemble all the characters into one big string
	l[21] = chars_to_nexus_fmt(charmatrix, namewidth, require_question_marks, ambiguity_check)
	l[22] = ";"
	l[23] = ""
	l[24] = "END;"
	
	
	
	write.table(l, nexus_outfn, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	headf(nexus_outfn)
	tailf(nexus_outfn)
	
	cat("Character matrix output to NEXUS: '", nexus_outfn, "'\n")
	}







nexus_data_to_chardf <- function(nexusdata)
	{
	# Convert character data read in by APE's read.data.nexus
	# into a simple dataframe table (rows: taxa, columns: characters)
	# 
	# This table can then be written out to NEXUS or TNT
	
	# 
	chardf1 = as.data.frame(nexusdata)
	chardf2 = as.data.frame(t(chardf1))

	# Get the number of characters (columns)
	numchars = dim(chardf2)[2]

	# make the row names a column in the table
	tmp_rownames = rownames(chardf2)
	
	# blank them out
	rownames(chardf2) = NULL
	
	# add a column with the rownames
	chardf3 = cbind(tmp_rownames, chardf2)
	
	# make the column headers
	col_headers1 = paste("char", 1:numchars, sep="")
	col_headers2 = c("rownames", col_headers1)
	
	# put the column headers into the chardf
	names(chardf3) = col_headers2
	
	#
	cat("\n\nPreview of character dataframe (first 8 rows & 15 columns:\n\n")
	print(chardf3[1:8, 1:15])
	
	cat("\nTotal chardf dimensions: ", dim(chardf3)[1], " rows & ", dim(chardf3)[2], " columns.\n\n", sep="")
	
	return(chardf3)
	}







charnames_to_tnt_text <- function(inarray, quotemark="")
	{
	tmpstr = c()
	for (i in 1:length(inarray) )
		{
		tmptmp = paste("{", i-1, " ", quotemark, inarray[i], quotemark, ";\n", sep="")
		tmpstr = paste(tmpstr, tmptmp, sep="")
		}
	return(tmpstr)
	}


chars_to_tnt_fmt <- function(charmatrix, namewidth=12)
	{
	require(gdata)	# for "trim" function
	
	outstr = ""
	width = ncol(charmatrix)

	for (i in 1:nrow(charmatrix))
		{
		
		tmpstr1 = as.character(charmatrix$rownames[i])
		#len = nchar(tmpstr1, type="chars")
		#numspaces_to_add = namewidth - len
		#spaces_list = rep(" ", numspaces_to_add)
		tmpstr2 = array_to_text(c(tmpstr1, "	"), spacer="")
		
		
		chars = charmatrix[i,2:width]
		tmpstr3 = array_to_text(inarray=chars, spacer="")

		# check for "&" chars: convert 1&2 to (1 2)
		if (grepl("&", tmpstr3))
			{
			for (j in 1:(width-1))
				{
				if (grepl("&", chars[j]))
					{
					char = strsplit(as.character(chars[j]), "&")
					newchars_str = array_to_text(char[[1]], spacer=" ")
					newchar = paste("[", newchars_str, "]", sep="")
					chars[j] = newchar
					cat("Ambiguous character detected & converted...", charmatrix[i,1], ", char #", j, ": ", newchar, "\n")
					}
				}
			tmpstr3 = array_to_text(inarray=chars, spacer="")
			}		
		tmpline = paste("", tmpstr2, tmpstr3, "\n", sep="") 
		outstr = paste(outstr, tmpline, sep="")
		}
	outstr = trim(outstr)
	#outstr = paste("", outstr, sep="")
	return(outstr)
	}




charmatrix_to_TNT_file <- function(charmatrix, tnt_outfn, quotemark="")
	{
	# initialize array of strings
	# l = lines
	l = c()
	
	l[1] = "xread"
	l[2] = paste(ncol(charmatrix)-1, nrow(charmatrix)) 

	tmpstr = chars_to_tnt_fmt(charmatrix)
	
	l[3] = tmpstr
	l[4] = ";"

	l[5] = ""
	l[6] = "cnames"
	l[7] = charnames_to_tnt_text( names(charmatrix)[2:ncol(charmatrix)] )
	l[8] = ";"
	l[9] = ""
	l[10] = ""
	l[11] = "proc /;"
	l[12] = "comments 0"
	l[13] = ";"
	l[14] = ""
	l[15] = ""
	
	
	write.table(l, tnt_outfn, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	headf(tnt_outfn)
	tailf(tnt_outfn)
	
	cat("Character matrix output to TNT: '", tnt_outfn, "'\n")
	}




scale_2values_axis <- function(observed_nums, observed_axis_inbreaks, desired_axis)
	{
	xlims_v1_min = min(observed_nums)
	xlims_v1_max = max(observed_nums)
	xlims_min_inbreaks = min(observed_axis_inbreaks)
	xlims_max_inbreaks = max(observed_axis_inbreaks)
	xlims_min_desired = min(desired_axis)
	xlims_max_desired = max(desired_axis)

	x = c(xlims_v1_min, xlims_v1_max)
	y = c(xlims_min_inbreaks, xlims_max_inbreaks)
	reg_result = lm(y~x)
	
	slope = reg_result$coefficients[2]
	intercept = reg_result$coefficients[1]
	
	xlims_touse_min = slope * xlims_min_desired + intercept
	xlims_touse_max = slope * xlims_max_desired + intercept
	
	xlims_touse = c(xlims_touse_min, xlims_touse_max)
	
	return(xlims_touse)
	}


make_multihist <- function(l, desired_xlims, colors, numbreaks, tmptitle, tmp_xlabel=NULL, tmp_ylabel=NULL)
	{
	
	for (i in 1:length(l))
		{
		if (i==1)
			{
			hist(l[[i]], breaks=numbreaks, col=colors[i], border="white", xlim=desired_axis, main=tmptitle, add=FALSE, xlab=tmp_xlabel, ylab=tmp_ylabel)
			} else {
			hist(l[[i]], breaks=numbreaks, col=colors[i], border="white", xlim=desired_axis, add=TRUE, ylab=tmp_ylabel)
			}

		#col=colors, border=c("white"), xlim=xlims_touse, plot.it=FALSE)
		#plot(tmphist, col=colors[i], border="white", xlim=desired_axis)

		#tmphist = hist(l[[i]], breaks=numbreaks, plot=FALSE)
		#plot(tmphist, col=colors[i], border="white", add=TRUE)
		}
	
	#border=c("white"),
	
	#offsets = min(tmphist$breaks)
	#barwidth = tmphist$breaks[1]-tmphist$breaks[0]
	
	#xrange = c(0, max(tmphist$breaks))
	
	
	#barplot(tmphist$out, offset=offsets), width=barwidth, xlim=c(min(xrange), max(xrange)), ylim=c(0, 1.2*max(tmphist$out)))
	#, beside=TRUE, col=colors, xlim=desired_axis, ylim=c(0, 1.2*max(tmphist$out)), bty="o", xaxt="n")
	
	#axis(1, at=pretty(), label=pretty(desired_axis), tick=TRUE, tcl=ticklength, pos=0)	
	}


# Limit each numeric column in a dataframe to e.g. 4 significant digits
signif_digits_df <- function(dtf, numsig=4, printout=FALSE)
	{
	dtf_classes = cls.df(dtf, printout)
	
	for (i in 1:length(dtf_classes$cls_col_list))
		{
		cls_col = dtf_classes$cls_col_list[i]
		if (cls_col == "numeric")
			{
			# Convert numeric column to one with signif digits
			colname = dtf_classes$dtf_names[i]
			cmdstr = paste("dtf$", colname, " = signif(dtf$", colname, ", digits=", numsig, ")", sep="")
			eval(parse(text = cmdstr))				
			}
		}
	return(dtf)
	}

setup = '
char=";"
'

# remove ending character from each line of a file
remove_ending_character_fn <- function(fn, outfn, char=";")
	{
	# Read file to tmplines
	tmplines = scan(fn, what="character", sep="\n")
	
	# Remove the last character
	for (i in 1:length(tmplines))
		{
		tmplines[i] = remove_ending_character(tmplines[i], char)
		}
	
	write_table_good_no_header(tmplines, outfn)
	headf(outfn)
	return(outfn)
	}

# remove ending character, if it equals "char"
remove_ending_character <- function(tmpline, char=";")
	{
	tmpchars = strsplit(tmpline, split="")[[1]]
	
	indices_that_match = match_list1_in_list2_return_all_indices(char, tmpchars)
	if (sum(indices_that_match) > 0)
		{
		last_match = indices_that_match[length(indices_that_match)]
	
		if (tmpchars[last_match] == tmpchars[length(tmpchars)])
			{
			newline = substr(tmpline, start=1, stop=nchar(tmpline)-1)
			} else {
			newline = tmpline
			}
		} else {
		newline = tmpline
		}
	
	return(newline)
	}

stripquotes_spaces <- function(tmpstr)
	{
	require(gdata)
	
	# Remove quotes
	tmpstr = gsub(pattern='\\"', replacement="", x=tmpstr)

	# Remove quotes
	tmpstr = gsub(pattern="'", replacement="", x=tmpstr)


	# Remove \
	tmpstr = gsub(pattern='\\\\', replacement="", x=tmpstr)


	# trim whitespace
	tmpstr = trim(tmpstr)
	
	return(tmpstr)
	}


#stripquotes_spaces_df <- function(tmpdf)
#	{
#	tmpdf = apply(tmpdf, 2, stripquotes_spaces)
#	return(tmpdf)
#	}
	
read_table_good <- function(fn, tmpskip=0)
	{
	# Read in the data, store in variable d
	# This has all of Nick's preferred options
	dtf = read.table(fn, header=TRUE, skip=tmpskip, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
	return(dtf)
	}

read_table_good_w_rownames <- function(fn, tmpskip=0)
	{
	# Read in the data, store in variable d
	# This has all of Nick's preferred options
	dtf = read.table(fn, header=TRUE, skip=tmpskip, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE, row.names=1)
	return(dtf)
	}


read_table_good_no_header <- function(fn, tmpskip=0)
	{
	# Read in the data, store in variable d
	# This has all of Nick's preferred options
	dtf = read.table(fn, header=FALSE, skip=tmpskip, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE)
	return(dtf)
	}

write_table_good <- function(dtf, outfn, ...)
	{
	write.table(x=dtf, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=FALSE, col.names=TRUE, ...)
	}

write_table_good_no_header <- function(dtf, outfn, sepval="\n")
	{
	write.table(dtf, file=outfn, append=FALSE, quote=FALSE, sep=sepval, row.names=FALSE, col.names=FALSE)
	}


write_table_good_with_rownames <- function(dtf, outfn)
	{
	write.table(dtf, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=TRUE, col.names=TRUE)
	}

write_table_good_no_colnames_or_rownames <- function(dtf, outfn)
	{
	write.table(dtf, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=FALSE, col.names=FALSE)
	}


write_line_good <- function(dtf, outfn)
	{
	write.table(dtf, file=outfn, append=TRUE, quote=FALSE, sep="	", row.names=FALSE, col.names=TRUE)
	}

write_lines_good <- function(dtf, outfn, sepval="\n", tmpappend=FALSE)
	{
	write.table(dtf, file=outfn, append=tmpappend, quote=FALSE, sep=sepval, row.names=FALSE, col.names=FALSE)
	}

write_table_good_w_rownames <- function(dtf, outfn, ...)
	{
	write.table(dtf, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=TRUE, col.names=TRUE, ...)
	}



write_table_good_w_rownames_yesappend <- function(dtf, outfn, ...)
	{
	write.table(dtf, file=outfn, append=TRUE, quote=FALSE, sep="	", row.names=TRUE, col.names=TRUE, ...)
	}


write_table_good_w_rownames_no_colnames <- function(dtf, outfn, ...)
	{
	write.table(dtf, file=outfn, append=FALSE, quote=FALSE, sep="	", row.names=TRUE, col.names=FALSE)
	}





# Confusion matrix
make_confusion_matrix <- function(truth, predicted)
	{
		
	#########################################
	# Make a confusion matrix
	#########################################
	#
	# Following:
	# NASA Remote Sensing Tutorial
	# http://rst.gsfc.nasa.gov/Sect13/Sect13_3.html
	# http://rst.gsfc.nasa.gov/Sect13/originals/Fig13_7.jpg
	# 
	# Rows are "truth" (from ground/high-res image interpretation)
	# Cols are "predicted" (from Landsat)
	#
	
	total_num_cells = length(predicted)
	
	unique_classes = sort(unique(c(truth, predicted)))
	
	confusion_matrix = matrix(data=NA, nrow=length(unique_classes)+3, ncol=length(unique_classes)+4)
	
	numclasses = length(unique_classes)
	
	for (i in 1:numclasses)
		{
		for (j in 1:numclasses)
			{
			truth_val = unique_classes[i]
			pred_val = unique_classes[j]
			
			truth_matching_TF = truth == truth_val
			pred_matching_TF = predicted == pred_val
			
			cell_val = sum(truth_matching_TF + pred_matching_TF == 2)
			
			confusion_matrix[i,j] = cell_val
			}
		}
	confusion_matrix
	
	
	# Do the totals -- totals of columns
	for (i in 1:numclasses)
		{
		confusion_matrix[numclasses+1,i] = sum(confusion_matrix[,i], na.rm=TRUE)
		}
	confusion_matrix
	
	
	# Do the totals -- totals of rows
	for (i in 1:(numclasses+1))
		{
		confusion_matrix[i, numclasses+1] = sum(confusion_matrix[i,], na.rm=TRUE)
		}
	confusion_matrix
	
	
	# Do the commission and omission errors
	# omissions = errors of omission = pixels missed that should have been found (false negatives)
	# 1-ommission error = Producers's Accuracy
	#
	# comissions = errors of comission = pixels put in the category that should NOT have been (false positives)
	# 1-commission error = User's Accuracy
	#	
	# omission errors -- across rows
	# 1-ommission error = Producers's Accuracy
	omissions_ttl = 0
	for (i in 1:numclasses)
		{
		total = confusion_matrix[i, numclasses+1]
		number_of_hits = confusion_matrix[i,i]
		number_of_misses = total - number_of_hits
		omissions_ttl = omissions_ttl + number_of_misses
		confusion_matrix[i, numclasses+2] = number_of_misses / total * 100
		confusion_matrix[i, numclasses+3] = (1 - number_of_misses / total) * 100
		}
	confusion_matrix[numclasses+1, numclasses+2] = omissions_ttl/total_num_cells * 100
	confusion_matrix[numclasses+1, numclasses+3] = 100 - omissions_ttl/total_num_cells * 100
	
	
	# Comission errors -- down columns
	# 1-commission error = User's Accuracy
	comissions_ttl = 0
	for (i in 1:numclasses)
		{
		total = confusion_matrix[numclasses+1, i]
		number_of_hits = confusion_matrix[i,i]
		number_of_misses = total - number_of_hits
		comissions_ttl = comissions_ttl + number_of_misses
		confusion_matrix[numclasses+2, i] = number_of_misses / total * 100
		confusion_matrix[numclasses+3, i] = (1 - number_of_misses / total) * 100
		}
	confusion_matrix[numclasses+2, numclasses+1] = comissions_ttl/total_num_cells * 100
	confusion_matrix[numclasses+3, numclasses+1] = 100 - comissions_ttl/total_num_cells * 100
	
	
	# Mapping accuracy -- across rows
	for (i in 1:numclasses)
		{
		total = confusion_matrix[i, numclasses+1]
		number_of_hits = confusion_matrix[i,i]
		confusion_matrix[i, numclasses+4] = number_of_hits / total * 100
		}
	confusion_matrix
	
	
	# Overall accuracy (total number correct / total)
	correct_ttl = 0
	for (i in 1:numclasses)
		{
		correct_ttl = correct_ttl + confusion_matrix[i,i]
		}
	total = confusion_matrix[numclasses+1, numclasses+1]
	mapping_accuracy = correct_ttl / total * 100
	confusion_matrix[nrow(confusion_matrix), ncol(confusion_matrix)] = mapping_accuracy
	
	confusion_matrix = adf(confusion_matrix)
	
	names1 = paste("pred_", as.character(unique_classes), sep="")
	names2 = paste("true_", as.character(unique_classes), sep="")
	names(confusion_matrix) = c(names1, "total", "omission", "users_accuracy", "mapping_accuracy")
	row.names(confusion_matrix) = c(names2, "total", "comission", "producers_accuracy")
	confusion_matrix

	return(confusion_matrix)
	}

# Get best guess for each row
# (e.g., prediction with the highest probability)
get_highest_val <- function(preds)
	{
	tmpnames = names(adf(preds))
	
	maxvals = apply(X=preds, MARGIN=1, FUN=max)
	
	colnums = 1:length(tmpnames)
	
	colnames_matrix = maxvals
	
	for (i in 1:nrow(preds))
		{
		TF = preds[i, ] == maxvals[i]
		colnames_matrix[i] = tmpnames[TF]
		}
	return(colnames_matrix)
	}


calc_kappa_from_confusion_matrix <- function(confusion_matrix, numclasses)
	{
	#
	# Kappa is basically: 
	# ( observed agreement - chance agreement ) /
	# (                  1 - chance agreement )
	# Following:
	# http://kogs-www.informatik.uni-hamburg.de/PROJECTS/censis/satsim/node10.html
	#
	# Better: Congalton 1991
	#
	# Congalton, R. G. 1991. A review of assessing the accuracy
	# of classifications of remotely sensed data. Remote Sensing
	# of Environment 37:35-46.
	# http://dx.doi.org/10.1016/0034-4257(91)90048-B
	# 
	# Interpretation:
	# http://en.wikipedia.org/wiki/Cohen%27s_kappa
	#
	# Nonetheless, magnitude guidelines have appeared in the literature. Perhaps the
	# first was Landis and Koch,[11] who characterized values < 0 as indicating no
	# agreement and 0–.20 as slight, .21–.40 as fair, .41–.60 as moderate, .61–.80
	# as substantial, and .81–1 as almost perfect agreement. This set of guidelines
	# is however by no means universally accepted; Landis and Koch supplied no evidence
	# to support it, basing it instead on personal opinion. It has been noted that these
	# guidelines may be more harmful than helpful.[2] Fleiss's[12]:218 equally arbitrary
	# guidelines characterize kappas over .75 as excellent, .40 to .75 as fair to good,
	# and below .40 as poor.
	# 
	# A very clear presentation:
	# http://www.geocomputation.org/1999/044/gc_044.htm
	#
	# For variance of Kappa, see especially Congalton 1999, p. 106
	# http://www.amazon.com/Assessing-Accuracy-Remotely-Sensed-Data/dp/1420055127#reader_1420055127
	#
	# Get the total # of pixels
	N = confusion_matrix[numclasses+1,numclasses+1]
	
	# Get p0
	p0 = 0
	for (i in 1:numclasses)
		{
		p0 = p0 + confusion_matrix[i,i]
		}
	p0 = p0
	
	# Get pz
	pz = 0
	for (i in 1:numclasses)
		{
		#sum1
		sum1 = 0
		sum2 = 0
		for (j in 1:numclasses)
			{
			sum1 = sum1 + confusion_matrix[i,j]
			sum2 = sum2 + confusion_matrix[j,i]
			}
		pz = pz + (sum1 * sum2)
		}
	#pz = pz * (1/ (N^2))	
	
	
	kappaval = ((N * p0) - pz) / ((N^2) - pz)
	return(kappaval)
	}




calc_kappavar_from_confusion_matrix <- function(confusion_matrix, numclasses)
	{
	# Following:
	# http://kogs-www.informatik.uni-hamburg.de/PROJECTS/censis/satsim/node10.html
	
	# Better: Congalton 1991
	#
	# Congalton, R. G. 1991. A review of assessing the accuracy
	# of classifications of remotely sensed data. Remote Sensing
	# of Environment 37:35-46.
	# http://dx.doi.org/10.1016/0034-4257(91)90048-B
	# 
	# Interpretation:
	# http://en.wikipedia.org/wiki/Cohen%27s_kappa
	#
	# Nonetheless, magnitude guidelines have appeared in the literature. Perhaps the
	# first was Landis and Koch,[11] who characterized values < 0 as indicating no
	# agreement and 0–.20 as slight, .21–.40 as fair, .41–.60 as moderate, .61–.80
	# as substantial, and .81–1 as almost perfect agreement. This set of guidelines
	# is however by no means universally accepted; Landis and Koch supplied no evidence
	# to support it, basing it instead on personal opinion. It has been noted that these
	# guidelines may be more harmful than helpful.[2] Fleiss's[12]:218 equally arbitrary
	# guidelines characterize kappas over .75 as excellent, .40 to .75 as fair to good,
	# and below .40 as poor.
	# 
	# A very clear presentation:
	# http://www.geocomputation.org/1999/044/gc_044.htm
	#
	# For variance of Kappa, see especially Congalton 1999, p. 106
	# http://www.amazon.com/Assessing-Accuracy-Remotely-Sensed-Data/dp/1420055127#reader_1420055127
	#
	# Get the total # of pixels
	n = confusion_matrix[numclasses+1,numclasses+1]
	
	# Get theta1
	theta1 = 0
	for (i in 1:numclasses)
		{
		theta1 = theta1 + confusion_matrix[i,i]
		}
	theta1 = 1/n * theta1


	# Get theta2
	theta2 = 0
	for (i in 1:numclasses)
		{
		#sum1
		sum1 = 0
		sum2 = 0
		for (j in 1:numclasses)
			{
			sum1 = sum1 + confusion_matrix[i,j]
			sum2 = sum2 + confusion_matrix[j,i]
			}
		theta2 = theta2 + (sum1 * sum2)
		}
	theta2 = 1/(n^2) * theta2


	# Get theta3
	theta3 = 0
	for (i in 1:numclasses)
		{
		nii = confusion_matrix[i,i]
		#sum1
		sum1 = 0
		sum2 = 0
		for (j in 1:numclasses)
			{
			sum1 = sum1 + confusion_matrix[i,j]
			sum2 = sum2 + confusion_matrix[j,i]
			}
		theta3 = theta3 + nii * (sum1 + sum2)
		}
	theta3 = 1/(n^2) * theta3
	

	# Get theta4
	theta4 = 0
	isum = 0
	for (i in 1:numclasses)
		{
		nii = confusion_matrix[i,i]
		#sum1
		
		jsum = 0
		for (j in 1:numclasses)
			{
			sum1 = 0
			for (ii in 1:numclasses)
				{
				sum1 = sum1 + confusion_matrix[j,ii]
				sum2 = sum2 + confusion_matrix[ii,i]
				}
			nij = confusion_matrix[i,j]
			jsum = jsum + (nii * (sum1 + sum2)^2)
			}
		isum = isum + jsum
		}
	theta4 = 1/(n^3) * isum
		
	
	# Super-formula in multiple parts
	part1 = theta1*(1-theta1) / ((1-theta2)^2)
	
	part2num = 2 * (1-theta1) * ((2*theta1*theta2) - theta3)
	part2den = (1-theta2)^3
	
	part3num = ((1-theta1)^2) * (theta4 - 4*(theta2^2))
	part3den = (1-theta2)^4
	
	varK = (1/n) * (part1 + (part2num / part2den) + (part3num / part3den))
	
	return(varK)
	}




is.not.NA <- function(x)
	{
	return(is.na(x) == FALSE)
	}
# 
# is.not.na <- function(x)
# 	{
# 	return(is.na(x) == FALSE)
# 	}


# remove NA rows from list
na.remove <- function(x)
	{
	outx = x[is.not.NA(x)]
	return(outx)
	}
	
allneg <- function(inlist)
	{
	x = (inlist <= 0)
	sumneg = sum(x, na.rm=TRUE)
	lenneg = length(na.remove(x))
	if (lenneg == sumneg)
		{
		return(TRUE)
		}
	else
		{
		return(FALSE)
		}
	}

allpos <- function(inlist)
	{
	x = (inlist >= 0)
	sumpos = sum(x, na.rm=TRUE)
	lenpos = length(na.remove(x))
	if (lenpos == sumpos)
		{
		return(TRUE)
		}
	else
		{
		return(FALSE)
		}
	}


remove_NA_rows <- function(dtf, colname)
	{
	cat("\n")
	cat("remove_NA_rows: initial #rows = ", nrow(dtf), "\n")

	cmdstr = paste("rows_to_keep = is.not.NA(dtf$", colname, ")", sep="")
	#print(cmdstr)
	eval(parse(text = cmdstr))
	
	#print(dim(df))
	#print(length(rows_to_keep))
	#print(rows_to_keep)
	
	outdf = dtf[rows_to_keep, ]

	cat("remove_NA_rows: ending #rows = ", nrow(outdf), "\n")
	return(outdf)
	}




# check if all elements match
do_all_match <- function(arr1, arr2)
	{
	matches = (arr1 == arr2)
	if (sum(matches) == length(matches))
		{
		return(TRUE)
		} else {
		return(FALSE)
		}
	}



# split strings on whitespace
strsplit_on_tabs_remove_whitespace <- function(tmpline)
	{
	# split on 1 or more whitespaces
	temp = strsplit(tmpline, "[\t]")
	
	# get the list
	list_of_strs = temp[[1]]
	
	# remove any leading/trailing ""
	list_of_strs = list_of_strs[list_of_strs != ""]
	
	for (i in 1:length(list_of_strs))
		{
		tmpstr = list_of_strs[i]
		list_of_strs[i] =  gsub(" ", "", tmpstr)
		}
		
	return(list_of_strs)
	}



# catch "integer(0)" e.g. from non-matching grep
is.numzero <- function(x)
	{
	if (length(x) == 0)
		{
		return(TRUE)
		}
	else
		{
		return(FALSE)
		}
	}

is.inf <- function(x)
	{
	return(x == Inf)
	}


is.not.inf <- function(x)
	{
	return(is.inf(x) == FALSE)
	}


# return_items_not_NA <- function(x)
# 	{
# 	y = x[is.not.NA(x)]
# 	return(y)
# 	}

plot_point <- function(x, y, pt_title, color)
	{
	points(x, y, type="p", pch=19, cex=2, col=color)
	text(x, y, pt_title, col=color, pos=4, bg="white")

	}


# print out the last 5 rows of a data frame
foot <- function(dtf)
	{
	nrows = nrow(dtf)
	print(dtf[((nrows-5) : nrows), ])
	}
	
# print out the header, footer, and basic info on a dataframe (df)
sumdf <- function(dtf)
	{
	print(dim(dtf))
	print(head(dtf))
	foot(dtf)
	
	# print the class of each column
	for (col in 1:ncol(dtf))
		{
		cat("Col#", col, " class = ", class(dtf[1,col]), "\n", sep="")
		}
	
	}


	
make_cdf <- function(list_rel_probs)
	{
	if (sum(list_rel_probs) == 0)
		{
		print("make_cdf() ERROR: sum(list_rel_probs) is ZERO")
		print(list_rel_probs)
		return(NA)
		}
	
	if (class(list_rel_probs) != "data.frame")
		{
		dim(list_rel_probs) = c(length(list_rel_probs), 1)
		list_rel_probs = data.frame(list_rel_probs)
		}
	cdf_rel_probs = rep(0, nrow(list_rel_probs))
	dim(cdf_rel_probs) = c(length(cdf_rel_probs), 1)
	cdf_rel_probs = data.frame(cdf_rel_probs)
	cdf_rel_probs = cbind(cdf_rel_probs, cdf_rel_probs)
	names(cdf_rel_probs) = c("bottom", "top")
	
	for (i in 1:nrow(list_rel_probs))
		{
		if (i==1)
			{
			cdf_rel_probs$bottom[i] = 0
			cdf_rel_probs$top[i] = list_rel_probs[i, 1]
			}
		else
			{
			cdf_rel_probs$bottom[i] = cdf_rel_probs$top[i-1]
			cdf_rel_probs$top[i] = cdf_rel_probs$top[i-1] + list_rel_probs[i, 1]
			}
		}
	
	cdf_rel_probs_new = cdf_rel_probs / sum(list_rel_probs)
	return(cdf_rel_probs_new)
	}


# empirical p-values (for empirical PDF, cdf, percentiles, quantiles)
empirical_pval_fast <- function(x, observed_val)
	{
	num_lt_observed = sum(x <= observed_val)
	
	pval = num_lt_observed / length(x)
	
	return(pval)
	}


# empirical p-values (for empirical PDF, cdf, percentiles, quantiles)
empirical_pval_slow <- function(x, observed_val)
	{
	freqtab = f.freq(x)
	
	pval = freqtab$cumul_fraction[freqtab$x == observed_val]
	
	return(pval)
	}




f.freq <- function(x)
	{
	freqtab <- data.frame(table(x))
	freqtab$fraction <- freqtab$Freq / length(x)
	freqtab$cumul_fraction[1] <- freqtab$fraction[1]
	for(i in 2:length(freqtab[,1]))
		{
		freqtab$cumul_fraction[i] = freqtab$cumul_fraction[i-1] + freqtab$fraction[i]
		}
	return(freqtab)
	}



# Convert unique items to characters
unique_items_to_chars <- function(tempy, y, namestr)
	{
	col_w_name = match(namestr, names(tempy), nomatch=NaN)
	tmpname = names(y)[col_w_name]
	cmdstr = paste("charstates = unique(tempy$", tmpname, ")", sep="")
	eval(parse(text = cmdstr))	
	charstates = charstates[charstates != "?"]
	charindexes = seq(1, length(charstates))
	cat("\n")
	cat(namestr, " - ", "# of states: ", length(charstates), "\n", sep="")
	cat("code	state\n")
	# Display the characters
	for (i in 1:length(charstates))
		{
		#y$General_Field[tempy$General_Field == charstates[i]] = mesquite_character_states[charindexes[i]]
		cat(mesquite_character_states[charindexes[i]], ":	", charstates[i], "\n", sep="")
		cmdstr = paste("y$", tmpname, "[tempy$", tmpname, " == charstates[i]] = mesquite_character_states[charindexes[i]]", sep="")
		eval(parse(text = cmdstr))
		}
	return(y)
	}


# Retrieve chars
# Note: only works if there is a 1-to-1 mapping of codes to text
# 
# inchars_outchars = retrieve_char_text(tempy, y, namestr)
# inchars = inchars_outchars[[1]]
# outchars = inchars_outchars[[2]]
retrieve_char_text <- function(tempy, y, namestr)
	{
	tmpname = namestr
	
	# Get numerical characters
	cmdstr = paste("charstates = unique(tempy$", tmpname, ")", sep="")
	eval(parse(text = cmdstr))	

	# Get text states
	cmdstr = paste("textstates = unique(y$", tmpname, ")", sep="")
	eval(parse(text = cmdstr))	
	
	inchars = c()
	outchars = c()
	for (i in 1:length(charstates))
		{
		charstate = charstates[i]
		
		charstate_index = match(charstate, charstates)
		textstate = textstates[charstate_index]
		
		inchars = c(inchars, charstate)
		outchars = c(outchars, textstate)
		
		}
		
	inchars_outchars = list(inchars, outchars)
	return(inchars_outchars)
	}



# Count the number of matching characters
count_chars <- function(X, char=",")
	{
	# Count the number of matching characters in X,
	# where X is a string or list of strings
	numchars_that_match = count.fields(textConnection(X), sep=char)
	
	return(numchars_that_match)
	}


print_chars <- function(y, namestr)
	{
	col_w_name = match(namestr, names(y), nomatch=NaN)
	tmpname = names(y)[col_w_name]
	cmdstr = paste("charstates = unique(y$", tmpname, ")", sep="")
	eval(parse(text = cmdstr))	
	cat("\n")
	cat(namestr, " - ", "# of states: ", length(charstates), "\n", sep="")
	cat("List of the character states:\n")
	# Display the characters
	for (i in 1:length(charstates))
		{
		cat(charstates[i], "\n", sep="")
		}	
	}

print_char_key <- function(inchars, outchars)
	{
	cat("\n")
	cat("Printing codes for characters:\n")
	for (i in 1:length(inchars))
		{
		cat(outchars[i], ": ", inchars[i], "\n", sep="")
		}
	}

print_chars_catlist <- function(y, namestr)
	{
	mesquite_character_states = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "K", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "a", "b", "c", "d", "e", "f", "g", "h", "l", "m", "n", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")

	
	col_w_name = match(namestr, names(y), nomatch=NaN)
	tmpname = names(y)[col_w_name]
	cmdstr = paste("charstates = unique(y$", tmpname, ")", sep="")
	eval(parse(text = cmdstr))	
	cat("\n")
	cat(namestr, " - ", "# of states: ", length(charstates), "\n", sep="")
	cat("\n")
	cat("List of the character states:\n")


	# Display the characters -- text input
	cat("inchars = c(")
	for (i in 1:(length(charstates)-1))
		{
		cat('"', charstates[i], '", \n', sep="")
		}	
	cat('"', charstates[length(charstates)], '")\n', sep="")
	cat("\n")

	# Display the characters -- suggested coded output
	cat("outchars = c(")
	charcount = 0
	for (i in 1:(length(charstates)-1))
		{
		if (charstates[i] != "?")
			{
			charcount = charcount + 1
			cat('"', mesquite_character_states[charcount], '", \n', sep="")
			}
		else
			{
			cat('"?"', ', \n', sep="")			
			}
		}	
	if (charstates[length(charstates)] != "?")
		{
		charcount = charcount + 1
		cat('"', mesquite_character_states[charcount], '")\n', sep="")
		}
	else
		{
		cat('"?"', ')\n', sep="")			
		}
	cat("\n")


	}


convert_items_to_items <- function(y, namestr, name1list, name2list)
	{
	for (i in 1:length(name1list))
		{
		name1 = name1list[i]
		name2 = name2list[i]
		y = convert_item_to_item(y, namestr, name1, name2)
		}
	return(y)
	}


convert_item_to_item <- function(y, namestr, name1, name2)
	{
	col_w_name = match(namestr, names(y), nomatch=NaN)
	tmpname = names(y)[col_w_name]
	cmdstr = paste("y$", tmpname, "[y$", tmpname, ' == "', name1, '"] = "', name2, '"', sep="")
	print(cmdstr)
	eval(parse(text = cmdstr))	
	return(y)
	}

convert_grepl_item_to_item <- function(y, namestr, name1, name2)
	{
	col_w_name = match(namestr, names(y), nomatch=NaN)
	tmpname = names(y)[col_w_name]
	cmdstr = paste("rows_to_change = grepl(name1, y$", tmpname, ")", sep="")
	eval(parse(text = cmdstr))
	cmdstr = paste("y$", tmpname, "[rows_to_change] = '", name2, "'", sep="")
	print(cmdstr)
	eval(parse(text = cmdstr))	
	return(y)
	}



# Extract rows with indices matching rows_to_keep
gather_rows_with_indices_old <- function(tmptable, rows_to_keep)
	{
	# make empty matrix
	newtable = matrix(data=NA, nrow = length(rows_to_keep), ncol=ncol(tmptable))
	names(newtable) = names(tmptable)
	for (i in 1:nrow(newtable))
		{
		# insert the new row 
		newtable = tmptable[rows_to_keep[i], ]
		}
	return(newtable)
	}



# Extract rows with indices matching rows_to_keep
# Fixing a bug that was caused by row numbers being characters instead of
# numbers
gather_rows_with_indices <- function(tmptable, rows_to_keep)
	{
	rows_to_keep = as.numeric(rows_to_keep)
	
	# make empty matrix
	newtable = matrix(data=NA, nrow = length(rows_to_keep), ncol=ncol(tmptable))
	names(newtable) = names(tmptable)
	for (i in 1:nrow(newtable))
		{
		# insert the new row 
		newtable = tmptable[rows_to_keep[i], ]
		}
	return(newtable)
	}


# remove rows from temp_table, based on col matching list_to_remove
remove_rows_with_list_to_remove <- function(temp_table, col, list_to_remove)
	{
	# get TRUE/FALSE for which rows in species col are undefined, according to list_to_remove
	truefalse_rows_match = match_list1_in_list2(col, list_to_remove)
	cat("remove_rows_with_list_to_remove(): initial table #rows = ", nrow(temp_table), "\n")
	
	cat("Removing ", (sum(truefalse_rows_match)), " rows with c(", list2str(c(list_to_remove)), ") in them.\n", sep="")
	temp_table = temp_table[truefalse_rows_match == FALSE, ]
	cat("remove_rows_with_list_to_remove(): revised table #rows = ", nrow(temp_table), "\n")
		
	return(temp_table)
	}	



# This should be a faster list2 string
# This should be a faster list2 string
list2str_fast <- function(list1, spacer="")
	{
	# convert to character
	# split into list of lists of characters
	# merge lists with unlist
	# paste with collapse argument
	tmpstr = paste(unlist(strsplit(as.character(list1), split="")), collapse=spacer)
	return(tmpstr)
	}


list2str_fast_nosplit <- function(list1, spacer="")
	{
	# convert to character
	# split into list of lists of characters
	# merge lists with unlist
	# paste with collapse argument
	tmpstr = paste(unlist(as.character(list1)), collapse=spacer)
	return(tmpstr)
	}



list2str_unlist <- function(list1, spacer="	")
	{
	outlist = c("startlist")
	for (i in 1:length(list1))
		{
		listitem = list1[i]
		if (length(listitem) > 1)
			{
			outlist = append(outlist, list2str(unlist(listitem), spacer="	"))
			}
		else
			{
			addstr = as.character(list1[i])
			#tmpstr = paste(tmpstr, addstr, sep=spacer)
			outlist = append(outlist, addstr)
			}
		}
	tmpstr = list2str(outlist, spacer="	")
	return(tmpstr)
	}


return_unlist <- function(list1, spacer="	")
	{
	outlist = c("startlist")
	for (i in 1:length(list1))
		{
		listitem = list1[i]
		if (length(listitem) > 1)
			{
			outlist = append(outlist, unlist(listitem), spacer="	")
			}
		else
			{
			addstr = as.character(list1[i])
			#tmpstr = paste(tmpstr, addstr, sep=spacer)
			outlist = append(outlist, addstr)
			}
		}
	
	return(outlist)
	}


# unlist_dtf_cols <- function(dtf, printflag=FALSE)
# 	{
# 	# Sometimes cbind makes each column a list, this can screw up use/searching of
# 	#  the column later on.  
# 	# Unlist each column...
# 	for (i in 1:ncol(dtf))
# 		{
# 		tmpstr = paste("unlisting col: ", names(dtf)[i], "...", sep="")
# 		prflag(tmpstr, printflag=printflag)		
# 		
# 		# catch a possible error from unlisting
# 		# the number of rows needs to stay the same!!
# 		tmpcol = unlist(dtf[, i])
# 		if (length(tmpcol) != length(dtf[, i]))
# 			{
# 			tmpstr = paste("...failed! unlist(col) length=", length(tmpcol), "; nrow(dtf) = ", nrow(dtf), sep="")
# 			prflag(tmpstr, printflag=printflag)
# 			} 
# 		else
# 			{
# 			dtf[, i] = tmpcol
# 			tmpstr = paste(" ", " ", sep="")
# 			prflag(tmpstr, printflag=printflag)
# 			}
# 		}
# 	
# 	#dtf2 = adf(dtf)
# 	
# 	return(dtf)
# 	}




char_descriptions_to_NEXUS <- function(namestr_list, inchars_list, outchars_list)
	{
	# Create l, an empty list of lines...
	l = list()
	
	#l[(h=h+1)] = "	CHARSTATELABELS"
	
	tmpcharstr_list = c()
	for (i in 1:length(namestr_list))
		{
		tmpname = namestr_list[i]

		# filter out ? and ''
		inchars = inchars_list[[i]]
		outchars = outchars_list[[i]]

		chars_to_keep1 = outchars != ""
		chars_to_keep2 = outchars != "?"
		
		# remove the character states which are "" or "?" in the coded section
		chars_to_keep = ((chars_to_keep2 + chars_to_keep1) == 2)
		inchars2 = inchars[chars_to_keep]

		oldstr = "'"
		newstr = ""	
		tmp_charstates_list_noapost = gsub(oldstr, newstr, inchars2)

		oldstr = " "
		newstr = "_"	
		tmp_charstates_list_nospace = gsub(oldstr, newstr, tmp_charstates_list_noapost)

		
		#tmp_charstates_str = list2str(inchars_list[[i]], spacer='" "')
		tmp_charstates_str = list2str(tmp_charstates_list_nospace, spacer="' '")
		
		tmp_charstates_str2 = paste("'", tmp_charstates_str, "'", sep="")
		
		tmpcharstr = paste(i, tmpname, "/ ", tmp_charstates_str2)
		
		tmpcharstr_list = c(tmpcharstr_list, tmpcharstr)
		}
	
	bigstr_chardescs = list2str(tmpcharstr_list, spacer=", ")
	charstatelabels = paste("		", bigstr_chardescs, " ;", sep="")
	
	
	return(charstatelabels)
	}


# Add character state labels to a preexisting NEXUS file...
add_charstatelabels_to_NEXUS <- function(charstatelabels, infn, outfn)
	{
	lines = scan(file=infn, what="character", sep="\n")
	
	# Go through lines
	line_with_charstatelabels = NULL
	
	for (i in 1:length(lines))
		{
		# Get line
		line = lines[i]
		
		# Check if line matches CHARSTATELABELS
		if (grepl("CHARSTATELABELS", line))
			{
			line_with_charstatelabels = i
			}
		else
			{
			blah=TRUE
			}
		}
	
	
	# Assuming you found line_with_charstatelabels...
	lines[line_with_charstatelabels+1] = charstatelabels
	
	write.table(lines, file=outfn, quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	return(lines)
	
	}




count_rows_NOT_containing_substr_in_col <- function(dtf, col, tmp_substr)
	{
	list_to_kill = grepl(tmp_substr, col, fixed=TRUE)
	list_to_keep = (list_to_kill == FALSE)
	print(sum(list_to_keep))
	return(sum(list_to_keep))
	}

count_rows_containing_substr_in_col <- function(dtf, col, tmp_substr)
	{
	list_to_keep = grepl(tmp_substr, col, fixed=TRUE)
	print(sum(list_to_keep))
	return(sum(list_to_keep))
	}



remove_rows_containing_substr_in_col <- function(dtf, col, tmp_substr)
	{
	list_to_kill = grepl(tmp_substr, col, fixed=TRUE)
	list_to_keep = (list_to_kill == FALSE)
	return(dtf[list_to_keep, ])
	}


retain_rows_containing_substr_in_col <- function(dtf, col, tmp_substr)
	{
	list_to_keep = grepl(tmp_substr, col, fixed=TRUE)
	return(dtf[list_to_keep, ])
	}


extract_rows_with_col_matching_findthis_str <- function(temp_table, col, findthis)
	{
	# Extract North-America-specific modern extinctions
	matches = find_str_inlist(findthis, tmplist=col)
	cat("# of rows matching '", findthis, "': ", sum(matches), "\n", sep="")
	
	temp_table = temp_table[matches, ]
	return(temp_table)
	}


# 2014-01-30_NJM:
# moving to BioGeoBEARS readwrite
# for process_DIVA_output

# Find a string inside another string...
find_instring <- function(findthis, string)
	{
	result = grep(findthis, string)
	if (length(result) > 0)
		{
		#print("match")
		match = TRUE
		}
	else
		{
		match = FALSE
		}
	return(match)
	}

# 2014-01-30_NJM:
# moving to BioGeoBEARS readwrite
# for process_DIVA_output
find_str_inlist <- function(findthis, tmplist)
	{
	# make empty list
	outlist = rep(NA, length(tmplist))
	for (i in 1:length(tmplist))
		{
		string = tmplist[i]
		outlist[i] = find_instring(findthis, string)
		
		}
	return(outlist)
	}


match_list1_in_list2_return_all_indices <- function(list1, list2)
	{
	# returns indices of list2 that match something in list1
	matchlist = list2 %in% list1
	
	matching_indices = seq(from=1, to=length(list2), by=1)[matchlist]
	
	return(matching_indices)
	}

# 
# 
# # return matching TRUE/FALSE values
# # list1 (.e.g. a big list) TRUE if it is found in list2 (e.g. a smaller list)
# match_list1_in_list2 <- function(list1, list2)
# 	{
# 	matchlist = list1 %in% list2
# 	return(matchlist)
# 	}


# return matching TRUE/FALSE values
# list1 (.e.g. a big list) TRUE if it is found in list2 (e.g. a smaller list)
match_item_in_list2 <- function(item, list2, return_indices=TRUE)
	{
	matchlist_TF = list2 %in% item
	
	if (return_indices == TRUE)
		{
		matching_indices = seq(from=1, to=length(list2), by=1)[matchlist_TF]
		return(matching_indices)
		}

	if (return_indices == FALSE)
		{
		return(matchlist_TF)
		}
	}

compare_2cols_against_another_2cols_find_matching_pairs_in_list2 <- function(list1, list2)
	{
	list1_paste = paste(list1[,1], ",", list1[,2], sep="")
	list2_pasteA = paste(list2[,1], ",", list2[,2], sep="")
	list2_pasteB = paste(list2[,1], ",", list2[,2], sep="")
	
	list2_matches1 = match_list1_in_list2(list2_pasteA, list1_paste)
	list2_matches2 = match_list1_in_list2(list2_pasteB, list1_paste)
	
	list2_matches = list2_matches1 + list2_matches2 > 0
	return(list2_matches)
	}



# NOTE!!! THESE MATCH FUNCTIONS JUST RETURN THE *FIRST* MATCH, *NOT* ALL MATCHES
# (argh)
# return indices in 2nd list matching the first list
# It WILL return one match for each item in the list, though...
# get_indices_where_list1_occurs_in_list2 <- function(list1, list2)
# 	{
# 	match_indices = match(list1, list2)
# 	return(match_indices)
# 	}


get_indices_where_string_in_list1_occurs_in_list2 <- function(list1, list2)
	{
	match_indices = unlist(lapply(list1, get_index_in_list2_matching_string, list2))
	match_indices[is.na(match_indices) == TRUE] = ""
	return(match_indices)
	}

# returns the index of the FIRST hit
get_index_in_list2_matching_string <- function(tmpstr, list2)
	{
	matches_TF = lapply(tmpstr, grepl, list2)[[1]]
	first_hit_index = match(TRUE, matches_TF)
	return(first_hit_index)
	}


# # return indices in 2nd list matching the first list
# get_indices_where_list1_occurs_in_list2_noNA <- function(list1, list2)
# 	{
# 	match_indices = match(list1, list2)
# 	match_indices = return_items_not_NA(match_indices)
# 	return(match_indices)
# 	}

# this is f-ed up!!
get_real_indices_where_list1_occurs_in_list2 <- function(list1, list2)
	{
	print("get_real_indices_where_list1_occurs_in_list2()")
	print("WARNING: THIS FUNCTION RETURNS THE INDICES OF NOT NA, NOT THE INDICES OF THE MATCHES")
	indices = 1:length(list1)
	
	matches_or_NAs = get_indices_where_list1_occurs_in_list2(list1, list2)
	not_NA = is.not.NA(matches_or_NAs)
	indices_to_return = indices[not_NA]
	return(indices_to_return)
	}


# place items with string matching list2 in new column of size list2
#place_items_matching_str_in_list1_in_newcol <- function()



# colmax, max.col

colmax <- function(dtf)
	{
	maxvals = apply(dtf, 2, max, na.rm=TRUE)
	return(maxvals)
	}

colmin <- function(dtf)
	{
	minvals = apply(dtf, 2, min, na.rm=TRUE)
	return(minvals)
	}

normalize_by_colmax <- function(dtf)
	{
	maxvals = colmax(dtf)
	# this outputs the results into columns, instead of corresponding rows 
	x = apply(dtf, 1, "/", maxvals)
	
	# so transpose
	normalized_df = t(x)
	normalized_df = as.data.frame(normalized_df)
	names(normalized_df) = names(dtf)
	
	return(normalized_df)
	}




##############################################
# Statistical tests
##############################################


# Welch's t-test (comparing means of two samples with different 
# sizes and possibly different sample variances)
setup = '
samp1mean = mean_rfd_genetree_to_genetree
samp1std = std_rfd_genetree_to_genetree
samp2mean = mean_rfd_genetree_to_sptree
samp2std = std_rfd_genetree_to_sptree
'
welchs_ttest_using_means <- function(samp1mean, samp2mean, samp1std, samp2std, n1, n2)
	{
	#print(mean_x1)
	#print(mean_x2)
	mean_x1 = samp1mean
	mean_x2 = samp2mean
	
	std_x1 = samp1std
	std_x2 = samp2std
	
	# Welch's t-test:
	# http://en.wikipedia.org/wiki/Welch%27s_t_test
	
	t_top = mean_x1 - mean_x2
	s_x1_x2 = ((std_x1^2/n1) + (std_x2^2/n2))
	t_bot = sqrt(s_x1_x2)
	
	t_welch = t_top / t_bot
	
	# Calculate degrees of freedom
	df_top = ((std_x1^2/n1) + (std_x2^2/n2))^2
	
	v1 = n1 - 1
	v2 = n2 - 1
	df_bot1 = std_x1^4 / (n1^2 * v1)
	df_bot2 = std_x2^4 / (n2^2 * v2)
	
	degfree = df_top / (df_bot1 + df_bot2)
	
	
	# one-tailed test, is observed t greater than null t?
	# (i.e., is null t less than observed t)
	#
	# (i.e., are between-treecloud distances (samp2) same as null (within-cloud distances)? )
	# 
	prob_obs_t_GT_null = pt(t_welch, degfree, lower.tail=TRUE)   # i.e., between dists >0
	
	#prob_between_dists_x1_is_greater_than_x2 = pt(t_welch, degfree)
	#cat("t_top: ", t_top, "; prob_null_between_dists_x1_not_greater_than_x2: ", prob_null_between_dists_x1_not_greater_than_x2, "; prob_between_dists_x1_is_greater_than_x2: ", prob_between_dists_x1_is_greater_than_x2, "\n", sep="")

	cat("samp1=", samp1mean, "+/-", samp1std, "	samp2=", samp2mean, "+/-", samp2std, "	t_welch=", t_welch, "	df = ", degfree, ", p(null < obs t)=", prob_obs_t_GT_null, " =P(samp2=>samp1); ", 1-prob_obs_t_GT_null, "= P(samp2<samp1) \n", sep="")

	
	tmprow = c(samp1mean, samp2mean, t_top, std_x1, std_x2, t_welch, degfree, prob_obs_t_GT_null, 1-prob_obs_t_GT_null)
	
	return(tmprow)
	}



min.f1f2 <- function(x, mu1, mu2, sd1, sd2)
	{
    f1 <- dnorm(x, mean=mu1, sd=sd1)
    f2 <- dnorm(x, mean=mu2, sd=sd2)
    pmin(f1, f2)
	}




##############################################
# Matrix functions
##############################################
# symmetric flip across the diagonal of a matrix
flipdiag <- function(mat)
	{
	newmat = mat
	for (i in ncol(mat))
		{
		for (j in nrow(mat))
			{
			newmat[j,i] = mat[i,j]
			}
		}
	return(newmat)
	}


nondiagTF <- function(m)
	{
	# Returns TF for non-diagonals
	diag(m) = "blah"
	
	TF = m != "blah"
	
	return(TF)
	}

nondiag <- function(m)
	{
	# Returns non-diagonals
	nondiag_m = m[nondiagTF(m)]
	
	return(nondiag_m)
	}



# Compare lower and upper diagonals
compare_lower_upper_diag <- function(x)
	{
	xtri = x[lower.tri(x)]
	
	ytmp = flipdiag(x)
	ytri = ytmp[lower.tri(ytmp)]
	
	(cor(cbind(xtri, ytri)))
	dif = c(mean(xtri))-c(mean(ytri))
	prstr = paste(mean(c(xtri)), mean(c(ytri)), dif, sep=" ")
	print(prstr)
	
	plot(xtri, ytri, type="p", pch=".")
	
	return(dif)
	}

# Compare row/col
compare_row_col <- function(x, rownum)
	{
	xnums = c(x[rownum, ])
	ynums = c(x[, rownum])
	
	(cor(cbind(xnums, ynums)))
	dif = c(mean(xnums))-c(mean(ynums))
	prstr = paste(mean(c(xnums)), mean(c(ynums)), dif, sep=" ")
	print(prstr)
	
	plot(xnums, ynums, type="p", pch=".")
	
	return(dif)
	}


# Find rows contains
# (note that this 
findrows_containing <- function(x, thing_to_find)
	{
	rows_to_keep = 1:nrow(x)
	for (i in 1:nrow(x))
		{
		if(sum(x[i, ] == thing_to_find, na.rm = TRUE) > 0)
			{
			rows_to_keep[i] = TRUE
			}
		else
			{
			rows_to_keep[i] = FALSE
			}
		}
	return(rows_to_keep)
	}


# Find rows containing no NAs
findrows_containing_noNA <- function(x)
	{
	rows_to_keep = is.na(rowSums(x))
	return(rows_to_keep == FALSE)
	}


# Find rows containing no NAs
findrows_containing_noINF <- function(x)
	{
	rows_to_keep = is.inf(rowSums(x))
	return(rows_to_keep == FALSE)
	}

# Find rows containing no NAs, no Infs
return_rows_containing_noNA_noInf <- function(x)
	{
	rows_to_keep = 1:nrow(x)
	rows_to_keep = findrows_containing_noNA(x)
	y = x[rows_to_keep, ]
	rows_to_keep = findrows_containing_noINF(y)	
	z = y[rows_to_keep, ]
	return(z)
	
	}


# Linear color scheme
defaults = '
nums = rnorm(1000, mean=500, sd=200)
numcolors = length(nums)
minval = 300
maxval = 700
mincolor = "yellow"
maxcolor = "red"
scaling=1

nums=avg_percent_diff
minval=0.35
maxval=0.52
'

colorscheme <- function(nums = rnorm(1000, mean=500, sd=200), numcolors = length(nums), minval = 300, maxval = 700, mincolor = "red", maxcolor = "yellow", scaling=1)
	{
	# for smoothColors:
	# smoothColorscalculates a sequence of colors. If two color names in the arguments
	# are separated by a number, that number of interpolated colors will be inserted between
	# the two color endpoints. 
	require(plotrix)
	
	# make a color table
	tmpdf1 = cbind(1:length(nums), nums)
	
	# sort it
	tmpdf2 = tmpdf1[order(nums), ]
	
	numvals_below_min = sum(nums <= minval)
	numvals_above_max = sum(nums >= maxval)
	minval_index = numvals_below_min
	maxval_index = (length(nums) - numvals_above_max) 
	numvals_that_change = maxval_index - minval_index - 2

	vals_below_min = 1:numvals_below_min
	vals_above_max = (length(nums) - numvals_above_max) : length(nums)
	
	
	# put in the changing colors
	changing_colors = smoothColors(mincolor, numvals_that_change, maxcolor)
	
	
	# rank scaling
	if (scaling == "rank")
		{
		changing_colors_final = changing_colors
		}
	else
		{
		# scaling is an exponent
		ttmpdf3 = tmpdf2[(minval_index+1) : (maxval_index), ]
		ttmpdf4 = ttmpdf3
		
		numer = ttmpdf3[, 2]-min(ttmpdf3[, 2])
		denom = max(ttmpdf3[, 2]) - min(ttmpdf3[, 2])
		color_indices = 1+round(numvals_that_change * (numer/denom)^scaling  )
		
		changing_colors_final = changing_colors[color_indices]
		
		#plot(color_indices, rep(1, length(color_indices)), col=changing_colors_final)
		#title("Colors that are changing")
		#xy = dotPlot(nums)
		#xy
		}	

	color_list = c( rep(mincolor, numvals_below_min), changing_colors_final, rep(maxcolor, numvals_above_max) )
	print(length(color_list))
	
	# Get the xy coordinates of the dotPlot
	xylist = get_dotPlot_coords(nums)

	# Plot color scheme
	# set new params but save old params
	old_params = par(mfrow=c(3, 1))
	
	# Points ranked by value
	plot(1:length(nums), rep(1, length(nums)), col=color_list)
	title("Colors in order")
	
	# Stacked-dot plot (histogram-like)
	plot(xylist[,1], xylist[,2], col="black", pch=1)
	title("Stacked-dot plot of input numbers")
	
	plot(xylist[,1], xylist[,2], col=color_list, pch=19)
	title("Stacked-dot plot of input numbers, colored by approximate assigned color")
	
	tmpdf3 = cbind(tmpdf2, color_list)
	tmpdf4 = tmpdf3[order(as.numeric(tmpdf3[,1])), ]
	tmpdf4 = adf(tmpdf4)
	names(tmpdf4) = c("orig_order", "nums", "color_list")
	
	par(old_params)
	
	return(tmpdf4)
	}




get_dotPlot_coords <- function (x, group, xlim, ylim, col, xlab, ylab, pch, cex, breaks, stacked = TRUE, ...) 
	{
	xlist = NULL
	ylist = NULL
	
    DB = FALSE
    pch.size = "O"
    grouped = TRUE
    parList = list(...)
    if (missing(xlab)) 
        xlab = deparse(substitute(x))
    x = x[!is.na(x)]
    if (missing(xlim)) 
        xlim = range(x)
    if (missing(ylim)) 
        ylim = c(0, 1)
    if (missing(ylab)) 
        ylab = ""
    if (missing(cex)) 
        cex = 1
    if (missing(group)) 
        group = rep(1, length(x))
    if (length(unique(group)) == 1) 
        grouped = FALSE
    if (missing(pch) || length(unique(group)) > length(pch)) 
        pch = 1:length(unique(group))
    if (missing(col) || length(unique(group)) > length(col)) 
        col = 1:length(unique(group))
    if (missing(breaks))
    	{
        plot(1, 1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, cex = cex, xlab = xlab, ylab = ylab)
        slotSizeX = strwidth(pch.size, units = "user", cex = cex)
        if (DB) 
            print(paste("slotSizeX:", slotSizeX))
        span = diff(range(x))
        temp1 = ppoints(2 * ceiling(span/slotSizeX))
        temp2 = numeric(length(temp1) + 2)
        temp2[2:(length(temp1) + 1)] = temp1
        temp2[1] = temp2[1] - 1.01 * diff(c(temp1[1], temp1[2]))
        temp2[length(temp2)] = rev(temp1)[1] + 1.01 * diff(c(temp1[1], 
            temp1[2]))
        temp2 = temp2 * span + min(x)
        temp = min(x) + ppoints(span/slotSizeX) * span
        breaks = numeric(length(temp) + 2)
        breaks[2:(length(temp) + 1)] = temp
        breaks[1] = temp[1] - diff(c(temp[1], temp[2])) * 1.001
        breaks[length(breaks)] = rev(temp)[1] + diff(c(temp[1], 
            temp[2])) * 1.001
        breaks = temp2
	    }
    slotSizeY = strheight(pch.size, units = "user", cex = cex)
    if (DB) 
        print(paste("slotSizeY:", slotSizeY))
    span = diff(ylim)
    temp1 = ppoints(2 * ceiling(span/slotSizeY))
    temp2 = numeric(length(temp1) + 2)
    temp2[2:(length(temp1) + 1)] = temp1
    temp2[1] = temp2[1] - 1.01 * diff(c(temp1[1], temp1[2]))
    temp2[length(temp2)] = rev(temp1)[1] + 1.01 * diff(c(temp1[1], 
        temp1[2]))
    yVec = temp2 * span + min(ylim)
    if (yVec[1] < 0) 
        yVec = yVec + abs(yVec[1])
    else yVec = yVec - yVec[1]
    if (DB) 
        print(paste("temp2:", temp2))
    if (DB) 
        print(paste("breaks:", breaks))
    histObj = hist(x, breaks = breaks, right = FALSE, plot = FALSE)
    hMids = histObj$mids
    hCounts = histObj$counts
    hMids = histObj$mids
    mat = matrix(NA, nrow = length(x), ncol = length(hMids))
    colMat = mat
    groupmat = mat
    numVec = 1:nrow(mat)
    cutOff = 1
    groupList = vector(mode = "list", length = length(unique(group)))
    for (k in unique(group))
    	{
        histObj = hist(x[group == k], breaks = breaks, plot = FALSE)
        hMids = histObj$mids
        hCounts = histObj$counts
        hMids = histObj$mids
        for (i in seq(along = hMids))
        	{
            value = pch[k]
            colValue = col[k]
            from = 0
            from = numVec[is.na(mat[, i])][1]
            to = from
            if (hCounts[i] == 0) 
                value = NA
            if (hCounts[i] >= 1) 
                to = to + hCounts[i] - 1
            if (to > cutOff) 
                cutOff = to
            if (DB)
            	{
                print(paste("from:", from))
                print(paste("to:", to))
                print(paste("i:", i))
                print(paste("value:", value))
    	        }
            mat[from:to, i] = value
            colMat[from:to, i] = colValue
	        }
        groupList[[k]] = groupmat
		}
    if (grouped && !stacked)
    	{
        groupIndex = unique(group)
        par(mfrow = c(length(groupIndex), 1))
        #for (i in groupIndex) dotPlot(x[group == i], xlim = xlim, breaks = breaks, cex = cex, xlab = xlab, ylab = ylab, col = col, pch = pch, ...)
	    }
    else
    	{
        mat = mat[1:cutOff, ]
        if (!is.matrix(mat)) 
            mat = matrix(mat, nrow = 1)
        if (DB) 
            print(mat)
        #plot(1, 1, xlim = xlim, ylim = ylim, type = "n", cex = cex, xlab = xlab, ylab = ylab, ...)
        for (i in 1:nrow(mat))
        	{
            x = hMids[!is.na(mat[i, ])]
            y = rep(i * 0.3, times = length(x))
            y = rep(yVec[i], times = length(x))
            col = colMat[i, !is.na(mat[i, ])]
            pch = mat[i, !is.na(mat[i, ])]
            #points(x, y, col = col, pch = pch, cex = cex)
            xlist = c(xlist, x)
            ylist = c(ylist, y)
            
        	}
	    }
	xylist = cbind(xlist, ylist)
	xylist = xylist[order(xlist), ]
	
	return(xylist)
    #if (DB) 
    #    print(hMids)
    #invisible(mat)
}




##############################################
# Plotting functions
##############################################
# Modified from mrds package, http://www.inside-r.org/packages/cran/mrds/docs/histline , e.g.
# Need "mrds:::" as the author left @export out the .R file...
# mrds:::histline(height=histvals$density, breaks=histvals$breaks, lineonly=FALSE, outline=TRUE, fill=FALSE)
# mrds:::histline(height=histvals$density, breaks=histvals$breaks, lineonly=TRUE, outline=TRUE, fill=FALSE)

# TURN THIS BACK ON AFTER UPGRADING TO OPTIMX 2013

junk='
mrds:::histline
function (height, breaks, lineonly = FALSE, outline = FALSE, 
    fill = FALSE, ylim = range(height), xlab = "x", ylab = "y", 
    det.plot = FALSE, ...) 
{
    n = length(height)
    if (length(breaks) != (n + 1)) 
        stop("breaks must be 1 longer than height")
    if (outline) {
        y = c(0, rep(height, times = rep(2, n)), 0)
        x = rep(breaks, times = rep(2, (n + 1)))
    }
    else {
        y = rep(0, 4 * n)
        x = rep(0, 4 * n + 2)
        for (i in 1:n) {
            y[((i - 1) * 4 + 1):(i * 4)] = c(0, rep(height[i], 
                2), 0)
            x[((i - 1) * 4 + 1):(i * 4)] = c(rep(breaks[i], 2), 
                rep(breaks[i + 1], 2))
        }
        x = x[1:(4 * n)]
    }
    if (lineonly) {
        if (!fill) 
            lines(x, y, ...)
        else polygon(x, y, ...)
    }
    else {
        if (!fill) 
            plot(x, y, type = "l", ylim = ylim, xlab = xlab, 
                ylab = ylab, ...)
        else {
            if (det.plot) {
                plot(x, y, type = "n", ylim = ylim, xlab = xlab, 
                  ylab = ylab, yaxp = c(0, 1, 5))
            }
            else {
                plot(x, y, type = "n", ylim = ylim, xlab = xlab, 
                  ylab = ylab)
            }
            polygon(x, y, ...)
        }
    }
}
'

# Square the x and y dimensions
squareplot <- function(x, y)
	{
	# Scatter plot with axes drawn on the same scale
	# I'd like to produce some scatter plots where N units on the X axis are > equal to N units
	# on the Y axis (as measured with a ruler, on screen or paper).
	# x <- sample(10:200,40)
	# y <- sample(20:100,40)
	# windows(width = max(x),height = max(y))
	# plot(x,y)
	
	# try:
	plot(x, y, asp = 1)

	# or, better:
	# library(MASS)
	# eqscplot(x,y)

	# or
	# library(lattice)
	# xyplot(y ~ x, aspect = "iso")
	}



# Make some subplots
subplots <- function(nrows=1, ncols=2)
	{
	q = quartz(width=7*ncols, height=7*nrows)
	par(mfcol = c(nrows, ncols))
	return(q)
	}


# Pairs plot with histograms
panel.hist <- function(x, ...)
	{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
	}

pairs_with_hist <- function(datamat, ...)
	{
	pairs(datamat, diag.panel=panel.hist, ...)
	return()
	}



setup='
xvals = UK_input[,2]
yvals = UK_input[,1]
nbins=10
'
meansplot = function(xvals, yvals, nbins=10, tmppch=1)
	{
	
	bins_to_plot = seq(min(xvals), max(xvals), (max(xvals)-min(xvals))/nbins )
	mean_xvals = NULL
	mean_yvals = NULL
	for (i in 1:(length(bins_to_plot)-1) )
		{
		binstart = bins_to_plot[i]
		binend = bins_to_plot[i+1]
		
		TF1 = xvals >= binstart
		TF2 = xvals < binend
		TF = TF1 + TF2 == 2
		mean_xvals = c(mean_xvals, mean(xvals[TF]))
		mean_yvals = c(mean_yvals, mean(yvals[TF]))
		}
	lm1 = linear_regression_plot(mean_xvals, mean_yvals, tmppch=tmppch)
	return(lm1)
	}



setup = '
x = ipd$rating
y = ipd$num_plays
xlabel="len_sec"
ylabel="num_plays"
'

# If slope1=TRUE, subtract 1:1 slope from the line and test for differences

linear_regression_plot <- function(x, y, xlabel="x", ylabel="y", tmppch=".", pointscol="black", tmplinecol="black", tmplty=1, tmplwd=1, plottext=TRUE, legend_title="", textcol="black", legend_x="topleft", legend_y=NULL, xlim=minmax_pretty(x), ylim=minmax_pretty(y), increment_fraction=NULL, legend_cex=1, axis_cex=1, slope1=FALSE, intercept_in_legend=TRUE, add_to_plot=FALSE, printall=TRUE, ...)
	{
	# Make a linear regression plot
	model1 = lm(y~x)
	slope = model1$coefficients[2]
	intercept = model1$coefficients[1]
	
	if (printall)
		{
		print("Summary of model 1 (standard regression, y predicted by x)")
		print(summary(model1))
		}
	



	# Set up Legend
	if (is.null(increment_fraction))
		{
		increment_fraction = 0.05
		}

	# Legend positions
	if (legend_x == "topleft")
		{
		legend_x = min(xlim)
		legend_y = 1 * max(ylim)
		}
	
	increment = increment_fraction * (max(ylim) - min(ylim))
	
	
	# Plot the full plot, or just add points
	if (add_to_plot == FALSE)
		{
		plot(x, y, xlab=xlabel, ylab=ylabel, pch=tmppch, xlim=xlim, ylim=ylim, col=pointscol, cex.axis=axis_cex, ...)
		} else {
		points(x, y, pch=tmppch, xlim=xlim, ylim=ylim, col=pointscol, ...)
		}
	tmpx1 = min(x, na.rm=TRUE)
	tmpx2 = max(x, na.rm=TRUE)
	tmpy1 = slope*tmpx1 + intercept
	tmpy2 = slope*tmpx2 + intercept
	
	segments(x0=tmpx1, y0=tmpy1, x1=tmpx2, y1=tmpy2, col=tmplinecol, lty=tmplty, lwd=tmplwd)
	
	
	if (plottext)
		{
		# Subtract 1:1 slope if desired
		# (for Patrick Shih cyanobacteria analysis)
		if (slope1 == TRUE)
			{
			# Subtract 1:1 from the y values
			ynew = y-x
			model2 = lm(ynew~x)
			pval = summary(model2)$coefficients[2,4]

			print("Summary of model 2 (standard regression, (y-x) predicted by x; this means we've substracted out the 1:1 expected relationship.)")
			print(summary(model2))
			} else {
			pval = summary(model1)$coefficients[2,4]
			} # END if (slope1 == TRUE)
		
		# Start ypos for legend text at 0
		ypos = 0
		# Plot a legend title if desired
		if (legend_title != "")
			{
			text(x=legend_x, y=legend_y - ypos*increment, pos=4, labels=legend_title, cex=legend_cex, col=textcol)
			} # END if (legend_title != "")
	
		
		if (intercept_in_legend==TRUE)
			{		
			R2 = summary(model1)$r.squared
			slope = summary(model1)$coefficients[2,1]
			slopeSE = summary(model1)$coefficients[2,2]
			slope95 = 1.96*slopeSE

			intercept = summary(model1)$coefficients[1,1]
			interceptSE = summary(model1)$coefficients[1,2]
			intercept95 = 1.96*interceptSE
		
			R2txt = bquote(italic(R)^2 == .(format(R2, digits=3)))
			text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=R2txt, cex=legend_cex, col=textcol)
			slopetxt = bquote(italic(m) == .(format(slope, digits=3)) %+-% .(format(slope95, digits=3)))
			text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=slopetxt, cex=legend_cex, col=textcol)
			#intercepttxt = paste("intercept=", format(intercept, digits=3), " +/- ", format(intercept95, digits=3), sep="")
			intercepttxt = bquote(italic(b) == .(format(intercept, digits=3)) %+-% .(format(intercept95, digits=3)))
			text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=intercepttxt, cex=legend_cex, col=textcol)
		
			if (slope1 == TRUE)
				{
				#pvaltxt = paste("p = ", format(pval, digits=3), " (null: slope is 1:1)", sep="")
				pvaltxt = bquote(italic(p) == .(format(pval, digits=3))~" (null: slope is 1:1)")
				text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=pvaltxt, cex=legend_cex, col=textcol)
				} else {
				#pvaltxt = paste("p = ", format(pval, digits=3), sep="")
				pvaltxt = bquote(italic(p) == .(format(pval, digits=3)))
				text(x=legend_x, y=legend_y - (ypos=ypos+1)*increment, pos=4, labels=pvaltxt, cex=legend_cex, col=textcol)
				} # END if (slope1 == TRUE)

			} # END if (intercept_in_legend==TRUE)
		
# 		if (intercept_in_legend==TRUE)
# 			{
# 			txt_to_plot = paste(R2txt, slopetxt, intercepttxt, pvaltxt, sep="\n")
# 			} else {
# 			txt_to_plot = paste(R2txt, slopetxt, pvaltxt, sep="\n")
# 			}
		
		# http://stackoverflow.com/questions/3761410/how-can-i-plot-my-r-squared-value-on-my-scatterplot-using-r
		# bty suppresses box
		# print(legend_x)
		# print(legend_y)
		# Discounted
		#legend(x=legend_x, y=legend_y, bty="n", legend=txt_to_plot, cex=0.9)

		}
	
	if (printall == TRUE)	
		{
		print(model1)
		print(summary(model1))
		}
	
	
	
	return(model1)
	}




extract_lm_stats <- function(tmplm)
	{
	R2 = summary(tmplm)$r.squared
	
	# Get the slp information
	slp = summary(tmplm)$coefficients[2,1]
	slpSE = summary(tmplm)$coefficients[2,2]
	slp_2sigma = 1.96*slpSE
	slp0025 = slp + 1.96*slpSE
	slp0975 = slp - 1.96*slpSE
	slp_pval = summary(tmplm)$coefficients[2,4]
	slp_tval = summary(tmplm)$coefficients[2,3]
	slp_sig = pvals_to_sigstars(slp_pval)
	
	
	
	# Get the intercent (int) information
	int = summary(tmplm)$coefficients[1,1]
	intSE = summary(tmplm)$coefficients[1,2]
	int_2sigma = 1.96*intSE
	int0025 = int + 1.96*intSE
	int0975 = int - 1.96*intSE
	int_pval = summary(tmplm)$coefficients[1,4]
	int_tval = summary(tmplm)$coefficients[1,3]
	int_sig = pvals_to_sigstars(int_pval)
	
	tmprow = c(R2, slp, slpSE, slp_2sigma, slp0025, slp0975, slp_pval, slp_tval, slp_sig, int, intSE, int_2sigma, int0025, int0975, int_pval, int_tval, int_sig)
	
	tmpcolnames = c("R2", "slp", "slpSE", "slp_2sigma", "slp0025", "slp0975", "slp_pval", "slp_tval", "slp_sig", "int", "intSE", "int_2sigma", "int0025", "int0975", "int_pval", "int_tval", "int_sig")
	
	return(tmprow)
	}




extract_glm_coefs <- function(glmobj)
	{
	#R2 = summary(glmobj)$r.squared
	
	# Get the slp information
	coefs = summary(glmobj)$coefficients
	names_coefs = rownames(coefs)
	names_params = colnames(coefs)
	coefs_list = rep(NULL, times=nrow(coefs)*ncol(coefs))
	names_list = rep(NULL, times=nrow(coefs)*ncol(coefs))
	

	conf_intervals = confint(glmobj)
	lower025 = rep(NULL, times=nrow(coefs))
	upper975 = rep(NULL, times=nrow(coefs))
	names_lower025 = paste("lower025_", names_coefs, sep="")
	names_upper975 = paste("upper975_", names_coefs, sep="")

	
	onum = 0
	for (n in 1:length(names_coefs))
		{
		
		for (m in 1:length(names_params))
			{
			onum = onum+1
			coefs_list[onum] = coefs[n,m]
			
			# Fix names
			tmpname = paste(names_coefs[n], "_", names_params[m], sep="")
			tmpname0 = gsub("(Intercept)", "Intercept", tmpname, fixed=TRUE)
			tmpname1 = gsub("Std. Error", "stdErr", tmpname0)
			tmpname2 = gsub("z value", "zval", tmpname1)
			tmpname3 = gsub(" ", "_", tmpname2)
			tmpname4 = gsub("\\.", "_", tmpname3)
			tmpname5 = gsub("Pr(>|z|)", "pvalGTz", tmpname4, fixed=TRUE)
			names_list[onum] = tmpname5
			}

		# Don't calculate using the normal distribution, instead use confint()
		# lower025[n] = coefs[n,1] - (1.96 * coefs[n,2])
		# upper975[n] = coefs[n,1] + (1.96 * coefs[n,2])
		lower025[n] = conf_intervals[n,1]
		upper975[n] = conf_intervals[n,2]

		}
	coefs_list
	names_list
	
	
	tmprow = c(coefs_list, lower025, upper975)
	tmpcolnames = c(names_list, names_lower025, names_upper975)
	
	names(tmprow) = tmpcolnames
	
	glmcoefs = tmprow
	return(glmcoefs)
	}


extract_glm_stats <- function(glmobj)
	{
	a = summary(glmobj)$deviance
	b = summary(glmobj)$aic
	c = summary(glmobj)$df.residual
	d = summary(glmobj)$null.deviance
	e = summary(glmobj)$df.null
	
	tmprow = c(a, b, c, d, e)
	tmpnames = c("deviance", "aic", "df.residual", "null.deviance", "df.null")
	names(tmprow) = tmpnames
	
	tmpstats = tmprow
	return(tmpstats)
	}	


# Extract the approximate 95% confidence intervals
# as brunch() doesn't have a confint method
confint_brunch <- function(brunchMod)
	{
	coefs = summary(brunchMod)$coefficients
	nrows = nrow(coefs)
	brunch_CI = NULL

	#check_for_NA
	if (is.na(coefs[1]))
		{
		return(NA)
		}

	
	for (j in 1:nrows)
		{
		lower025 = coefs[j,"Estimate"] - 1.96 * coefs[1,"Std. Error"]
		upper975 = coefs[j,"Estimate"] + 1.96 * coefs[1,"Std. Error"]
		CIs = c(lower025, upper975)
		
		CIs2 = matrix(CIs, nrow=1, ncol=2)
		brunch_CI = rbind(brunch_CI, CIs2)
		}
	
	rownames(brunch_CI) = rownames(coefs)
	colnames(brunch_CI) = c("2.5 %", "97.5 %")
	brunch_CI
	
	return(brunch_CI)
	}



# Extract the approximate 95% confidence intervals
# as compar.gee() doesn't have a confint method
confint_gee <- function(geeobj)
	{
	coefs = extract_gee_coefmat(geeobj)
	nrows = nrow(coefs)
	brunch_CI = NULL

	#check_for_NA
	if (is.na(coefs[1]))
		{
		return(NA)
		}

	
	for (j in 1:nrows)
		{
		lower025 = coefs[j,"Estimate"] - 1.96 * coefs[1,"S.E."]
		upper975 = coefs[j,"Estimate"] + 1.96 * coefs[1,"S.E."]
		CIs = c(lower025, upper975)
		
		CIs2 = matrix(CIs, nrow=1, ncol=2)
		brunch_CI = rbind(brunch_CI, CIs2)
		}
	
	rownames(brunch_CI) = rownames(coefs)
	colnames(brunch_CI) = c("2.5 %", "97.5 %")
	brunch_CI
	
	return(brunch_CI)
	}




extract_brunch_coefs <- function(brunchMod)
	{
	require(caper)	# For brunch
	#######################################################
	# caper's brunch algorithm
	# caper.pdf
	# "The 'brunch' algorithm calculates contrasts for models that include binary categorical variables.
	# Contrasts are identified and calculated for all variables in the model for a set of nodes where each
	# side can be unequivocally attributed to one or other of the categories. Unlike 'crunch', nested contrasts
	# are not calculated and each row of data at the tips is used only once. This follows Burt (1989):
	# contrasts whose paths do not meet or cross at any point will be phylogenetically independent."
	#######################################################

	#R2 = summary(brunchMod)$r.squared
	
	# Get the slp information
	coefs = summary(brunchMod)$coefficients
	
	#check_for_NA
	if (is.na(coefs[1]))
		{
		return(NA)
		}
	
	names_coefs = rownames(coefs)
	names_params = colnames(coefs)
	coefs_list = rep(NULL, times=nrow(coefs)*ncol(coefs))
	names_list = rep(NULL, times=nrow(coefs)*ncol(coefs))
	

	conf_intervals = confint_brunch(brunchMod)
	lower025 = rep(NULL, times=nrow(coefs))
	upper975 = rep(NULL, times=nrow(coefs))
	names_lower025 = paste("lower025_", names_coefs, sep="")
	names_upper975 = paste("upper975_", names_coefs, sep="")

	
	onum = 0
	for (n in 1:length(names_coefs))
		{
		
		for (m in 1:length(names_params))
			{
			onum = onum+1
			coefs_list[onum] = coefs[n,m]
			
			# Fix names
			tmpname = paste(names_coefs[n], "_", names_params[m], sep="")
			tmpname0 = gsub("(Intercept)", "Intercept", tmpname, fixed=TRUE)
			tmpname1 = gsub("Std. Error", "stdErr", tmpname0)
			tmpname2 = gsub("z value", "zval", tmpname1)
			tmpname3 = gsub(" ", "_", tmpname2)
			tmpname4 = gsub("\\.", "_", tmpname3)
			tmpname5 = gsub("Pr(>|z|)", "pvalGTz", tmpname4, fixed=TRUE)
			tmpname6 = gsub("Pr(>|t|)", "pvalGTt", tmpname5, fixed=TRUE)
			names_list[onum] = tmpname6
			}

		# Don't calculate using the normal distribution, instead use confint()
		# lower025[n] = coefs[n,1] - (1.96 * coefs[n,2])
		# upper975[n] = coefs[n,1] + (1.96 * coefs[n,2])
		lower025[n] = conf_intervals[n,1]
		upper975[n] = conf_intervals[n,2]

		}
	coefs_list
	names_list
	
	
	tmprow = c(coefs_list, lower025, upper975)
	tmpcolnames = c(names_list, names_lower025, names_upper975)
	
	names(tmprow) = tmpcolnames
	
	glmcoefs = tmprow
	return(glmcoefs)
	}






extract_brunch_stats <- function(brunchMod)
	{
	#check_for_NA
	coefs = summary(brunchMod)$coefficients
	if (is.na(coefs[1]))
		{
		return(NA)
		}
	
	a = summary(brunchMod)$sigma
	b = summary(brunchMod)$df
	c = summary(brunchMod)$r.squared
	d = summary(brunchMod)$adj.r.squared
	e = summary(brunchMod)$fstatistic
	f = summary(brunchMod)$cov.unscaled
	
	tmprow = c(a, b, c, d, e, f)
	tmpnames = c("sigma", "df1", "df2", "df3", "r.squared", "adj.r.squared", "Fstat_value", "Fstat_numdf", "Fstat_dendf", "cov.unscaled")
	names(tmprow) = tmpnames
	
	tmpstats = tmprow
	return(tmpstats)
	}	



extract_gee_stats <- function(geeobj)
	{
	#check_for_NA
	coefs = extract_gee_coefmat(geeobj)
	if (is.na(coefs[1]))
		{
		return(NA)
		}

	# 
	# nobs - the number of observations.
	# QIC - the quasilikelihood information criterion as defined by Pan (2001).
	# coefficients - the estimated coefficients (or regression parameters).
	# scale	- the scale (or dispersion parameter).
	# W	- the variance-covariance matrix of the estimated coefficients.
	# dfP - the phylogenetic degrees of freedom (see Paradis and Claude for details on this).
	# 

	
	a = geeobj$nobs
	b = geeobj$dfP
	c = geeobj$scale
	d = geeobj$QIC
	
	tmprow = c(a, b, c, d)
	tmpnames = c("nobs", "dfP", "scale", "QIC")
	names(tmprow) = tmpnames
	
	tmpstats = tmprow
	return(tmpstats)
	}	





extract_gee_coefmat <- function(geeobj)
	{
	#R2 = summary(geeobj)$r.squared
	
	# modified from print.compar.gee
    NAs <- is.na(geeobj$coef)
    coefmat <- geeobj$coef[!NAs]
    cnames <- names(coefmat)
    coefmat <- matrix(rep(coefmat, 4), ncol = 4)
    dimnames(coefmat) <- list(cnames, c("Estimate", "S.E.", "t", 
        "Pr(T > |t|)"))
    df <- geeobj$dfP - dim(coefmat)[1]
    coefmat[, 2] <- sqrt(diag(geeobj$W))
    coefmat[, 3] <- coefmat[, 1]/coefmat[, 2]
    
    # pt = prob. t-distribution
    if (df < 0)
    	{
        warning("not enough degrees of freedom to compute P-values.")
        coefmat[, 4] <- NA
	    } else {
		coefmat[, 4] <- 2 * (1 - pt(abs(coefmat[, 3]), df))
		}
	coefmat
	
	return(coefmat)
	}

extract_gee_coefs <- function(geeobj)
	{
	coefmat = extract_gee_coefmat(geeobj)
	coefmat

	coefs = coefmat
	names_coefs = rownames(coefs)
	names_params = colnames(coefs)
	coefs_list = rep(NULL, times=nrow(coefs)*ncol(coefs))
	names_list = rep(NULL, times=nrow(coefs)*ncol(coefs))
	

	conf_intervals = confint_gee(geeobj)
	lower025 = rep(NULL, times=nrow(coefs))
	upper975 = rep(NULL, times=nrow(coefs))
	names_lower025 = paste("lower025_", names_coefs, sep="")
	names_upper975 = paste("upper975_", names_coefs, sep="")

	
	onum = 0
	for (n in 1:length(names_coefs))
		{
		
		for (m in 1:length(names_params))
			{
			onum = onum+1
			coefs_list[onum] = coefs[n,m]
			
			# Fix names
			tmpname = paste(names_coefs[n], "_", names_params[m], sep="")
			tmpname0 = gsub("(Intercept)", "Intercept", tmpname, fixed=TRUE)
			tmpname1 = gsub("S.E.", "stdErr", tmpname0)
			tmpname2 = gsub("Intercept t", "Int_tval", tmpname1)
			tmpname3 = gsub("Pr(T > |t|)", "pvalGTt", tmpname2, fixed=TRUE)
			tmpname4 = gsub(" ", "_", tmpname3)
			tmpname5 = gsub("\\.", "_", tmpname4)
			names_list[onum] = tmpname5
			}

		# Don't calculate using the normal distribution, instead use confint()
		# lower025[n] = coefs[n,1] - (1.96 * coefs[n,2])
		# upper975[n] = coefs[n,1] + (1.96 * coefs[n,2])
		lower025[n] = conf_intervals[n,1]
		upper975[n] = conf_intervals[n,2]

		}
	coefs_list
	names_list
	
	
	tmprow = c(coefs_list, lower025, upper975)
	tmpcolnames = c(names_list, names_lower025, names_upper975)
	
	names(tmprow) = tmpcolnames
	
	geecoefs = tmprow
	return(geecoefs)
	}




lmstats_to_string_for_plot <- function(tmprow)
	{
	R2txt = paste("R2 = ", format(as.numeric(tmprow[1]), digits=3), sep="")
	slopetxt = paste("m=", format(as.numeric(tmprow[2]), digits=3), " +/- ", format(as.numeric(tmprow[4]), digits=3), sep="")
	pvaltxt = paste("p = ", format(as.numeric(tmprow[7]), digits=3), sep="")
	
	txt_to_plot = paste(R2txt, slopetxt, pvaltxt, sep="\n")
	
	return(txt_to_plot)
	}


pvals_to_sigstars <- function(pval, pval1star=0.05, pval2star=0.01, pval3star=0.001)
	{
	sigstar = " "
	
	if (pval < 0.05)
		{
		sigstar = "*"
		}
	if (pval < 0.01)
		{
		sigstar = "**"
		}
	if (pval < 0.001)
		{
		sigstar = "***"
		}
	
	return(sigstar)
	}


# Linear colorscale
# Take some data, assign a color to each point, linearly on the colorscale

setup='
# 
datapoints = exp(tmp_points_to_plot$response) / (1 + exp(tmp_points_to_plot$response))
colorscale = colorlist
datapoints=tmp_datapoints
colorscale=tmp_colorscale
rescale=FALSE
minscale=""
maxscale=""
minscale=0
maxscale=1
minscale=cellStats(marsgrid4, "min")
maxscale=cellStats(marsgrid4, "max")
'
linear_colorscale <- function(datapoints, colorscale, rescale=FALSE, minscale="", maxscale="", minout="", maxout="")
	{
	length_colorscale = length(colorscale)
	minval = min(datapoints)
	maxval = max(datapoints)
		
	# for linear scaling from 1 to numpoints
	if (rescale == TRUE)
		{
		if (minscale == "")
			{
			colorscale_indices = ceiling( ( ((datapoints - minval) / (maxval-minval)) * (length_colorscale-1) ))
			} else {
			#minscale = min(datapoints)
			#maxscale = max(datapoints)
			# no rescaling, just raw values on the 0-1 scale
			colorscale_indices = ceiling(  ( ((datapoints - minval) / (maxval-minval)) * (maxscale - minscale) + minscale ) * (length_colorscale-1) )
			}
		} else {
		if (minscale == "")
			{
			# no rescaling, just raw values on the 0-1 scale
			colorscale_indices = ceiling( datapoints * (length_colorscale-1) )
			
			# Correction, if the points are negative, make that the new zero
			if (minval < 0)
				{
				datapoints_new = datapoints - minval
				datapoints_new = datapoints_new * maxval
				colorscale_indices = ceiling( datapoints_new * (length_colorscale-1) )
				}
			} else {
			#minscale = min(datapoints)
			#maxscale = max(datapoints)
			# no rescaling, just raw values on the minscale-maxscale scale
			colorscale_indices = ceiling( datapoints * (length_colorscale-1) )
			colorscale_indices[datapoints >= maxscale] = length_colorscale
			colorscale_indices[datapoints <= minscale] = 1
			
			}		
		}
	
	colors_to_plot_per_datapoint = colorscale[colorscale_indices]
	return(colors_to_plot_per_datapoint)
	}

setup='
datapoints=tmp_datapoints
colorscale=tmp_colorscale
minin=tmpmin
maxin=tmpmax
minin=""
maxin=""
minout=""
maxout=""
'
linear_colorscale2 <- function(datapoints, colorscale, rescale=FALSE, minin="", maxin="", minout="", maxout="")
	{
	length_colorscale = length(colorscale)
	
	if (minin == "")
		{
		minin = min(datapoints)
		}
	if (maxin == "")
		{
		maxin = max(datapoints)
		}
	
	datapoints[datapoints < minin] = minin
	datapoints[datapoints > maxin] = maxin
	
	# minout / maxout -- the proportions of the colorscale that you want to occupy...
	if (minout == "")
		{
		minout = 0
		}
	if (maxout == "")
		{
		maxout = 1
		}
		
	colorscale_indices = ceiling(   ( ((datapoints - minin) / (maxin-minin)) * (maxout-minout) + minout) * (length_colorscale-1)        )

	colors_to_plot_per_datapoint = colorscale[colorscale_indices]
	return(colors_to_plot_per_datapoint)
	}


make_grid_from_coords <- function(xcoords, ycoords)
	{
	xcoords1 = rep(xcoords, length(ycoords))
	ycoords1 = rep(ycoords, length(xcoords))
	xy = cbind(xcoords1, ycoords1)
	#plot(xy)
		
	# Convert to SpatialPoints
	S = SpatialPoints(xy)
	class(S)
	proj4string(S) = envproj
	
	# Convert to SpatialPixels
	G = S
	gridded(G) = TRUE
	#plot(G)
	class(G)
	
	return(G)
	}















# =========================================================
# Google-related functions
# =========================================================
# Copy the corrected google numbers to each matching quote 
# (one could also split the results amongst the matching quotes; see next function)
copy_fixed_quantchars_to_each_matching <- function(quantchars, quotes_list_no_extra_spaces_to_google, fixed_goog, cols_in_fixed_goog, unique_quotes_to_google, divide_by_frequency = TRUE)
	{
	quantchar_copy_nums_to_repeated_quotes = NULL
	q = NULL
	
	# the full list of quotes (from whatever source; full quantchars dataset, or subset)
	rownames = quantchars$rownames
	quotetxt = quotes_list_no_extra_spaces_to_google
	
	numcols = ncol(quantchars) - 1
	tmprow = rep(NA, numcols)
	
	# store the accumulating rows here
	tmpmat = NULL

	names_cols_in_fixed_goog = names(fixed_goog)[cols_in_fixed_goog]
	
	# default (no spreading out of Google Hits across multiple hits of quote)
	divide_by_factor = 1
	for (i in 1:nrow(quantchars))
		{
		tmp_quotetxt_in_bigmatrix_list = c(quotetxt[i])
		tmp_quotetxt_in_bigmatrix = quotetxt[i]
		matching_unique_quotetxt_row = get_indices_where_list1_occurs_in_list2(tmp_quotetxt_in_bigmatrix_list, unique_quotes_to_google)
		#print(matching_unique_quotetxt_row)
		
		# get number of other guys showing that match in the big list
		num_in_biglist_matching = sum(quotes_list_no_extra_spaces_to_google == tmp_quotetxt_in_bigmatrix)
		print(num_in_biglist_matching)
		
		# if you do want to change the numbers by dividing equally amongst multiple hits of the quote,
		# do so here...
		if (divide_by_frequency == TRUE)
			{
			divide_by_factor = num_in_biglist_matching
			} else {
			pass = "blah"
			}

		tmprow = as.numeric(fixed_goog[matching_unique_quotetxt_row, cols_in_fixed_goog]) / divide_by_factor
		numchars = nchar(tmp_quotetxt_in_bigmatrix)
		numwords = length(strsplit(tmp_quotetxt_in_bigmatrix, " ")[[1]])
		
		tmprow2 = c(as.data.frame(numchars), as.data.frame(numwords), tmprow)
		tmpmat = rbind(tmpmat, tmprow2)
		}
	tmpmat2 = cbind(as.data.frame(rownames), as.data.frame(quotetxt), as.data.frame(tmpmat, row.names=1:nrow(tmpmat)))

	outdf = as.data.frame(tmpmat2, row.names=1:nrow(tmpmat2))
	names(outdf) = c("rownames", "quotetxt", "numchars", "numwords", names_cols_in_fixed_goog)

	quantchar_copy_nums_to_repeated_quotes = outdf
	return(quantchar_copy_nums_to_repeated_quotes)

	}


# If a list has [[1]], [[2]], or is a dataframe, kill that with this
get_items_from_nested_list <- function(nested_list)
	{
	tmpout = c()
	for (i in 1:length(nested_list))
		{
		item = nested_list[[i]]
		tmpout = c(tmpout, item)
		}
	return(tmpout)
	}

# If a list has [[1]], [[2]], or is a dataframe, kill that with this
get_items_from_nonnested_list <- function(nested_list)
	{
	tmpout = c()
	for (i in 1:length(nested_list))
		{
		item = nested_list[i]
		tmpout = c(tmpout, item)
		}
	return(tmpout)
	}




#######################################################
# Compiler switcheroo
#######################################################
junk='
switch_compiler_to_new()
switch_compiler_to_default()
junk='


# Change the compiler to e.g. g++46
switch_compiler_to_new <- function(old_compiler_symlink="/usr/bin/g++", new_compiler_fn="/my_gcc/bin/g++46", new_R_Makevars="/Users/nickm/.R/Makevars46")
	{
	
	# Replace the symlink to the compiler
	cmdstr = paste("rm ", old_compiler_symlink, "; ln -s ", new_compiler_fn, " ", old_compiler_symlink, sep="")
	system(cmdstr)
	
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	# Check the compiler
	cmdstr2 = "g++ -v"
	system(cmdstr2)
	
	return(cmdstr)
	}


# Change the compiler back to e.g. default g++42
switch_compiler_to_default <- function(old_compiler_symlink="/usr/bin/g++", new_compiler_fn="/usr/bin/llvm-g++-4.2", new_R_Makevars="/Users/nickm/.R/Makevars_default")
	{
	#######################################################
	# NOTE: ln -s works like this:
	# 
	# ln -s THING_TO_LINK_TO THING_THAT_POINTS_TO_SOMETHING_ELSE
	# 
	#######################################################
	
	# old_compiler_symlink="/usr/bin/g++"; new_compiler_fn="/usr/bin/llvm-g++-4.2"; new_R_Makevars="/Users/nickm/.R/Makevars_default"
	
	# Replace the symlink to the compiler
	cmdstr = paste("rm ", old_compiler_symlink, "; ln -s ", new_compiler_fn, " ", old_compiler_symlink, sep="")
	system(cmdstr)
	
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	# Check the compiler
	cmdstr2 = "g++ -v"
	system(cmdstr2)

	return(cmdstr)
	}


# Change the gfortran to e.g. 4.6
switch_gfortran_to_new <- function(old_compiler_symlink="/usr/local/bin/gfortran", new_compiler_fn="/usr/local/gfortran/bin/gfortran")
	{
	
	# Replace the symlink to the compiler
	cmdstr = paste("rm ", old_compiler_symlink, "; ln -s ", new_compiler_fn, " ", old_compiler_symlink, sep="")
	system(cmdstr)
		
	# Check the compiler
	cmdstr2 = "gfortran -v"
	system(cmdstr2)
	
	return(cmdstr)
	}


# Change the gfortran back to e.g. default gfortran42
switch_gfortran_to_default <- function(old_compiler_symlink="/usr/local/bin/gfortran", new_compiler_fn="/usr/bin/gfortran-4.2")
	{
	
	# Replace the symlink to the compiler
	cmdstr = paste("rm ", old_compiler_symlink, "; ln -s ", new_compiler_fn, " ", old_compiler_symlink, sep="")
	system(cmdstr)
	
	
	# Check the compiler
	cmdstr2 = "gfortran -v"
	system(cmdstr2)

	return(cmdstr)
	}



# Change default Makevars for install
new_Makevars <- function(new_R_Makevars="/Users/nickm/.R/Makevars_phyRmcmc")
	{
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	return(cmdstr)
	}

# Change default Makevars for install
old_Makevars <- function(new_R_Makevars="/Users/nickm/.R/Makevars_default")
	{
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	return(cmdstr)
	}


# Check gcc version
check_gcc <- function(gccstr="/usr/bin/gcc")
	{
	cmdstr = paste(gccstr, " -v", sep="")
	system(cmdstr)
	}

# Check g++ version
check_gpp <- function(gccstr="/usr/bin/g++")
	{
	cmdstr = paste(gccstr, " -v", sep="")
	system(cmdstr)
	}

# Check gfortran version
# use "which gfortran" to find out where it is
check_gfortran <- function(gccstr="/usr/bin/local/gfortran")
	{
	cmdstr = paste("which gfortran", sep="")
	system(cmdstr)

	cmdstr = paste(gccstr, " -v", sep="")
	system(cmdstr)
	}

