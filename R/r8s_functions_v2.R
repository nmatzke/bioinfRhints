# R functions for interacting with r8s (Sanderson 200x)

# Load with:
#sourcedir = '/Dropbox/_njm/'
#source3 = 'r8s_functions_v1.R'
#source(paste(sourcedir, source3, sep=""))

# It requires:
# R functions for dealing with trees
#sourcedir = '/Dropbox/_njm/'
#source3 = '_R_tree_functions_v1.R'
#source(paste(sourcedir, source3, sep=""))

# Generic R functions - for checkwd
sourcedir = '/drives/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))




# Run LF or NPRS algorithm in r8s
# Assumes 1 single calibration point
# Returns the name of the log file
defaults = '
tr
calibration_node_tip_specifiers

r8s_method="LF"
addl_cmd=""
calibration_age=1.0
nsites=100
tmpwd=wd
r8s_nexus_fn="r8s_nexus_fn.nex"
r8s_logfile_fn=""
r8s_path="/Applications/r8s"

tr = tmptree3
r8s_method="LF"
addl_cmd=""
calibration_age=1.0
nsites=100
tmpwd=outdir
r8s_nexus_fn="r8s_nexus_fn.nex"
r8s_logfile_fn=""
r8s_path="/Applications/r8s"

'

run_r8s_1calib <- function(tr, calibration_node_tip_specifiers, r8s_method="LF", addl_cmd="", calibration_age=1.0, nsites=10000, tmpwd=getwd(), r8s_nexus_fn="r8s_nexus_fn.nex", r8s_logfile_fn="", r8s_path="/Applications/r8s")
	{

	# checkwd to see if wd ends in slash, which it should
	# checkwd is in _genericR_v1.R
	tmpwd = checkwd(tmpwd)
	
	#print(tmpwd)
	tmp_r8s_nexus_fn = "tmp_r8s_nexus_fn.nex"
	tmp_r8s_nexus_log_fn = "tmp_r8s_nexus_fn.nex.log"
	
	# Figure out name of (permanent) log file
	if (r8s_logfile_fn == "")
		{
		r8s_logfile_fn = paste(r8s_nexus_fn, ".log", sep="")	
		}
	
	
	
	# If check_for_small_branches==TRUE, collapse those branches to 0
	
	
	# Output cmd
	outcmd = paste(r8s_path, " -v -f ", paste(tmpwd, tmp_r8s_nexus_fn, sep=""), " > ", paste(tmpwd, tmp_r8s_nexus_log_fn, sep=""), sep="")

	print(outcmd)
	
	# Backup the old files, if they exist
	fn_bkup(tmpwd, tmp_r8s_nexus_fn)
	fn_bkup(tmpwd, tmp_r8s_nexus_log_fn)
	fn_bkup(tmpwd, r8s_nexus_fn)
	fn_bkup(tmpwd, r8s_logfile_fn)
	
	
	# a list/matrix of the strings to output to the text file
	s = c()

	# count up the line numbers, starting with 0
	h = 0

		
	# script header
	s[(h=h+1)] = "#NEXUS"
	s[(h=h+1)] = "[ Running r8s via Langley-Fitch, from R.             ]"
	s[(h=h+1)] = "[  (using an R function by Nick Matzke, Feb. 2011)   ]"
	s[(h=h+1)] = "[ To run:                                            ] "
	s[(h=h+1)] = paste("[ ", 	outcmd, "     ] ", sep="")
	s[(h=h+1)] = ""
	s[(h=h+1)] = "begin trees;"
	s[(h=h+1)] = "[Note: This tree was output from R.]"
	s[(h=h+1)] = paste("tree Rtree = ", write.tree(tr), sep="")
	s[(h=h+1)] = "end;"
	s[(h=h+1)] = ""
	s[(h=h+1)] = ""
	s[(h=h+1)] = "begin r8s;"
	#
	s[(h=h+1)] = paste("MRCA calib1 ", list2str(calibration_node_tip_specifiers), ";", sep="")
	ultrametric_txt = "no"
	if (is.ultrametric(tr))
		{
		ultrametric_txt = "yes"
		}
	s[(h=h+1)] = paste("blformat lengths=persite ultrametric=", ultrametric_txt, " round=yes nsites=", nsites, ";", sep="")

	s[(h=h+1)] = "collapse;"
	# collapse (removes all zero-length branches from all trees
	# and converts them to hard polytomies)


	s[(h=h+1)] = "[ lengths=persite means that the branch length is in						]"
	s[(h=h+1)] = "[     changes per site not total number of changes							]"
	s[(h=h+1)] = "[ ultrametric=no means that the input tree is not ultrametric				]"
	s[(h=h+1)] = ""
	s[(h=h+1)] = "[ Input some age constraints; a point for Clade1, a ]"
	s[(h=h+1)] = "[ fixage taxon=Clade1 age=150 sets the age of node Clade1 to 150			]"
	s[(h=h+1)] = paste("fixage taxon=calib1 age=", calibration_age, ";", sep="")
	s[(h=h+1)] = ""
	s[(h=h+1)] = "[ constrain taxon=node2 min_age=200 max_age=300 							]"
	s[(h=h+1)] = "[    forces node two to be between 200 and 300. These times are 			]"
	s[(h=h+1)] = "[    measured backwards from the present, so that min_age=200 means			]"
	s[(h=h+1)] = "[    that this divergence happened at least 200 (million?) years ago.		]"
	s[(h=h+1)] = "[    I believe the units are relative, depending on what you input.			]"
	s[(h=h+1)] = "[ constrain taxon=rt min_age=3.5 max_age=1.0; ]"
	s[(h=h+1)] = ""
	s[(h=h+1)] = "[ divtime method=LF starts the fitting algorithm using the Langley-Fitch	]"
	s[(h=h+1)] = "[    method which deduces node times using maximum likelihood of the 		]"
	s[(h=h+1)] = "[    branch lengths assuming a constant rate of substitution				]"
	s[(h=h+1)] = addl_cmd
	s[(h=h+1)] = paste("divtime method=", r8s_method, ";", sep="")
	s[(h=h+1)] = "describe plot=chronogram;"
	s[(h=h+1)] = "describe plot=chrono_description;"
	s[(h=h+1)] = "showage;"
	s[(h=h+1)] = "q;"
	s[(h=h+1)] = "end;"
	s[(h=h+1)] = ""

	# write the strings to a text file
	write.table(s, file=paste(tmpwd, tmp_r8s_nexus_fn, sep=""), quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)

	# Run r8s on the system
	system(outcmd)
	
	# Copy the temporary NEXUS file to the desired fn
	cmdstr = paste("cp ", tmp_r8s_nexus_fn, " ", r8s_nexus_fn)
	system(cmdstr)
	cmdstr = paste("cp ", tmp_r8s_nexus_log_fn, " ", r8s_logfile_fn)
	system(cmdstr)
			
	# Return name of logfile
	return(r8s_logfile_fn)
	}


extract_tree_from_r8slog <- function(logfn="r8s_nexus_fn.nex.log", delimiter=" = ", printall=TRUE)
	{
	# Read a text file into a list of strings
	lines = scan(logfn, what="character", sep="\n", blank.lines.skip = FALSE)
	
	if (printall == TRUE)
		{
		print(paste("\nextract_tree_from_r8slog(", logfn, ", delimiter='", delimiter, "'):\n", sep=""))
		print(lines)
		}
	
	word_to_find = "tree"
	wordnum = 1
	delimiter = " "
	nexusstr = extract_lines_with_word_at_wordnum(lines, word_to_find, wordnum, delimiter, printflag=TRUE)
	
	newickstr = extract_newickstr_from_nexusstr(nexusstr)
	
	utr = read.tree(file="", text=newickstr)
	return(utr)
	}


# This one sucks, compared to extract_rates_from_r8slog2
extract_rates_from_r8slog <- function(logfn="r8s_nexus_fn.nex.log")
	{
	# Read a text file into a list of strings
	lines = scan(logfn, what="character", sep="\n", blank.lines.skip = FALSE)
	
	rates_data = c()
	n = 0
	hyphenscount = 0
	save_data = FALSE
	 
	for (tmpline in lines)
		{
		if (tmpline == "")
			{
			next
			}
		words = strsplit_whitespace(tmpline)
		#print(words)
		if ((words[1] == "Node") & (words[2] == "Fix") & (words[3] == "[Mod]"))
			{
			headers = words
			save_data = TRUE
			next
			}
		if (save_data == TRUE)
			{
			# skip "--------------"
			if (substr(tmpline, start=1, stop=5) == "-----")
				{
				hyphenscount = hyphenscount + 1
				if (hyphenscount > 1)
					{
					save_data = FALSE
					next
					}
				else
					{
					next
					}
				}
			# save the data
			tmprow = strsplit_on_tabs_remove_whitespace(tmpline)
			rates_data = rbind(rates_data, tmprow)
			}
		}
	return(rates_data)
	}


extract_rates_from_r8slog2 <- function(logfn="r8s_nexus_fn.nex.log")
	{
	# Read a text file into a list of strings
	lines = scan(logfn, what="character", sep="\n", blank.lines.skip = FALSE)
	
	rates_data = c()
	n = 0
	hyphenscount = 0
	save_data = FALSE
	 
	for (tmpline in lines)
		{
		if (tmpline == "")
			{
			next
			}
		words = strsplit_whitespace(tmpline)
		#print(words)
		if ((words[1] == "Node") & (words[2] == "Fix") & (words[3] == "[Mod]"))
			{
			rates_data[[(n=n+1)]] = tmpline
			save_data = TRUE
			next
			}
		if (save_data == TRUE)
			{
			# skip "--------------"
			if (substr(tmpline, start=1, stop=5) == "-----")
				{
				hyphenscount = hyphenscount + 1
				if (hyphenscount > 1)
					{
					save_data = FALSE
					next
					}
				else
					{
					next
					}
				}
			# save the data
			rates_data[[(n=n+1)]] = tmpline
			}
		}
	
	fn = "tmp_rates_data.txt"
	write.table(rates_data, file=fn, quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)
	
	rates_dtf = read_table_good(fn)
	names(rates_dtf) = c("Node", "Fix..Mod.", "Min_Max", "Age", "X0", "Estimated", "Local", "X1", "X2")
	
	rates_dtf$X0 = NULL 
	rates_dtf$X1 = NULL 
	rates_dtf$X2 = NULL 
	
	return(rates_dtf)
	}


get_summary_rate_variation <- function(rates_dtf)
	{
	cat(" \n")
	cat("get_summary_rate_variation(rates_dtf): printing averages of branches\n")
	summary_rate_variation = c()
	summary_rate_variation$mean = mean(rates_dtf$Local, na.rm=TRUE)
	summary_rate_variation$sd = sd(rates_dtf$Local, na.rm=TRUE)
	summary_rate_variation$min = min(rates_dtf$Local, na.rm=TRUE)
	summary_rate_variation$max = max(rates_dtf$Local, na.rm=TRUE)
	summary_rate_variation$range = summary_rate_variation$max - summary_rate_variation$min
	summary_rate_variation$ratio = summary_rate_variation$max / summary_rate_variation$min
	
	#print(summary_rate_variation)
	return(summary_rate_variation)
	}


#run_chrono_r8s_lotsa_ways()
#	{
#	
#	}

