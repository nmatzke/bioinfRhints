#!/usr/bin/env Rscript
runtxt='
cd /GitHub/bioinfRhints/terminal/
Rscript --vanilla ex1.R 1 2
' # END runtxt

# Read in arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0)
	{
  stop("At least one input argument must be supplied.n", call.=FALSE)
	} else if (length(args)==1) {
  # default output file
  args[2] = 1
	} # END if (length(args)==0)

args = as.numeric(args)
sumval = args[1] + args[2]
txt = paste0(args[1], " + ", args[2], " = ", sumval)
cat(txt)
cat("\n")
