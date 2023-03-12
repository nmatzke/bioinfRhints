

#########################################################
# Loading an IQtree-derived tree, with bootstraps
#########################################################

library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches


nexfn = "/GitHub/bioinfRhints/flag/530_MotA_rename/530_sequences_Alignment_contree_reRootLadder_gIDs_protFirst_containsEQ.nexus"

# Ape's read.nexus fails:
tr = ape::read.nexus(nexfn)

# Error in FUN(X[[i]], ...) : 
# numbers of left and right parentheses in Newick string not equal

# Phytools readNexus fails:
tr = phytools::readNexus(nexfn, format="standard")

# But reading as raxml works:
tr = phytools::readNexus(nexfn, format="raxml")

# But nodelabels don't work
tr$node.label

# Get the 
sourceall("/GitHub/BEASTmasteR/R/")
file=nexfn; digits=4
printflag=FALSE


stats_in_brackets = extractBEASTstats_orig(file=nexfn, digits=4, printflag=FALSE) 

# Tree with statistics attached
tr2stats = read.beast_original(file=nexfn, digits=4) 
names(tr2stats)
head(tr2stats$nodenum)
head(tr2stats$bs)
tail(tr2stats$nodenum)
tail(tr2stats$bs)


# Stable of statistics
trtable = read.beast.table_original(file=nexfn, digits=4) 
