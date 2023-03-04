

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