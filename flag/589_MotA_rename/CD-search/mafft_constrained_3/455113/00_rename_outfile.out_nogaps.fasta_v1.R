library(ape)	# for read/write NEXUS
library(BioGeoBEARS)	# for list2str, moref
library(gdata)	# for trim
library(phytools) # for midpoint.root
library(seqinr)				# for read.fasta

source("/GitHub/bioinfRhints/Rsrc/protein_bioinf_v1.R")
source("~/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R")
#source('/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R')

#wd = "/drives/GDrive/__GDrive_projects/2016-06-16_venerid3/02_TNT/"
wd = "/GitHub/bioinfRhints/flag/589_MotA_rename/CD-search/mafft_constrained_3/455113/"
setwd(wd)


# Source of good names:
alnfn = "motA_hmmalign589_full_tr2_orderMiddle.fasta"
aln = read.fasta(alnfn)
fullnames = fullnames_from_readFasta(aln)

# Alignment to rename
alnfn_to_rename = "outfile.out_nogaps.fasta"
aln_to_rename = read.fasta(alnfn_to_rename)
aln2 = relabel_fastaAln_wGID_first(aln_to_rename, fullnames, split="\\|")
outfn = gsub(pattern=".fasta", replacement="_rename.fasta", x=alnfn_to_rename)
write.fasta(sequences=aln2, names=fullnames_from_readFasta(aln2), file.out=outfn)

