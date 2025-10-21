# Installation
# (hit 'n' for no when asked to update a zillion other packages)
install_txt='
library(BiocManager)
BiocManager::install(c("annotate", "Biostrings", "genomes"))
install.packages("rentrez") # for entrez_fetch
'



# Setup
library(seqinr)
library(BioGeoBEARS)	# for cls.df
library(gdata)			# for trim
library(XML)
#library(parallel)
#library(S4Vectors)
#library(stats4)
#library(IRanges)
#library(XVector)
#library(Biostrings)
#library(BiocGenerics)
library(rentrez) # Fetching full records: entrez_fetch()

library(openxlsx)
library(ape)

#library(annotate)		# for Bioconductor's 
						#"genbank" (efetch by 
						# gene Accession #s)
# BLAST and download related sequences


# HELP:
# http://www.ncbi.nlm.nih.gov/books/NBK3837/figure/EntrezHelp.F4/?report=objectonly
# http://news.open-bio.org/news/2009/06/ncbi-einfo-biopython/
# http://www.ncbi.nlm.nih.gov/books/NBK179288/
# http://biopython-documentation-chinese-translate.googlecode.com/git/_build/epub/ch8.html
# http://www.ncbi.nlm.nih.gov/Class/NAWBIS/Modules/InfoHubs/Exercises/infohubs_qa_taxon.html
# http://www.ncbi.nlm.nih.gov/refseq/rsg/about/

# Source this function:
# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html
#source("/drives/GDrive/__GDrive_projects/2017-08-08_Matt_Baker_flagellum_speaker/01_BLAST/blastsequences_v4.R")

# Code to use:
source("~/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")


#wd = "/drives/GDrive/__GDrive_projects/2017-08-08_Matt_Baker_flagellum_speaker/01_BLAST/startdata/"
wd = "~/GitHub/bioinfRhints/ex/download_geneIDs/"
setwd(wd)


#######################################################
# Read in Masafumi's CP sequence IDs from Excel
#######################################################
library(openxlsx)

xlsfn = "Lenti_seqs.xlsx"
xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet="Utricularia", startRow=1, colNames=TRUE)

# Things we want to download
IDs = xls$matK

# Subset xls table to just the rows for which we have matKs

# Convert NA to "/"
xls2 = xls
naTA = is.na(xls2$matK)
xls2$matK[naTA == TRUE] = "/"
TF = xls2$matK != "/"
xls2 = xls2[TF,]
xls2$matK
xls2$Species

#######################################################
# Download a list of sequence IDs, via a loop to label things
#######################################################

list_of_download_fns = NULL

cat(paste0("\nDownloading ", length(xls2$matK), " sequence IDs using rentrez::entrez_fetch()\n"))
for (i in 1:length(xls2$matK))
#for (i in 1:5) # test on just 5
	{
	cat(i, ",", sep="")
	seq_ids = c(xls2$matK[i])
	downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")
	outfn = paste0(xls2$Species[i], "_", seq_ids[1], ".fasta")
	list_of_download_fns = c(list_of_download_fns, outfn)
	write(x=downloaded_seqs, file=outfn, sep="")
	} # END for (i in 1:length(xls2$matK))


# Check Utricularia gibba OR619942.1
downloaded_seqs = entrez_fetch(db="nuccore", id=c("OR619942.1"), rettype="fasta_cds_na")
# Write out as a once off
write(x=downloaded_seqs, file="Utricularia_gibba_OR619942.1.fasta", sep="")


# Look at example download
moref("Utricularia_tenuicaulis_NC_058517.1.fasta")
list.files()
list_of_download_fns

list_of_fastas = list.files()
TF1 = grepl(pattern="Utricularia_", x=list_of_fastas)
TF2 = grepl(pattern=".fasta", x=list_of_fastas)
TF = (TF1 + TF2) == 2
list_of_fastas = list_of_fastas[TF]
list_of_fastas

# Pull out the matKs
matK_seqs_txt = list()
list_i = 0
matK_seqs_names = NULL
for (i in 1:length(list_of_fastas))
	{
	fn = list_of_fastas[i]
	#seqs = ape::read.FASTA(file=fn, type="DNA") # THE 'type' IS VERY IMPORTANT HERE
	seqs = seqinr::read.fasta(file=fn, seqtype="DNA")
	#seqnames = names(seqs)
	
	# Just names doesnt work for seqinr seqs, which uses "attributes", specifically "Annot", to hold the full name
	# to get 1 attribute for the 1st sequence in the list:
	attr(x=seqs[[1]], which="Annot")
	
	# sapply applies/runs a function over each item in a list:
	seqnames = sapply(X=seqs, FUN=attr, which="Annot")

	
	# Search for gene=matK in the sequence headers
	TF1 = grepl(pattern="gene=matK", x=seqnames)
	TF2 = grepl(pattern="protein=maturase K", x=seqnames) # alternative way to identify matKs
	TF = (TF1 + TF2) >= 1
	
	# Error trap
	if (sum(TF) < 1)
		{
		stop_txt = paste0("Stop: no matK found for i=", i, ", fn=", fn)
		stop(stop_txt)
		}
	matching_names = seqnames[TF]
	print(matching_names)
	
	# Extract the matK
	nums = 1:length(TF)
	num_for_matK = nums[TF]
	#seqs # from APE, sequences in a "binary" format to save space
	#seqstxt = as.character(seqs) # sequences as letters
	#seqstxt = sapply(X=as.character(seqs), FUN=paste0, collapse="") # sequences as strings (desired)
	#matKtxt = seqstxt[[num_for_matK]]
	matKtxt = seqs[[num_for_matK]]
	matK_seqs_txt[[(list_i=list_i+1)]] = matKtxt
	matK_seqs_names = c(matK_seqs_names, matching_names)
	}
matK_seqs_names = unname(matK_seqs_names) # take off the names on the names
matK_seqs_names

# Remove first character (which is a >)
substr(x=matK_seqs_names[1], start=2, stop=nchar(matK_seqs_names[1]))
matK_seqs_names = sapply(X=matK_seqs_names, FUN=substr, start=2, stop=nchar(matK_seqs_names))

names(matK_seqs_txt) = matK_seqs_names
#matK_seqs_txt

#matK_seqs = as.DNAbin(matK_seqs_txt)
matK_seqs = matK_seqs_txt
matK_seqs

seqinr::write.fasta(sequences=matK_seqs, names=matK_seqs_names, file.out="all_Utricularia_matK_origLabels.fasta")

system("open all_Utricularia_matK_origLabels.fasta")



# Put in species names as sequence labels
spnames = gsub(pattern=".fasta", replacement="_matK", x=list_of_fastas)
spnames

seqinr::write.fasta(sequences=matK_seqs, names=spnames, file.out="all_Utricularia_matK.fasta")
system("open all_Utricularia_matK.fasta")


# Database fields (db)
# Table 1
# – Valid values of &retmode and &rettype for EFetch (null = empty string)
# https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/

# protein
#downloaded_seqs = entrez_fetch(db="protein", id=seq_ids, rettype="fasta")

# get FASTA file of CDS (coding sequences) previously identified by genome project
# rettype = return type = type of file to return
# CDS nucleotide FASTA
downloaded_seqs = entrez_fetch(db="nuccore", id=seq_ids, rettype="fasta_cds_na")

outfn = "CDS_nucleotide_FASTA.fasta"
write(x=downloaded_seqs, file=outfn, sep="")
opd()






#######################################################
# Get sequences by BLAST
#######################################################

install_cmds='
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("rBLAST", repos = "https://mhahsler.r-universe.dev")
' # END install_cmds

library(Biostrings)

library(rBLAST)

