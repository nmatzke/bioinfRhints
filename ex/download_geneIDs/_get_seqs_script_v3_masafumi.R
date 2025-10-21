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



#######################################################
# Download a list of sequence IDs
#######################################################
seq_ids = c("")

downloaded_seqs = entrez_fetch(db="protein", id=seq_ids, rettype="fasta")

outfn = "62_downloaded_seqs.fasta"
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

