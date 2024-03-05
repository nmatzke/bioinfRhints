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
source("/GitHub/bioinfRhints/ex/download_geneIDs/blastsequences_v4.R")

#wd = "/drives/GDrive/__GDrive_projects/2017-08-08_Matt_Baker_flagellum_speaker/01_BLAST/startdata/"
wd = "/GitHub/bioinfRhints/ex/download_geneIDs/"




seq_ids = c("QEI18926.1","WP_002721884.1","RLG45601.1","RMG20999.1","DGBBEKCF_02484","QEI19744.1","CAG37394.1","ABD81419.1","MBU2985819.1","ABD82476.1","MBU2985330.1","SCA58681.1","WP_015920816.1","ABD81321.1","MBU2985725.1","QEI19442.1","UEO04809.1","UEO06636.1","ABD81792.1","MBU2984437.1","QEI18804.1","ABD82818.1","MBU2986144.1","QEI20720.1","UEO04322.1","WP_002720852.1","ABD82214.1","MBU2985168.1","UEO05991.1","UEO03725.1","QEI19895.1","CAG37138.1","DGBBEKCF_01299","UEO07578.1","UEO05650.1","QEI18324.1","CAG37703.1","WP_015921325.1","MBU2986356.1","DGBBEKCF_01621","ABD79617.1","DGBBEKCF_00127","DGBBEKCF_00463","QEI19039.1","UEO02151.1","QEI18407.1","DGBBEKCF_00016","PKPEBJJI_00297","ABD82072.1","MBU2984729.1","CAG37702.1","PKPEBJJI_00564","QEI18325.1","ABD79618.1","MBU2986357.1","MBU2984728.1","ABD82073.1","DGBBEKCF_01764","QEI18408.1","DGBBEKCF_01496","MBU2987451.1","PKPEBJJI_00344")

downloaded_seqs = entrez_fetch(db="protein", id=seq_ids, rettype="fasta")

outfn = "62_downloaded_seqs.fasta"
write(x=downloaded_seqs, file=outfn, sep="")
opd()


# I WILL HAVE TO MANUALLY ADD THESE, FROM PROKKA: DGBBEKCF..., PKPEBJJI...

