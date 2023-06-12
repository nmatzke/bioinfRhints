



#######################################################
# 2023-03-13 Lab meeting
# Scripts: manipulating alignments, sequence names / tip names, etc.
#          Extracting gene order from a list of downloaded genomes
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R

# Set working directory
#wd = "/GitHub/bioinfRhints/flag/gene_order_example/" # example
wd = "~/Downloads/Full_genomes/"	 # full dataset
setwd(wd)

# Get the list of genome directories that have been saved somewhere

# DOWNLOAD FROM:
#
# OneDrive > Caroline Puente-Lelievre > Flagellum > data > genomes https://uoa-my.sharepoint.com/personal/cpue388_uoa_auckland_ac_nz/_layouts/15/onedrive.aspx?e=5%3Ae28ad7a76b0041e784fca85609803fd4&at=9&FolderCTID=0x0120001084CF6CCC61BD4F9CD6A8F986921623&id=%2Fpersonal%2Fcpue388%5Fuoa%5Fauckland%5Fac%5Fnz%2FDocuments%2FFlagellum%2Fdata%2Fgenomes&view=0


# Run with:
runcmds='
cd /GitHub/bioinfRhints/flag/assemble_all_genome_feature_tables/
Rscript assemble_all_genome_feature_tables_v1.R
'


genome_dirs = list.files(path="genomes", pattern=NULL, recursive=FALSE)
genome_dirs

prot_feature_tables_all_fn = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/2023-06-12_prot_feature_tables_all_v1.txt"

local_unzipped_location = paste0("~/Downloads/", lastword(prot_feature_tables_all_fn))


unzipTF = TRUE
list_of_protIDs_in_genome = list()
cat("\nMerging protein feature tables (gene order tables) from ", length(genome_dirs), " genome directories...\n")
for (i in 1:length(genome_dirs))
	{
	if (i == 1)
		{
		cat("\ngenome_dirs #", i, ",", sep="")
		} else {
		cat(i, ",", sep="")
		}
	genome_dir = genome_dirs[i]
	# Access that directory, look for protein in list
	gene_order_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_feature_table.txt")
	if ((file.exists(gene_order_table_fn) == FALSE) && (unzipTF == TRUE))
		{
		gene_order_table_zipfn = slashslash(paste0(wd, "/", "genomes/", genome_dir, "/", genome_dir, "_feature_table.txt.gz"))
		cmdstr = paste0("gunzip ", gene_order_table_zipfn)
		system(cmdstr)
		}

	prot_feature_table_df = read_gene_order_table(gene_order_table_fn=gene_order_table_fn)
	dim(prot_feature_table_df)

	# Also get the species name from the 
	# GCA_000005845.2_ASM584v2_assembly_report.txt
	assembly_report_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_assembly_report.txt")
	spname = get_spname_from_assembly_report(assembly_report_table_fn)

	# Add some stuff to the front of the line
	id = genome_dir
	prot_feature_table_df = cbind(id, prot_feature_table_df)
	prot_feature_table_df = cbind(spname, prot_feature_table_df)
	
	#prot_feature_tables_all_df = rbind(prot_feature_tables_all_df, prot_feature_table_df)

	if (i == 1)
		{
		write.table(prot_feature_table_df, file=prot_feature_tables_all_fn, sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE)
		} else {
		write.table(prot_feature_table_df, file=prot_feature_tables_all_fn, sep="\t", row.names=FALSE, append=TRUE, quote=FALSE, col.names=FALSE)
		}
	} # END for (i in 1:length(genome_dirs))


file.copy(prot_feature_tables_all_fn, to=local_unzipped_location, overwrite=TRUE)


junk='
prot_feature_tables_all_df = read.table(prot_feature_tables_all_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
dim(prot_feature_tables_all_df)

head(prot_feature_tables_all_df)

#prot_feature_tables_all_df[3832:3836,]
dim(prot_feature_tables_all_df)
'

# system(paste0("open ", prot_feature_tables_all_fn))