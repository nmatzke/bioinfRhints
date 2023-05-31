#######################################################
# Make a table linking fasta-file names to protein IDs, genome IDs, species names, phyla, etc.
#######################################################

library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(openxlsx)				# for read.xlsx (replaces gdata::read.xls and xlsx::read.xlsx)
library(stringr) 			# for regmatches

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R

# Set working directory
# wd = "/GitHub/bioinfRhints/flag/gene_order_example/" # example
wd = "/GitHub/bioinfRhints/flag/protIDs_to_genomes/"
setwd(wd)

protID = "AAC74960.1"  # Escherichia coli str. K12 substr. MG1655


# Get the list of downloaded genomes
genomes_path = "~/Downloads/Full_genomes/genomes"	 # full dataset
genome_dirs = list.files(path=genomes_path, pattern=NULL, recursive=FALSE)
genome_dirs = paste0(genomes_path, "/", genome_dirs)
head(genome_dirs)

# Name & load the Excel file containing metadata (e.g. genome filenames, species name etc.)
xlsfn = "species_list_08052023.xlsx"
xls = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet=1)

# MotA sequence data
#alnfn = "/GitHub/bioinfRhints/flag/gene_order_example/motA_alignment_refined.fasta"
alnfn = "motA_hmmalign589_full_unaligned.fasta"
aln = read.fasta(alnfn)
aln_annotations = fullnames_from_readFasta(aln)
gids = names(aln_annotations)
head(aln_annotations)

# Get the genome taxon
xlstaxa = rep("", times=nrow(xls))
for (i in 1:nrow(xls))
	{
	tmp = paste(xls$Genus[i], xls$Species[i], xls$Strain[i], sep=" ", collapse="")
	xlstaxa[i] = stringr::str_squish(tmp)
	}

# Get the genome director for one protID
protID = "AAC74960.1"  # Escherichia coli str. K12 substr. MG1655
genome_dir = get_genome_name_from_protID(protID, genome_dirs, aln_annotations, xls, xlstaxa)
genome_dir



# Get genome IDs for everything in a FASTA file
protIDs = gids
matches = sapply(X=protIDs, FUN=get_genome_name_from_protID2, genome_dirs=genome_dirs, aln_annotations=aln_annotations, xls=xls, xlstaxa=xlstaxa, returnwhat="genome_dir", max.distance=0.13)
