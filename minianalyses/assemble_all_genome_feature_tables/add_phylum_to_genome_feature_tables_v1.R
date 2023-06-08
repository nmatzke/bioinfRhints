
#######################################################
# Combine feature table information with e.g. Phylum information
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches
library(openxlsx)			# for openxlsx::read.xlsx

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R

# Protein feature tables (from all genomes of interest; see: assemble_all_genome_feature_tables_v1.R
prot_feature_tables_all_fn = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/2023-06-02_prot_feature_tables_all_v1.txt"
# Takes a few seconds
prot_feature_tables_all_df = read.table(prot_feature_tables_all_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
head(prot_feature_tables_all_df)

# Genomes information table
xlsfn = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/species_list_08052023.xlsx"
xlsx = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet=1)
head(xlsx)
sort(unique(xlsx$Phylum))

# Replace some terms
old_terms = c(
"Patescibacteria (CPR)",
"Patescibacteria group",
"Deinococcus-Thermus",
"Candidatus Omnitrophica")

new_terms = c(
"Patescibacteria",
"Patescibacteria",
"Deino_Thermus",
"Omnitrophica")

for (i in 1:length(old_terms))
	{
	TF = xlsx$Phylum == old_terms[i]
	xlsx$Phylum[TF] = new_terms[i]
	}
sort(unique(xlsx$Phylum))


TF = xlsx$Phylum == "Archaea"
sort(unique(xlsx$Class[TF]))
old_terms = c("Asgard group",
"Euryarchaeota",
"TACK group")
new_terms = c("ARCHAEA_Asgard",
"ARCHAEA_Euryarchaeota",
"ARCHAEA_TACK")

for (i in 1:length(old_terms))
	{
	TF = xlsx$Class == old_terms[i]
	xlsx$Phylum[TF] = new_terms[i]
	}
sort(unique(xlsx$Phylum))

