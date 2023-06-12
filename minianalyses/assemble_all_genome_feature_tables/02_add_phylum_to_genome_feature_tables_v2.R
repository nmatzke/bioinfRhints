
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
prot_feature_tables_all_fn = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/2023-06-12_prot_feature_tables_all_v1.txt"
# Takes a few seconds
prot_feature_tables_all_df = read.table(prot_feature_tables_all_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
head(prot_feature_tables_all_df)

# Genomes information table
xlsfn = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/species_list_10062023_NJM.xlsx"
xlsx = openxlsx::read.xlsx(xlsxFile=xlsfn, sheet=1)
head(xlsx)
sort(unique(xlsx$Phylum))

# Genomes to spnames file
genomes_to_spnames_fn = gsub(pattern=".xlsx", replacement="+group+spname_v1.txt", x=xlsfn)



#######################################################
# Now, create ideal labels for "phylum" and species/strain
#######################################################

# "Group" column - phylum or class, depending
group = xlsx$Phylum


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
	group[TF] = new_terms[i]
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
	group[TF] = new_terms[i]
	}
sort(unique(group))


# Replace some terms
sort(unique(xlsx$Class))

old_terms = c(
"Alphaproteobacteria",
"Betaproteobacteria",
"Gammaproteobacteria",
"Deltaproteobacteria",
"Epsilonproteobacteria",
"Zetaproteobacteria")

new_terms = c(
"Alpha-proteo",
"Beta-proteo",
"Gamma-proteo",
"Delta-proteo",
"Epsilon-proteo",
"Zeta-proteo")

for (i in 1:length(old_terms))
	{
	TF = xlsx$Class == old_terms[i]
	group[TF] = new_terms[i]
	}

sort(unique(group))



# Create "spname", a short & program-universal species/strain name
# ESPECIALLY, REMOVE THESE CHARACTERS, WHICH MESS UP NEWICK FILES:
# ( ) ; , \ / [ ]
strain_name = paste0(xlsx$Genus, " ", xlsx$Species, " ", xlsx$Strain, sep="")
spname = paste0(xlsx$Genus, " ", xlsx$Species, " ", xlsx$Strain)
spname = clean_gb_spname(spname)
sort(unique(spname))


#######################################################
# Check that all genomes in xlsx have matches in prot_feature_tables_all_df
#######################################################

head(prot_feature_tables_all_df)

spnames = sort(unique(prot_feature_tables_all_df$spname))
genome_ids = sort(unique(prot_feature_tables_all_df$assembly))

TF = spnames %in% strain_name
sum(TF)
length(TF)

TF = strain_name %in% spnames
sum(TF)
length(TF)

TF = xlsx$GenBank.ID %in% genome_ids
sum(TF)
length(TF)
sort(xlsx$GenBank.ID[TF == FALSE])

TF = genome_ids %in% xlsx$GenBank.ID
sum(TF)
length(TF)
sort(genome_ids[TF == FALSE])


##############################################
# Newish genomes, 06 2023:
#
# GCA_000009985.1 # Magnetospirillum magneticum AMB-1 (a-proteobacteria)
# GCA_000016645.1	# Flavobacterium johnsoniae UW101 (CFB group bacteria)
# GCA_003258315.1 # Bradymonas sediminis (d-proteobacteria)
# GCA_003812925.1 # Klebsiella oxytoca (enterobacteria)
# GCA_024459955.1	# Anaerovorax odorimutans (firmicutes)
# GCF_000348725.1	# Pseudobdellovibrio exovorus JSS (bacteria)
# GCF_000691605.1 # Bdellovibrio bacteriovorus (bacteria)
##############################################

xlsx2 = cbind(xlsx, group, spname, strain_name)

write.table(xlsx2, file=genomes_to_spnames_fn, sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE)

moref(genomes_to_spnames_fn)






