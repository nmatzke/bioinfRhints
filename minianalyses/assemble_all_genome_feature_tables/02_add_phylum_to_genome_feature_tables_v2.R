
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

spname = paste0(xlsx$Genus, " ", xlsx$Species, " ", xlsx$Strain)
spname

spname = gsub(pattern="serovar", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="unclassified", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="uncultured", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Candidatus", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="chromosome ", replacement="chrom", x=spname, ignore.case=TRUE)
spname = gsub(pattern="sp.", replacement="sp", x=spname, ignore.case=TRUE)
spname = gsub(pattern="WGS isolate:", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="plasmid ", replacement="pl_", x=spname, ignore.case=TRUE)

spname = stringr::str_trim(stringr::str_squish(spname))
sum(grepl("Wolfebacteria Wolfebacteria", x=spname))
sum(grepl("Wolfebacteria_Wolfebacteria", x=spname))


spname = gsub(pattern="Korarchaeota Candidatus", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Chlamydiales bacterium Chlamydiales bacterium", replacement="Chlamydiales bacterium", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Candidatus Prometheoarchaeum Candidatus Prometheoarchaeum", replacement="Prometheoarchaeum", x=spname, ignore.case=TRUE)

spname = gsub(pattern="of Drosophila melanogaster isolate wMel", replacement="Dmelanogaster", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Beckwithbacteria Beckwithbacteria", replacement="Beckwithbacteria", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Campbellbacteria Campbellbacteria", replacement="Campbellbacteria", x=spname, ignore.case=TRUE)

spname = gsub(pattern="Moranbacteria Moranbacteria", replacement="Moranbacteria", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Uhrbacteria Uhrbacteria", replacement="Uhrbacteria", x=spname, ignore.case=TRUE)

spname = gsub(pattern="Wolfebacteria Wolfebacteria", replacement="Wolfebacteria", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Woesebacteria Woesebacteria", replacement="Woesebacteria", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Wolfebacteria_Wolfebacteria", replacement="Wolfebacteria", x=spname, ignore.case=TRUE)
spname = gsub(pattern="Woesebacteria_Woesebacteria", replacement="Woesebacteria", x=spname, ignore.case=TRUE)

sum(grepl("Wolfebacteria Wolfebacteria", x=spname))
sum(grepl("Wolfebacteria_Wolfebacteria", x=spname))

spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)




spname = gsub(pattern=";", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="\\(", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="\\)", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="\\[", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="\\]", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="=", replacement="EQ", x=spname, ignore.case=TRUE)
spname = gsub(pattern="\\\\", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="str. ", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="\\.", replacement="", x=spname, ignore.case=TRUE)
spname = gsub(pattern="/", replacement="-", x=spname, ignore.case=TRUE)
spname = gsub(pattern=":", replacement="", x=spname, ignore.case=TRUE)

# Remove extra spaces
spname = stringr::str_trim(stringr::str_squish(spname))

# Add underscores
spname = gsub(pattern=" ", replacement="_", x=spname, ignore.case=TRUE)
sum(grepl("Wolfebacteria_Wolfebacteria", x=spname))

#spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
#spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
#spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)
#spname = gsub(pattern="", replacement="", x=spname, ignore.case=TRUE)

sort(unique(spname))


#######################################################
# Check that all genomes in xlsx have matches in prot_feature_tables_all_df
#######################################################

head(prot_feature_tables_all_df)

spnames = sort(unique(prot_feature_tables_all_df$spname))
genome_ids = sort(unique(prot_feature_tables_all_df$assembly))

strain_name = paste0(xlsx$Genus, " ", xlsx$Species, " ", xlsx$Strain, sep="")
strain_name


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


xlsx2 = cbind(xlsx, group, spname)

write.table(xlsx2, file=genomes_to_spnames_fn, sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE)

moref(genomes_to_spnames_fn)






