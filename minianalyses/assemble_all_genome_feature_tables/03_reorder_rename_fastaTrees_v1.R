
#######################################################
# TASK 03:
# Produce various renamings & reorderings
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches
library(openxlsx)			# for openxlsx::read.xlsx

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R
sourceall("/GitHub/bioinfRhints/R/BEASTmasteR/") 
# for remove_equals_from_tips()
# for read.beast.table_original


wd = "/GitHub/bioinfRhints/minianalyses/assemble_all_genome_feature_tables/"
setwd(wd)

prot_feature_tables_all_fn = "~/Downloads/2023-06-12_prot_feature_tables_all_v1.txt"
genomes_to_spnames_fn = "species_list_10062023_NJM+group+spname_v1.txt"

# Takes a few seconds
prot_feature_tables_all_df = read.table(prot_feature_tables_all_fn, header=TRUE, comment.char="%", quote="\"", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
head(prot_feature_tables_all_df)

genomes_to_spnames_df = read.table(genomes_to_spnames_fn, header=TRUE, comment.char="%", quote="", sep="\t", fill=TRUE, stringsAsFactors=FALSE)
genomes_to_spnames_df[1:10,]
dim(genomes_to_spnames_df)



# Order an alignment against a tree
trfn = "motA_hmmalign589_full_tr2_order_domainsOrder_cutNonHom_ed1_MiddleConstrained3.fasta.treefile_FigTree.nexus"

alnfn = "motA_hmmalign589_full_tr2_order_domainsOrder_cutNonHom_ed1_MiddleConstrained3.fasta"

# Read tree, remove "'"
tr = read.nexus(trfn)
names(tr)

# Get the node labels (e.g. bootstraps)
#tr = phytools::readNexus(trfn, format="raxml")
# Extract the internal node labels
nodelabels_table = extractBEASTstats_orig(file=trfn, digits=4, printflag=FALSE) 
bootstraps = nodelabels_table$bs

# Add to tree object
if (length(bootstraps) == (tr$Nnode-1))
	{
	# Add a root node label
	bootstraps = c(0, bootstraps)
	}

if (length(bootstraps) == tr$Nnode)
	{
	tr$node.label = bootstraps
	} else {
	txt = "WARNING: the length of bootstraps does not equal the number of internal nodes, tr$Nnode"
	warning(txt)
	}


tr$tip.label = gsub(pattern="'", replacement="", x=tr$tip.label)
tip_gids = tr$tip.label
head(tip_gids)

aln = seqinr::read.fasta(alnfn, seqtype="AA")
fullnames = fullnames_from_readFasta(aln)
seq_gids = names(aln)
head(fullnames)
head(seq_gids)

# Check quickly if all tree tipname Genbank IDs (gids) are in sequence alignment file
TF = tip_gids %in% seq_gids
sum(TF)
length(TF)

# Check quickly if all sequence alignment file Genbank IDs (gids) are in tree tipnames
TF = seq_gids %in% tip_gids
sum(TF)
length(TF)

#######################################################
# Put alignment in tip-order
#######################################################
sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R
aln_to_tip_order = match(x=tip_gids, table=seq_gids)
aln2 = aln[aln_to_tip_order]
fullnames = fullnames_from_readFasta(aln2)
list_of_strings = extract_protein_name_info(fullnames)
abbr_protnames = classify_MotAfam_labels(list_of_strings)

# Check the outputs
sort(unique(abbr_protnames))

head(fullnames)
head(tr$tip.label)

TF = names(aln2) == tip_gids
sum(TF)
length(TF)

head(names(aln2))
head(tr$tip.label)





#######################################################
# Check for any tip IDs not found in the prot_features_tables_all_df
#######################################################

# Find tip_gids in prot_features_table
tip_gids_in_df_matches = match(tip_gids, table=prot_feature_tables_all_df$product_accession)
tip_gids_in_df_matches

nonmatches = tip_gids[is.na(tip_gids_in_df_matches)]
tmplist = c()

if (length(nonmatches) > 0)
	{
	for (i in 1:length(nonmatches))
		{
		num = grep(pattern=nonmatches[i], x=fullnames)
		cat("\n")
		cat(fullnames[num])
		tmplist = c(tmplist, fullnames[num])
		}
	cat("\n")
	}

sort(unique(extract_last_brackets(tmplist)))



#######################################################
# Convert alignment and tree to full tip labels
#######################################################
sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R

# Add new labels, in tree tip order
group_first_tipnames = NULL
seqlength_first_tipnames = NULL
protname_first_tipnames = NULL

cat("\nMaking names for ", length(tip_gids), " tip protein IDs / gids: ")
for (i in 1:length(tip_gids))
	{
	cat(i, ",", sep="")
	
	prot_feature_tables_all_matchnum = match(x=tip_gids[i], table=prot_feature_tables_all_df$product_accession)
	prot_feature_tables_all_df[prot_feature_tables_all_matchnum,]

	assembly = prot_feature_tables_all_df$assembly[prot_feature_tables_all_matchnum]
	group = get_group_from_genome_ID(assembly=assembly, genomes_to_spnames_df=genomes_to_spnames_df, printwarnings=FALSE)	
	spname = get_spname_from_genome_ID(assembly=assembly, genomes_to_spnames_df=genomes_to_spnames_df, printwarnings=FALSE)	
	
	group_first_tipname = paste0(group, "|", abbr_protnames[i], "|", spname, "|", tip_gids[i])
	seqlength = prot_feature_tables_all_df$product_length[prot_feature_tables_all_matchnum]
	seqlength_first_tipname = paste0(seqlength, "|", abbr_protnames[i], "|", spname, "|", tip_gids[i])
	protname_first_tipname = paste0(abbr_protnames[i], "|", group, "|", spname, "|", tip_gids[i])

	group_first_tipnames = c(group_first_tipnames, group_first_tipname)
	seqlength_first_tipnames = c(seqlength_first_tipnames, seqlength_first_tipname)
	protname_first_tipnames = c(protname_first_tipnames, protname_first_tipname)
	} # END for (i in 1:length(tip_gids))
cat("...done.\n")

head(group_first_tipnames)
head(seqlength_first_tipnames)
head(protname_first_tipnames)



# Group first (phylum or subphylum)
tr_group = tr
tr_group$tip.label = group_first_tipnames
outfn = paste0(all_but_suffix(trfn), "_groupFirst.nexus")
write.nexus(tr_group, file=outfn)

# Write alignment in reverse order if desired (ie to match ape plot)
#revnums = rev(1:length(aln2))
aln2fn = paste0(all_but_suffix(alnfn), "_groupFirst.fasta")
write.fasta(sequences=aln2, names=group_first_tipnames, file.out=aln2fn)


# Protein sequence length first
tr_seqlength = tr
tr_seqlength$tip.label = seqlength_first_tipnames
outfn = paste0(all_but_suffix(trfn), "_seqLength.nexus")
write.nexus(tr_seqlength, file=outfn)

aln2fn = paste0(all_but_suffix(alnfn), "_seqLength.fasta")
write.fasta(sequences=aln2, names=seqlength_first_tipnames, file.out=aln2fn)


# Protein abbreviation first
tr_protFirst = tr
tr_protFirst$tip.label = protname_first_tipnames
outfn = paste0(all_but_suffix(trfn), "_protFirst.nexus")
write.nexus(tr_protFirst, file=outfn)

aln2fn = paste0(all_but_suffix(alnfn), "_protFirst.fasta")
write.fasta(sequences=aln2, names=protname_first_tipnames, file.out=aln2fn)






