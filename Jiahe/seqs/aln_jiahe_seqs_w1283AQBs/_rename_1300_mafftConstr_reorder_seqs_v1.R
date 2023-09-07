
#######################################################
# Setup
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches
library(seqinr) 			# for read.fasta
library(openxlsx)			# for openxlsx::read.xlsx

sourceall("/GitHub/bioinfRhints/Rsrc/") # for protein_bioinf_v1.R

#wd = "~/Downloads/iqtree_genomes_wLitLinks/"
#wd = "~/GitHub/bioinfRhints/flag/AQB_classification/"
#wd = "~/GitHub/bioinfRhints/flag/AQB_classification/"
wd = "~/GitHub/bioinfRhints/Jiahe/seqs/aln_jiahe_seqs_w1283AQBs/"
main_wd = "~/GitHub/bioinfRhints/Jiahe/seqs/aln_jiahe_seqs_w1283AQBs/"
setwd(main_wd)



#######################################################
# Seqs to remove
#######################################################

# This one sucks, manually delete:
# 'Gamma|_||ABY96489.1|integral_membrane_sensor_signal_transduction_histidine_kinase [Pseudomonas putida GB-1]'
gid_duds = c("ABY96489.1")


#######################################################
# Open the PowerSource notes
#######################################################
wd = "~/GitHub/bioinfRhints/flag/AQB_classification/powerSource/"
setwd(wd)

powerSource_fn = "powerSource_notes_was_flag_OTUs_v3.xlsx"
powerSource_df = openxlsx::read.xlsx(xlsxFile=powerSource_fn, sheet="powerSource2023")
head(powerSource_df)

keepTF = !is.na(powerSource_df$MotA_tr_tipname)
powerSource_df[!keepTF,]

powerSource_df = powerSource_df[keepTF,]
head(powerSource_df[,1:15])
tail(powerSource_df[,1:15])


#######################################################
# Read translation table
# generated at: ~/Downloads/z_genomes_wLitLinks_processing/_cmds_unzip_cat_HMMER_v1.txt
#######################################################
wd = "~/GitHub/bioinfRhints/flag/AQB_classification/"
setwd(wd)

translation_seqsfn = "1444seqs.fasta"
infn = "1444seqs_table.txt"
translate_df = read.table(file=infn, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
head(translate_df)
tail(translate_df)

# Remove gids from desc
for (i in 1:nrow(translate_df))
	{
	text_to_remove = paste0(translate_df$gids[i], " ")
	translate_df$desc[i] = gsub(pattern=text_to_remove, replacement="", x=translate_df$desc[i])
	}
head(translate_df$desc)

# Alignment used for tree
#alnfn = "1208_seqs_merged.clustal.fasta"
setwd(main_wd)
alnfn = "1283_AQBs_plusJiahe_aln_mafftMiddleConstrained2.fasta"
aln = seqinr::read.fasta(alnfn, seqtype="AA")
fullnames = fullnames_from_readFasta(aln)
aln_gids = names(aln)
for (i in 1:length(aln_gids))
	{
	#words = strsplit(aln_gids[i], split="\\|")[[1]]
	#aln_gids[i] = words[2]
	}
head(aln_gids)
length(aln_gids)
length(unique(aln_gids))
head(rev(sort(table(aln_gids))))
duplicate_gids = names(rev(sort(table(aln_gids))))[rev(sort(table(aln_gids))) > 1]
duplicate_gids

nums_to_drop = rep(0, times=length(duplicate_gids))
for (i in 1:length(duplicate_gids))
	{
	TF = aln_gids == duplicate_gids[i]
	nums = (1:length(aln_gids))[TF]
	nums_to_drop[i] = nums[2]
	}
nums_to_keep_TF = ((1:length(aln_gids)) %in% nums_to_drop) == FALSE
nums_to_keep = (1:length(aln_gids))[nums_to_keep_TF]
length(nums_to_keep)

aln = aln[nums_to_keep]
aln_gids = aln_gids[nums_to_keep]


# Remove "dud" sequences
nums_to_delete = match(x=gid_duds, table=aln_gids)
nums_to_keep = (1:length(aln_gids))[nums_to_keep_TF]
aln = aln[nums_to_keep]
aln_gids = aln_gids[nums_to_keep]


#######################################################
# Load tree
#######################################################
# Original tree, with node labels (bootstraps)
trfn = "1283_AQBs_plusJiahe_aln_mafftMiddleConstrained2.fasta.treefile"
tr = read.tree(trfn)

# Remove "dud" sequences
tr = drop.tip(phy=tr, tip=gid_duds)
tr2 = phytools::midpoint.root(tr)
tr3 = ladderize(phy=tr2, right=TRUE)
orig_trtable = prt(tr3)
tail(orig_trtable$label)
tail(tr3$node.label)

# Load FigTree-rooted tree
#trfn = "1283_AQBs_hmmcore2_ABY96489div100_midroot.newick"
trfn = "1283_AQBs_plusJiahe_aln_mafftMiddleConstrained2.newick"
tr = read.tree(trfn)
tr$tip.label = gsub(pattern="'", replacement="", x=tr$tip.label)

# Remove "dud" sequences
tr = drop.tip(phy=tr, tip=gid_duds)
new_trtable = prt(tr)
tr3 = read.tree(file="", text=write.tree(tr, file=""))

# Match new table node labels to old table nodelabels
new_to_old = match(x=new_trtable$tipnames, table=orig_trtable$tipnames)
head(new_to_old)
tail(new_to_old)
head(new_trtable$tipnames)
head(orig_trtable$tipnames[new_to_old])

tail(new_trtable$tipnames)
tail(orig_trtable$tipnames[new_to_old])
internal_nodenums = (length(tr3$tip.label)+1):(length(tr3$tip.label)+tr$Nnode)
new_trtable$label[internal_nodenums] = orig_trtable$label[internal_nodenums]
tr3$node.label = orig_trtable$label[internal_nodenums]

plot(tr3, show.tip.label=FALSE)
tr3


# Parse tip labels
tipnames = tr3$tip.label
head(tipnames)
tipnames = gsub(pattern="'", replacement="", x=tipnames)
tipnames3 = tipnames

for (i in 1:length(tipnames))
	{
	words = strsplit(x=tipnames[i], split="\\|")[[1]]
	#tipnames3[i] = words[1]
	tipnames3[i] = words[6]
	
	if (is.na(tipnames3[i]))
		{
		print(tipnames[i])
		}
	}

head(tipnames3)
tail(tipnames3)



# Match to table
matches1 = match(x=tipnames3, table=translate_df$gids)
matches2 = match(x=translate_df$gids, table=tipnames3)
TF = is.na(matches1)
sum(TF)

tipnames3[is.na(matches1)]


# Put the seqs into tr3 tip order
head(tipnames3)
head(translate_df$gids[matches1])
translate_df3 = translate_df[matches1,]
head(translate_df3)

short_desc3 = classify_MotAfam_labels(list_of_strings=translate_df3$desc)
head(short_desc3)
rev(sort(table(short_desc3)))

cbind(short_desc3, translate_df3$desc)
rev(sort(table(short_desc3)))

translate_df3 = cbind(translate_df3, short_desc3)

#######################################################
# Rename seqs and tree tips in new files
#######################################################
xlsfn = "species_list_10071623_NJMg.xlsx"
xls = openxlsx::read.xlsx(xlsfn)

# Naming group & by previous paper
rev(sort(table(xls$group)))

TF = xls$group == "Alpha"
xls$group[TF] = "Alpha_FB21"

TF = xls$group == "Alpha1"
xls$group[TF] = "Alpha1_LO7"

TF = xls$group == "Alpha2"
xls$group[TF] = "Alpha2_LO7"

TF = xls$group == "Beta"
xls$group[TF] = "Beta_FB21"

TF = xls$group == "Delta"
xls$group[TF] = "Delta_FB21"

TF = xls$group == "Epsilon"
xls$group[TF] = "Epsilon_FB21"

TF = xls$group == "Gamma"
xls$group[TF] = "Gamma_FB21"

TF = xls$group == "non-enteric Gamma"
xls$group[TF] = "nonEntGamma_L07"

TF = xls$group == "enteric Gamma"
xls$group[TF] = "EntGamma_L07"

TF = xls$group == "non-labeled"
xls$group[TF] = "nonLab_FB21"

# Gamma by taxonomy
groupTax = rep("", times=nrow(xls))
TF = grepl(pattern="Alpha", x=xls$Class); sum(TF)
groupTax[TF] = "Alpha"
TF = grepl(pattern="Beta", x=xls$Class); sum(TF)
groupTax[TF] = "Beta"
TF = grepl(pattern="Delta", x=xls$Class); sum(TF)
groupTax[TF] = "Delta"
TF = grepl(pattern="Gamma", x=xls$Class); sum(TF)
groupTax[TF] = "Gamma"
TF1 = grepl(pattern="Entero", x=xls$Order); sum(TF1)
TF2 = grepl(pattern="Gamma", x=xls$Class); sum(TF2)
TF = (TF1 + TF2) == 2; sum(TF)
groupTax[TF] = "Entero"
TF = grepl(pattern="Epsilon", x=xls$Class); sum(TF)
groupTax[TF] = "Epsilon"
TF = grepl(pattern="Zeta", x=xls$Class); sum(TF)
groupTax[TF] = "Zeta"

unassigned_TF = groupTax == ""
groupTax[unassigned_TF] = xls$Phylum[unassigned_TF]
groupTax = gsub(pattern="Candidatus ", replacement="", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="\\(CPR\\)", replacement="", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="Desulfobacterota_I", replacement="Desulfobacterota", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="Thermodesulfobacteriota", replacement="Desulfobacterota", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="Deinococcota", replacement="Deinococcus-Thermus", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="Firmicutes_C", replacement="Firmicutes", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="Patesci group", replacement="Patesci", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="Spirochaetota", replacement="Spirochaetes", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="bacteria", replacement="", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="bacteriota", replacement="", x=groupTax, ignore.case=FALSE)
groupTax = gsub(pattern="bacterota", replacement="", x=groupTax, ignore.case=FALSE)

groupTax = gsub(pattern="  ", replacement=" ", x=groupTax)
groupTax = gsub(pattern="  ", replacement=" ", x=groupTax)
groupTax = gdata::trim(groupTax)
table(groupTax)
rev(sort(table(groupTax)))

table(xls$group)
rev(sort(table(xls$group)))

# Primary and secondary systems
PrimSec = rep("", times=nrow(xls))
TF = stringr::str_ends(string=xls$spname_in_lit, pattern=" 1_2"); sum(TF, na.rm=TRUE)
PrimSec[TF] = "1_2_FB21"


#TF = stringr::str_ends(string=xls$spname_in_lit, pattern=" 1"); sum(TF, na.rm=TRUE)
#PrimSec[TF] = "1_FB21"

#TF = stringr::str_ends(string=xls$spname_in_lit, pattern="Shewanella oneidensis MR 1")
#PrimSec[TF] = ""

#TF = stringr::str_ends(string=xls$spname_in_lit, pattern=" 2"); sum(TF, na.rm=TRUE)
#PrimSec[TF] = "2_FB21"

# Flag-1, 2, 3a, 3b, 4, 5
entflag = rep("", times=nrow(xls))

xls$"flag-1"[is.na(xls$"flag-1")] = ""
TF1 = ((xls$"flag-1" != "-") + (xls$"flag-1" != "")) == 2

xls$"flag-2"[is.na(xls$"flag-2")] = ""
TF2 = ((xls$"flag-2" != "-") + (xls$"flag-2" != "")) == 2

xls$"flag-3a"[is.na(xls$"flag-3a")] = ""
TF3a = ((xls$"flag-3a" != "-") + (xls$"flag-3a" != "")) == 2

xls$"flag-3b"[is.na(xls$"flag-3b")] = ""
TF3b = ((xls$"flag-3b" != "-") + (xls$"flag-3b" != "")) == 2

xls$"flag-4"[is.na(xls$"flag-4")] = ""
TF4 = ((xls$"flag-4" != "-") + (xls$"flag-4" != "")) == 2

xls$"flag-5"[is.na(xls$"flag-5")] = ""
TF5 = ((xls$"flag-5" != "-") + (xls$"flag-5" != "")) == 2

entflag_matrix = matrix(data=c("1","2","3a","3b","4","5"), nrow=nrow(xls), ncol=6, byrow=TRUE)
head(entflag_matrix)
entflag_matrix[,1][TF1==FALSE] = ""
entflag_matrix[,2][TF2==FALSE] = ""
entflag_matrix[,3][TF3a==FALSE] = ""
entflag_matrix[,4][TF3b==FALSE] = ""
entflag_matrix[,5][TF4==FALSE] = ""
entflag_matrix[,6][TF5==FALSE] = ""
head(entflag_matrix)
entflag_matrix[TF1,]
# Convert to text
entflag = c(apply(X=entflag_matrix, MARGIN=c(1), FUN=paste0, collapse=""))
entflag[entflag != ""]


PrimSec_entflag_matrix = cbind(PrimSec, entflag_matrix)
PrimSec_entflag = c(apply(X=PrimSec_entflag_matrix, MARGIN=c(1), FUN=paste0, collapse=""))
PrimSec_entflag[PrimSec_entflag != ""]

species_xls2 = cbind(xls, PrimSec_entflag, groupTax)


# Rename the tipnames
length(tipnames3)
dim(translate_df3)
length(short_desc3)

# Get genome name for each tip
head(tipnames3)

# GenBank genome matches
matches_to_xls1 = match(x=translate_df3$genome_id, table=species_xls2$GenBank.ID)
sum(is.na(matches_to_xls1))
translate_df3[is.na(matches_to_xls1),]

# RefSeq genome matches
matches_to_species_xls2 = match(x=translate_df3$genome_id[is.na(matches_to_xls1)], table=species_xls2$RefSeq)
matches_to_species_xls2

matches_to_xls1[is.na(matches_to_xls1)] = matches_to_species_xls2
translate_df3[is.na(matches_to_xls1),] # All match


head(translate_df3$genome_id)
head(species_xls2$GenBank.ID[matches_to_xls1])


rev(sort(table(species_xls2$groupTax)))
rev(sort(table(species_xls2$group)))

# Experimental powersource in name
head(powerSource_df[,1:15])
powerSource_into_tipnames = match(x=tipnames3, table=powerSource_df$MotA_tr_tipname)
head(tipnames3)
expPower = powerSource_df$powerSource[powerSource_into_tipnames]

tipnames3_new = paste0(species_xls2$groupTax[matches_to_xls1], "|", expPower, "|sp", species_xls2$group[matches_to_xls1], "|", species_xls2$PrimSec_entflag[matches_to_xls1], "|", tipnames3, "|", short_desc3, " [", translate_df3$taxon, "]")
tipnames3_new = gsub(pattern="\\|NA\\|", replacement="|_|", x=tipnames3_new)

sum(grepl(pattern="=", x=tipnames3_new))
sum(grepl(pattern=";", x=tipnames3_new))
sum(grepl(pattern=",", x=tipnames3_new))


tipnames3_new = gsub(pattern=" = ", replacement="eq", x=tipnames3_new)
tipnames3_new = gsub(pattern="=", replacement="eq", x=tipnames3_new)
tipnames3_new = paste0("'", tipnames3_new, "'")

head(tipnames3_new)
sum(grepl(pattern="=", x=tipnames3_new))
sum(grepl(pattern=";", x=tipnames3_new))
sum(grepl(pattern=",", x=tipnames3_new))




#######################################################
# Check for duplicate tips
#######################################################
tipnames_that_are_duplicated = names(rev(sort(table(tipnames3_new))))
TF = rev(sort(table(tipnames3_new))) > 1
tipnames_that_are_duplicated = tipnames_that_are_duplicated[TF]
tipnames_that_are_duplicated

tipnums_to_drop = rep(0, length(tipnames_that_are_duplicated))
for (i in 1:length(tipnames_that_are_duplicated))
	{
	TF = tipnames3_new %in% tipnames_that_are_duplicated[i]
	tipnums_to_drop[i] = (1:length(tipnames3_new))[TF][2]
	}
tipnums_to_drop

if (!is.na(tipnums_to_drop) == TRUE)
	{
	tr3_uniq = drop.tip(tr3, tip=tipnums_to_drop)
	tipnames3_new_uniq = tipnames3_new[-tipnums_to_drop]
	tipnames3_uniq = tipnames3[-tipnums_to_drop]
	} else {
	tr3_uniq = tr3
	tipnames3_new_uniq = tipnames3_new
	tipnames3_uniq = tipnames3
	}
tr3_groupFirst = tr3_uniq
tr3_groupFirst$tip.label = tipnames3_new_uniq

length(tipnames3_new_uniq)
length(unique(tipnames3_new_uniq))



# Reshuffle alignment
prefix = "groupTax_"
prefix2 = "groupTaxRev_"

outtrfn = paste0(prefix, trfn)
write.tree(phy=tr3_groupFirst, file=outtrfn)


pdffn = paste0(prefix, trfn, ".pdf")
pdf(file=pdffn, width=18, height=96)

ape::plot.phylo(ladderize(tr3_groupFirst, right=FALSE), cex=0.5)
#plot.phylo(tr3_groupFirst, show.tip.label=FALSE)
add.scale.bar()
title("IQtree LG+G+F on 1283 MotA homologs, orig+flag1-5+litLinks etc., mafft-constrained")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)






# Reshuffle alignment
prefix = "groupTax_"
prefix2 = "groupTaxRev_"

outtrfn = paste0(prefix, trfn)
write.tree(phy=tr3_groupFirst, file=outtrfn)


pdffn = paste0(prefix, trfn, "+bs.pdf")
pdf(file=pdffn, width=18, height=96)

ape::plot.phylo(ladderize(tr3_groupFirst, right=FALSE), cex=0.5)

internal_nodenums = (length(tr3_groupFirst$tip.label)+1):(length(tr3_groupFirst$tip.label)+tr$Nnode)
edge_beneath_each_node = new_trtable$parent_br[internal_nodenums]
#nodelabels(text=tr3_groupFirst$node.label, node=internal_nodenums)
edgelabels(text=tr3_groupFirst$node.label, edge=edge_beneath_each_node, frame="none", bg="none", adj=c(0.5,0.0))
#plot.phylo(tr3_groupFirst, show.tip.label=FALSE)
add.scale.bar()
title("IQtree LG+G+F on 1283 MotA homologs, orig+flag1-5+litLinks etc., mafft-constrained")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)








match_aln_to_tipnames3 = match(x=tipnames3_uniq, table=aln_gids)
head(tipnames3_uniq)
head(aln_gids[match_aln_to_tipnames3])



aln3 = aln[match_aln_to_tipnames3]
head(names(aln3))
head(tipnames3_uniq)
head(tipnames3_new_uniq)



names(aln3) = aln_gids[match_aln_to_tipnames3]
outalnfn = paste0(prefix, alnfn)
seqinr::write.fasta(sequences=aln3, names=tipnames3_new_uniq, file.out=outalnfn)
outalnfn = paste0(prefix2, alnfn)
seqinr::write.fasta(sequences=rev(aln3), names=rev(tipnames3_new_uniq), file.out=outalnfn)

# Lengths of each sequence (non-gapped)
num_aas <- function(seq)
	{
	sum(seq != "-")
	}
len = sapply(X=aln3, FUN=num_aas)
hist(len)

#######################################################
# Open Excel file in tip order, confer matching rows to new tree
#######################################################
tipnames_xlsfn = "groupTax_1282_mafftConstr_2023-08-07_edit.xlsx"
tipnames_df = openxlsx::read.xlsx(tipnames_xlsfn)
head(tipnames_df)

convert_xlsx_to_tipnames3 = match(x=tipnames3_uniq, table=tipnames_df$tipnames3_uniq)
head(tipnames3_uniq)
head(tipnames_df$tipnames3_uniq[convert_xlsx_to_tipnames3])
sum(is.na(convert_xlsx_to_tipnames3))
xlsnew = tipnames_df[convert_xlsx_to_tipnames3,]

# Insert the info you have for the new tipnames
convert_translate_df3_to_tipnames3 = match(x=tipnames3_uniq, table=translate_df3$gids)
sum(is.na(convert_translate_df3_to_tipnames3))

head(cbind(tipnames3_uniq, translate_df3$gids[convert_translate_df3_to_tipnames3]))

header_names = names(translate_df3)
for (i in 1:length(header_names))
	{
	txt = paste0("xlsnew$", header_names[i], " = translate_df3$", header_names[i], "[convert_translate_df3_to_tipnames3]")
	eval(parse(text=txt))
	}
head(xlsnew)	
xlsnew[12,]

newTF = is.na(xlsnew$tipnames3_uniq)
sum(newTF)
new = rep("", times=nrow(xlsnew))
new[newTF] = "y"
xlsnew$tipnames3_uniq = tipnames3_uniq
xlsnew$tipnames3_new_uniq = tipnames3_new_uniq

# Add the "new" column
if ("new" %in% names(xlsnew))
	{
	xlsnew$new = new
	} else {
	xlsnew = cbind(new, xlsnew)
	}

# Fill in the genome parts of the "new" rows
matches1 = match(x=xlsnew$genome_id[newTF], table=species_xls2$GenBank.ID)
matches2 = match(x=xlsnew$genome_id[newTF][is.na(matches1)], table=species_xls2$RefSeq)
matches1[is.na(matches1)] = matches2
matches1

header_names = names(species_xls2)
for (i in 1:length(header_names))
	{
	txt = paste0("xlsnew$'", header_names[i], "'[newTF] = species_xls2$'", header_names[i], "'[matches1]")
	eval(parse(text=txt))
	}


# Add the literature power source info
head(powerSource_df[,1:15])
tail(powerSource_df[,1:15])

powerSource_gids = powerSource_df$MotA_tr_tipname
litPower_entries = powerSource_df$powerSource
xls_to_powerSource_df_matches = match(x=powerSource_gids, table=xlsnew$gids)

# Filter for non-matching gids (none as of 2023-08-07)
powerSource_gids[is.na(xls_to_powerSource_df_matches)]
powerSource_gids[!is.na(xls_to_powerSource_df_matches)]
powerSource_gids = powerSource_gids[!is.na(xls_to_powerSource_df_matches)]
litPower_entries = litPower_entries[!is.na(xls_to_powerSource_df_matches)]
xls_to_powerSource_df_matches = xls_to_powerSource_df_matches[!is.na(xls_to_powerSource_df_matches)]

litPower_entries
xlsnew$gids[xls_to_powerSource_df_matches]

litPower = rep(NA, times=nrow(xlsnew))
litPower[xls_to_powerSource_df_matches] = litPower_entries

# Length of each protein


# Split table, insert litPower & length columns
cols_start1 = 1
cols_end1 = (1:ncol(xlsnew))[names(xlsnew) == "lat_polar"]
cols_start2 = (1:ncol(xlsnew))[names(xlsnew) == "tipnames3_uniq"]
cols_end2 = ncol(xlsnew)
xlsnew = cbind(xlsnew[,cols_start1:cols_end1], len, litPower, xlsnew[,cols_start2:cols_end2])


# Write to file
out_tipnames_xlsfn = "groupTax_1282_mafftConstr_2023-08-07.xlsx"
openxlsx::write.xlsx(x=xlsnew, file=out_tipnames_xlsfn)
#out_tipnames_xlsfn = "groupTax_1282_hmmcore_2023-08-05_WORKED.xlsx"
#openxlsx::write.xlsx(x=xlsnew, file=out_tipnames_xlsfn)

system(paste0("open ", out_tipnames_xlsfn))







#######################################################
# MAKE A BIG PLOT
#######################################################
edited_tipinfo_xlsfn = "groupTax_1282_mafftConstr_2023-08-07_edit.xlsx"
xlsnew = openxlsx::read.xlsx(edited_tipinfo_xlsfn)


tr_to_plot = ladderize(tr3_groupFirst, right=FALSE)
tree_age = get_max_height_tree(tr_to_plot)
ntips = length(tr_to_plot$tip.label)


internal_nodenums = (length(tr_to_plot$tip.label)+1):(length(tr_to_plot$tip.label)+tr_to_plot$Nnode)
new_trtable_to_plot = prt(tr_to_plot)
edge_beneath_each_node = new_trtable_to_plot$parent_br[internal_nodenums]
#nodelabels(text=tr3_groupFirst$node.label, node=internal_nodenums)

#plot.phylo(tr3_groupFirst, show.tip.label=FALSE)






# Plot with labels and boxes
pdffn = paste0(prefix, trfn, ".pdf")
pdf(file=pdffn, width=18, height=96)

#nf <- layout(layout_boxes)
#layout.show(nf)
#points(x=1:10, y=1:10, pch="*", cex=10)
#points(x=1:10, y=1:10, pch="*", cex=10)

ape::plot.phylo(tr_to_plot, cex=0.4, show.tip.label=TRUE, label.offset=4, align.tip.label=TRUE)

add.scale.bar()
title("IQtree LG+G+F on 1282 MotA homologs, orig+flag1-5+litLinks etc., mafft-constrained")

edgelabels(text=tr_to_plot$node.label, edge=edge_beneath_each_node, frame="none", bg="none", adj=c(0.5,0.0), cex=0.5)


points(x=tree_age+0, y=ntips)
points(x=tree_age+1, y=ntips)
points(x=tree_age+2, y=ntips)
points(x=tree_age+3, y=ntips)
points(x=tree_age+4, y=ntips)

# Plot data
xlsnums = 1:nrow(xlsnew)

# Plot power source
xval = 0.1
tipnums = xlsnums[!is.na(xlsnew$litPower)]
ys = ntips - (tipnums - 1)
xs = rep(tree_age+xval, times=length(tipnums))
labs = xlsnew$litPower[tipnums]
labs = gsub(pattern="H+_hi_pH_N+", replacement="H>pH>Na", x=labs)
labs = gsub(pattern="\\+", replacement="", x=labs)
labs = gsub(pattern="PROB", replacement="?", x=labs)
labs = gsub(pattern="maybe", replacement="?", x=labs)
labs = gsub(pattern="suggested", replacement="?", x=labs)
labs = gsub(pattern=" or ", replacement="/", x=labs)
cols = rep("gray50", times=length(tipnums))
cols[labs=="H"] = "yellow3"
cols[labs=="H?"] = "yellow3"
cols[labs=="H/Na"] = "orange2"
cols[labs=="Na/H"] = "orange2"
cols[labs=="Na"] = "red"
cols[labs=="Na?"] = "red"
cols[labs=="Na/K"] = "purple"
cols[labs=="Mg2Ca2St2"] = "purple"
text(xs, ys, labels=labs, col=cols, cex=0.65)


# Plot motor type
xval = 0.6
tipnums = xlsnums[!is.na(xlsnew$"flag1-5")]
ys = ntips - (tipnums - 1)
xs = rep(tree_age+xval, times=length(tipnums))
labs = xlsnew$"flag1-5"[tipnums]
cols = rep("gray50", times=length(tipnums))
cols[labs=="ExbB"] = "pink"
cols[labs=="TolQ"] = "grey"
cols[grepl(pattern="Agl", x=labs)] = "gold2"
cols[grepl(pattern="glide", x=labs)] = "gold2"
cols[labs=="MotC"] = "darkblue"
cols[labs=="F5"] = "darkblue"
cols[labs=="MotP"] = "purple"
cols[labs=="PomA"] = "red3"
cols[labs=="F1"] = "coral"
cols[labs=="F2"] = "green"
cols[labs=="F3a"] = "lightblue"
cols[labs=="F3b"] = "blue"
cols[labs=="F4"] = "darkolivegreen"
cols[grepl(pattern="lateral", x=labs)] = "yellow3"
#cols[labs=="F5"] = "pink"
text(xs, ys, labels=labs, col=cols, cex=0.65)


# Make transparent boxes
fraction_treeheight_to_left = 0.1
xleft_orig = 0.0
xleft = xleft_orig - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft_orig - (3/5*fraction_treeheight_to_left * tree_age)
xright = tree_age


key = "TolQ"
tipnums = xlsnums[xlsnew$L2 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)

tmpcol = col2rgb(cols[labs==key])[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=cols[labs==key])


key = "ExbB"
tipnums = xlsnums[xlsnew$L2 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)

tmpcol = col2rgb(cols[labs==key])[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=)



key = "PomA"
tipnums = xlsnums[xlsnew$L2 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)

tmpcol = col2rgb(cols[labs==key])[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=cols[labs==key])



key = "Bacillus_alkaline_MotP"
tipnums = xlsnums[xlsnew$L2 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)

colortxt = "blue"
tmpcol = col2rgb(colortxt)[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=1, col=colortxt)




key = "MotC"
tipnums = xlsnums[xlsnew$L2 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)

tmpcol = col2rgb(cols[labs==key])[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=cols[labs==key])


key = "AglR"
tipnums = xlsnums[xlsnew$L3 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)

tmpcol = col2rgb(cols[labs==key])[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=cols[labs==key])


key = "LafT"
tipnums = xlsnums[xlsnew$L3 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)


#tmpcol = col2rgb(cols[labs=="lateral"])[,1]
colortxt = "darkgreen"
tmpcol = col2rgb(colortxt)[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=colortxt)






key = "entero"
tipnums = xlsnums[xlsnew$L3 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)


#tmpcol = col2rgb(cols[labs=="lateral"])[,1]
colortxt = "seagreen1"
tmpcol = col2rgb(colortxt)[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=colortxt)



key = "BetaGamma"
tipnums = xlsnums[xlsnew$L3 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)


#tmpcol = col2rgb(cols[labs=="lateral"])[,1]
colortxt = "green1"
tmpcol = col2rgb(colortxt)[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=0.8, col=colortxt)




key = "Beta"
tipnums = xlsnums[xlsnew$L3 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)


#tmpcol = col2rgb(cols[labs=="lateral"])[,1]
colortxt = "green2"
tmpcol = col2rgb(colortxt)[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=1, col=colortxt)





key = "Alpha"
tipnums = xlsnums[xlsnew$L3 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)


#tmpcol = col2rgb(cols[labs=="lateral"])[,1]
colortxt = "green3"
tmpcol = col2rgb(colortxt)[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=1, col=colortxt)








key = "DUF3450"
tipnums = xlsnums[xlsnew$L2 == key]
ybottom = ntips - max(tipnums, na.rm=TRUE) + 0.5
ytop = ntips - min(tipnums, na.rm=TRUE) + 1.5
node_for_clade = getMRCA(phy=tr_to_plot, tip=tipnums[!is.na(tipnums)])
age_of_clade = new_trtable_to_plot$time_bp[node_for_clade]
xleft = xleft_orig + tree_age * ((tree_age-age_of_clade)/tree_age) - (fraction_treeheight_to_left * tree_age)
xleft_group_label = xleft + (2/5*fraction_treeheight_to_left * tree_age)


#tmpcol = col2rgb(cols[labs=="lateral"])[,1]
colortxt = "gray70"
tmpcol = col2rgb(colortxt)[,1]
color = rgb(red=tmpcol["red"], green=tmpcol["green"], blue=tmpcol["blue"], alpha=100, maxColorValue=255)
rect(xleft, ybottom, xright, ytop, col=color)

text(x=xleft_group_label, y=mean(c(ytop, ybottom)), label=key, srt=90, pos=NULL, cex=2, col=colortxt)








dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)




