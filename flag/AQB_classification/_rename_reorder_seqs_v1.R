
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

wd = "~/Downloads/iqtree_genomes_wLitLinks/"
setwd(wd)



#######################################################
# Read translation table
# generated at: ~/Downloads/z_genomes_wLitLinks_processing/_cmds_unzip_cat_HMMER_v1.txt
#######################################################
translation_seqsfn = "AQBs_orig+flag1-5_wLitLinks.fasta"
infn = "AQBs_orig+flag1-5_wLitLinks_table.txt"
translate_df = read.table(file=infn, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
head(translate_df)
tail(translate_df)

# Alignment used for tree
alnfn = "1214_seqs_merged_mafftMiddleConstrained2.fasta"
aln = seqinr::read.fasta(alnfn, seqtype="AA")
fullnames = fullnames_from_readFasta(aln)
aln_gids = names(aln)
for (i in 1:length(aln_gids))
	{
	words = strsplit(aln_gids[i], split="\\|")[[1]]
	aln_gids[i] = words[2]
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


# Tree
trfn = "1214_seqs_merged_mafftMiddleConstrained2.fasta.treefile"
tr = read.tree(trfn)
tr2 = phytools::midpoint.root(tr)
tr3 = ladderize(phy=tr2, right=TRUE)
tr3 = read.tree(file="", text=write.tree(tr3, file=""))
plot(tr3, show.tip.label=FALSE)
tr3


# Parse tip labels
tipnames = tr3$tip.label
head(tipnames)
tipnames3 = tipnames

for (i in 1:length(tipnames))
	{
	words = strsplit(x=tipnames[i], split="\\|")[[1]]
	tipnames3[i] = words[2]
	}

tipnames3


# Match to table
matches1 = match(x=tipnames3, table=translate_df$gids)
matches2 = match(x=translate_df$gids, table=tipnames3)
TF = is.na(matches1)
sum(TF)


# Put the seqs into tr3 tip order
head(tipnames3)
head(translate_df$gids[matches1])
translate_df3 = translate_df[matches1,]
head(translate_df3)

short_desc3 = classify_MotAfam_labels(list_of_strings=translate_df3$desc)
head(short_desc3)
rev(sort(table(short_desc3)))


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

xls2 = cbind(xls, PrimSec_entflag, groupTax)


# Rename the tipnames
length(tipnames3)
dim(translate_df3)
length(short_desc3)

# Get genome name for each tip
head(tipnames3)

matches_to_xls1 = match(x=translate_df3$genome_id, table=xls2$GenBank.ID)
sum(is.na(matches_to_xls1))
translate_df3[is.na(matches_to_xls1),] # All match

head(translate_df3$genome_id)
head(xls2$GenBank.ID[matches_to_xls1])


rev(sort(table(xls2$groupTax)))
rev(sort(table(xls2$group)))
tipnames3_new = paste0(xls2$groupTax[matches_to_xls1], "|", xls2$group[matches_to_xls1], "|", xls2$PrimSec_entflag[matches_to_xls1], "|", tipnames3, "|", short_desc3, " [", translate_df3$taxon, "]")
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

tr3_uniq = drop.tip(tr3, tip=tipnums_to_drop)
tipnames3_new_uniq = tipnames3_new[-tipnums_to_drop]
tipnames3_uniq = tipnames3[-tipnums_to_drop]

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
title("IQtree LG+G+F on 1187 MotA homologs, orig+flag1-5")

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

