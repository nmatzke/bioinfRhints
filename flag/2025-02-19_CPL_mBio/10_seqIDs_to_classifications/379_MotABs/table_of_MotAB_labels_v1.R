#######################################################
# Make tables of label names for MotA ids, and MotB ids
#######################################################

# 4.2NCBI identical protein retrieval
# https://bioconductor.riken.jp/packages/3.21/bioc/vignettes/ginmappeR/inst/doc/ginmappeR.html
library(ape)
library(BioGeoBEARS)
library(reutils) # for efetch
#library(rentrez) # for entrez_fetch
library(ginmappeR)
library(queryup)

sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/379_MotABs/"
setwd(wd)

motB_seqs_w_MotA_prefix_fn = "motB_filtered.fasta"
motB_seqs_w_MotA_prefix = read_FASTA_safe(motB_seqs_w_MotA_prefix_fn, type="AA")
names(motB_seqs_w_MotA_prefix)
motB_seqnames = all_but_prefixes(names(motB_seqs_w_MotA_prefix), split="_")
motB_seqnames

motA_seqIDs_from_FASTA = get_firstwords(names(motB_seqs_w_MotA_prefix), split="_")
motB_seqIDs_from_FASTA = get_firstwords(motB_seqnames, split="_")
motAB_table_df = as.data.frame(cbind(motA_seqIDs_from_FASTA, motB_seqIDs_from_FASTA), stringsAsFactors=FALSE)
names(motAB_table_df) = c("motA", "motB")
head(motAB_table_df)
write.table(x=motAB_table_df, file="379_motAB_table_df.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



#######################################################
# Get the original labels from GenBank
#######################################################
# Request too large!
# As_from_genbank = rentrez::entrez_fetch(db="protein", id=motA_seqIDs_from_FASTA, rettype="fasta")
# Bs_from_genbank = rentrez::entrez_fetch(db="protein", id=motB_seqIDs_from_FASTA, rettype="fasta")
As_from_genbank = get_fasta_seqs_from_seqids(seqids=motA_seqIDs_from_FASTA, tmpfn="As_from_genbank.fasta", db="protein", type="AA")
Bs_from_genbank = get_fasta_seqs_from_seqids(seqids=motB_seqIDs_from_FASTA, tmpfn="Bs_from_genbank.fasta", db="protein", type="AA")

As_from_genbank_labels = full_seqnames_to_just_label(list_of_strings=names(As_from_genbank), split=" ")
Bs_from_genbank_labels = full_seqnames_to_just_label(list_of_strings=names(Bs_from_genbank), split=" ")

head(As_from_genbank_labels)
head(Bs_from_genbank_labels)




#######################################################
# MotA Beast2 MCC tree
#######################################################
trfn = "A379tree.nexus"
tr = read.nexus(trfn)
seqnames = tr$tip.label
head(seqnames)
seqidsA = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqidsA[1:10])

#######################################################
# MotB Beast2 MCC tree
#######################################################
trfn = "B379tree.nexus"
trB = read.nexus(trfn)
trB$tip.label = gsub(pattern="'", replacement="", x=trB$tip.label)
seqnames = trB$tip.label
head(seqnames)
seqidsB = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqidsB[1:10])

#######################################################
# MotAB Beast2 MCC tree
#######################################################
trfn = "AB379tree.nexus"
trAB = read.nexus(trfn)
trA$tip.label = gsub(pattern="'", replacement="", x=trA$tip.label)
seqnames = trA$tip.label
head(seqnames)
seqidsAB = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqidsAB[1:10])


#######################################################
# Read in metadata from saved tables
#######################################################

NOT_NEEDED='
AQBs_379_df = read.table(file="379_AQBs_uniProt_info_table_v1.txt", header=TRUE, sep="\t", quote="", row.names=NULL, fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE)

BRDs_379_df = read.table(file="379_BRDs_uniProt_info_table_v1.txt", header=TRUE, sep="\t", quote="", row.names=NULL, fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE)
' # NOT_NEEDED

tip_specifers_FIT = c("QDU25289", "ALF47869")
tip_specifers_TGI5 = c("QDU25289", "QXB11670")
tip_specifers_TGI4 = c("ADC89566", "ALF47869")
tip_specifers_GIT = c("ADV63710", "AFE07995")
tip_specifers_CCD2star = c("ADV63710", "QDU94732")
tip_specifers_CCD3star = c("ABG59284", "ACF43024")
tip_specifers_CCD2 = c("AKL98121", "AJJ55217")
tip_specifers_CCD3 = c("QAT17590", "AFE07995")

tip_specifers_nums = match_grepl(patterns=tip_specifers_FIT, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_FIT = tmptr$tip.label
length(tips_in_FIT) # 107
seqids_in_FIT = get_firstwords(OTUs=tips_in_FIT, split="_")

tip_specifers_nums = match_grepl(patterns=tip_specifers_TGI5, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_TGI5 = tmptr$tip.label
length(tips_in_TGI5) # 107
seqids_in_TGI5 = get_firstwords(OTUs=tips_in_TGI5, split="_")

tip_specifers_nums = match_grepl(patterns=tip_specifers_TGI4, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_TGI4 = tmptr$tip.label
length(tips_in_TGI4) # 107
seqids_in_TGI4 = get_firstwords(OTUs=tips_in_TGI4, split="_")

tip_specifers_nums = match_grepl(patterns=tip_specifers_GIT, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_GIT = tmptr$tip.label
length(tips_in_GIT) # 107
seqids_in_GIT = get_firstwords(OTUs=tips_in_GIT, split="_")

tip_specifers_nums = match_grepl(patterns=tip_specifers_CCD2, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_CCD2 = tmptr$tip.label
length(tips_in_CCD2) # 107
seqids_in_CCD2 = get_firstwords(OTUs=tips_in_CCD2, split="_")

tip_specifers_nums = match_grepl(patterns=tip_specifers_CCD3, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_CCD3 = tmptr$tip.label
length(tips_in_CCD3) # 107
seqids_in_CCD3 = get_firstwords(OTUs=tips_in_CCD3, split="_")

tip_specifers_nums = match_grepl(patterns=tip_specifers_CCD2star, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_CCD2star = tmptr$tip.label
length(tips_in_CCD2star) # 107
seqids_in_CCD2star = get_firstwords(OTUs=tips_in_CCD2star, split="_")

tip_specifers_nums = match_grepl(patterns=tip_specifers_CCD3star, x=trAB$tip.label)
node_for_clade = getMRCA(phy=trAB, tip=trAB$tip.label[tip_specifers_nums$matchnums])
tmptr = extract.clade(phy=trAB, node=node_for_clade)
tips_in_CCD3star = tmptr$tip.label
length(tips_in_CCD3star) # 107
seqids_in_CCD3star = get_firstwords(OTUs=tips_in_CCD3star, split="_")


As_from_genbank_labels_orig = As_from_genbank_labels
Bs_from_genbank_labels_orig = Bs_from_genbank_labels

As_from_genbank_labels = classify_MotAfam_labels(list_of_strings=As_from_genbank_labels)
Bs_from_genbank_labels = classify_MotBfam_labels(list_of_strings=Bs_from_genbank_labels)

rev(sort(table(As_from_genbank_labels)))
rev(sort(table(Bs_from_genbank_labels)))

#######################################################
# Use the table() function to get the most common genbank labels
#######################################################
# MotA labels on the AB tree
matchnums = match(x=seqids_in_FIT, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
FIT_labels_A = match_labels

matchnums = match(x=seqids_in_TGI5, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
TGI5_labels_A = match_labels

matchnums = match(x=seqids_in_TGI4, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
TGI4_labels_A = match_labels

matchnums = match(x=seqids_in_GIT, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
GIT_labels_A = match_labels

matchnums = match(x=seqids_in_CCD2, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
CCD2_labels_A = match_labels

matchnums = match(x=seqids_in_CCD3, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
CCD3_labels_A = match_labels

matchnums = match(x=seqids_in_CCD2star, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
CCD2star_labels_A = match_labels

matchnums = match(x=seqids_in_CCD3star, table=motA_seqIDs_from_FASTA)
match_labels = As_from_genbank_labels[matchnums]
CCD3star_labels_A = match_labels


# MotB labels on the AB tree
matchnums = match(x=motAB_table_df[match(seqids_in_FIT, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
FIT_labels_B = match_labels

matchnums = match(x=motAB_table_df[match(seqids_in_TGI5, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
TGI5_labels_B = match_labels

matchnums = match(x=motAB_table_df[match(seqids_in_TGI4, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
TGI4_labels_B = match_labels

matchnums = match(x=motAB_table_df[match(seqids_in_GIT, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
GIT_labels_B = match_labels

matchnums = match(x=motAB_table_df[match(seqids_in_CCD2, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
CCD2_labels_B = match_labels

matchnums = match(x=motAB_table_df[match(seqids_in_CCD3, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
CCD3_labels_B = match_labels

matchnums = match(x=motAB_table_df[match(seqids_in_CCD2star, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
CCD2star_labels_B = match_labels

matchnums = match(x=motAB_table_df[match(seqids_in_CCD3star, table= motAB_table_df$motA),2], table=motB_seqIDs_from_FASTA)
match_labels = Bs_from_genbank_labels[matchnums]
CCD3star_labels_B = match_labels

c1 = names(rev(sort(table(As_from_genbank_labels)))[1:5])
c2 = names(rev(sort(table(Bs_from_genbank_labels)))[1:5])
c3 = names(rev(sort(table(FIT_labels_A)))[1:5])
c4 = names(rev(sort(table(FIT_labels_B)))[1:5])
c5 = names(rev(sort(table(TGI5_labels_A)))[1:5])
c6 = names(rev(sort(table(TGI5_labels_B)))[1:5])
c7 = names(rev(sort(table(TGI4_labels_A)))[1:5])
c8 = names(rev(sort(table(TGI4_labels_B)))[1:5])
c9 = names(rev(sort(table(GIT_labels_A)))[1:5])
c10 = names(rev(sort(table(GIT_labels_B)))[1:5])
c11 = names(rev(sort(table(CCD2_labels_A)))[1:5])
c12 = names(rev(sort(table(CCD2_labels_B)))[1:5])
c13 = names(rev(sort(table(CCD3_labels_A)))[1:5])
c14 = names(rev(sort(table(CCD3_labels_B)))[1:5])
c15 = names(rev(sort(table(CCD2star_labels_A)))[1:5])
c16 = names(rev(sort(table(CCD2star_labels_B)))[1:5])
c17 = names(rev(sort(table(CCD3star_labels_A)))[1:5])
c18 = names(rev(sort(table(CCD3star_labels_B)))[1:5])

top_labels_A = cbind(c1, c3, c5, c7, c9, c11, c13, c15, c17)
top_labels_B = cbind(c2, c4, c6, c8, c10, c12, c14, c16, c18)
top_labels = cbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18)
top_labels


c1 = rev(sort(table(As_from_genbank_labels)))[1:5]
c2 = rev(sort(table(Bs_from_genbank_labels)))[1:5]
c3 = rev(sort(table(FIT_labels_A)))[1:5]
c4 = rev(sort(table(FIT_labels_B)))[1:5]
c5 = rev(sort(table(TGI5_labels_A)))[1:5]
c6 = rev(sort(table(TGI5_labels_B)))[1:5]
c7 = rev(sort(table(TGI4_labels_A)))[1:5]
c8 = rev(sort(table(TGI4_labels_B)))[1:5]
c9 = rev(sort(table(GIT_labels_A)))[1:5]
c10 = rev(sort(table(GIT_labels_B)))[1:5]
c11 = rev(sort(table(CCD2_labels_A)))[1:5]
c12 = rev(sort(table(CCD2_labels_B)))[1:5]
c13 = rev(sort(table(CCD3_labels_A)))[1:5]
c14 = rev(sort(table(CCD3_labels_B)))[1:5]
c15 = rev(sort(table(CCD2star_labels_A)))[1:5]
c16 = rev(sort(table(CCD2star_labels_B)))[1:5]
c17 = rev(sort(table(CCD3star_labels_A)))[1:5]
c18 = rev(sort(table(CCD3star_labels_B)))[1:5]

top_label_counts_A = cbind(c1, c3, c5, c7, c9, c11, c13, c15, c17)
top_label_counts_B = cbind(c2, c4, c6, c8, c10, c12, c14, c16, c18)
top_label_counts = cbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18)
top_label_counts

# MotAs
tmpmat = matrix(data=NA, nrow=5, ncol=18)
nums = seq(1, 18, by=2)
tmpmat[,nums] = top_label_counts_A
nums = seq(2, 18, by=2)
tmpmat[,nums] = top_labels_A

top_hits_As_df = as.data.frame(tmpmat, stringsAsFactors=FALSE)
names(top_hits_As_df) = c("c1", "all", "c2", "FIT", "c3", "TGI5", "c4", "TGI4", "c5", "GIT", "c6", "CCD2", "c7", "CCD3", "c8", "CCD2star", "c9", "CCD3star")
top_hits_As_df

# MotBs
tmpmat = matrix(data=NA, nrow=5, ncol=18)
nums = seq(1, 18, by=2)
tmpmat[,nums] = top_label_counts_B
nums = seq(2, 18, by=2)
tmpmat[,nums] = top_labels_B

top_hits_Bs_df = as.data.frame(tmpmat, stringsAsFactors=FALSE)
names(top_hits_Bs_df) = c("c1", "all", "c2", "FIT", "c3", "TGI5", "c4", "TGI4", "c5", "GIT", "c6", "CCD2", "c7", "CCD3", "c8", "CCD2star", "c9", "CCD3star")
top_hits_Bs_df


write.table(x=top_hits_As_df, file="top_hits_As_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

write.table(x=top_hits_Bs_df, file="top_hits_Bs_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



motAB_table_df[motAB_table_df$motB=="AVQ29523",]

issue_to_fix_TF = motAB_table_df$motA==motAB_table_df$motB
motAB_table_df[issue_to_fix_TF,]
sum(issue_to_fix_TF)

