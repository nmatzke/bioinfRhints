library(ape)
library(BioGeoBEARS)
library(ginmappeR)

sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/"
setwd(wd)

# Compare Jiahe MotA labels and Caroline MotA labels
Jiahe_MotAs_fn = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/Jiahe_2024-08_439_MotHomologs/439A/439_SubunitA_AAs_raw.fasta"
Jiahe_MotAs = as.character(read_FASTA_safe(Jiahe_MotAs_fn, type="AA"))
Jiahe_MotA_names = names(Jiahe_MotAs)
Jiahe_MotA_seqids = get_leading_seqids_from_name(strings=Jiahe_MotA_names, split="_")
Jiahe_MotA_seqids

Caroline_MotAs_fn = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/A379tree.nexus"
Caroline_MotAs_tr = read.nexus(Caroline_MotAs_fn)
Caroline_MotA_names = Caroline_MotAs_tr$tip.label
Caroline_MotA_seqids = get_leading_seqids_from_name(strings=Caroline_MotA_names, split="_")
Caroline_MotA_seqids

TF = Caroline_MotA_seqids %in% Jiahe_MotA_seqids
sum(TF)

# 321 of Caroline's MotAs match to Jiahe's; these are before GldLM additiosn

# These are extra Caroline additions; get the MotB seqnums in alignment, then get
# the MotBs, then BLAST those
MotA_seqIDs_to_get_to_379 = Caroline_MotA_seqids[TF == FALSE]

catn(MotA_seqIDs_to_get_to_379)

ex='
AAC07752
AAC73831
ABQ60495
ACL94252
ADB16663
ADE86105
ADV60824
ADV62354
ADV63187
ADV63710
ADY58078
ADY58376
AEI88772
AEM21069
AEM22241
AEO47025
AGR80198
AIK95558
AJJ31614
AJJ31712
AJJ32465
AJJ56099
AKJ65464
ALJ56692
AVM45096
AVM45415
AVQ28424
AVQ29260
AWV88602
BAC58952
BAC59320
BAC61768
BAC62899
CCG56832
EKT86397
EKT86693
EKT87439
EKT88184
EKT88755
QDU26963
QEO41881
QEO43084
QJE97381
QKY94846
QNQ12173
QOY85398
QOY87643
QOY88067
QOY88866
QQS88170
QQS88496
QSR84237
QTD50766
QTD53836
QUV79767
QVL32579
QXB11428
UPH48256
'

nums = (1:length(TF))[TF == FALSE]
nums











# Try and get the MotBs for the above MotAs here:

xlsfn = "/GitHub/bioinfRhints/flag/get_MotBs/groupTax_1282_mafftConstr_2023-08-07_edit_wGeneOrder_BestMotBs.xlsx"

xls = openxlsx::read.xlsx(xlsfn)

rownums = match_grepl(patterns=MotA_seqIDs_to_get_to_379, x=xls$tipnames3_uniq, return_counts=FALSE)
rownums

xls$MotB_bestMatch_gid[rownums]

cbind(MotA_seqIDs_to_get_to_379, xls$MotB_bestMatch_gid[rownums])

# Adjacent proteins, putative MotBs
xls$acc0[rownums][is.na(xls$MotB_bestMatch_gid[rownums])]
xls$acc1[rownums][is.na(xls$MotB_bestMatch_gid[rownums])]


catn(xls$acc0[rownums][is.na(xls$MotB_bestMatch_gid[rownums])])

# QSR84237.1 is a MotA
# the corresponding MotB position is blank


seqids_to_get = xls$acc1[rownums][is.na(xls$MotB_bestMatch_gid[rownums])]
seqids_to_get[is.na(seqids_to_get)==FALSE]
seqids_to_get[is.na(seqids_to_get)]

seqids_to_download = seqids_to_get[is.na(seqids_to_get)==FALSE]

fasta_entries <- reutils::efetch(uid=seqids_to_download, db="protein", rettype = "fasta", retmode = "text")
write(reutils::content(fasta_entries, "text"), file="putative_MotBs_by_gene_order_to_blast.fasta")






# 
motB399_fn = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/motB399_2.fasta"
motB399 = as.character(read_FASTA_safe(motB399_fn, type="AA"))
names(motB399)


TF = Jiahe_MotA_seqids %in% Caroline_MotA_seqids
sum(TF)



fasta_fn = "motB399_2.fasta"
seqs = read_FASTA_safe(fasta_fn, type="AA")

seqs = as.character(seqs)
seqs = sapply(X=seqs, FUN=paste0, collapse="")
seqs

longest_fragments = NULL
for (i in 1:length(seqs))
	{
	words = strsplit(seqs[[i]], split="-")[[1]]
	words = words[words != ""]
	lengths = sapply(X=words, FUN=nchar)
	TF = lengths == max(lengths)
	longest_word = words[TF][1]
	longest_fragments[[i]] = longest_word
	}

names(longest_fragments) = names(seqs)
longest_fragments

#######################################################
# Now match to the aligned MotBs (BDRs)
#######################################################
unaligned_motBs_fasta_fn = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/Jiahe_2024-08_439_MotHomologs/439B/439_SubunitB_AAs_raw.fasta"

unaligned_motBs_w_MotA_labels = read_FASTA_safe(unaligned_motBs_fasta_fn, type="AA")

unaligned_motBs_w_MotA_labels = as.character(unaligned_motBs_w_MotA_labels)
unaligned_motBs_w_MotA_labels = sapply(X=unaligned_motBs_w_MotA_labels, FUN=paste0, collapse="")
unaligned_motBs_w_MotA_labels

MotA_names_for_MotBs = names(unaligned_motBs_w_MotA_labels)

# Go through the 399 MotBs with original names; link to the MotBs with MotA names
i = 1
matchnums_table = NULL
for (i in 1:length(longest_fragments))
	{
	pattern = toupper(longest_fragments[[i]])
	matchnums = match_grepl(pattern=toupper(longest_fragments[[i]]), x=unaligned_motBs_w_MotA_labels, return_counts=TRUE)
	matchnums_table = rbind(matchnums_table, matchnums)
	}
matchnums_table

