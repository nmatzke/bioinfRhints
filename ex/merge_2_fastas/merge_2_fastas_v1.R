library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches

# Script at:
# https://github.com/nmatzke/bioinfRhints/tree/main/ex/merge_2_fastas



wd = "/GitHub/bioinfRhints/ex/merge_2_fastas/"
setwd(wd)

fn1 = "motA_alignment_refined.fasta"
fn2 = "motA_jhmmer11052023.fasta"

seqs1 = read.fasta(fn1)
seqs2 = read.fasta(fn2)

names1 = attr(seqs1, "name")
names2 = attr(seqs2, "name")

# Remove the /100-132
names2_subset = names2
for (i in 1:length(names2_subset))
	{
	tmpwords = strsplit(names2_subset[i], split="/")[[1]]
	names2_subset[i] = tmpwords[1]
	}


allnames = c(names1, names2_subset)

length(names1)
length(names2_subset)
length(allnames)

length(unique(allnames))
uniq_allnames = unique(allnames)

merged_seqs = list()
j = 0
h = 0
jmmer_seqs = list()
jmmer_seqnames = NULL
jmmer_seqIDs = NULL
for (i in 1:length(uniq_allnames))
	{
	cat(i, sep=",")
	
	TF = names1 %in% uniq_allnames[i]
	nums = (1:length(names1))[TF]
	num = nums[1]
	
	if (sum(TF) > 0)
		{
		tmpseq = seqs1[[num]]
		} else {
		TF = names2_subset %in% uniq_allnames[i]
		nums = (1:length(names2_subset))[TF]
		num = nums[1]
		tmpseq = seqs2[[num]]
		jmmer_seqs[[(h=h+1)]] = seqs2[[num]]
		jmmer_seqnames = c(jmmer_seqnames, attr(seqs2[[num]], "Annot"))
		jmmer_seqIDs = c(jmmer_seqIDs, names2_subset[num])
		}
	
	merged_seqs[[(j=j+1)]] = tmpseq	
	} # END for (i in 1:length(uniq_allnames))

merged_seqs[[1]]
merged_seqs[[2]]
merged_seqs[[572]]
merged_seqs[[592]]
jmmer_seqs


attr(jmmer_seqs, "name")
sapply(X=jmmer_seqs, FUN=attr, "name")


# Write out the sequences
merged_seqs_names = sapply(X=merged_seqs, FUN=attr, "Annot")

# remove ">"
merged_seqs_names = gsub(pattern=">", replacement="", x=merged_seqs_names)



# WRITE THE MERGED FILE TO FASTA
write.fasta(merged_seqs, names=merged_seqs_names, file.out="592seqs_merged.fasta")


# Seqs to get originals of:
jmmer_seqs_names = sapply(X=jmmer_seqs, FUN=attr, "Annot")
cat(jmmer_seqs_names, sep="\n")

cat(jmmer_seqIDs, sep="\n")




# jmmer_seqIDs that were added:
txt='
QEI18926.1
WP_002721884.1
RLG45601.1
RMG20999.1
DGBBEKCF_02484
QEI19744.1
CAG37394.1
ABD81419.1
MBU2985819.1
ABD82476.1
MBU2985330.1
SCA58681.1
WP_015920816.1
ABD81321.1
MBU2985725.1
QEI19442.1
UEO04809.1
UEO06636.1
ABD81792.1
MBU2984437.1
QEI18804.1
ABD82818.1
MBU2986144.1
QEI20720.1
UEO04322.1
WP_002720852.1
ABD82214.1
MBU2985168.1
UEO05991.1
UEO03725.1
QEI19895.1
CAG37138.1
DGBBEKCF_01299
UEO07578.1
UEO05650.1
QEI18324.1
CAG37703.1
WP_015921325.1
MBU2986356.1
DGBBEKCF_01621
ABD79617.1
DGBBEKCF_00127
DGBBEKCF_00463
QEI19039.1
UEO02151.1
QEI18407.1
DGBBEKCF_00016
PKPEBJJI_00297
ABD82072.1
MBU2984729.1
CAG37702.1
PKPEBJJI_00564
QEI18325.1
ABD79618.1
MBU2986357.1
MBU2984728.1
ABD82073.1
DGBBEKCF_01764
QEI18408.1
DGBBEKCF_01496
MBU2987451.1
PKPEBJJI_00344
'




