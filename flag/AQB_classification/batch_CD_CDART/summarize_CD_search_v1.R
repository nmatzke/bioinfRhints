#######################################################
# Search the concise hits of a NCBI CD (Conserved Domain) batch search
# 
# Identify input proteins that don't match any domains of interest
#
# (automated CDART, basically)
#######################################################

# Batch CDART
# Automated CDART
# Batch CDD
# Batch conserved domain
# Automated CDD
# Automated conserved domain

# NCBI Batch Web CD-Search Tool:
# https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
#
# Help page on the same:
# https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#SearchMethodBatchQuery

# 1. Get gids from fasta file:
#
cmds="
cd /GitHub/bioinfRhints/flag/AQB_classification/batch_CD_CDART/
grep -o -E "^>\w+" 1283_AQBs_noQs.fasta | tr -d ">"
grep -o '^>\w+' 1283_AQBs_noQs.fasta

# Extract the 5th element in a "|" delimited FASTA header, with ">" in front
awk 'BEGIN{FS="|"}{if(/^>/){print ">"$5}else{print $0}}' 1283_AQBs_noQs.fasta | head

# Extract the 5th element in a "|" delimited FASTA header, without ">" in front
awk 'BEGIN{FS="|"}{if(/^>/){print ""$5}else{print $0}}' 1283_AQBs_noQs.fasta | head

# Extract the 5th element in a "|" delimited FASTA header, without ">" in front; DON'T print sequence
awk 'BEGIN{FS="|"}{if(/^>/){print ""$5}}' 1283_AQBs_noQs.fasta | head

# Extract 5th item to gids file, output to text
cd /GitHub/bioinfRhints/flag/AQB_classification/batch_CD_CDART/
awk 'BEGIN{FS="|"}{if(/^>/){print ""$5}}' 1283_AQBs_noQs.fasta > 1283_AQBs_gids.txt
more 1283_AQBs_gids.txt

wc -l 1283_AQBs_gids.txt
# 1282 - 600 = 682
head -50 1283_AQBs_gids.txt > 1283_AQBs_gids_pt0.txt
head -600 1283_AQBs_gids.txt > 1283_AQBs_gids_pt1.txt
tail -682 1283_AQBs_gids.txt > 1283_AQBs_gids_pt2.txt

# DOESN'T WORK: Merge files at end (but this duplicates headers)
# cat 1283_AQBs_gids_pt1_hitdata.txt 1283_AQBs_gids_pt2_hitdata.txt > 1283_AQBs_gids_hitdata.txt 
wc -l 1283_AQBs_gids_hitdata.txt


wc -l 1283_AQBs_gids_pt2_hitdata.txt
# 698

tail -690 1283_AQBs_gids_pt2_hitdata.txt | more
tail -690 1283_AQBs_gids_pt2_hitdata.txt > 1283_AQBs_gids_pt2_hitdata_noHeader.txt

cat 1283_AQBs_gids_pt1_hitdata.txt 1283_AQBs_gids_pt2_hitdata_noHeader.txt > 1283_AQBs_gids_hitdataB.txt 
head 1283_AQBs_gids_hitdataB.txt 
wc -l 1283_AQBs_gids_hitdataB.txt 
# 1334
" # END cmds

# 
# 1. Put in a fasta file of proteins (1000 max, or gids as above)
#
# 1a. Extract 
#
# 2. For the results Download, click:
#    * Domain Hits, Data Mode: Concise, click "Superfamily Only"
#		 * Align details: BLAST Text
#    * Click Download
#
# 3. (Download Full versions etc. as you like for visualization)

library(gdata)				# for trim
library(BioGeoBEARS)	# for sourceall
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(stringr)			# for str_squish, 
library(openxlsx)			# for read.xlsx
sourceall("/GitHub/bioinfRhints/Rsrc/")

#wd = "/GitHub/bioinfRhints/flag/589_MotA_rename/CD-search/"
wd = "/GitHub/bioinfRhints/flag/AQB_classification/batch_CD_CDART/"
setwd(wd)

search_results_fn = "1283_AQBs_gids_pt1_hitdata.txt"
tdf1 = read_concise_CDbatch_results(search_results_fn)
search_results_fn = "1283_AQBs_gids_pt2_hitdata.txt"
tdf2 = read_concise_CDbatch_results(search_results_fn)
tdf = rbind(tdf1, tdf2)

search_results_fn = "1283_AQBs_gids_hitdata.txt"
write.table(tdf, file=search_results_fn, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
head(tdf)

# Most common hits
rev(sort(table(tdf$Short.name)))

# Proteins with more than one CD hit:
TF = table(tdf$Query) > 1
multidomain_hits = rev(sort(table(tdf$Query)[TF]))
multidomain_hits

# Assemble proteins with multiple domain hits
domain_hits_per_protein_df = multidomain_prots_to_single_line(tdf)
head(domain_hits_per_protein_df)

domain_hits_per_protein_df[domain_hits_per_protein_df$numdomains > 1,]


#######################################################
# Which proteins DO NOT have any domains of interest
#######################################################

domains_of_interest = c(
"MotA", 
"TolQ", 
"ExbB")

TF = get_domainsTF_in_CDbatch_results(domain_hits_per_protein_df, domains_of_interest=domains_of_interest)

proteins_not_matching_df = domain_hits_per_protein_df[TF==FALSE,]

dim(proteins_not_matching_df)


IDs_to_remove = get_first_words(proteins_not_matching_df$Query)
cat(IDs_to_remove, sep="\n")


#######################################################
# Check for interesting domains
#######################################################
homologous_hits_per_protein_df = domain_hits_per_protein_df[TF==TRUE,]
interesting = names(rev(sort(table(homologous_hits_per_protein_df$Short.name))))
interesting = interesting[-c(1:4,7)]
interesting

TF = domain_hits_per_protein_df$Short.name %in% interesting
interesting_names = domain_hits_per_protein_df$Query[TF]
interesting_gids = get_first_words(domain_hits_per_protein_df$Query[TF])
interesting_gids

# Proteins to remove as non-homologous
names_to_remove = proteins_not_matching_df$Query
gids_to_remove = get_first_words(proteins_not_matching_df$Query)
both = c(names_to_remove, interesting_names)
both_gids = c(gids_to_remove, interesting_gids)


# Load alignment to reorder
alnfn = "groupTax_1283_AQBs_mafftMiddleConstrained2.fasta"
aln = read.fasta(alnfn)
gids = unname(sapply(X=aln, FUN=attr, which="name"))
gids = gsub(pattern=">", replacement="", x=gids)
gids
fullnames = sapply(X=aln, FUN=attr, which="Annot")
fullnames = gsub(pattern=">", replacement="", x=fullnames)



interesting_gids_in_aln_nums = sort(find_gid_positions_in_nameslist(tmpgids=interesting_gids, nameslist=fullnames))

gids_to_remove_in_aln_nums = sort(find_gid_positions_in_nameslist(tmpgids=gids_to_remove, nameslist=fullnames))

fullnames[gids_to_remove_in_aln_nums]

nums_to_move = c(interesting_gids_in_aln_nums, gids_to_remove_in_aln_nums)
other_nums_TF = (1:length(fullnames)) %in% nums_to_move
nums_to_keep = (1:length(fullnames))[other_nums_TF == FALSE]
nums_to_keep

new_order_nums = c(interesting_gids_in_aln_nums, gids_to_remove_in_aln_nums, nums_to_keep)
length(new_order_nums)

aln2 = list()
fullnames2 = rep("", times=length(new_order_nums))
for (i in 1:length(new_order_nums))
	{
	num = new_order_nums[i]
	aln2[[i]] = aln[[num]]
	fullnames2[i] = fullnames[num]
	}

outfn = gsub(pattern=".fasta", replacement="_domainsOrder.fasta", x=alnfn)
write.fasta(aln2, names=fullnames2, file.out=outfn)



# Additional weird ones to cut;

# QDU31046.1 hypothetical protein ETAA8_61990 [Anatilimnocola aggregata]

# Positions 548-809 -- realign outside of this


#######################################################
# Read CD hit table, identify domains to cut
#######################################################
wd = "/GitHub/bioinfRhints/flag/AQB_classification/batch_CD_CDART/"
setwd(wd)

library(openxlsx)			# for read.xlsx

cdhits_fn = "1283_AQBs_gids_hitdata_v1.xlsx"
xlsfn = "/GitHub/bioinfRhints/flag/AQB_classification/groupTax_1282_mafftConstr_2023-08-07_edit.xlsx"

cdxls = read.xlsx(cdhits_fn)
xls = read.xlsx(xlsfn)

gids_in_xls = xls$tipnames3_uniq
head(rev(sort(table(gids_in_xls))))
hist(xls$len, breaks=50)
# A single-core domain protein can be up to 350 aas

numhits_by_gid = rev(sort(table(cdxls$Query)))
max(numhits_by_gid)
min(numhits_by_gid)
head(numhits_by_gid)

# 4 domains
cdxls[cdxls$Query == "QNQ08565.1", ]
# 3 domains
cdxls[cdxls$Query == "QJE97904.1", ]
cdxls[cdxls$Query == "ABE50669.1", ]
cdxls[cdxls$Query == "AAZ96310.1", ]
cdxls[cdxls$Query == "AAU91782.1", ]

# ...manually edit 3-4 domain hits entries, to be just 2 domains
cdhits_fn = "1283_AQBs_gids_hitdata_v2merge34.xlsx"
cdxls = read.xlsx(cdhits_fn)
numhits_by_gid = rev(sort(table(cdxls$Query)))
max(numhits_by_gid)
min(numhits_by_gid)
head(numhits_by_gid)

source("/GitHub/str2phy/Rsrc/str2phy_v1.R")
max_core_length = 225
proteins_over_this_will_be_split = max_core_length * 1.5
proteins_over_this_will_be_split
core_domains = c("MotA_ExbB superfamily", "TolQ superfamily", "MotA superfamily")
cdhits_df = process_CDhits(cdxls, xls, max_core_length=max_core_length, core_domains=core_domains)


head(cdhits_df)
tail(cdhits_df)
dim(cdhits_df)


cdhits_outfn = gsub(pattern=".xlsx", replacement="_processed.xlsx", x=cdhits_fn)
write.xlsx(x=cdhits_df, file=cdhits_outfn)
cmdtxt = paste0("open ", cdhits_outfn)
cmdtxt
system(cmdtxt)







