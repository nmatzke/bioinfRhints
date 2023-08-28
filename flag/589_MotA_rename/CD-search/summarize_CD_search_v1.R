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

# NCBI Batch Web CD-Search:
# https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
#
# Help page on the same:
# https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#SearchMethodBatchQuery

# 1. Put in a fasta file of proteins (1000 max)
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
sourceall("/GitHub/bioinfRhints/Rsrc/")

wd = "/GitHub/bioinfRhints/flag/589_MotA_rename/CD-search/"
setwd(wd)

search_results_fn = "589seqs_hitdata_concise.txt"
tdf = read_concise_CDbatch_results(search_results_fn)
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
alnfn = "motA_hmmalign589_full_tr2_order.fasta"
aln = read.fasta(alnfn)
gids = unname(sapply(X=aln, FUN=attr, which="name"))
gids = gsub(pattern=">", replacement="", x=gids)
gids
fullnames = sapply(X=aln, FUN=attr, which="Annot")
fullnames = gsub(pattern=">", replacement="", x=fullnames)



interesting_gids_in_aln_nums = sort(find_gid_positions_in_nameslist(tmpgids=interesting_gids, nameslist=fullnames))

gids_to_remove_in_aln_nums = sort(find_gid_positions_in_nameslist(tmpgids=gids_to_remove, nameslist=fullnames))

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

