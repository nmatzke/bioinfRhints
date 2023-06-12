#######################################################
# Scripts: manipulating alignments, sequence names / tip names, etc.
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches

sourceall("/GitHub/bioinfRhints/Rsrc/") 
sourceall("/GitHub/bioinfRhints/R/BEASTmasteR/") 
# for remove_equals_from_tips()
# for read.beast.table_original


#wd = "/GitHub/bioinfRhints/flag/530_MotA_rename/"
wd = "/GitHub/bioinfRhints/flag/589_MotA_rename_MiddleConstrained/"
setwd(wd)

#alnfn = "motA_alignment_refined.fasta"
alnfn = "589_seqs_mafftMiddleConstrained2.fasta"
aln = read.fasta(alnfn)

#trfn = "530_sequences_Alignment_contree_reRootLadder_gIDs.newick"
#nexfn = "530_sequences_Alignment_contree_reRootLadder_gIDs.nexus"
#trfn = "motA589_IQtree_midpoint_FigTree.newick"
#nexfn = "motA589_IQtree_midpoint_FigTree.nexus"
trfn = "589_seqs_mafftMiddleConstrained2.fasta.contree"

tr = read.tree(trfn)
tipnames = tr$tip.label
#tr2 = read.nexus(nexfn)



# "aln" is an R list, each element is a sequence, with a name and some other attributes:
attributes(aln[[1]])


# Get the full names of all the sequences
# (lapply = "list apply" = apply the function "attr" to each element in the list "aln")
# (unlist turns the list of names into a vector)
fullnames = unlist(lapply(X=aln, FUN=attr, which="Annot"))
fullnames = gsub(pattern=">", replacement="", x=fullnames)
gid = unlist(lapply(X=aln, FUN=attr, which="name")) # gid = Genbank IDs


TF1 = gid %in% tipnames
TF2 = tipnames %in% gid

matches = match(tipnames, table=gid)
gid[matches][1:5]
tipnames[1:5]
# Match!

# Reorder alignment
tmplist = list()
for (i in 1:length(aln))
	{
	tmplist[[i]] = aln[[matches[i]]]
	}

alnfn_out = gsub(pattern=".fasta", replacement="_tipOrder.fasta", x=alnfn)
write.fasta(sequences=tmplist, names=gid[matches], file.out=alnfn_out)








# Extract the species (between brackets)
# Examples:
tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
regmatches(tmptxt, gregexpr("\\[.+?\\]", tmptxt))

tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium] [hello]"
regmatches(tmptxt, gregexpr("\\[.+?\\]", tmptxt))

# Function to extract last item in [square brackets like this] in a string
# ...repeat for each item in a list
list_of_strings = fullnames; replace_spaces=TRUE

# Run the function
species_names = extract_last_brackets(list_of_strings=fullnames)
species_names

seqlengths = get_seqlengths(aln)
hist(seqlengths, breaks=50)




#######################################################
# Parse the names for protein name info
#######################################################

# Test
protein_name = extract_protein_name_info(">>ACF14374.1 MotA/TolQ/ExbB proton channel [Chloroherpeton thalassium ATCC 35110]")
short_protname = classify_MotAfam_labels(list_of_strings=protein_name)
short_protname

# Run on everything
fullnames = unlist(lapply(X=aln, FUN=attr, which="Annot"))
fullnames = gsub(pattern=">", replacement="", x=fullnames)

protein_name = extract_protein_name_info(fullnames)
short_protname = classify_MotAfam_labels(list_of_strings=protein_name)
sort(table(protein_name))




list_of_strings = protein_name
sum(grepl(pattern="flagellar stator protein MotA", x=list_of_strings, ignore.case=TRUE))
protein_name = extract_protein_name_info(fullnames)

short_protname = classify_MotAfam_labels(list_of_strings=protein_name)
short_protname

sort(table(short_protname))

# Remove "=" from tipnames before renaming
species_names_wSpaces = gsub(pattern="=", replacement="EQ", x=species_names_wSpaces) # remove "="
species_names = gsub(pattern="=", replacement="EQ", x=species_names) # remove "="


newnames1_df = cbind(species_names, short_protname, gid, seqlengths)
newnames2_df = cbind(short_protname, species_names, gid, seqlengths)
newnames3_df = cbind(seqlengths, short_protname, species_names, gid)
newnames1 = apply(X=newnames1_df, MARGIN=1, paste0, collapse="|")
newnames2 = apply(X=newnames2_df, MARGIN=1, paste0, collapse="|")
newnames3 = apply(X=newnames3_df, MARGIN=1, paste0, collapse="|")


nexfn2 = gsub(pattern=".nexus", replacement="_noEQ.nexus", x=nexfn)
remove_equals_from_tips(nexfn=nexfn, outfn=nexfn2, format="raxml")


# Write the new names into the NEXUS file; this keeps the bootstraps etc.
outfn = gsub(pattern="_noEQ.nexus", replacement="_newNames1.nexus", x=nexfn2)
tmpstrs = readLines(nexfn)

for (i in 1:length(tmpstrs))
	{
	tmpstr = tmpstrs[i]
#	tmpstr = gsub(pattern="'", replacement="", x=tmpstr)
	for (j in 1:length(gid))
		{
		tmpstr = gsub(pattern=gid[j], replacement=newnames1[j], x=tmpstr)
		}
	tmpstrs[i] = tmpstr
	}
writeLines(text=tmpstrs, con=outfn)
outfn1 = gsub(pattern=".nexus", replacement="_nameFirst.nexus", x=nexfn)
writeLines(text=tmpstrs, con=outfn1)


# Write the new names into the NEXUS file; this keeps the bootstraps etc.
outfn = gsub(pattern=".nexus", replacement="_newNames2.nexus", x=nexfn)
tmpstrs = readLines(nexfn)

for (i in 1:length(tmpstrs))
	{
	tmpstr = tmpstrs[i]
	#tmpstr = gsub(pattern="'", replacement="", x=tmpstr)
	
	for (j in 1:length(gid))
		{
		tmpstr = gsub(pattern=gid[j], replacement=newnames2[j], x=tmpstr)
		}
	tmpstrs[i] = tmpstr
	}
writeLines(text=tmpstrs, con=outfn)
outfn2 = gsub(pattern=".nexus", replacement="_protFirst.nexus", x=nexfn)
writeLines(text=tmpstrs, con=outfn2)



# Write the new names into the NEXUS file; this keeps the bootstraps etc.
outfn = gsub(pattern=".nexus", replacement="_newNames3.nexus", x=nexfn)
tmpstrs = readLines(nexfn)

for (i in 1:length(tmpstrs))
	{
	tmpstr = tmpstrs[i]
	#tmpstr = gsub(pattern="'", replacement="", x=tmpstr)
	for (j in 1:length(gid))
		{
		tmpstr = gsub(pattern=gid[j], replacement=newnames3[j], x=tmpstr)
		}
	tmpstrs[i] = tmpstr
	}
writeLines(text=tmpstrs, con=outfn)
outfn3 = gsub(pattern=".nexus", replacement="_lengthFirst.nexus", x=nexfn)
writeLines(text=tmpstrs, con=outfn3)


sort(table(short_protname))






#######################################################
# Coloring in the branches by size
#######################################################

# Classify the sizes
size_classes = rep(0, times=length(seqlengths))
categories = c(0, 100, 200, 300, 400, 500, 600, 700, 800, (max(seqlengths)+1))
charval = c(0, 1, 2, 3, 4, 5, 6, 7, 8)

for (i in 1:(length(categories)-1))
	{
	minval = categories[i]
	maxval = categories[i+1]
	
	TF1 = seqlengths >= minval
	TF2 = seqlengths < maxval
	TF = (TF1 + TF2) == 2
	size_classes[TF] = charval[i]
	}

hist(size_classes)

names(size_classes) = newnames3

library(ape)
library(phytools)

# Ladderize will flip up/down, so it looks like FigTree

tr2 = phytools::readNexus(outfn3, format="standard")
checktree(tr2)


if (is.binary(tr2) == FALSE)
	{
	tr2 = multi2di(tr2)
	
	brlens = tr2$edge.length
	TF1 = is.na(brlens)
	TF2 = is.nan(brlens)
	brlens[TF1] = 0.0
	brlens[TF2] = 0.0
	TF3 = brlens <= 0.0
	brlens[TF3] = min(brlens > 0.0)
	tr2$edge.length = brlens
	}


if (has.singles(tr2) == TRUE)
	{
	tr2 = collapse.singles(tr2)
	brlens = tr2$edge.length
	TF1 = is.na(brlens)
	TF2 = is.nan(brlens)
	brlens[TF1] = 0.0
	brlens[TF2] = 0.0
	TF3 = brlens <= 0.0
	brlens[TF3] = min(brlens > 0.0)
	tr2$edge.length = brlens
	}
checktree(tr2)



TF = sort(tr2$tip.label) == sort(newnames3)

namevals = cbind(sort(tr2$tip.label), sort(newnames3))
namevals[TF==FALSE,]
sum(TF == FALSE)



#stats_table = read.beast.table_original(file=outfn3, digits=4) 


tr2 = ladderize(tr2, right=FALSE)
tr2$tip.label = gsub(pattern="'", replacement="", x=tr2$tip.label)
#plot(tr2, show.tip.label=FALSE)


tip_rownum_in_sorted = rep(0, times=length(size_classes))
char_rownum_in_charnames = rep(0, times=length(size_classes))
for (i in 1:nrow(namevals))
	{
	tip_rownum_in_sorted[i] = match(namevals[i,1], table=tr2$tip.label)
	char_rownum_in_charnames[i] = match(namevals[i,2], table=newnames3)
	}
orders = cbind(tip_rownum_in_sorted, char_rownum_in_charnames)
tiporder = order(orders[,1])
char_into_tip_order = orders[tiporder,2]

size_classes = size_classes[char_into_tip_order]

sort(names(size_classes))[1:5]
sort(tr2$tip.label)[1:5]

# Orders match?
names(size_classes)[1:5]
tr2$tip.label[1:5]



# Create an ordered character
params_matrix = matrix(data=0, nrow=length(unique(size_classes)), ncol=length(unique(size_classes)))
for (i in 1:nrow(params_matrix))
	{
	for (j in 1:ncol(params_matrix))
		{
		if (j == i+1)
			{
			params_matrix[i,j] = 1
			}
		if (i == j+1)
			{
			params_matrix[i,j] = 1
			}
		}
	}
params_matrix	

res = ace(x=unname(size_classes), phy=tr2, type="discrete", model=params_matrix)


#pdffn = "530_sequences_Alignment_contree_reRootLadder_gIDs_protFirst_sizeColors.pdf"
pdffn = gsub(pattern=".nexus", replacement="_sizeColors.pdf", x=outfn2)
pdf(file=pdffn, width=12, height=48)

titletxt = "Colored plot of sequence lengths\n(red=low, yellow=high"
plot.phylo(tr2, type="phylogram", label.offset=0.0, cex=0.45, font=1, align.tip.label=TRUE)#, tip.color="white")
title(titletxt)

colors <- heat.colors(n=length(unique(size_classes)))
tiplabels(pch=22, bg=colors[as.numeric(size_classes)], cex=2, adj=0.5)
res$lik.anc[res$lik.anc < 0] = 0.0
nodelabels(pie=res$lik.anc, piecol=colors, cex=0.4)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)







#######################################################
# Order alignment by size...
#######################################################
alnfn = "motA_hmmalign589_full.fasta"
aln = read.fasta(alnfn)

# "aln" is an R list, each element is a sequence, with a name and some other attributes:
attributes(aln[[1]])

# Get the full names of all the sequences
# (lapply = "list apply" = apply the function "attr" to each element in the list "aln")
# (unlist turns the list of names into a vector)
fullnames = unlist(lapply(X=aln, FUN=attr, which="Annot"))
gids = unlist(lapply(X=aln, FUN=attr, which="name")) # gid = Genbank IDs

# Fix GIDs with "_"
underscores_in_gids_TF = grepl(pattern="_", x=gids)
gids_w_underscores = gids[underscores_in_gids_TF]
gids_wo_underscores = gsub(pattern="_", replacement="", x=gids_w_underscores)

tipnames = tr2$tip.label
tipnames = gsub(pattern="'", replacement="", x=tipnames)
tipnames

# Fix tipnames with "_" in gid:
for (j in 1:length(gids_w_underscores))
	{
	matches_alignment_TF = grepl(pattern=gids_w_underscores[j], x=gids)
	gids[matches_alignment_TF] = gsub(pattern=gids_w_underscores[j], replacement=gids_wo_underscores[j], x=gids[matches_alignment_TF])
	
	matches_tips_TF = grepl(pattern=gids_w_underscores[j], x=tipnames)
	tipnames[matches_tips_TF] = gsub(pattern=gids_w_underscores[j], replacement=gids_wo_underscores[j], x=tipnames[matches_tips_TF])
	}


tipnames2 = rep("", times=length(tipnames))
for (i in 1:length(tipnames))
	{
	words = strsplit(tipnames[i], split="_")[[1]]
	trgid = words[1]
	tipnames2[i] = trgid
	}

tmporder1 = match(gids, table=tipnames2)
tmporder1

# Are there any unmatched?
cat(unname(gids[is.na(tmporder1)]), sep="\r")


tmporder2 = match(tipnames2, table=gids)
tmporder2

gids[is.na(tmporder1)]
tipnames2[is.na(tmporder2)]


alnfn_out = gsub(pattern=".fasta", replacement="_tr2_order.fasta", x=alnfn)
aln2 = aln[tmporder2]
names2 = fullnames[tmporder2]

write.fasta(sequences=aln2, names=names2, file.out=alnfn_out)

