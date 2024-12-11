#######################################################
# Taking BayesTraits outputs, and mapping on a phylogeny
# by Nick Matzke for Jamiema Philip
# flagellum presence/absence
#######################################################

library(ape)
library(BioGeoBEARS)

# Set working directory
wd = "/GitHub/bioinfRhints/ex/BayesTraits_on_phylo/"
setwd(wd)


###########################################################
# Load tree
###########################################################
trfn = "pruned_bacterial_rooted_tree.nex"
tr = ape::read.nexus(trfn)

# Basic tree info
numtips = length(tr$tip.label)
numnodes = tr$Nnode
total_numnodes = numtips + numnodes
numtips
numnodes
total_numnodes

# APE's node numbers (may not be the same as BayesTrait nodenums)
tipnums = 1:numtips
nodenums = (numtips+1):(numtips+numnodes)
###########################################################


# Load the MCMC samples (manually cut to the MCMC samples table)
# into a data.frame
###########################################################
mcmc_fn = "mcmc_samples.txt"
mcmc_df = read.table(file=mcmc_fn, header=TRUE, sep="\t", quote="", strip.white=TRUE)

head(mcmc_df)
###########################################################


###########################################################
# Plot the tree
pdffn = "pruned_bacterial_rooted_tree.pdf"
pdf(file=pdffn, width=12, height=18)

plot(tr)
title("Backbone tree of bacteria from [INSERT]")
add.scale.bar(x=0, y=numtips)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


# Plot the tree, with APE R package node numbers
pdffn = "pruned_bacterial_rooted_tree_wNodeNums.pdf"
pdf(file=pdffn, width=12, height=18)

plot(tr)
title("Backbone tree of bacteria from [INSERT],\nwith APE node numbers")
add.scale.bar(x=0, y=numtips)

tiplabels(text=tipnums, tip=tipnums)
nodelabels(text=nodenums, node=nodenums)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)
###########################################################

# Since the root node is 75 in both APE and BayesTraits, maybe the 
# node numbering is the same

# We'll plot the last MCMC sample
dim(mcmc_df)
# 90 153

last_sample = mcmc_df[90,]
last_sample

# The first 7 items are other info; then we have 
# 73 P(N) = probability (nonmotile)
# 73 P(M) = probability (motile)

# Get these into a table
probs = c(last_sample[1,8:length(last_sample)])
N = probs[seq(from=2, to=length(probs), by=2)]
M = probs[seq(from=1, to=length(probs), by=2)]
probs_matrix = cbind(N,M)
probs_df = as.data.frame(probs_matrix, stringsAsFactors=FALSE)
row.names(probs_df) = nodenums
cls.df(probs_df)
probs_df = unlist_df(probs_df)
cls.df(probs_df)

# Add the tip probabilities
tip_fn = "discrete_meta.txt"
tip_states = read.table(file=tip_fn, header=FALSE, sep="\t", quote="", strip.white=TRUE)
names(tip_states) = c("name", "state")
tip_states

tip_probs = matrix(data=0.0, nrow=nrow(tip_states), ncol=2)
for (i in 1:nrow(tip_states))
	{
	if (tip_states$state[i] == "N")
		{
		tip_probs[i,] = c(1.0, 0.0)
		}
	if (tip_states$state[i] == "M")
		{
		tip_probs[i,] = c(0.0, 1.0)
		}
	if (tip_states$state[i] == "-")
		{
		tip_probs[i,] = c(0.5, 0.5)
		}
	}

# Fix the row names
tip_probs_df = as.data.frame(tip_probs, stringsAsFactors=FALSE)
row.names(tip_probs_df) = tip_states$name
names(tip_probs_df) = c("N", "M")
tip_probs_df

# Put the tips in order of the tree tip labels
tiplabels = tr$tip.label
ordering = match(tiplabels, table=row.names(tip_probs_df))

# Double-check
head(tip_probs_df[ordering,])
head(tr$tip.label)

# Reorder
tip_probs_df = tip_probs_df[ordering,]

# Merge
all_probs_df = rbind(tip_probs_df, probs_df)
head(all_probs_df)
dim(all_probs_df)



# Plot onto tree
colors = c("green3", "blue")

###########################################################
# Plot the tree
pdffn = "pruned_bacterial_rooted_tree_w1_mcmc.pdf"
pdf(file=pdffn, width=7, height=10)

plot(tr, label.offset=0.03, cex=0.7)
title("Backbone tree of bacteria from [INSERT]\nwith 1 MCMC's ancestral states")
add.scale.bar(x=0, y=numtips)

tiplabels(tip=tipnums, pie=tip_probs_df, piecol=colors, cex=0.5)
nodelabels(node=nodenums, pie=probs_df, piecol=colors, cex=0.5)


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)



