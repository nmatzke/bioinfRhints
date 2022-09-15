
#######################################################
# Example: converting to/from some trees in R
#######################################################


#######################################################
# Getting some standard Newick/NEXUS trees for reference
#######################################################
library(ape)	# for read/write NEXUS
library(BioGeoBEARS)	# for list2str, moref

# Directory with some example trees
wd = "/drives/Dropbox/_njm/__packages/TNTR_setup/example_data/"
setwd(wd)

trfn = "Psychotria_5.2.newick"
tr = read.tree(trfn)
plot(tr)

nexfn_translateTRUE = "Psychotria_5.2_translateTRUE.nex"
write.nexus(tr, file=nexfn_translateTRUE, translate=TRUE)
nexfn = "Psychotria_5.2.nex"
write.nexus(tr, file=nexfn, translate=FALSE)

moref(nexfn_translateTRUE)
moref(nexfn)





# This file has a problem
library(ape)

nexfn = "aquickiewStats_notags3.nex"
tr = read.nexus(file=nexfn)
# Oops, error!
# Error in start:end : argument of length 0
# 
# This broke due to these:
# begin ;
# ...and...
# end ;


# Fix it:
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')

# begin ; -> begin;
# end ; -> end;
infn = "aquickiewStats_notags3.nex"
outfn = "aquickiewStats_notags3_fixed.nex"

# Run UNIX 'sed' command to fix
rsed(oldtxt=" ;", newtxt=";", infn=infn, outfn=outfn, print_warnings=TRUE)
 
tr = read.nexus(file=outfn)


tr = read_TNT_nexus(file=infn, correction="sed_nofix") 
plot(tr)


tr = read_TNT_nexus(file=infn, correction="gsub_nofix") 
plot(tr)



#######################################################
# OK, try reading multiple TNT-NEXUS trees, with tags
#######################################################
# Source TNTR package
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')

# Directory with some example trees
wd = "/drives/Dropbox/_njm/__packages/TNTR_setup/example_data/"
setwd(wd)

# TNT tree filename
nexfn = "aquickiewStats.tre"

# read_TNT_nexus() deals with the " ;" problem, and
# the trailing slashes if there are TNT labels
trs = read_TNT_nexus(file=nexfn, correction="sed_nofix")
tr = trs[[1]]

# Update prt to deal with parsimony branchlengths
source('/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R')
trtable = prt(tr)
trtable$daughter_nds


#######################################################
# Read a tree in TNT format -- no branchlengths
#######################################################
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')

# Read a tnt tree with just topology
tntfn = "input_chars_best_trees.tnttrees"
tr = tntfile2R(tntfn=tntfn, brlens=FALSE)
plot(tr)
trtable = prt(tr, printflag=FALSE)

# Read a tree with branchlengths
# (if brlens=FALSE, branchlengths will be read as node labels)
tntfn = "input_chars_best_tree_Oth_brlen.tnttree"
tr = tntfile2R(tntfn=tntfn, brlens=TRUE)
plot(tr)
trtable = prt(tr, printflag=FALSE)
trtable


# Read a tree with Bootstrap supports
# (if brlens=FALSE, branchlengths will be read as node labels)
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
tntfn = "single_tree_w_Bootstraps.tnttree"
tr = tntfile2R(tntfn=tntfn, brlens=FALSE, branchlabels="=")
plot(tr)
trtable = prt(tr, printflag=FALSE)
trtable


# Read a tree with branchlengths, to the labels
# (if brlens=FALSE, branchlengths will be read as node labels)
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
tntfn = "input_chars_best_tree_Oth_brlen.tnttree"
tr = tntfile2R(tntfn=tntfn, brlens=FALSE, branchlabels="=")
plot(tr)
trtable = prt(tr, printflag=FALSE)
trtable

# Get the labels for the branchlengths below the tips
tr_out = tr
cat(tr_out$tip.label, sep="\n")
sep_tipnames_branchlabels(tip_labels=tr_out$tip.label, branchlabels="=")





# Write an APE phylo tree to Newick string, then to tnt string
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
tntfn = "single_tree_w_Bootstraps.tnttree"
tr = tntfile2R(tntfn=tntfn, brlens=FALSE, branchlabels="=")
tr

# Write to TNT string
tr$edge.length = NULL
newickstr = write.tree(tr, file="")
tntstr = newickstr2tntstr(newickstr)

# Write TNT string to file
outfn = "tree.tnt"
tntstr2file(tntstr, file=outfn)
moref(outfn)

# Check the reading of the file
tr = tntfile2R(tntfn=outfn, brlens=FALSE, branchlabels="=")
tr


# Write an APE phylo object to TNT file
# single tree
outfn = "tree.tnt"
write_tree_TNT(tr, file=outfn)
moref(outfn)

# Read it in tnt
tnt_cmds = '
cd /drives/Dropbox/_njm/__packages/TNTR_setup/example_data/
tnt;
proc tree_data.tnt;
proc tree.tnt;
'


# Write an APE phylo object to TNT file
# multiple trees
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
tntfn = "treefile"
trs = tntfile2R(tntfn=tntfn, brlens=FALSE, branchlabels="=")
outfn = "trees.tnt"
write_tree_TNT(trs, file=outfn)
moref(outfn)

# Read them in tnt
tnt_cmds = '
cd /drives/Dropbox/_njm/__packages/TNTR_setup/example_data/
tnt;
proc treefile_data.tnt;
proc trees.tnt;
'


# Read multiple TNT trees
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
tntfn = "treefile"
trs = tntfile2R(tntfn=tntfn, brlens=FALSE, branchlabels="=")
class(trs)


# Read multiple trees, by first converting to Newick
source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
tntfn = "treefile"
trfn = tntfile2newick(tntfn=tntfn, brlens=FALSE)
moref(trfn)

trs = read.tree(trfn)
par(mfrow=c(2,2))
for (i in 1:4)
	{
	plot(trs[[i]])
	}

# Read multiple trees, by first converting to NEXUS
tntfn = "treefile"
nexfn = tntfile2newick(tntfn=tntfn, brlens=FALSE, options="nexus")
moref(nexfn)

trs = read.nexus(nexfn)
par(mfrow=c(2,2))
for (i in 1:4)
	{
	plot(trs[[i]])
	}


#######################################################
# Read the results of "auto.run"
#######################################################
library(gdata)	# for trim

source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
source('/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R')

wd = "/drives/GDrive/__GDrive_projects/2015-08-01_AFA_cladogram/_02_TNT/2015-08-25_runs/allchars_MP/"
setwd(wd)

# Working directory
settings = NULL
settings$wd = "."
# Filenames
settings$script_fn = "auto.run"
settings$MP_topologies_fn = "MP_topologies.tnt"
settings$indata_fn = "tmpdata.tnt"
#settings$intrees_fn = "tmptrees.tnt"
settings$nodenums_fn = "auto_nodenums.tnt"
settings$synapos_fn = "auto_synapos.tnt"
settings$branchlengths_fn = "auto_branchlengths.tnt"
settings$strictcon_fn = "auto_strict_consensus.tnt"
settings$bremer_abs_fn = "auto_BremerAbsolute.tnt"
settings$bremer_rel_fn = "auto_BremerRelative.tnt"
settings$bootstraps_fn = "auto_Bootstraps.tnt"

# Tree with branch lengths (trbl)
trbl = tntfile2R(settings$branchlengths_fn, brlens=TRUE)

# Tree table
trtable = read_auto_results(settings=settings)
head(trtable)

# Extract synapomorphies
logfn = "auto_logfile.txt"
synapos = extract_synapos(logfn, trtable=trtable)

head(synapos)
tail(synapos)



# Extract stats
logfn = "auto_logfile.txt"

# Get tree-wide stats
treestats = autoget_treestats(logfn)
treestats

# Get per-character stats
charstats = autoget_charstats(logfn)
charstats


# Compare TNT vs. manual calculation
# (NOTE: all of these numbers assume that autapomorphies were left in; 
#  excluding them will change the statistics somewhat)
sum(cscores)
# 1263
treestats$TL
# 1263

CIsum = calc_CI(cscores=sum(charstats$cscores), minsteps=sum(charstats$minsteps))
CIsum
# 0.4196358
treestats$CI
# 0.42

RIsum = calc_RI(cscores=sum(charstats$cscores), minsteps=sum(charstats$minsteps), maxsteps=sum(charstats$maxsteps))
RIsum
# 0.7864219
treestats$RI
# 0.786

calc_RI(cscores=1263, minsteps=530, maxsteps=3962)
# 0.7864219
(3962 - 1263) / (3962 - 530)
# 0.7864219

RCIsum = calc_RCI(cscores=sum(charstats$cscores), minsteps=sum(charstats$minsteps), maxsteps=sum(charstats$maxsteps))
RCIsum
# 0.3300108
treestats$RCI
# 0.3300108

write.table(trtable[,-8], file="trtable.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(synapos, file="synapos.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(treestats, file="treestats.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(charstats, file="charstats.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



#######################################################
# Make some nice plots
#######################################################
ltr = ladderize(trbl, right=FALSE)
ltr_table = prt(ltr, printflag=FALSE, get_tipnames=TRUE)
reorder_trtable = match(trtable$tipnames, table=ltr_table$tipnames)
reorder_trtable

# Stick in TNT columns, ordering to match ladderized tree
cols_to_get = c("TNT_nodenums", "absBremer", "relBremer", "freqBoot", "gcBoot", "TNT_synapos")
col_positions = match(cols_to_get, table=names(trtable))
col_positions = col_positions[!is.na(col_positions)]
ltr_table = cbind(ltr_table, trtable[reorder_trtable,][,col_positions])
ltr_table

# Reorder the columns so it looks nicer
cols_to_get = c("x", "node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label", "TNT_nodenums", "absBremer", "relBremer", "freqBoot", "gcBoot", "TNT_synapos", "tipnames")
col_positions = match(cols_to_get, table=names(ltr_table))
col_positions = col_positions[!is.na(col_positions)]
ltr_table = ltr_table[, col_positions]
head(ltr_table)


#######################################################
# Make the PDF
#######################################################
pdffn = "auto_trplots_v1.pdf"
pdf(pdffn, width=9, height=12)

# Settings
tip_hadj = -0.1
tipcex = 0.7
edgecex = 0.7
scalebar_cex = 0.7
titletxt = "MP tree"

ntips = length(ltr$tip.label)


# Branch lengths
subtitle_txt = "branch lengths"
plot(ltr, show.tip.label=FALSE)
add.scale.bar(x=10, y=round(0.9*ntips), length=5, cex=scalebar_cex)
txt = paste0(titletxt, "\n(", subtitle_txt, ")")
title(txt)

# Add tip labels
tiplabels(text=ltr_table$label[1:ntips], tip=1:ntips, adj=tip_hadj, cex=tipcex, frame="none", bg="none", font=1)

# Edge labels
PPtxt = ltr_table$edge.length
PPtxt[PPtxt==0] = ""
edges = ltr_table$parent_br
edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=edgecex)



# Number of unambiguous synapomorphies
subtitle_txt = "# synapomorphies"
plot(ltr, show.tip.label=FALSE)
add.scale.bar(x=10, y=round(0.9*ntips), length=5, cex=scalebar_cex)
txt = paste0(titletxt, "\n(", subtitle_txt, ")")
title(txt)

# Add tip labels
tiplabels(text=ltr_table$label[1:ntips], tip=1:ntips, adj=tip_hadj, cex=tipcex, frame="none", bg="none", font=1)

# Edge labels
TNT_synapos = ltr_table$TNT_synapos
matches = gregexpr(pattern="\\|", text=TNT_synapos) 
num_synapos = sapply(X=matches, FUN=length)+1
num_synapos[ sapply(X=matches, FUN=isNegOne) ] = 1
num_synapos[TNT_synapos==""] = 0
num_synapos

PPtxt = num_synapos
PPtxt[PPtxt==0] = ""
edges = ltr_table$parent_br
edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=edgecex)


# Absolute Bremer support (decay index)
subtitle_txt = "absolute Bremer support/decay index"
plot(ltr, show.tip.label=FALSE)
add.scale.bar(x=10, y=round(0.9*ntips), length=5, cex=scalebar_cex)
txt = paste0(titletxt, "\n(", subtitle_txt, ")")
title(txt)

# Add tip labels
tiplabels(text=ltr_table$label[1:ntips], tip=1:ntips, adj=tip_hadj, cex=tipcex, frame="none", bg="none", font=1)

# Edge labels
PPtxt = ltr_table$absBremer
PPtxt[PPtxt==0] = ""
edges = ltr_table$parent_br
edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=edgecex)



# Relative Bremer support (decay index)
subtitle_txt = "relative Bremer support/decay index"
plot(ltr, show.tip.label=FALSE)
add.scale.bar(x=10, y=round(0.9*ntips), length=5, cex=scalebar_cex)
txt = paste0(titletxt, "\n(", subtitle_txt, ")")
title(txt)

# Add tip labels
tiplabels(text=ltr_table$label[1:ntips], tip=1:ntips, adj=tip_hadj, cex=tipcex, frame="none", bg="none", font=1)

# Edge labels
PPtxt = ltr_table$relBremer
PPtxt[PPtxt==0] = ""
edges = ltr_table$parent_br
edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=edgecex)




# Bootstrap frequency
subtitle_txt = "bootstrap support/frequency"
plot(ltr, show.tip.label=FALSE)
add.scale.bar(x=10, y=round(0.9*ntips), length=5, cex=scalebar_cex)
txt = paste0(titletxt, "\n(", subtitle_txt, ")")
title(txt)

# Add tip labels
tiplabels(text=ltr_table$label[1:ntips], tip=1:ntips, adj=tip_hadj, cex=tipcex, frame="none", bg="none", font=1)

# Edge labels
PPtxt = ltr_table$freqBoot
PPtxt[PPtxt==0] = ""
edges = ltr_table$parent_br
edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=edgecex)



# Bootstrap gc
subtitle_txt = "bootstrap support/gc"
plot(ltr, show.tip.label=FALSE)
add.scale.bar(x=10, y=round(0.9*ntips), length=5, cex=scalebar_cex)
txt = paste0(titletxt, "\n(", subtitle_txt, ")")
title(txt)

# Add tip labels
tiplabels(text=ltr_table$label[1:ntips], tip=1:ntips, adj=tip_hadj, cex=tipcex, frame="none", bg="none", font=1)

# Edge labels
PPtxt = ltr_table$gcBoot
PPtxt[PPtxt==0] = ""
edges = ltr_table$parent_br
edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=edgecex)




dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


