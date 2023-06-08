
#######################################################
# Read the results of "auto.run"
#######################################################


#######################################################
# Getting some standard Newick/NEXUS trees for reference
#######################################################
library(ape)	# for read/write NEXUS
library(BioGeoBEARS)	# for list2str, moref
library(gdata)	# for trim
library(phytools) # for midpoint.root
library(seqinr)				# for read.fasta

source("/GitHub/bioinfRhints/Rsrc/protein_bioinf_v1.R")
source("~/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R")
#source('/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R')

#wd = "/drives/GDrive/__GDrive_projects/2016-06-16_venerid3/02_TNT/"
wd = "/Users/nmat471/HD/GitHub/bioinfRhints/flag/589_MotA_rename/CD-search/mafft_constrained_3/02_tnt/"
setwd(wd)


# Source of good names:
#alnfn = "motA_hmmalign589_full_tr2_orderMiddle.fasta"
alnfn = "outfile.out_nogaps_rename.fasta"
aln = read.fasta(alnfn)
fullnames = fullnames_from_readFasta(aln)
head(fullnames)
TF = grepl(pattern="\\(", x=fullnames)
sum(TF)


# Working directory
settings = NULL
settings$wd = "."
# Filenames
settings$script_fn = "auto.run"
settings$MP_topologies_fn = "MP_topologies.tnt"
settings$indata_fn = "tmpdata.tnt"
#settings$intrees_fn = "tmptrees.tnt"
settings$nodenums_fn = "auto_nodenums.tnt"
settings$synapos_fn = "auto_synapos2.tnt"
settings$branchlengths_fn = "auto_branchlengths.tnt"
settings$strictcon_fn = "auto_strict_consensus.tnt"
settings$bremer_abs_fn = "auto_BremerAbsolute.tnt"
settings$bremer_rel_fn = "auto_BremerRelative.tnt"
settings$bootstraps_fn = "auto_Bootstraps.tnt"

# Tree with branch lengths (trbl)
trbl = tntfile2R(settings$branchlengths_fn, brlens=TRUE)
trbl$root.edge = NULL
ltr = ladderize(trbl, right=FALSE)
plot(ltr, show.tip.label=FALSE)

pdffn = "midpoint_tree_wBrlens.pdf"
pdf(file=pdffn, width=18, height=72)
plot(ltr, show.tip.label=TRUE, cex=0.8)
add.scale.bar(x=1, y=length(ltr$tip.label))

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)



write.nexus(ltr, file="MPstrict_w_brlens.nexus")
write.tree(ltr, file="MPstrict_w_brlens.newick")

# No branch lengths
ltr2 = ltr
ltr2$edge.length = NULL
write.nexus(ltr2, file="MPstrict_wo_brlens.nexus")
write.tree(ltr2, file="MPstrict_wo_brlens.newick")



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
sum(charstats$cscores)
# 712
treestats$TL
# 712

CIsum = calc_CI(cscores=sum(charstats$cscores), minsteps=sum(charstats$minsteps))
CIsum
# 0.6095506
treestats$CI
# 0.6095506

RIsum = calc_RI(cscores=sum(charstats$cscores), minsteps=sum(charstats$minsteps), maxsteps=sum(charstats$maxsteps))
RIsum
0.687991
treestats$RI
# 0.688

calc_RI(cscores=712, minsteps=434, maxsteps=1325)
# 0.687991
(1325 - 712) / (1325 - 434)
# 0.687991

RCIsum = calc_RCI(cscores=sum(charstats$cscores), minsteps=sum(charstats$minsteps), maxsteps=sum(charstats$maxsteps))
RCIsum
# 0.4193653
treestats$RCI
# 0.4193653

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
pdf(pdffn, width=18, height=72)

# Settings
tip_hadj = -0.1
tipcex = 0.7
edgecex = 0.7
scalebar_cex = 0.7
titletxt = "MP tree"

ntips = length(ltr$tip.label)


# Branch lengths
subtitle_txt = "branch lengths"
ltr$root.edge = NULL
ape::plot.phylo(ltr, show.tip.label=FALSE)
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











#######################################################
# Output trees, with good labels, midpoint or not
#######################################################
ltr_relabeled = relabel_tree_tips_wGID_first(tr=ltr, fullnames=fullnames)
internal_nodenums = (length(ltr_relabeled$tip.label)+1):(length(ltr_relabeled$tip.label)+ltr_relabeled$Nnode)
#trbl = phytools::midpoint.root(trbl)

columns_to_save = c("TNT_synapos",
"absBremer",
"relBremer",
"freqBoot",
"gcBoot")


for (i in 1:length(columns_to_save))
	{
	replacement = paste0("_", columns_to_save[i], ".newick")
	outfn = gsub(pattern=".tnt", replacement=replacement, x=settings$branchlengths_fn)
	
	ltr_relabeled_w_nodeLabels = ltr_relabeled
	
	
	head(ltr_relabeled_w_nodeLabels$node.label)
	tail(ltr_relabeled_w_nodeLabels$node.label)
	
	cmdtxt = paste0("ltr_relabeled_w_nodeLabels$node.label[internal_nodenums] = ltr_table$", columns_to_save[i], "[internal_nodenums]")
	eval(parse(text=cmdtxt))

	head(ltr_relabeled_w_nodeLabels$node.label)
	tail(ltr_relabeled_w_nodeLabels$node.label)

	
	write.tree(ltr_relabeled_w_nodeLabels, file=outfn)


	# Output midpoint rooting as well
	replacement = paste0("_", columns_to_save[i], "_midroot.newick")
	outfn = gsub(pattern=".tnt", replacement=replacement, x=settings$branchlengths_fn)
	write.tree(phytools::midpoint.root(ltr_relabeled_w_nodeLabels), file=outfn)	
	}




