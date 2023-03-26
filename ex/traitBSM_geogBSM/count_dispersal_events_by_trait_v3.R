# Load the package (after installation, see above).
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
                  # seems to sometimes fail on simple problems (2-3 parameters)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)

#######################################################
# 2018-10-10 update: I have been putting the 
# updates on CRAN/GitHub
# You should use:
# rexpokit version 0.26.6 from CRAN
# cladoRcpp version 0.15 from CRAN
# BioGeoBEARS version 1.1 from GitHub, install with:
# library(devtools)
# devtools::install_github(repo="nmatzke/BioGeoBEARS")
#######################################################
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

# Double-check your working directory with getwd()
getwd()

#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory 
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# 
# Each computer operating system might install BioGeoBEARS in a different place, 
# depending on your OS and settings. 
# 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)




#######################################################
# Load bodysize-as-discrete-trait analysis
#######################################################


#wd = "C:/Users/chloe/OneDrive/Documents/BGB/M0BSTrait/"
wd = "/GitHub/bioinfRhints/ex/traitBSM_geogBSM/M0BSTrait/"
setwd(wd)


trfn = "Cetacea_Trait_Safe_MAP_sub223B.newick"

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
pdffn = "tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = read.tree(trfn)
tr
plot(tr, cex=0.5)
title("Cetacea_Trait_Safe_MAP_sub223B.newick")
axisPhylo() # plots timescale

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

geogfn = "Trait_Safe.txt"

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges_bodysize = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges_bodysize

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges_bodysize@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
# 1. Set the maximum range size to 1 area by setting {{max_range_size=1}}
max_range_size = 1


# Load whale Mk bodysize trait ML result
resfn = "whales_Mk_bodysize_M0_unconstrained_v1.Rdata"
load(resfn)
resMk = res
rn(resMk)


# Loads to: RES_clado_events_tables
load(file="RES_clado_events_tables.Rdata")
# Loads to: RES_ana_events_tables
load(file="RES_ana_events_tables.Rdata")

trait_RES_clado_events_tables = RES_clado_events_tables
trait_RES_ana_events_tables = RES_ana_events_tables

length(trait_RES_clado_events_tables)
length(trait_RES_ana_events_tables)

dim(trait_RES_clado_events_tables[[1]])
dim(trait_RES_ana_events_tables[[1]])


#######################################################
########################################################
# Load geog results
#######################################################
#######################################################
#wd = "C:/Users/chloe/OneDrive/Documents/BGB/M0BSTrait/"
wd = "/GitHub/bioinfRhints/ex/traitBSM_geogBSM/M0Trait_geogOnly/"
setwd(wd)
max_range_size = 8


#tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()
trfn = "Cetacea_Trait_Safe_MAP_sub223B.newick"
geogfn = "Geographic_Trait_Safe_data.txt"

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Load the geog ML results
resfn = "BAYAREALIKEj_inf_geogOnly_wPreviousML.Rdata"
load(resfn)
resBAYAREALIKEj = res
rn(resBAYAREALIKEj)


# Load the BSM results

# Loads to: geog_RES_clado_events_tables
load(file="geog_RES_clado_events_tables.Rdata")
# Loads to: geog_RES_ana_events_tables
load(file="geog_RES_ana_events_tables.Rdata")

# Extract BSM output
clado_events_tables = geog_RES_clado_events_tables
ana_events_tables = geog_RES_ana_events_tables

length(clado_events_tables)
length(ana_events_tables)

dim(clado_events_tables[[1]])
dim(ana_events_tables[[1]])


#######################################################
# Simulate the source areas
#######################################################
areanames = names(tipranges@df)
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res=resBAYAREALIKEj, clado_events_tables=clado_events_tables_wTrait, ana_events_tables=ana_events_tables_wTrait, areanames=areanames)




#######################################################
# Count dispersal events by trait
# (NOTE: assumes non-stratified inputs
#######################################################

# Goal: make a table of dispersal events, with time, direction, etc.
# For each, check each trait stochastic map to determine what the body size was
# Output to text file (to avoid the slower-and-slower problem with appending lines 
#   to large tables)

# Actually, to keep this from exploding in size, we will just line up
# geog BSM #1 with trait BSM #1, and append to geog_BSM #1 the trait

# Copy the events tables
clado_events_tables_wTrait = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables_wTrait = BSMs_w_sourceAreas$ana_events_tables

i = 1
for (i in 1:length(clado_events_tables_wTrait))
{
if (i == 1)
	{
	cat("\nProcessing ", length(clado_events_tables_wTrait), " BSMs. Doing BSM #", i, ",", sep="")
	} else {
	cat(i, ",", sep="")
	}
	

clado_geog = clado_events_tables_wTrait[[i]]
ana_geog = ana_events_tables_wTrait[[i]]

# Add a traits column
clado_traits = rep(0, times=nrow(clado_geog))
ana_traits = rep(0, times=nrow(ana_geog))

# Bodysize trait BSM
trait_ana_BSM = trait_RES_ana_events_tables[[i]]
trait_clado_BSM = trait_RES_clado_events_tables[[i]]

# Go through the anagenetic events table
j=1

for (j in 1:nrow(ana_geog))
	{
	nodenum_at_top_of_branch = ana_geog$nodenum_at_top_of_branch[j]
	
	# Time-distance from highest tip (time=0 mya)
	abs_event_time = ana_geog$abs_event_time[j]
	
	# Time-distance from bottom of branch
	event_time = ana_geog$event_time[j]
	
	# Find the same timepoint on the bodysize traits table
	nodeTF = trait_clado_BSM$node == nodenum_at_top_of_branch
	
	# Anagenetic dispersal/extinction events
	if (trait_clado_BSM$anagenetic_events_txt_below_node[nodeTF] == "none")
		{
		ana_traits[j] = trait_clado_BSM$sampled_states_AT_nodes[nodeTF]
		} else {
		ana_events_on_this_branch_TF = trait_ana_BSM$nodenum_at_top_of_branch == nodenum_at_top_of_branch
		ana_events_on_this_branch = trait_ana_BSM[ana_events_on_this_branch_TF,]
		
		# Find which timespan the d/e event fit in
		# Start: 0 (bottom of branch, counting towards present
		# Stop: total branch length
		event_time_starts = c(0, ana_events_on_this_branch$event_time)
		event_time_stops = c(ana_events_on_this_branch$event_time, ana_geog$brlen[j])
		TF1 = event_time > event_time_starts
		TF2 = event_time <= event_time_stops
		TF = (TF1 + TF2) == 2
		
		matching_ana_event = ana_events_on_this_branch[TF,]
		current_rangetxt = matching_ana_event$current_rangetxt
		current_rangenum_1based = matching_ana_event$current_rangenum_1based
		
		ana_traits[j] = current_rangenum_1based
		}	# END if (trait_clado_BSM$anagenetic_events_txt_below_node[nodeTF] == "none")
	} # END for (j in 1:nrow(ana_geog))

bodysize = ana_traits
ana_geog = cbind(ana_geog, bodysize)
ana_events_tables_wTrait[[i]] = ana_geog
	
# Cladogenetic events
j=1
for (j in 1:nrow(clado_geog))
	{
	nodenum_at_top_of_branch = clado_geog$node[j]
	
	# Time-distance from highest tip (time=0 mya)
	abs_event_time = clado_geog$time_bp[j]
	
	# Find the same timepoint on the bodysize traits table
	nodeTF = trait_clado_BSM$node == nodenum_at_top_of_branch
	
	clado_row = trait_clado_BSM[nodeTF,]
	
	trait_at_node = clado_row$sampled_states_AT_nodes
	clado_traits[j] = trait_at_node
	} # END for (j in 1:nrow(clado_geog))
	
bodysize = clado_traits
clado_geog = cbind(clado_geog, bodysize)
clado_events_tables_wTrait[[i]] = clado_geog


} # END for (i in 1:length(clado_events_tables_wTrait))

clado_events_tables_wTrait[[1]]
ana_events_tables_wTrait[[1]]


tail(ana_events_tables_wTrait[[1]])
tail(ana_events_tables_wTrait[[2]])


wd = "/GitHub/bioinfRhints/ex/traitBSM_geogBSM/"
setwd(wd)

save(clado_events_tables_wTrait, file="clado_events_tables_wTrait.Rdata")
save(ana_events_tables_wTrait, file="ana_events_tables_wTrait.Rdata")

save(clado_events_tables_wTrait, file="clado_events_tables_wTrait_WORKED.Rdata")
save(ana_events_tables_wTrait, file="ana_events_tables_wTrait_WORKED.Rdata")


# Loads to clado_events_tables_wTrait
load("clado_events_tables_wTrait_WORKED.Rdata")
# Loads to ana_events_tables_wTrait
load("ana_events_tables_wTrait_WORKED.Rdata")



#######################################################
# Simulate the source areas of dispersal events
#######################################################

clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables


# Subset by bodysize trait: 1 vs. 2
clado_events_tables_wTrait_state1 = clado_events_tables_wTrait
clado_events_tables_wTrait_state2 = clado_events_tables_wTrait
ana_events_tables_wTrait_state1 = ana_events_tables_wTrait
ana_events_tables_wTrait_state2 = ana_events_tables_wTrait

for (i in 1:length(clado_events_tables_wTrait))
	{
	tmptable = clado_events_tables_wTrait[[i]]
	state1_TF = tmptable$bodysize == 1
	state2_TF = tmptable$bodysize == 2

	clado_events_tables_wTrait_state1[[i]] = tmptable[state1_TF,]
	clado_events_tables_wTrait_state2[[i]] = tmptable[state2_TF,]
	}

for (i in 1:length(ana_events_tables_wTrait))
	{
	tmptable = ana_events_tables_wTrait[[i]]
	state1_TF = tmptable$bodysize == 1
	state2_TF = tmptable$bodysize == 2

	ana_events_tables_wTrait_state1[[i]] = tmptable[state1_TF,]
	ana_events_tables_wTrait_state2[[i]] = tmptable[state2_TF,]
	}


timeperiods = c(0.0, 10.0)
areanames = names(tipranges@df)

counts_all = count_ana_clado_events(clado_events_tables=clado_events_tables_wTrait, ana_events_tables=ana_events_tables_wTrait, areanames=areanames, actual_names=areanames, timeperiod=timeperiods) 

counts_bodysize1 = count_ana_clado_events(clado_events_tables=clado_events_tables_wTrait_state1, ana_events_tables=ana_events_tables_wTrait_state1, areanames=areanames, actual_names=areanames, timeperiod=timeperiods) 

counts_bodysize2 = count_ana_clado_events(clado_events_tables=clado_events_tables_wTrait_state2, ana_events_tables=ana_events_tables_wTrait_state2, areanames=areanames, actual_names=areanames, timeperiod=timeperiods) 


counts_all$summary_counts_BSMs
counts_bodysize1$summary_counts_BSMs
counts_bodysize2$summary_counts_BSMs



# Count dispersals by trait, by timebin

timeperiods = c(20.0, 21.0)
areanames = names(tipranges@df)

counts_all = count_ana_clado_events(clado_events_tables=clado_events_tables_wTrait, ana_events_tables=ana_events_tables_wTrait, areanames=areanames, actual_names=areanames, timeperiod=timeperiods) 

counts_bodysize1 = count_ana_clado_events(clado_events_tables=clado_events_tables_wTrait_state1, ana_events_tables=ana_events_tables_wTrait_state1, areanames=areanames, actual_names=areanames, timeperiod=timeperiods) 

clado_events_tables=clado_events_tables_wTrait_state2; ana_events_tables=ana_events_tables_wTrait_state2; areanames=areanames; actual_names=areanames; timeperiod=timeperiods
counts_bodysize2 = count_ana_clado_events(clado_events_tables=clado_events_tables_wTrait_state2, ana_events_tables=ana_events_tables_wTrait_state2, areanames=areanames, actual_names=areanames, timeperiod=timeperiods) 



#######################################################
# Convert the BSMs to phytools format to add up branchlengths
#######################################################
#wd = "C:/Users/chloe/OneDrive/Documents/BGB/M0BSTrait/"
wd = "/GitHub/bioinfRhints/ex/traitBSM_geogBSM/M0BSTrait/"
setwd(wd)


#######################################################
# Checking that BSMs_to_phytools works
#######################################################
trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)

res=resMk; clado_events_tables=trait_RES_clado_events_tables; ana_events_tables=trait_RES_ana_events_tables

phytools_SM = BSM_to_phytools_SM(res, clado_events_table=trait_RES_clado_events_tables[[1]], ana_events_table=trait_RES_ana_events_tables[[1]])

phytools_SMs = BSMs_to_phytools_SMs(res=resMk, clado_events_tables=trait_RES_clado_events_tables, ana_events_tables=trait_RES_ana_events_tables)






#######################################################
# Adding up the branchlengths in each state, for BSM # i=1
#######################################################

clado_events_table = trait_RES_clado_events_tables[[1]]
ana_events_table = trait_RES_ana_events_tables[[1]]

timewidths_by_state = count_brlen_in_each_state(timeperiods=c(0,1), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state)
tail(timewidths_by_state)

colSums(timewidths_by_state)
sum(colSums(timewidths_by_state))


timewidths_by_state = count_brlen_in_each_state(timeperiods=c(0,100), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state)
tail(timewidths_by_state)

colSums(timewidths_by_state)
sum(colSums(timewidths_by_state))
sum(tr$edge.length)

# Check that it works in a division
timewidths_by_state0 = count_brlen_in_each_state(timeperiods=c(0,20), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state0)
tail(timewidths_by_state0)

colSums(timewidths_by_state0)
sum(colSums(timewidths_by_state0))


timewidths_by_state1 = count_brlen_in_each_state(timeperiods=c(20,100), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state1)
tail(timewidths_by_state1)

colSums(timewidths_by_state1)
sum(colSums(timewidths_by_state1))

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
sum(colSums(timewidths_by_state0)) + sum(colSums(timewidths_by_state1))

sum(tr$edge.length)

# WORKS!











#######################################################
# Adding up the branchlengths in each state, for BSM # i=2
#######################################################

clado_events_table = trait_RES_clado_events_tables[[2]]
ana_events_table = trait_RES_ana_events_tables[[2]]

timewidths_by_state = count_brlen_in_each_state(timeperiods=c(0,1), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state)
tail(timewidths_by_state)

colSums(timewidths_by_state)
sum(colSums(timewidths_by_state))


timewidths_by_state = count_brlen_in_each_state(timeperiods=c(0,100), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state)
tail(timewidths_by_state)

colSums(timewidths_by_state)
sum(colSums(timewidths_by_state))
sum(tr$edge.length)

# Check that it works in a division
timewidths_by_state0 = count_brlen_in_each_state(timeperiods=c(0,20), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state0)
tail(timewidths_by_state0)

colSums(timewidths_by_state0)
sum(colSums(timewidths_by_state0))


timewidths_by_state1 = count_brlen_in_each_state(timeperiods=c(20,100), res, trtable, clado_events_table, ana_events_table)
head(timewidths_by_state1)
tail(timewidths_by_state1)

colSums(timewidths_by_state1)
sum(colSums(timewidths_by_state1))

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
sum(colSums(timewidths_by_state0)) + sum(colSums(timewidths_by_state1))

sum(tr$edge.length)

# WORKS!






#######################################################
# Check for all 50 BSMs
#######################################################

for (i in 1:length(trait_RES_clado_events_tables))
	{
	print(i)
	clado_events_table = trait_RES_clado_events_tables[[i]]
	ana_events_table = trait_RES_ana_events_tables[[i]]

	timewidths_by_state = count_brlen_in_each_state(timeperiods=c(0,1), res, trtable, clado_events_table, ana_events_table)
	head(timewidths_by_state)
	tail(timewidths_by_state)

	colSums(timewidths_by_state)
	sum(colSums(timewidths_by_state))


	timewidths_by_state = count_brlen_in_each_state(timeperiods=c(0,100), res, trtable, clado_events_table, ana_events_table)
	head(timewidths_by_state)
	tail(timewidths_by_state)

	colSums(timewidths_by_state)
	sum(colSums(timewidths_by_state))
	sum(tr$edge.length)

	# Check that it works in a division
	timewidths_by_state0 = count_brlen_in_each_state(timeperiods=c(0,20), res, trtable, clado_events_table, ana_events_table)
	head(timewidths_by_state0)
	tail(timewidths_by_state0)

	colSums(timewidths_by_state0)
	sum(colSums(timewidths_by_state0))


	timewidths_by_state1 = count_brlen_in_each_state(timeperiods=c(20,100), res, trtable, clado_events_table, ana_events_table)
	head(timewidths_by_state1)
	tail(timewidths_by_state1)

	colSums(timewidths_by_state1)
	sum(colSums(timewidths_by_state1))

	sum1 = sum(colSums(timewidths_by_state0)) + sum(colSums(timewidths_by_state1))

	sum2 = sum(tr$edge.length)
	
	if (sum1 != sum2)
		{
		print("ERROR! sum1 != sum2. Printing...")
		print(sum1)
		print(sum2)
		} else {
		print("Sums match")
		}
	}
# WORKS!





