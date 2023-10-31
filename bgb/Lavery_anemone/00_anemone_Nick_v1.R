
# Load the package (after installation, see above).
library(optimx)         # You need to have some version of optimx available
# as it is a BioGeoBEARS dependency; however, if you
# don't want to use optimx, and use optim() (from R core) 
# you can set:
# BioGeoBEARS_run_object$use_optimx = FALSE
# ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)

#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files"

# Set your working directory for output files
# default here is your home directory ("~")
# Change this as you like
#wd = np("~")
#wd = np("~/Documents/programs/R packages/BioGeoBears")

# NJM
# wd = np("~/Documents/desktop stuff/anemones-Jasper Vining/anemone BioGeoBears")   ###anemone
wd = np("/GitHub/bioinfRhints/bgb/Lavery_anemone/")   ###anemone
setwd(wd)

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

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# Example files are given below. To run your own data,
# make the below lines point to your own files, e.g.
# trfn = "/mydata/frogs/frogBGB/tree.newick"
# geogfn = "/mydata/frogs/frogBGB/geog.data"

#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################
# This is the example Newick file for Hawaiian Psychotria
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
#trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))

#trfn = "test_newick.txt"
#trfn = "pruned_tree_NI_newick.txt"     #################
trfn = "anemone_tree-newick.txt"     ###anemone

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)
tr


pdffn = "anemone_tree.pdf"
pdf(file=pdffn, width=11, height=11)

plot(tr, cex=0.8)
title("anemone phylogeny, Lavery lab")                    ###anemone
axisPhylo() # plots timescale
mtext(text="Millions of years ago (Mega-annum, Ma)", side=1, line=2)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


#######################################################
# Geography file
# Notes:
# 1. This is a PHLYIP-formatted file. This means that in the 
#    first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
# 2. This is the same format used for C++ LAGRANGE geography files.
# 3. All names in the geography file must match names in the phylogeny file.
# 4. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
# 5. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################

# This is the example geography file for Hawaiian Psychotria
# (from Ree & Smith 2008)
#geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))
#geogfn = "geog_data_test.txt"
#geogfn = "geog_data_NI.txt"         #################
geogfn = "anemone_geog_data.txt"         ###anemone

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
#max_range_size = 3         ###
#max_range_size = 8         ###anemone
max_range_size = 6         ###anemone
#max_range_size = 4         ###anemone
#max_range_size = 2         ###anemone

####################################################
####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size.
# If you have more than about 500-600 states, the calculations will get REALLY slow,
# since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
# will just sit there and crunch, and never get through the calculation of the first
# likelihood.
# 
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
# To check the number of states for a given number of ranges, try:
#numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
#numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=12, maxareas=6, include_null_range=TRUE)    ###anemone
#numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)

# Large numbers of areas have problems:
#numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)

# ...unless you limit the max_range_size:
#numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)
####################################################
####################################################

####################################################
#TIME-STRATIFIED DISPERSAL 
#NOTE #2: To make these analyses work, you have to make sure 
#you give the BioGeoBEARS_run_object the correct 
#path+filename of each file you want it to see.  You also 
#have to uncomment some/all of these lines in the example 
#script, and edit the file names to match whatever your 
#filenames are: 
  
#  http://phylo.wikidot.com/biogeobears#toc7 
#================================ 
  
  # Set up the stratified part 
#BioGeoBEARS_run_object$timesfn = "speriods.txt" 
#BioGeoBEARS_run_object$dispersal_multipliers_fn =  "manual_dispersal_multipliers.txt" 
#BioGeoBEARS_run_object$dispersal_multipliers_fn ="dispersal_matrix_test.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt" 

# Divide the tree up by strata 
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE) 
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object,make_master_table=TRUE, plot_pieces=FALSE) 



#######################################################
#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################
#######################################################
# NOTE: The BioGeoBEARS "DEC" model is identical with 
# the Lagrange DEC model, and should return identical
# ML estimates of parameters, and the same 
# log-likelihoods, for the same datasets.
#
# Ancestral state probabilities at nodes will be slightly 
# different, since BioGeoBEARS is reporting the 
# ancestral state probabilities under the global ML
# model, and Lagrange is reporting ancestral state
# probabilities after re-optimizing the likelihood
# after fixing the state at each node. These will 
# be similar, but not identical. See Matzke (2014),
# Systematic Biology, for discussion.
#
# Also see Matzke (2014) for presentation of the 
# DEC+J model.
#######################################################
#######################################################

#######################################################
#######################################################

#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$timesfn = "anemone1_timeperiods.txt"      ####anemone

#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn ="dispersal_matrix_NI_H0.txt"      #################
BioGeoBEARS_run_object$dispersal_multipliers_fn ="anemone_dispersal_matrix_A.txt"      ###anemone

#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
#BioGeoBEARS_run_object$distsfn = "anemone distance matrix A.txt"     ###anemone###########################################
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
#BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$num_cores_to_use = 4     ###anemone
library(parallel)     ###anemone
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)


# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)     ###anemone
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table     ###anemone

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = FALSE
#resfn = "TestH0_DEC_M0_unconstrained_v1.Rdata"       ############
resfn = "Test1_DEC.Rdata"       ###anemone
if (runslow)
	{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
	} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
	}

#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"       ############
BioGeoBEARS_run_object$timesfn = "anemone1_timeperiods.txt"      ####anemone

#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_matrix_NI_H0.txt"       ############
BioGeoBEARS_run_object$dispersal_multipliers_fn ="anemone_dispersal_matrix_A.txt"      ###anemone
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"

#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
#BioGeoBEARS_run_object$distsfn = "anemone distance matrix A.txt"     ###anemone
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)     ###anemone
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table     ###anemone

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#resfn = "Psychotria_DEC+J_M0_unconstrained_v1.Rdata"
#resfn = "TestH0_DEC+J_M0_unconstrained_v1.Rdata"        ###########
resfn = "Test1_DEC+J.Rdata"        ####anemone
runslow = FALSE
if (runslow)
	{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
	} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
	}

#######################################################
# PDF plots
#######################################################
#pdffn = "TestH0_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"        ###########
pdffn = "Test1_DEC_vs_DEC+J.pdf"        ####anemone
pdf(pdffn, width=12, height=12)

#######################################################
# Plot ancestral states - DEC
#######################################################
#analysis_titletxt ="BioGeoBEARS DEC on TestH0_unconstrained"      ###########
analysis_titletxt ="BioGeoBEARS DEC on Anemones"      ###########

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Anemones"     #######

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
#######################################################
# DIVALIKE AND DIVALIKE+J ANALYSIS
#######################################################
#######################################################
# NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
# Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's
# processes the way DEC does, but only allowing the 
# processes DIVA allows (widespread vicariance: yes; subset
# sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
# DIVALIKE is a likelihood interpretation of parsimony
# DIVA, and it is "like DIVA" -- similar to, but not
# identical to, parsimony DIVA.
#
# I thus now call the model "DIVALIKE", and you should also. ;-)
#######################################################
#######################################################

#######################################################
# Run DIVALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"        ####anemone
BioGeoBEARS_run_object$timesfn = "anemone1_timeperiods.txt"        ####anemone
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_matrix_NI_H0.txt"    ##########
BioGeoBEARS_run_object$dispersal_multipliers_fn ="anemone_dispersal_matrix_A.txt"      ###anemone

#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"

#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
#BioGeoBEARS_run_object$distsfn = "anemone distance matrix A.txt"     ###anemone
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4     ###anemone
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)     ###anemone
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table     ###anemone

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = FALSE
#resfn = "Psychotria_DIVALIKE_M0_unconstrained_v1.Rdata" 
#resfn = "TestH0_DIVALIKE_M0_unconstrained_v1.Rdata"        ########
resfn = "Test1_DIVALIKE.Rdata"        ###anemone
if (runslow)
	{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
	} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
	}

#######################################################
# Run DIVALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"       ########
BioGeoBEARS_run_object$timesfn = "anemone1_timeperiods.txt"       ###anemone
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_matrix_NI_H0.txt"      #####
BioGeoBEARS_run_object$dispersal_multipliers_fn ="anemone_dispersal_matrix_A.txt"      ###anemone

#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"

#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
#BioGeoBEARS_run_object$distsfn = "anemone distance matrix A.txt"     ###anemone
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4      ###anemone
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)     ###anemone
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table     ###anemone

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#resfn = "Psychotria_DIVALIKE+J_M0_unconstrained_v1.Rdata"
#resfn = "TestH0_DIVALIKE+J_M0_unconstrained_v1.Rdata"       ########
resfn = "Test1_DIVALIKE+J.Rdata"       ###anemone
runslow = FALSE
if (runslow)
	{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
	} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
	}

#######################################################
# PDF plots
#######################################################
#pdffn = "Psychotria_DIVALIKE_vs_DIVALIKE+J_M0_unconstrained_v1.pdf"
pdffn = "Test1_DIVALIKE_vs_DIVALIKE+J.pdf"       ###anemone
pdf(pdffn, width=12, height=12)

#######################################################
# Plot ancestral states - DIVALIKE
#######################################################
#analysis_titletxt ="BioGeoBEARS DIVALIKE on Psychotria M0_unconstrained"
analysis_titletxt ="BioGeoBEARS DIVALIKE on Anemones"    #####

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DIVALIKE+J
#######################################################
#analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Psychotria M0_unconstrained"
analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Anemones"     ######

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
#######################################################
# BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
#######################################################
#######################################################
# NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
# not identical with the full Bayesian model implemented 
# in the "BayArea" program of Landis et al. (2013). 
#
# Instead, this is a simplified likelihood interpretation
# of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
# "d" and "e" work like they do in the DEC model of Lagrange 
# (and BioGeoBEARS), and then BayArea's cladogenesis assumption
# (which is that nothing in particular happens at cladogenesis) is 
# replicated by BioGeoBEARS.
#
# This leaves out 3 important things that are in BayArea:
# 1. Distance dependence (you can add this with a distances 
#    matrix + the "x" parameter in BioGeoBEARS, however)
# 2. A correction for disallowing "e" events that drive
#    a species extinct (a null geographic range)
# 3. The neat Bayesian sampling of histories, which allows
#    analyses on large numbers of areas.
#
# The main purpose of having a "BAYAREALIKE" model is 
# to test the importance of the cladogenesis model on 
# particular datasets. Does it help or hurt the data 
# likelihood if there is no special cladogenesis process?
# 
# BAYAREALIKE is a likelihood interpretation of BayArea,
# and it is "like BayArea" -- similar to, but not
# identical to, Bayesian BayArea.
# I thus now call the model "BAYAREALIKE", and you should also. ;-)
#######################################################
#######################################################

#######################################################
# Run BAYAREALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"      ###anemone
BioGeoBEARS_run_object$timesfn = "anemone1_timeperiods.txt"      ###anemone
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_matrix_NI_H0.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn ="anemone_dispersal_matrix_A.txt"      ###anemone

#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"

#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
#BioGeoBEARS_run_object$distsfn = "anemone distance matrix A.txt"     ###anemone

# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)     ###anemone

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = FALSE
#resfn = "Psychotria_BAYAREALIKE_M0_unconstrained_v1.Rdata"
resfn = "TestH0_BAYAREALIKE_M0_unconstrained_v1.Rdata"
if (runslow)
	{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
	} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
	}

#######################################################
# Run BAYAREALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"     ###
BioGeoBEARS_run_object$timesfn = "anemone1_timeperiods.txt"     ###anemone
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_matrix_NI_H0.txt"     ####
BioGeoBEARS_run_object$dispersal_multipliers_fn ="anemone_dispersal_matrix_A.txt"      ###anemone

#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"

#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
#BioGeoBEARS_run_object$distsfn = "anemone distance matrix A.txt"     ###anemone

# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)     ###anemone
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table     ###anemone

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#resfn = "Psychotria_BAYAREALIKE+J_M0_unconstrained_v1.Rdata"
resfn = "Test1_BAYAREALIKE+J.Rdata"    #####
runslow = FALSE
if (runslow)
	{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
	} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
	}

#pdffn = "Psychotria_BAYAREALIKE_vs_BAYAREALIKE+J_M0_unconstrained_v1.pdf"
pdffn = "Test1_BAYAREALIKE_vs_BAYAREALIKE+J.pdf"      #######
pdf(pdffn, width=12, height=12)

#######################################################
# Plot ancestral states - BAYAREALIKE
#######################################################
#analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Psychotria M0_unconstrained"
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Anemones"      ######

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - BAYAREALIKE+J
#######################################################
#analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Psychotria M0_unconstrained"
analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Anemones"     ######

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=1.5, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
# 
# CALCULATE SUMMARY STATISTICS TO COMPARE
# DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
# 
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# REQUIRED READING:
#
# Practical advice / notes / basic principles on statistical model 
#    comparison in general, and in BioGeoBEARS:
# http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
#########################################################################
#########################################################################

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Statistics -- DEC vs. DEC+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# RESULTS: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
#########################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")

# Look at the results!!
restable
teststable

#######################################################
# Save the results tables for later -- check for e.g.
# convergence issues
#######################################################

# Loads to "restable"
save(restable, file="restable_v1.Rdata")
load(file="restable_v1.Rdata")

# Loads to "teststable"
save(teststable, file="teststable_v1.Rdata")
load(file="teststable_v1.Rdata")

# Also save to text files
write.table(restable, file="restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\t")

#######################################################
# Model weights of all six models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AIC")
restable_AICc_rellike
free_params = row.names(resDECj$output@params_table[resDECj$output@params_table$type=="free",])
names(restable_AICc_rellike) = c("LnL", "numparams", free_params, "AICc", "AICc_wt")

# Also save to text files
write.table(restable_AIC_rellike, file="restable_AIC_rellike.txt", quote=FALSE, sep="\t")
write.table(restable_AICc_rellike, file="restable_AICc_rellike.txt", quote=FALSE, sep="\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc_rellike), file="restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\t")

#######################################################
# Your results should be:
#######################################################

# > restable
# 
#                 LnL numparams            d            e         j
# DEC           -34.5         2 3.504546e-02 2.835632e-02 0.0000000
# DEC+J         -20.9         3 1.000000e-12 1.000000e-12 0.1142811
# DIVALIKE      -33.1         2 4.474416e-02 1.000000e-12 0.0000000
# DIVALIKE+J    -21.1         3 2.001000e-09 1.000000e-12 0.1157199
# BAYAREALIKE   -40.3         2 1.738085e-02 3.040188e-01 0.0000000
# BAYAREALIKE+J -21.6         3 1.000000e-12 1.000000e-12 0.1081158

#######################################################
# The p-value of the LRT (Likelihood Ratio Test) tells you whether or not you can reject the
# null hypothesis that DEC and DEC+J confer equal likelihoods on the data
#
# AIC and AICc model weights are also shown, giving a sense of the relative probability of the two models.
#
# (One could easily do model weights between all 6 models, but this is not done here; s
#  see Brian O'Meara's AIC webpage: http://www.brianomeara.info/tutorials/aic )
#######################################################

# > teststable
# 
#              alt        null LnLalt LnLnull DFalt DFnull DF Dstatistic    pval        test       tail  AIC1  AIC2 AICwt1  AICwt2 AICweight_ratio_model1 AICweight_ratio_model2
# 1          DEC+J         DEC -20.95  -34.54     3      2  1      27.19 1.8e-07 chi-squared one-tailed  47.9 73.08   1.00 3.4e-06                 294893                3.4e-06
# 11    DIVALIKE+J    DIVALIKE -21.09  -33.15     3      2  1      24.13 9.0e-07 chi-squared one-tailed 48.17  70.3   1.00 1.6e-05                  63797                1.6e-05
# 12 BAYAREALIKE+J BAYAREALIKE -21.09  -33.15     3      2  1      24.13 9.0e-07 chi-squared one-tailed 48.17  70.3   1.00 1.6e-05                  63797                1.6e-05









# NJM

# Make a plot of "per-area probability (probability of each each *area* at each node)" instead of a plot of the probability of each *range*
# Following instructions here: http://phylo.wikidot.com/example-biogeobears-scripts#per-area_probabilities

#######################################################
# Plot DEC/DEC+J M0 per-area probabilities -- with COLOR
# 
# (This assumes you have run a standard BioGeoBEARS analysis
#  from the main example script first)
#######################################################

#wd = "set_this_to_wherever_you_like"
#setwd(wd)

library(cladoRcpp)

pdffn = "per-area_probabilities_anemones_all6_models_v1_wColor.pdf"
pdf(pdffn, height=12, width=12)

areas = getareas_from_tipranges_object(tipranges)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

barwidth_proportion = 0.01
barheight_proportion = 0.013
split_reduce_x = 0.85    # fraction of sizes for the split probabilities
split_reduce_y = 0.85    # fraction of sizes for the split probabilities

# Manual offsets, if desired (commented out below)
offset_nodenums = NULL
offset_xvals = NULL
offset_yvals = NULL

root.edge = TRUE






# DEC
titletxt = "Anemones: probability of occupancy per *area*\nunder DEC model"

offset_nodenums = NULL
offset_xvals = NULL
offset_yvals = NULL

# Plot at nodes
probs_each_area = plot_per_area_probs(tr, res=resDEC, areas=areas, states_list_0based=states_list_0based, titletxt=titletxt, cols_each_area=TRUE, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
Sys.sleep(1)    # wait for this to finish

# Calculate per-area probabilities for corners
relprobs_matrix = resDEC$ML_marginal_prob_each_state_at_branch_bottom_below_node
probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)

# Add to left corners
#offset_nodenums = c(20, 21)
#offset_xvals = c(0, 0)
#offset_yvals = c(-1, -0.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="left", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)

# Add to right corners
#offset_nodenums = c(19, 17, 20, 15)
#offset_xvals = c(0, 0, 0, 0)
#offset_yvals = c(1.25, 2, 1, 1.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="right", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)


# DEC+J
titletxt = "Anemones: probability of occupancy per *area*\nunder DEC+J model"

offset_nodenums = NULL
offset_xvals = NULL
offset_yvals = NULL

# Plot at nodes
probs_each_area = plot_per_area_probs(tr, res=resDECj, areas=areas, states_list_0based=states_list_0based, titletxt=titletxt, cols_each_area=TRUE, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
Sys.sleep(1)    # wait for this to finish

# Calculate per-area probabilities for corners
relprobs_matrix = resDECj$ML_marginal_prob_each_state_at_branch_bottom_below_node
probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)

# Add to left corners
#offset_nodenums = c(20, 21)
#offset_xvals = c(0, 0)
#offset_yvals = c(-1, -0.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="left", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)

# Add to right corners
#offset_nodenums = c(19, 17, 20, 15)
#offset_xvals = c(0, 0, 0, 0)
#offset_yvals = c(1.25, 2, 1, 1.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="right", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)




# DIVALIKE
titletxt = "Anemones: probability of occupancy per *area*\nunder DIVALIKE model"

offset_nodenums = NULL
offset_xvals = NULL
offset_yvals = NULL

# Plot at nodes
probs_each_area = plot_per_area_probs(tr, res=resDIVALIKE, areas=areas, states_list_0based=states_list_0based, titletxt=titletxt, cols_each_area=TRUE, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
Sys.sleep(1)    # wait for this to finish

# Calculate per-area probabilities for corners
relprobs_matrix = resDIVALIKE$ML_marginal_prob_each_state_at_branch_bottom_below_node
probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)

# Add to left corners
#offset_nodenums = c(20, 21)
#offset_xvals = c(0, 0)
#offset_yvals = c(-1, -0.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="left", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)

# Add to right corners
#offset_nodenums = c(19, 17, 20, 15)
#offset_xvals = c(0, 0, 0, 0)
#offset_yvals = c(1.25, 2, 1, 1.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="right", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)



# DIVALIKE +J
titletxt = "Anemones: probability of occupancy per *area*\nunder DIVALIKE+J model"

offset_nodenums = NULL
offset_xvals = NULL
offset_yvals = NULL

# Plot at nodes
probs_each_area = plot_per_area_probs(tr, res=resDIVALIKEj, areas=areas, states_list_0based=states_list_0based, titletxt=titletxt, cols_each_area=TRUE, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
Sys.sleep(1)    # wait for this to finish

# Calculate per-area probabilities for corners
relprobs_matrix = resDIVALIKEj$ML_marginal_prob_each_state_at_branch_bottom_below_node
probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)

# Add to left corners
#offset_nodenums = c(20, 21)
#offset_xvals = c(0, 0)
#offset_yvals = c(-1, -0.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="left", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)

# Add to right corners
#offset_nodenums = c(19, 17, 20, 15)
#offset_xvals = c(0, 0, 0, 0)
#offset_yvals = c(1.25, 2, 1, 1.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="right", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)




# BAYAREALIKE
titletxt = "Anemones: probability of occupancy per *area*\nunder BAYAREALIKE model"

# Plot at nodes
# res=resBAYAREALIKE; areas=areas; states_list_0based=states_list_0based; titletxt=titletxt; cols_each_area=TRUE; barwidth_proportion=barwidth_proportion; barheight_proportion=barheight_proportion; offset_nodenums=offset_nodenums; offset_xvals=offset_xvals; offset_yvals=offset_yvals; root.edge=root.edge
border="grey75"
trcol="black"
plot_rangesizes=FALSE
states_list = states_list_0based

probs_each_area = plot_per_area_probs(tr, res=resBAYAREALIKE, areas=areas, states_list_0based=states_list_0based, titletxt=titletxt, cols_each_area=TRUE, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, trcol=trcol, plot_rangesizes=plot_rangesizes, border=border)
print("Waiting...")
Sys.sleep(2)    # wait for this to finish

# Calculate per-area probabilities for corners
relprobs_matrix = resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_bottom_below_node
probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)

# Add to left corners
#offset_nodenums = c(20, 21)
#offset_xvals = c(0, 0)
#offset_yvals = c(-1, -0.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="left", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
Sys.sleep(1)    # wait for this to finish

# Add to right corners
#offset_nodenums = c(19, 17, 20, 15)
#offset_xvals = c(0, 0, 0, 0)
#offset_yvals = c(1.25, 2, 1, 1.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="right", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
print("Waiting...")
Sys.sleep(1)    # wait for this to finish




# BAYAREALIKE+J
titletxt = "Anemones: probability of occupancy per *area*\nunder BAYAREALIKE+J model"

# Plot at nodes
border="grey75"
trcol="black"
plot_rangesizes=FALSE
states_list = states_list_0based

probs_each_area = plot_per_area_probs(tr, res=resBAYAREALIKEj, areas=areas, states_list_0based=states_list_0based, titletxt=titletxt, cols_each_area=TRUE, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, trcol=trcol, plot_rangesizes=plot_rangesizes, border=border)
print("Waiting...")
Sys.sleep(2)    # wait for this to finish

# Calculate per-area probabilities for corners
relprobs_matrix = resBAYAREALIKEj$ML_marginal_prob_each_state_at_branch_bottom_below_node
probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)

# Add to left corners
#offset_nodenums = c(20, 21)
#offset_xvals = c(0, 0)
#offset_yvals = c(-1, -0.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="left", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
Sys.sleep(1)    # wait for this to finish

# Add to right corners
#offset_nodenums = c(19, 17, 20, 15)
#offset_xvals = c(0, 0, 0, 0)
#offset_yvals = c(1.25, 2, 1, 1.25)
add_per_area_probs_to_corners(tr, probs_each_area, left_or_right="right", cols_each_area=TRUE, barwidth_proportion=split_reduce_x*barwidth_proportion, barheight_proportion=split_reduce_y*barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge, border=border)
print("Waiting...")
Sys.sleep(1)    # wait for this to finish



dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

