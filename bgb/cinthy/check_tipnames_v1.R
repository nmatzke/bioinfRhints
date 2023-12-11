library(BioGeoBEARS)

wd = "/GitHub/bioinfRhints/bgb/cinthy/"
setwd(wd)

trfn = "tr2_branches_fixed.newick"
tr = read.tree(trfn)

trtable = prt(tr)

#######################################################
# Make the tips match the geog
#######################################################

length(tr$tip.label)
# 106 tips


# Make geog_3areas.txt file in Excel...

geogfn = "geog_3areas.txt"

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Errors - identical names etc.
row.names(tipranges@df)

# Obviously different from:
cat(sort(tr$tip.label), sep="\n")
