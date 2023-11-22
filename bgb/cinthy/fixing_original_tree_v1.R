

library(BioGeoBEARS)

wd = "/GitHub/bioinfRhints/bgb/cinthy/"
setwd(wd)

trfn = "starbeast3_ANN_Epacrid_50_allruns_median.newick"
tr = read.tree(trfn)

trtable = prt(tr)

printall(trtable)


time_bp 131 18.471753 

18.471753-

should be:

18.471753-6.70155
11.7702

tip_dif1 = 72.94219-66.24064  
tip_dif1


Adjust all these nodes by:

multiplier1 = 11.7702 / 18.471753
0.6371999

39   3.3507758
40  12.5448828
43   6.9632225
41    2.576094
42    2.576094
44    8.157755
45    8.157755


tr$edge.length[39] = tr$edge.length[39] * multiplier1
tr$edge.length[40] = tr$edge.length[40] * multiplier1
tr$edge.length[41] = tr$edge.length[41] * multiplier1
tr$edge.length[42] = tr$edge.length[42] * multiplier1
tr$edge.length[43] = tr$edge.length[43] * multiplier1
tr$edge.length[44] = tr$edge.length[44] * multiplier1
tr$edge.length[45] = tr$edge.length[45] * multiplier1

plot(tr, show.tip.label=FALSE)

tip_dif2 = 6.939535-5.721052

69   69         69       15       tip       147   2.2366402      185              67.22114 5.721052    TRUE                      Styphelia_geniculata
70   70         70       15       tip       148   2.2366402      185              67.22114 5.721052    TRUE                          Styphelia_lucens

multiplier2 = (2.2366402-(6.939535-5.721052)) / 2.2366402

tr$edge.length[147] = tr$edge.length[147] * multiplier2
tr$edge.length[148] = tr$edge.length[148] * multiplier2


plot(tr, show.tip.label=FALSE)

is.ultrametric(tr)



trtable = prt(tr)

rev(sort(trtable$time_bp[1:length(tr$tip.label)]))


# Ultrametricize tree

tr2 = average_tr_tips(tr, fossils_older_than=0.5)

# Oops, fix 0-length branches

tr$edge.length[tr$edge.length < 0.001]

tr$edge.length[tr$edge.length < 0.001] = 0.01

tr2 = average_tr_tips(tr, fossils_older_than=0.5)

is.ultrametric(tr2)

tr2table = prt(tr2)

rev(sort(tr2table$time_bp[1:length(tr2$tip.label)]))

# Now it's ultrametric except for the fossils

write.tree(tr2, file="tr2_branches_fixed.newick")

