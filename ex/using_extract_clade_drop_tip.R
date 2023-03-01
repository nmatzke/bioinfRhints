library(ape)
library(BioGeoBEARS)

faketree_string = "((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);"

fulltree = read.tree(file="", text=faketree_string)
plot(fulltree)
add.scale.bar()
title("Full tree, before any subsetting")


trtable = prt(fulltree, printflag=FALSE, get_tipnames=TRUE)
trtable

# Look at trtable. Let's pretend "P_fauriei2,P_hathewayi_1,P_hawaiiensis_WaikamoiL1,P_mauiensis_Eke" is our "microevolution" clade - this is node 26
node_ancestral_to_micro_group = 26

microevo_tree = extract.clade(phy=fulltree, node=node_ancestral_to_micro_group)
plot(microevo_tree)
add.scale.bar()
title("Microevolution tree, extracted from fulltree")


# For the macroevolution tree, drop the tips you extracted for the micro tree
tips_to_drop = c("P_fauriei2","P_hathewayi_1","P_hawaiiensis_WaikamoiL1","P_mauiensis_Eke")
macroevo_tree = drop.tip(phy=fulltree, tip=tips_to_drop)

plot(macroevo_tree)
add.scale.bar()
title("Macroevolution tree, extracted from fulltree")



# Write out the tree to Newick files
write.tree(microevo_tree, file="microevo_tree.newick")
write.tree(macroevo_tree, file="macroevo_tree.newick")

