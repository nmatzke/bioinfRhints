
# Installation of Bioconductor
# (hit 'n' for no when asked to update a zillion other packages)
install_txt='
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
library(BiocManager)
BiocManager::install(c("BiocGenerics", "annotate", "Biostrings" , "genomes", "S4Vectors", "IRanges", "XVector", "zlibbioc"))


library(BiocManager)
BiocManager::install(c("ggtree"))
BiocManager::install(c("treeio"))
'


library(ggtree)
library(treeio) # for read.raxml
library(ape)
library(BioGeoBEARS)
wd = "/drives/GDrive/__GDrive_projects/2019-08-22_Craig_Millar_notho_fish/04_RAxML/89seqs_raxml3/"
setwd(wd)

raxml_tr_w_bootstraps_fn = "RAxML_bipartitionsBranchLabels.result"

raxml_tr_w_bootstraps = read.raxml(raxml_tr_w_bootstraps_fn)

trfn = "RAxML_bipartitionsBranchLabels.result.newick"
raxml2nwk(infile=raxml_tr_w_bootstraps_fn, outfile=trfn)





source('/drives/GDrive/__github/BEASTmasteR_USE_GITHUB_INSTEAD/R/tree_utils_v1.R', chdir = TRUE)
#source('/GitHub/BEASTmasteR/R/tree_utils_v1.R', chdir = TRUE)




nexfn = "RAxML_bipartitionsBranchLabels.result.nexus"
bootstraps = extractBEASTstats_orig(file=nexfn, digits=4, printflag=FALSE) 


# Read in and process the bootstrap tree
res = read_brackets_prt(file=nexfn, is_newick=TRUE)
