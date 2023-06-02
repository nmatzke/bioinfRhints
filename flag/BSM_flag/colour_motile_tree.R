library(ape)
library(BioGeoBEARS)
library("treeio")
library("ggtree")
library("ggplot2")



pdffn = "ATP_synthase_beta_v4.pdf"
pdf(file=pdffn, width=50, height=50)

try <- read.newick("./flag_partial/ATP_synthase_beta.nwk")

info = read.csv("./flag_partial/motile_annotations.txt", 
                sep = "\t",
                col.names = c("genome_name", "motility_number", "prediction", "motile"), 
                header = TRUE, 
                stringsAsFactors = FALSE)
names(info)[1] <- 'label'

try2 <- full_join(as.treedata(try), info, by='label')


groupInfo <- split(x = as.character(info[,1]),f = info[,4])
phylo <- groupOTU(try, groupInfo)
p1 <- ggtree(tr = phylo,layout = "circular",branch.length = "none",
             mapping = aes(color = group)) +
  scale_color_manual(values=c("black", "darkgreen", "red", "orange")) + 
  theme(legend.position="right")

plot(p1, cex=0.2)


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)
