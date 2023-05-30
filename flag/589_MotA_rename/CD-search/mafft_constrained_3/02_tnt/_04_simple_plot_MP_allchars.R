
#######################################################
# Plot of tip age, vs. parsimony node height
#######################################################
library(ape)	# for read/write NEXUS
library(BioGeoBEARS)	# for list2str, moref
library(gdata)	# for trim
library(Biograph) # for Date_as_year


source('/drives/Dropbox/_njm/__packages/TNTR_setup/tnt_R_utils_v1.R')
source('/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R')

wd = "/drives/GDrive/__GDrive_projects/2014-11-21_Randy_Irmis_autapomorphies/_02_tnt/"
setwd(wd)


pdffn = "simple_plot_MP_allchars.pdf"
pdf(file=pdffn, width=6, height=6)


trfn = "MPstrict_w_brlens.nexus"
tr = ladderize(read.nexus(trfn), right=FALSE)
plot(tr, show.tip.label=FALSE)
#add.scale.bar(x=10, y=60, length=10)
#axisPhylo2(minage=-170)
mtext(side=1, "steps above root", line=2.5)

atvals = c(0, 25, 50, 75, 100)
axis(side=1, at=atvals, labels=atvals)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)




#######################################################
# Copy to graphics
#######################################################
wd = "/drives/GDrive/__GDrive_projects/2014-11-21_Randy_Irmis_autapomorphies/_graphics/"
setwd(wd)


pdffn = "simple_plot_MP_allchars.pdf"
pdf(file=pdffn, width=4, height=4)

plot(tr, show.tip.label=FALSE)
#add.scale.bar(x=10, y=60, length=10)
#axisPhylo2(minage=-170)
mtext(side=1, "steps above root", line=2.5)

atvals = c(0, 25, 50, 75, 100)

axis(side=1, at=atvals, labels=atvals)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

