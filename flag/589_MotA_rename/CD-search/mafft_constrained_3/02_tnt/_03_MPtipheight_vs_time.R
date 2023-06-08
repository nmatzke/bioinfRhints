
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

#wd = "/drives/GDrive/__GDrive_projects/2015-08-01_AFA_cladogram/_02_TNT/2015-08-27_runs/allchars/"
#setwd(wd)

# Dates in decimal years before present
dates_fn = "/drives/GDrive/__GDrive_projects/2014-11-21_Randy_Irmis_autapomorphies/Muller_Reisz_2006_phylogeny_early_eureptiles/bills_decimal_years_ago.txt"
dates = read.table(dates_fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
time_height = (max(dates$decimal_year)) - dates$decimal_year
dates = cbind(dates, time_height)
dates

# Tree with branch lengths (trbl)
branchlengths_fn = "auto_branchlengths.tnt"
trbl = tntfile2R(branchlengths_fn, brlens=TRUE)

ltr = ladderize(trbl, right=TRUE)
plot(ltr)

ltr_table = prt(trbl, printflag=FALSE)
ltr_table
tipnums = 1:length(ltr$tip.label)

tipnames = ltr_table$label[tipnums]
tip_heights = ltr_table$node_ht[tipnums]

tip_ht_table1 = cbind(tipnames, tip_heights)

# Subtract the root height from the time heights of each set of tips
# (really just SEA)
rownums_in_dates = match(tip_ht_table1[,"tipnames"], table=dates$tipname)
time_height = dates$time_height[rownums_in_dates]
tip_ht_table1 = cbind(tip_ht_table1, time_height)
tip_ht_table1 = dfnums_to_numeric(adf2(tip_ht_table1))
tip_ht_table1


pdffn = "tip_height_vs_time.pdf"
pdf(file=pdffn, width=6, height=6)


source('/drives/Dropbox/_njm/_genericR_v1.R')
# Plot all points
x = tip_ht_table1$time_height
y = tip_ht_table1$tip_heights
model1 = linear_regression_plot(x, y, xlabel="time (years)", ylabel="character steps", tmppch=19, printall=TRUE, tmplinecol="black", tmplty=1, tmplwd=1, plottext=TRUE, legend_x="topleft", legend_y=NULL, xlim=minmax_pretty(x), ylim=minmax_pretty(y), slope1=FALSE, intercept_in_legend=FALSE, col="black")
title("Correlation of character changes with time")
x1 = x
y1 = y

# Plot just to set up
model3 = linear_regression_plot(x, y, xlabel="time (years)", ylabel="character steps", tmppch=".", printall=TRUE, tmplinecol="white", tmplty=1, tmplwd=1, plottext=FALSE, legend_x="topleft", legend_y=NULL, xlim=minmax_pretty(x), ylim=minmax_pretty(y), slope1=FALSE, intercept_in_legend=FALSE, col="white")


points(x=tip_ht_table1$time_height, y=tip_ht_table1$tip_heights, col="red")

slope = model1$coefficients[2]
intercept = model1$coefficients[1]
R2 = summary(model1)$r.squared
slope = summary(model1)$coefficients[2,1]
slopeSE = summary(model1)$coefficients[2,2]
slope95 = 1.96*slopeSE


intercept = summary(model1)$coefficients[1,1]
interceptSE = summary(model1)$coefficients[1,2]
intercept95 = 1.96*interceptSE
pval = summary(model1)$coefficients[2,4]
R2txt = paste("R2 = ", format(R2, digits=3), sep="")
slopetxt = paste("m=", format(slope, digits=3), " +/- ", format(slope95, digits=3), sep="")
pvaltxt = paste("p = ", format(pval, digits=3), sep="")
txt_to_plot = paste(R2txt, slopetxt, pvaltxt, sep="\n")
legend(x=0, y=70, bty="n", legend=txt_to_plot, cex=0.9, text.col="red")

tmpx1 = min(x1, na.rm=TRUE)
tmpx2 = max(x1[x1!=max(x1)], na.rm=TRUE)
tmpy1 = slope*tmpx1 + intercept
tmpy2 = slope*tmpx2 + intercept
segments(x0=tmpx1, y0=tmpy1, x1=tmpx2, y1=tmpy2, col="red", lty="solid", lwd=2)

title("Correlation of character changes with time")


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


