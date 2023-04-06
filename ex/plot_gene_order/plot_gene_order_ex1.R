
# Set working directory
wd = "/GitHub/bioinfRhints/ex/plot_gene_order/"
setwd(wd)

# Read in a tree
library(ape)
trstr = "(((chimp:1,human:1):1,gorilla:2):1,orang:3);"

tr = read.tree(file="", text=trstr)
plot(tr)
tr$tip.label = paste0(tr$tip.label, "                                 ", sep="")
plot(tr)
axisPhylo()

# Read in an Excel table
library(gdata)
xlsfn = "fakegene_positions.xlsx"
xls = read.xls(xlsfn)
head(xls)

# Plot a * at x=4, for gorilla
#points(x=4, y=3, pch="3", col="blue", cex=5)
#points(x=4, y=2, pch="2", col="blue", cex=5)
#points(x=4, y=1, pch="1", col="blue", cex=5)
#points(x=4, y=4, pch="4", col="blue", cex=5)


# Let's say all motA genes start at
# plot x coordinate xcoord
xcoord = 4
y_box_thickness = 0.1

# Let's scale # of nucleotides to x-axis
xunit_equals_this_many_nucleotides = 700

for (i in 1:length(tr$tip.label))
	{
	# Draw a box from xcoord to length of gene
	tipnum = i

	# Find row in Excel table matching that tip
	tipname = trim(tr$tip.label[tipnum])
	tipname

	# Excel row matching
	rownum = match(tipname, table=xls$tipname)
	rownum

	xleft = xcoord
	ybottom = tipnum-y_box_thickness
	xright = xcoord + (xls$MotA_end[tipnum] - xls$MotA_start[tipnum]) * (1/xunit_equals_this_many_nucleotides)
	ytop = tipnum+y_box_thickness

	rect(xleft, ybottom, xright, ytop, col="lightgreen")
	}
