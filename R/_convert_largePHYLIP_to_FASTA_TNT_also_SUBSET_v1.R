
#######################################################
# Convert a huge PHYLIP alignment to FASTA
# (and maybe do neighbor-joining)
#######################################################
library(phylotools)


fn = "/Volumes/tb/tb/r20180118/assemble_on_ref/r20180118_assembly1_outfiles/r20180118_assembly1.phy"
tmp = read.phylip(fn)

length(tmp)
class(tmp)
nchar(tmp[1])
nchar(tmp[2])
number_of_sites = nchar(tmp[2])
# > length(tmp)
# [1] 104
# > class(tmp)
# [1] "phy"
# > nchar(tmp[1])
# [1] 11
# > nchar(tmp[2])
# [1] 5172125



###############################################################
# Write each row to FASTA
###############################################################
outfn = gsub(pattern="\\.phy", replacement=".fasta", x=fn)
outfn

i = 2
cat("\n\nWriting ", length(tmp)-1, " PHYLIP sequences to FASTA, writing sequence #: ", sep="")
for (i in 2:length(tmp))
	{
	cat(paste0(i-1, " "), sep="")
	words = strsplit(tmp[i], split="\\s+")[[1]]
	length(words)
	
	seq_header = paste0(">", words[1])
	
	if (i == 2)
		{
		write(x=seq_header, file=outfn, sep="\n", append=FALSE)
		} else {
		write(x=seq_header, file=outfn, sep="\n", append=TRUE)
		}
	write(x=words[2], file=outfn, sep="\n", append=TRUE)
	}
cat("...done.\n\n")


###############################################################
# Write each row to TNT
###############################################################
tnt_outfn = gsub(pattern="\\.phy", replacement=".tnt", x=fn)
tnt_outfn

# Write the first lines
words = strsplit(tmp[1], split="\\s+")[[1]]
words
numcols = words[2]
num_species = words[1]
numcols
num_species

line1 = "xread"
line2 = paste(numcols, num_species, sep=" ")
line3 = "& [dna]"


seq_header = c(line1, line2, line3)
write(x=seq_header, file=tnt_outfn, sep="\n", append=FALSE)


i = 2
cat("\n\nWriting ", length(tmp)-1, " PHYLIP sequences to TNT, writing sequence #: ", sep="")
for (i in 2:length(tmp))
	{
	cat(paste0(i-1, " "), sep="")
	words = strsplit(tmp[i], split="\\s+")[[1]]
	length(words)
	
	tmpline = paste(words[1], words[2], sep="\t")
	write(x=tmpline, file=tnt_outfn, sep="\n", append=TRUE)
	}
cat("...done.\n\n")





###############################################################
###############################################################
# SUBSET ALIGNMENTS
###############################################################
###############################################################

###############################################################
# Write each row to PHYLIP - SUBSET
###############################################################
library(stringr)
numcols = number_of_sites
set.seed(12345)
# Take just 100,000 columns
num_columns_to_keep = 100000
random_samples = sample(x=1:numcols, size=num_columns_to_keep, replace=FALSE)
random_samples = sort(random_samples)
random_samples_fn = paste0("randomly_sampled_", num_columns_to_keep, "_columns_in_subset.txt")
write(x=random_samples, file=random_samples_fn)


tmp2 = tmp
header = strsplit(tmp2[1], split="\\s+")[[1]]
header[2] = as.character(length(random_samples))
outfn = gsub(pattern="\\.phy", replacement="_subset.phy", x=fn)
outfn
header_txt = paste0(header[1], "\t", header[2])
cat("\Writing PHYLIP header: ", header_txt, "\n...done.", sep="")
write(x=header_txt, file=outfn, sep="\n", append=FALSE)

i = 2
cat("\n\nWriting ", length(tmp)-1, " PHYLIP sequences to *SUBSET* PHYLIP, writing sequence #: ", sep="")
for (i in 2:length(tmp))
	{
	cat(paste0(i-1, " "), sep="")
	words = strsplit(tmp[i], split="\\s+")[[1]]
	seqname = words[1]
	tmpstring = words[2]
	subset_string = stringr::str_sub(string=tmpstring, start=random_samples, end=random_samples)
	newchars = paste0(subset_string, collapse="")
	newline = paste(seqname, newchars, sep="\t")
	tmp2[i] = newline
	write(x=newline, file=outfn, sep="\n", append=TRUE)
	}
cat("...done.\n\n")


###############################################################
# Write random subset of each row to FASTA
###############################################################
outfn = gsub(pattern="\\.phy", replacement="_subset.fasta", x=fn)
outfn

i = 2
cat("\n\nWriting ", length(tmp)-1, " PHYLIP sequences to FASTA -- AND SUBSETTING THEM -- writing sequence #: ", sep="")
for (i in 2:length(tmp))
	{
	cat(paste0(i-1, " "), sep="")
	words = strsplit(tmp[i], split="\\s+")[[1]]
	length(words)
	
	seq_header = paste0(">", words[1])
	
	if (i == 2)
		{
		write(x=seq_header, file=outfn, sep="\n", append=FALSE)
		} else {
		write(x=seq_header, file=outfn, sep="\n", append=TRUE)
		}
	
	tmpstring = words[2]
	subset_string = stringr::str_sub(string=tmpstring, start=random_samples, end=random_samples)
	newchars = paste0(subset_string, collapse="")
	write(x=newchars, file=outfn, sep="\n", append=TRUE)
	}
cat("...done.\n\n")




###############################################################
# Write each row to TNT - SUBSET
###############################################################
tnt_outfn = gsub(pattern="\\.phy", replacement="_subset.tnt", x=fn)
tnt_outfn

# Write the first lines
words = strsplit(tmp[1], split="\\s+")[[1]]
words
numcols = length(random_samples)
num_species = words[1]
numcols
num_species

line1 = "xread"
line2 = paste(numcols, num_species, sep=" ")
line3 = "& [dna]"


seq_header = c(line1, line2, line3)
write(x=seq_header, file=tnt_outfn, sep="\n", append=FALSE)


i = 2
cat("\n\nWriting ", length(tmp)-1, " PHYLIP sequences to TNT, writing sequence #: ", sep="")
for (i in 2:length(tmp))
	{
	cat(paste0(i-1, " "), sep="")
	words = strsplit(tmp[i], split="\\s+")[[1]]
	length(words)
	tmpstring = words[2]
	subset_string = stringr::str_sub(string=tmpstring, start=random_samples, end=random_samples)
	newchars = paste0(subset_string, collapse="")
	
	tmpline = paste(words[1], newchars, sep="\t")
	write(x=tmpline, file=tnt_outfn, sep="\n", append=TRUE)
	}
cat("...done.\n\n")










#######################################################
# MAFFT on subset:
# http://mafft.cbrc.jp/alignment/server/phylogeny.html
# takes 2 min to upload
#######################################################
# Uploaded, went with default (UPGMA, JC)




# WORKED
########################################################################################
Size = 103 sequences
Method = Alignment score → Rough distance → Average linkage (UPGMA) 
Alignment id = .180128100727343fs3QH04zhED

########################################################################################


Click tree, UPGMA rough tree

#########################################################################################
id = .171123154929651
________________________________________________________________________

Current load level of the server: load

Constructing a tree from alignment .1711231546mafftweb21280ntuHD0VLayj23mbO9SCoW.
This page will be reloaded 2.5 seconds later.
Reload now
Progress (in a new window)
You can bookmark 
https://mafft.cbrc.jp/alignment/server/spool/_nj.171123154929651.html
to retrieve the result when done.
#########################################################################################






#######################################################
# FASTTREE on subset
#######################################################
http://www.genome.jp/tools-bin/ete
http://www.microbesonline.org/fasttree/
Takes a 10+ minutes, just to upload!

2018-01-28
=========================
Phylogenetic analysis pipeline by ETE3
Your job ID is 1801281017105f19381e48111499f6002bc5a7e109e6fa5c12c3

Your job is still running.
Your result will be summarized at the following URL.
http://www.genome.jp/tools-bin/ete?id=1801281017105f19381e48111499f6002bc5a7e109e6fa5c12c3
=========================




#Download:
cd /Users/nickm/Downloads/fasttree
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

cd /Volumes/tb/tb/r20180118/analysis-fasttree 

./FastTree -nt r20180118_assembly1_v2_subset.fasta > r20180118_assembly1_v2_subset_fasttreeJCcat.newick


./FastTree -gtr -nt r20180118_assembly1_v2_subset.fasta > r20180118_assembly1_v2_subset_fasttreeGTR.newick
started 2017-09-20


#######################################################
# Trying MAFFT phylogeny with UPGMA here:
http://mafft.cbrc.jp/alignment/server/phylogeny.html
#######################################################
* upload takes a few minutes

* Running, 2017-09-20_4:34 pm
Job id: 1709201530s2486452aUZRijkgy8r 
Data size: 199 sequences × (2251552-2251552) nuc 
Command: input > output
Load level of the server: load

This frame will be reloaded 2.8 seconds later. 
Reload now 
Progress (in a new window)

Cancel

http://mafft.cbrc.jp/alignment/server/spool/_out1709201530s2486452aUZRijkgy8r.html?8


UPGMA tree:


id = 1709201658s24813059
________________________________________________________________________

Current load level of the server: load

Constructing a tree from alignment 1709201530s2486452aUZRijkgy8r.
This page will be reloaded 3.6 seconds later.
Reload now
Progress (in a new window)
You can bookmark 
http://mafft.cbrc.jp/alignment/server/spool/_nj1709201658s24813059.html
to retrieve the result when done.

...THIS ONLY GIVES UPGMA





#######################################################
# ExaML
#######################################################

install:
cd /Users/nickm/Downloads/ExaML-master/parser 
make -f Makefile.SSE3.gcc

cd ../examl/
make -f Makefile.SSE3.gcc
rm *.o


make -f Makefile.AVX.gcc
rm *.o


# Convert PHYLIP aligned data matrix to binary format:
cd /Users/nickm/Downloads/ExaML-master/parser

./parse-examl -s /Volumes/tb/tb/r20180118/assemble_on_ref/r20180118_assembly1_outfiles/r20180118_assembly1_subset.phy -m DNA -n r20180118_assembly1_subset_parser

==========================================
gappyness: 0.871604
Pattern compression: ON

Your alignment has 137424 unique patterns


Under CAT the memory required by ExaML for storing CLVs and tip vectors will be
902463408 bytes
881311 kiloBytes
860 MegaBytes
0 GigaBytes


Under GAMMA the memory required by ExaML for storing CLVs and tip vectors will be
3527811504 bytes
3445128 kiloBytes
3364 MegaBytes
3 GigaBytes

Please note that, these are just the memory requirements for doing likelihood calculations!
To be on the safe side, we recommend that you execute ExaML on a system with twice that memory.


Binary and compressed alignment file written to file r20180118_assembly1_v2_parser.binary

Parsing completed, exiting now ... 
==========================================

mv RAxML_info.r20180118_assembly1_v2_parser.binary /Volumes/tb/tb/r20180118/analysis-examl/ 
mv r20180118_assembly1_v2_parser.binary /Volumes/tb/tb/r20180118/analysis-examl/ 




Parsimonator:
https://sco.h-its.org/exelixis/web/software/parsimonator/index.html
https://github.com/stamatak/Parsimonator-1.0.2

unzip and:
cd /Users/nickm/Downloads/Parsimonator-1.0.2-master
make -f Makefile.SSE3.gcc

./parsimonator-SSE3
















#######################################################
# TNT
#######################################################
cd /Volumes/tb/tb/r20180118
mkdir analysis-tnt

# COPY IN:
# * tnt.command, tnt
# * the TNT-formatted sequence subset
# * converting_long_tipnames_to_TNT_v1.R
# * rename_taxa_TNT_v1.R
# * auto.run



#######################################################
# Subset to parsimony-informative characters!!
#######################################################
cd /Volumes/tb/tb/r20180118/analysis-tnt
./tnt
# Increase RAM
mxram 20000;
proc r20180118_assembly1.tnt;

# TNT - max. ram = 15000.00 Mbytes 
# No ram used yet 
# tnt*>proc r20180118_assembly1.tnt;
# Reading from r20180118_assembly1.tnt 
# Matrix (5172095x103, 16 states). Memory required for data: 6941.11 Mbytes 
# Matrix (100000x103, 16 states). Memory required for data: 134.36 Mbytes 
# OLD:
#Matrix (2251552x199, 16 states). Memory required for data: 3454.19 Mbytes      
#procedure  - out of ram (change "hold" and/or "mxram")                         
#proc &r20180118_assembly1_v2.fasta; 
#auto;


xinact;
# Save to parsinf-only
log parsinf.tnt;
xread -;
log /;

blocks;

# NOTE: Manually add "& [dna]" to line 4 of parsinf.txt

#######################################################
# END Subset to parsimony-informative characters!!
#######################################################

# 5145588 uninformative characters are inactive 
# 99535 uninformative characters are inactive
# OLD:
# 2240815 uninformative characters are inactive                                  
# xinact 9.00 secs. 

# This equals: 
# 26507 parsimony-informative sites

# Write the parsinf-characters to NEXUS
cd /Volumes/tb/tb/r20180118/analysis-tnt
./tnt
mxram 1000;
proc parsinf.tnt;
# Export parsinf characters to NEXUS:
export [ parsinf.nex;
# Saved data (Nexus format) to file parsinf.nex 
# Open in TextWrangler, convert {ACGT} to N

# Export parsinf characters to FASTA:
export ! parsinf.fasta;


# Cladogram of parsimony-informative characters
# NOTE: Manually add "& [dna]" to line 4 of parsinf.txt

cd /Volumes/tb/tb/r20180118/analysis-tnt
./tnt
mxram 1000;
proc parsinf.tnt;



xmult;       
# tnt*>proc parsinf.tnt;                                                         
# Reading from parsinf.tnt 
# Matrices will be read as alpha-numeric data 
# Matrix space: allocate for up to 8 states 
# 'Data saved from TNT'
# Matrix (26507x103, 8 states). Memory required for data:  12.64 Mbytes 
# tnt*>xmult;                                                                    
# Repl. Algor.     Tree        Score       Best Score   Time       Rearrangs.    
# 5     FUSE       5           ------      78058        0:01:54    18,716,352    
# Completed search.
# Total rearrangements examined: 18,716,352. 
# No target score defined.  Best score hit 1 times.                              
# Best score: 78058.  1 trees retained.

tsave * best_trees.tnt;
save *;
tsave / ;


# Runquickie chose this:
       ------------------------------------------------------------
        Search routine used: 
           a quick consensus estimation (Goloboff &Farris 2001),
           with 15 replications (each with default "xmult" but
           with 3 starting points instead of the default 5, and using 
           XSS --see Goloboff &Pol 2007).  The sectorial searches
           analyzing sectors of 60 or more taxa with a combined 
           strategy (5 starting points, 5 cycles of tree-drifting for
           each, fusing the results in 4 cycles).  Sectors selections:
           XSS dividing tree in 3 to 3 parts, 3 times; CSS
           and RSS with defaults.  For more details of CSS, RSS, 
           and tree-drifting, see Goloboff 1999; for details of
           XSS, see Goloboff &Pol 2007.  For details of "xmult",
           see documentation of TNT.
           Note: for consensus calculation, trees TBR-collapsed. 
        ------------------------------------------------------------ 

# Save tree with branchlengths
taxname =;
naked -;
ttag - ;
ttag = ;
tplot *0.;
blength *0.;
ttags /;
ttags;
tsave * auto_branchlengths.tnt;
save *;
tsave / ;


# Save topology with no branchlengths
tsave * auto_no_branchlengths.tnt;
save ;
tsave / ;


# Export tree to NEXUS
export - auto_branchlengths.nex;
export -> auto_branchlengths2.nex;
export > auto_branchlengths3.nex;
export < auto_branchlengths4.nex;

cd /Users/nickm/Downloads/ExaML-master/examl/
./examl -s r20180118_assembly1_v2_parser.binary -n examl_tree.newick -t auto_branchlengths4_v2.newick -m PSR



