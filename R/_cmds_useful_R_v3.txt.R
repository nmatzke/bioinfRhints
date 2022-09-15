###########################################################
###########################################################
# SETUP
###########################################################
###########################################################
# R, close all devices/windows/ quartz windows
graphics.off()

# Keep the source in loaded functions (so you can see comments)
options(keep.source=TRUE)

# Turn off all BS stringsAsFactors silliness
options(stringsAsFactors = FALSE)

options(width=200)



# Packages needed for life...
# required
install.packages(c("optimx", "FD", "ape", "phylobase", "rexpokit", "cladoRcpp"))##, "LaplacesDemon"))

# suggested
install.packages("xtable", "plotrix", "gdata")	

# for rexpokit
install.packages("diversitree", "expoRkit")	




# working install of BioGeoBEARS
install.packages(pkgs=c("BioGeoBEARS"), "/Library/Frameworks/R.framework/Resources/library/", repos="http://cran.cnr.Berkeley.edu", type='source', INSTALL_opts=c("--no-multiarch"))


#######################################################
# Install packages, KEEP COMMENTS FROM SOURCE
#######################################################
options("keep.source"=TRUE)
options("keep.source.pkgs"=TRUE)

install.packages(c("TreePar"), lib="/Library/Frameworks/R.framework/Resources/library/", contriburl=contrib.url(getOption("repos"), 'source'), type='source', dependencies=TRUE, configure.args=list(R_KEEP_PKG_SOURCE=TRUE))




options("keep.source"=TRUE)
options("keep.source.pkgs"=TRUE)

fn = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/SSE/Ranjard_etal_2014_competition/phyloland/"
install.packages(pkgs=fn, repos=NULL, type='source', configure.args=list(R_KEEP_PKG_SOURCE=TRUE))



#######################################################
# Working with tables
#######################################################

# Turn off strings as factors silliness
options(StringsAsFactors=FALSE)



library(BioGeoBEARS)

# Get a table
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))

geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))

tipranges = getranges_from_LagrangePHYLIP(geogfn)
tipranges

tmptable = tipranges@df

tmptable



# Is it a data.frame or matrix?
class(tmptable)

# The classes of the columns
cls.df(tmptable)


# Force all columns with some numbers to be
# class numeric instead of character or factor
dfnums_to_numeric(tmptable)

# Dimensions of a table
dim(tmptable)

# Sorting example
a = rownames(tmptable)
b = order(rownames(tmptable))
c = sort(rownames(tmptable))
d = rownames(tmptable)[order(rownames(tmptable))]

cbind(a,b,c,d)





https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file/731237#731237?newreg=fcd4a9430c7f4b62b7c77c366ae507e2
==================================================================================
         || visible in terminal ||   visible in file   || existing
  Syntax  ||  StdOut  |  StdErr  ||  StdOut  |  StdErr  ||   file   
==========++==========+==========++==========+==========++===========
    >     ||    no    |   yes    ||   yes    |    no    || overwrite
    >>    ||    no    |   yes    ||   yes    |    no    ||  append
          ||          |          ||          |          ||
   2>     ||   yes    |    no    ||    no    |   yes    || overwrite
   2>>    ||   yes    |    no    ||    no    |   yes    ||  append
          ||          |          ||          |          ||
   &>     ||    no    |    no    ||   yes    |   yes    || overwrite
   &>>    ||    no    |    no    ||   yes    |   yes    ||  append
          ||          |          ||          |          ||
 | tee    ||   yes    |   yes    ||   yes    |    no    || overwrite
 | tee -a ||   yes    |   yes    ||   yes    |    no    ||  append
          ||          |          ||          |          ||
 n.e. (*) ||   yes    |   yes    ||    no    |   yes    || overwrite
 n.e. (*) ||   yes    |   yes    ||    no    |   yes    ||  append
          ||          |          ||          |          ||
|& tee    ||   yes    |   yes    ||   yes    |   yes    || overwrite
|& tee -a ||   yes    |   yes    ||   yes    |   yes    ||  append
==================================================================================





# PROBLEM
system.file(package="roxygen2")
http://s270.codeinspot.com/q/1517951

better:
http://stackoverflow.com/questions/4616088/roxygen-how-to-set-a-default-parameter-including-backslash-to-functions

===============================================
I use Roxygen to generate Rd files of my packages under developement, but I have some problems with functions with default parameter set to '\n', e.g.:
  lineCount <- function(text, sep='\n') {
       ...
   }
Which purpose is to count new line ('\n') characters in a string. The problem is that R CMD check gives a warning about:
Codoc mismatches from documentation object 'lineCount':
lineCount
  Code: function(text, sep = "\n")
  Docs: function(text, sep = " ")
  Mismatches in argument default values:
    Name: 'sep' Code: "\n" Docs: " "
The problem seems to me that caused by writing to the Rd file (writing to standard LaTeX files via cat() always requires to double escape characters for some purpose, e.g.: \\newline - as I experienced). If I put an extra backslash to the separator, like:
  lineCount <- function(text, sep='\\n') {
       ...
   }
The problem still presists, as in the code it looks like '\\n', but in the docs (Rd files) it looks '\n'.
Is there an easy solution for my problem? May be an extra tag in Roxygen which could define how to write the function's params to the Rd file? Sorry if asked too obvious question, but I am lost after google-ing for a while.
History: http://permalink.gmane.org/gmane.comp.lang.r.roxygen/24
Source:StackExchange

 
Answers (1)
 cbeleites
I also ran into problems with too much escaped " and vanishing \t. I ended up changing the parse.formals function in roxygen's Rd2.R as follows:
  parse.formals <- function(partitum) {
    formals <- partitum$formals
    if (!is.null(formals)) {
      formals <- lapply(formals, trim)
      formals <- lapply(formals, paste, collapse=" ")
      name.defaults <- zip.c(names(formals), formals)
      args <-
        do.call(paste, c(Map(function(name.default) {
          name <- car(name.default)
          default <- cadr(name.default)
          if (! is.character (default)) {  # too much escaped. 
                                           # Not sure when escaping is needed. 
                                           # param = c ("x", "y", "z") works now
            default <- gsubfn("\"(.*)\"",
                              function(x)
                              sprintf("\"%s\"", gsub("\"", "\\\\\"", x)),
                              as.character(default))
          }
          default <- gsub ("\t", "\\\\t", default) # the tabs and newlines are already
          default <- gsub ("\n", "\\\\n", default) # tab and newline here.
          if (is.null.string(default))
            name
          else
            sprintf('%s=%s', name, default)
        },
                             name.defaults),
                         sep=', '))

      append.Rd(usageTag(parse.function.name(partitum), args))
    }
  }
===============================================
cd /Users/nickm/Library/R/2.15/library/roxygen2




SOLUTION TO FIX THE .Rd files in /man
"\|"
"\+"












cd /Dropbox/_njm/__packages
rm BioGeoBEARS_0.1.tar.gz
rm -R BioGeoBEARS.Rcheck
rm BioGeoBEARS/Read-and-delete-me
perl -e "s/require\(roxygen2\)/\ /g;" -pi $(find BioGeoBEARS -type f)
perl -e 's/("\\\\\|")/"\\\\\\\\\|"/g;' -pi $(find BioGeoBEARS/man -type f)
perl -e 's/("\\\\\+")/"\\\\\\\\\+"/g;' -pi $(find BioGeoBEARS/man -type f)




#######################################################
# R CMD build and check
#######################################################
# This fixes loose require(roxygen2) text:
# perl -e "s/require\(roxygen2\)/\ /g;" -pi $(find BioGeoBEARS -type f)
# 
# This fixes the problem that .Rd files have \ wherever the original function
# had \\:  NOT YET DONT WORK
# WORKS: changes "\+" --> "\\+"
# perl -e 's/("\\\\\|")/"\\\\\\\\\|"/g;' -pi $(find BioGeoBEARS/man -type f)
# WORKS: changes "\|" --> "\\|"
# perl -e 's/("\\\\\+")/"\\\\\\\\\+"/g;' -pi $(find BioGeoBEARS/man -type f)


cd /Dropbox/_njm/__packages
rm BioGeoBEARS_0.1.tar.gz
rm -R BioGeoBEARS.Rcheck
rm BioGeoBEARS/Read-and-delete-me
perl -e "s/require\(roxygen2\)/\ /g;" -pi $(find BioGeoBEARS -type f)
perl -e 's/("\\\\\|")/"\\\\\\\\\|"/g;' -pi $(find BioGeoBEARS/man -type f)
perl -e 's/("\\\\\+")/"\\\\\\\\\+"/g;' -pi $(find BioGeoBEARS/man -type f)

# remove the rcpp.skeleton automated insert:
# NO: Just change Rcpp.package.skeleton module=FALSE and files=FALSE
# rm /Dropbox/_njm/__packages/BioGeoBEARS/src/rcpp_module.cpp

R CMD build BioGeoBEARS
ls BioGeoBEARS*.gz
R CMD check --no-multiarch --as-cran BioGeoBEARS_0.1.tar.gz

mor


# Change back to standard R
cd /Library/Frameworks/R.framework/Versions
rm Current
ln -s 2.15 Current
ls -ls

# Or 
# R CMD check --no-multiarch --as-cran BioGeoBEARS



# R command line:
#install.packages(pkgs="/Dropbox/_njm/__packages/BioGeoBEARS_0.22.tar.gz", lib='/Library/Frameworks/R.framework/Resources/library/', repos=NULL, type='source', INSTALL_opts=c("--no-multiarch"))

# If it passes, FTP to CRAN, and send email to: CRAN@R-project.org
# http://cran.r-project.org/
# Check: ftp://cran.r-project.org/incoming/

# FTP instructions:
# http://www.hosting.com/support/ftp/ftp-from-mac-osx-terminal
cd /Dropbox/_njm/__packages/
ftp ftp://CRAN.R-project.org/incoming/
put BioGeoBEARS_0.1.tar.gz


















NOTES: When MAC is slow,
* MAKE SURE the file searching is not active!!


# R packages can be done with roxygen
# See directory:
#  /Dropbox/_njm/__packages/pseudoprime_setup/_compile_pseudoprime_v1.R 



# Set default time zone to UTC ALWAYS
Sys.setenv(TZ='UTC')

Sys.getenv(R_ARCH)


	# Get the number of taxa in each matching clade; i.e. number of commas + 1
	number_of_taxa_in_clade = 1 + count.fields(textConnection(tipnames_that_match), sep=",")




# Subplot formatting

# http://www.statmethods.net/advgraphs/layout.html
# The layout( ) function has the form layout(mat) where
# mat is a matrix object specifying the location of the N figures to plot.

# One figure in row 1 and two figures in row 2
attach(mtcars)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
hist(wt)
hist(mpg)
hist(disp)




#######################################################
# Latex from R
#######################################################


#######################################################
# Check the latex PDF file
#######################################################
sourcedir = '/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))
fn = "/Dropbox/_njm/__packages/rexpokit.Rcheck/rexpokit-manual.tex"
wd = "/Dropbox/_njm/__packages/rexpokit.Rcheck/"
setwd(wd)
runlatex(fn, wd)

# Error msg: I can't find file `phvr8t
# http://www.mail-archive.com/r-sig-mac@stat.math.ethz.ch/msg05185.html
# The problem is that BasicTeX installs only the Adobe 'times' font but not 'courier'. 
# For this reason I had to install package 'courier' (and 'helvetic') using 'TeX Live Utility.app' 

# THE ABOVE WORKS.  THE BELOW DOES NOT


# Installing new fonts:
# http://tug.org/fonts/fontinstall.html
# 
# Figure out where your TeXLive install is:
kpsewhich --var-value TEXMFLOCAL

# Result: /usr/local/texlive/2011basic/texmf-local
# Or probably: /usr/local/texlive/2012/texmf-local
# which is just another way of writing /usr/local/texlive/texmf-local

# Download .deb file from:
# http://packages.ubuntu.com/hardy/all/texlive-fonts-recommended/download

# Unpack as follow:
cd /Users/nickm/Downloads/
ar vx mypackage.deb
tar -xzvf data.tar.gz

cd /Users/nickm/Downloads/
ar vx texlive-fonts-recommended_2007-13ubuntu0.1_all.deb 
# Double-click data.tar.gz, or maybe: tar -xzvf data.tar.gz



# copy in from the other directories; -n so NO OVERWRITING happens
# http://stackoverflow.com/questions/9392735/linux-how-to-copy-but-not-overwrite
cd /usr/local/texlive/2011basic/texmf-local
sudo cp -n -R /Users/nickm/Downloads/data/usr/share/texmf-texlive/ .

# Update TeX database TeX about the new fonts
sudo -H mktexlsr

# Tell TeX about new font
sudo -H updmap-sys --enable Map=/Users/nickm/Downloads/data/usr/share/texmf-texlive/fonts/map/dvips/courier/ucr.map 

# ugh ERROR can't find map file...














#######################################################
# For uploading to CRAN
#######################################################
CRAN_cmds='
cd /Dropbox/_njm/__packages
rm rexpokit_0.2.tar.gz
rm -R rexpokit.Rcheck


# remove the rcpp.skeleton automated insert:
# NO: Just change Rcpp.package.skeleton module=FALSE and files=FALSE
# rm /Dropbox/_njm/__packages/rexpokit/src/rcpp_module.cpp

R CMD build rexpokit
ls rexpokit*.gz
R CMD check --no-multiarch --as-cran rexpokit_0.2.tar.gz

more /Dropbox/_njm/__packages/rexpokit.Rcheck/00install.out


# If it passes, FTP to CRAN, and send email to: CRAN@R-project.org
# http://cran.r-project.org/
# Check: ftp://cran.r-project.org/incoming/

# FTP instructions:
# http://www.hosting.com/support/ftp/ftp-from-mac-osx-terminal
cd /Dropbox/_njm/__packages/
ftp ftp://CRAN.R-project.org/incoming/
put rexpokit_0.2.tar.gz




install.packages(pkgs=c("/Dropbox/_njm/__packages/rexpokit_0.2.tar.gz"), "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source', INSTALL_opts=c("--no-multiarch"))
' # end CRAN_cmds
# Or 
# R CMD check --no-multiarch --as-cran rexpokit


#######################################################
# Check the latex PDF file
#######################################################
sourcedir = '/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))
fn = "/Dropbox/_njm/__packages/rexpokit.Rcheck/rexpokit-manual.tex"
wd = "/Dropbox/_njm/__packages/rexpokit.Rcheck/"
runlatex(fn, wd)


ftp://CRAN.R-project.org/incoming/ 
'

#######################################################
# R-Forge
#######################################################
# http://www.rubyrobot.org/tutorial/subversion-with-mac-os-x
#
# Package stored at:
# http://r-forge.r-project.org/projects/rexpokit/
#
# See the commits:
# https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=rexpokit 
# https://r-forge.r-project.org/scm/viewvc.php/pkg/R/rexpokit_v1.R?root=rexpokit&view=log
#
nmatzke@r-forge.r-project.org
1nepen



#######################################################
# Delete everything
#######################################################
# http://stackoverflow.com/questions/1461553/delete-all-files-from-svn-repository
cd /opt/local/Rpackages_svn_temp/

# Shallow checkout 
svn checkout --depth immediates svn+ssh://nmatzke@r-forge.r-project.org/svnroot/rexpokit/ copy_to_delete
cd copy_to_delete
svn rm *
svn ci -m "Deleting all"
svn status
svn update

#######################################################
# Upload new copy after deleting
#######################################################
cd /opt/local/Rpackages_svn_temp/
mkdir new_copy_to_copy_over_and_upload1

cd /opt/local/Rpackages_svn_temp/new_copy_to_copy_over_and_upload1

svn checkout svn+ssh://nmatzke@r-forge.r-project.org/svnroot/rexpokit/
cd rexpokit
cp -R /Dropbox/_njm/__packages/rexpokit/ .
svn add *
svn commit -m "Fresh upload"
svn status
svn update

#######################################################
# Check online
#######################################################
https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=rexpokit 

# Latest revision:
https://r-forge.r-project.org/scm/viewvc.php?root=rexpokit&view=rev







#######################################################
# For uploading to CRAN
#######################################################
CRAN_cmds='
cd /Dropbox/_njm/__packages
R CMD build rexpokit
ls rexpokit*.gz
R CMD check --as-cran rexpokit_0.1.tar.gz

more /Dropbox/_njm/__packages/rexpokit.Rcheck/00install.out

ftp://CRAN.R-project.org/incoming/ 
'

#######################################################
# R-Forge
#######################################################
# http://www.rubyrobot.org/tutorial/subversion-with-mac-os-x
#
# Package stored at:
# http://r-forge.r-project.org/projects/rexpokit/
#
# See the commits:
# https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=rexpokit 
# https://r-forge.r-project.org/scm/viewvc.php/pkg/R/rexpokit_v1.R?root=rexpokit&view=log
#
nmatzke@r-forge.r-project.org
1nepen

# Check out (wait an hour)
SVN_cmds='
svn checkout svn+ssh://nmatzke@r-forge.r-project.org/svnroot/rexpokit/

cp -R /Dropbox/_njm/__packages/rexpokit/ /opt/local/Rpackages/rexpokit/pkg/

cd /opt/local/Rpackages/rexpokit/pkg
svn add *
svn commit -m First commit by NJM; compiles locally on Intel Mac 10.7 with -arch x86 option; may need to turn on i386 to pass 
svn commit -m "Numerous edits to improve documentation and examples; I (NJM) probably won't make more edits." 
svn status
svn update

https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=rexpokit 

R CMD check
1nepen
'












#######################################################
# For uploading to CRAN
#######################################################
CRAN_cmds='
cd /Dropbox/_njm/__packages
rm modiscdata_0.12.tar.gz
rm -R modiscdata.Rcheck
rm modiscdata/Read-and-delete-me
R CMD build modiscdata
R CMD check --as-cran modiscdata_0.12.tar.gz



cd /Dropbox/_njm/__packages
rm modiscloud_0.12.tar.gz
rm -R modiscloud.Rcheck
rm modiscloud/Read-and-delete-me
R CMD build modiscloud
ls modiscloud*.gz
R CMD check --as-cran modiscloud_0.12.tar.gz

# CRAN update:
# If it passes, FTP to CRAN, and send email to: CRAN@R-project.org
# http://cran.r-project.org/
# Check: ftp://cran.r-project.org/incoming/

# FTP instructions:
# http://www.hosting.com/support/ftp/ftp-from-mac-osx-terminal
cd /Dropbox/_njm/__packages/
ftp ftp://CRAN.R-project.org/incoming/
put modiscloud_0.12.tar.gz
put modiscdata_0.12.tar.gz



#######################################################
# R-version switching, if desired
#######################################################
# Switch back to R-devel:
cd /Library/Frameworks/R.framework/Versions
ls -ls
rm Current
ln -s 3.0 Current
ls -ls
cd /Dropbox/_njm/__packages/


# Switch back to standard R:
cd /Library/Frameworks/R.framework/Versions
ls -ls
rm Current
ln -s 2.15 Current
ls -ls
cd /Dropbox/_njm/__packages/


# Doesn't work:
R_LIBS=
R CMD check -l /Library/Frameworks/R.framework/Versions/2.15/Resources/library --as-cran modiscloud_0.11.tar.gz
R CMD check -l . modiscdata_0.11.tar.gz --as-cran modiscloud_0.11.tar.gz


R CMD check --as-cran modiscdata_0.11.tar.gz modiscloud_0.11.tar.gz
R CMD check --as-cran modiscloud

more /Dropbox/_njm/__packages/modiscloud.Rcheck/00install.out


#######################################################
# R-Forge
#######################################################
# http://www.rubyrobot.org/tutorial/subversion-with-mac-os-x
#
# Package stored at:
# http://r-forge.r-project.org/projects/modiscloud/
#
# See the commits:
# https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=modiscloud 
# https://r-forge.r-project.org/scm/viewvc.php/pkg/R/modiscloud_v1.R?root=modiscloud&view=log
#
nmatzke@r-forge.r-project.org
1nepen



#######################################################
# Delete everything
#######################################################
# http://stackoverflow.com/questions/1461553/delete-all-files-from-svn-repository
cd /opt/local/Rpackages_svn_temp/

# Shallow checkout 
svn checkout --depth immediates svn+ssh://nmatzke@r-forge.r-project.org/svnroot/modiscloud/ copy_to_delete
cd copy_to_delete
svn rm *
svn ci -m "Deleting all"
svn status
svn update

#######################################################
# Upload new copy after deleting
#######################################################
cd /opt/local/Rpackages_svn_temp/
mkdir new_copy_to_copy_over_and_upload1

cd /opt/local/Rpackages_svn_temp/new_copy_to_copy_over_and_upload1

svn checkout svn+ssh://nmatzke@r-forge.r-project.org/svnroot/modiscloud/
cd modiscloud
cp -R /Dropbox/_njm/__packages/modiscloud/ .
svn add *
svn commit -m "Fresh upload"
svn status
svn update

#######################################################
# Check online
#######################################################
https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=modiscloud 

# Latest revision:
https://r-forge.r-project.org/scm/viewvc.php?root=modiscloud&view=rev
' # End CRAN commands














Use Rprof to run your code and see what is taking the time. 



#######################################################
# HPC stuff
#######################################################



Multicore parallel processing HPC high performance computing

# R!
> library(foreach)!
> library(doMC)!
> registerDoMC(cores=4)!
!
> system.time(foreach(i=1:10) %do% sum(runif(10000000)))!
!
user system elapsed !
4.796 0.448 5.245 !
!
> system.time(foreach(i=1:10) %dopar% sum(runif
(10000000)))!
!
user system elapsed !
4.332 0.609 1.459!

# R!
>library(multicore)!
>multicore:::detectCores()!
>options(cores = 8)!
>getOption('cores')!
>test <- lapply(1:10,function(x) rnorm(10000))!
!>
system.time(x <- lapply(test,function(x) loess.smooth
(x,x)))!
!
user system elapsed!
0.664 0.176 1.407!
!>
system.time(x <- mclapply(test,function(x) loess.smooth
(x,x)))!
!
user system elapsed!
0.008 0.008 0.351






#===================================================
# Starter source material
#===================================================
sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

sourcedir = '/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))

source2 = 'pbdbtools_v2.R'
source(paste(sourcedir, source2, sep=""))

source1 = 'pdb_utils_v6.R'
source(paste(sourcedir, source1, sep=""))

source8 = '_biogeog_sim_utils_v1.R'
source(paste(sourcedir, source8, sep=""))

# For latlong
source5 = '_phylogeostats_v1.R'
source(paste(sourcedir, source5, sep=""))
#===================================================


#===================================================
# Run with these commands
#===================================================
library(ape)
library(phylobase)
library(gdata) # needed for trim (whitespace trimming) function
library(lattice) # for histogram

# for mapping
library(dismo)
library(maptools)

# sp must be 0.9.60 or higher for rgdal

# Get rid of old sp
# remove.packages(c("sp"))

# better to install new sp manually, without dependencies
#install.packages(c("sp"), lib="/Library/Frameworks/R.framework/Resources/library/", contrib.url=contrib.url(getOption("repos"), 'source', type="source", dependencies=FALSE)

# 
library(sp)
library(rgdal) 	# for readOGR


# for loading ASCII grid
library(adehabitat)
library(maptools)

# For harmonic mean
# install.packages("psych")
library(psych)

# For variograms
library(gstat)

# For MARS
library(mda)

# Additional MARS functions for e.g. GLMs
# Source:
# http://onlinelibrary.wiley.com/doi/10.1111/j.1472-4642.2007.00340.x/suppinfo
# Jane Elith & John Leathwick (2007). "Predicting species distributions from 
# museum and herbarium records using multiresponse models fitted with multivariate
# adaptive regression splines." Diversity and Distributions, 13(3), 265-275.
# May 2007. DOI: 10.1111/j.1472-4642.2007.00340.x
sourcedir = '/Dropbox/_njm/'
source3 = 'mars.public.functions.3.1_njm1.R'
source(paste(sourcedir, source3, sep=""))
source3 = 'mmm.R'
source(paste(sourcedir, source3, sep=""))



##########################################
# Memory improvement (requires some of the above things for some reason)
##########################################
# based on comment in:
# https://r-forge.r-project.org/forum/message.php?msg_id=4025&group_id=294
#setOptions(chunksize = 1e+04, maxmemory = 1e+06)
setOptions(chunksize = 1e+05, maxmemory = 30e+06)

# http://r-sig-geo.2731867.n2.nabble.com/change-projection-of-large-raster-file-td6499121.html
# Doesn't seem to work: beginCluster(type="SOCK")

sessionInfo() 




GPU Graphics Processing Unit

GPU and R advantages:
>library(gputools)!
>matA <- matrix(runif(3*2), 3, 2)!
>matB <- matrix(runif(3*4), 3, 4)!
>gpuCrossprod(matA, matB) # Perform Matrix Cross-product
with a GPU!
!!>
numVectors <- 5!
>dimension <- 10!
>Vectors <- matrix(runif(numVectors*dimension),
>numVectors, dimension)!
>gpuDist(Vectors, "euclidean")!
>gpuDist(Vectors, "maximum")!
>gpuDist(Vectors, "manhattan")!
>gpuDist(Vectors, "minkowski", 4)!


bigmemory and other packages:
bigmemory: supports the creation, manipulation
and storage of large matrices.
S  bigalgebra: provides linear algebra functionality
with large matrices.
S  biganalytics: extends the functionality of
bigmemory.
S  bigtabulate: supports table(), split() and tapply()
like functionality for large matrices.
S  foreach + bigmemory: a winning combination for
massive data concurrent programming.


R profiling:
Pnmath
S  Profiling a program means determining how
much execution time a program spends in
various different sections of code.
S  We need to know where our code spends the
time to takes to compute our tasks.
S  R provides the tools for performance analysis.
O  The system.time function.
O  The Rprof for profiling R code.
O  The Rprofmem function for profiling memory usage.
S  In addition, the profr and proftools package
on CRAN can be used to visualize Rprof data.


Rprof("boot.out")!
##your code!
Rprof(NULL)!
!
##generates boot.out file!
!
Then run > R CMD Rprof boot.out!
!




R 2.15.2 contains improved 
library(parallel)
which includes mcmapply (multicore apply)

# Although, it doesn't seem to work on R.app (maybe)


library(parallel)
Ncores = parallel::detectCores()
parallel::mcmapply










#new_Makevars(new_R_Makevars="/Users/nickm/.R/Makevars_rexpokit")
old_Makevars()

# Check for 32-bit vs. 64-bit problems, i.e. i386 vs x86_64
# http://support.rstudio.org/help/discussions/problems/3120-i386x86_64-conflicts-in-libraries-on-mac-os-x-106-with-rstudio

Sys.getenv("R_ARCH")
# default returns:
# "/x86_64"

Sys.setenv(R_ARCH="/x86_64")
Sys.setenv(R_ARCH="/i386")
Sys.getenv("R_ARCH")


# This WORKS
# R.app is R64, so DON'T do multi-architecture
install.packages(packagename, lib="/Library/Frameworks/R.framework/Resources/library/", NULL, type='source', INSTALL_opts=c("--debug", "--no-multiarch"))

# (DON'T) Switch computer back to default Makevars
old_Makevars()







#######################################################
# 2012-05-08_ OpenBUGS on WINE on 10.7
#######################################################
http://www.davidbaumgold.com/tutorials/wine-mac/#part-0


openbugs on wine

no work:
/Users/nickm/.wine/drive_c

yes work:

Downloaded Winebottler from here:
http://winebottler.kronenberg.org/

...open OpenBUGS.exe in Winebottler

Installs to here:

/Users/nickm/Library/Application Support/Wine/prefixes/OpenBugs_WINE_v1/drive_c/Program Files/OpenBUGS/OpenBUGS321/

/Users/nickm/Library/Application\ Support/Wine/prefixes/OpenBugs_WINE_v1/drive_c/Program\ Files/OpenBUGS/OpenBUGS321/

open .

C:/Program Files/OpenBUGS/OpenBUGS321


...but -- characters are just [] [] [] [] [] !!!!!!!!!

#######################################
# RUNNING OPENBUGS ON A MAC -- SUCCESSFULLY
#######################################
Applications --> X11

cd /Users/nickm/Library/Application\ Support/Wine/prefixes/OpenBugs_WINE_v1/drive_c/Program\ Files/OpenBUGS/OpenBUGS321/

wine OpenBUGS.exe


######################
# Using OPENBUGS:
# After McCarthy 2007, pp. 250-251
######################
# Terminal:
cd /Users/nickm/Library/Application\ Support/Wine/prefixes/OpenBugs_WINE_v1/drive_c/Program\ Files/OpenBUGS/OpenBUGS321/

open .

(copy text files to this directory -- preferably odc files also, which are binaries I think)

# X11:
File --> Open --> text file
Inference --> Specification --> highlight "model" --> check model
highlight column headers or "list" for data --> load data
compile
highlight "list" for initial values --> load inits
or, gen inits to randomly generate them from the prio
compile

Model --> Update
Inference --> Samples --> type names of variables to sample
...or...
Inference --> Samples --> type * --> click trace
Click history, stats, density, etc.
Click update to run





To jiggle, wiggle, random points

jitter









multiline sed


# dput on an S4 object, or an S4 object within an S3 object, or something
# Source:
# http://stackoverflow.com/questions/3466599/dputting-an-s4-object
# can't get this modified version to insert commas appropriately, so we just
# have reformatted function here
dput2a <- function (x, file = "", control = c("keepNA", "keepInteger", "showAttributes"))
	{
    if (is.character(file))
    	{
        if (nzchar(file))
        	{
            file <- file(file, "wt")
            on.exit(close(file))
	        }
	    else
	    	{
	    	file <- stdout()
	    	}
	    }

    opts <- .deparseOpts(control)
    if (isS4(x)) # if #1
    	{
        cat("new(\"", class(x), "\"\n", file = file, sep = "")
        for (n in slotNames(x))
        	{
            cat("    ,", n, "= ", file = file)
            dput2a(slot(x, n), file = file, control = control)
        	}
        cat(")\n", file = file)
        invisible()
        } else if(length(grep('@',capture.output(str(x)))) > 0) # if #2
        {
        if(is.list(x))
        	{
			cat("list(\n", file = file, sep = "")
			for (i in 1:length(x))
				{
		        if(!is.null(names(x)))
		        	{
					n <- names(x)[i]
					if(n != '')
						{
						cat("    ,", n, "= ", file = file)
						}
			        }
          		dput2a(x[[i]], file = file, control = control)
          		}
			cat(")\n", file = file)
			invisible()
			} else {
			stop('S4 objects are only handled if they are contained within an S4 object or a list object')
			}
	    } else { #if #3
    	.Internal(dput(x, file, opts))
    	}
	}

# multiline sed; the example fixed dput2 output on S4 objects
multiline_sed <- function(fn, patternstr, newstr, outfn="sed_result.txt")
	{
	# R requires \\n here (for system() command)
	# UNIX/Terminal will just want \n
	
	#patternstr = ')\\nnew(\"Polygons\"'
	#newstr = '),new("Polygons"'
	
	# Do find/replace on mulitline (removes \n but who cares)
	
	# This sed works, although it removes the \n:
	#sed -n '1h;1!H;${;g;s|)\nnew(\"Polygons\"|),new("Polygons"|g;p;}' tmppoly2.txt > tmppoly2_edited.txt;

	k = 0
	cmds = NULL
	
	# Starter; from here:
	# http://austinmatzko.com/2008/04/26/sed-multi-line-search-and-replace/
	cmds[[(k=k+1)]] = "sed -n '1h;1!H;${;g;"
	
	# preface to starting string
	cmds[[(k=k+1)]] = "s|"
	
	# pattern string; user specifies escapes etc.
	#cmds[[(k=k+1)]] = ')\\nnew(\"Polygons\"'
	cmds[[(k=k+1)]] = patternstr
	
	# transition to replace string
	cmds[[(k=k+1)]] = "|"
	
	# replacement string; user specifies escapes etc.
	#cmds[[(k=k+1)]] = '),new("Polygons"'
	cmds[[(k=k+1)]] = newstr
	
	# End sed
	cmds[[(k=k+1)]] = "|g;p;}'"
	
	# input filename and output filename
	#outfn = gsub(pattern="\\.", replacement="_edited.", x=fn)
	#cmds[[(k=k+1)]] = paste(" ", fn, " > sed_result.txt;", sep="")
	cmds[[(k=k+1)]] = paste(" ", fn, " > ", outfn, ";", sep="")
	
	cmdstr = paste(cmds, collapse="", sep="")
	return(cmdstr)
	}







#or... for every 100 iterations

telliter <- 100
for( iter in 1:maxiter ) {
    #do some cool statistics
    if( iter %% telliter == 0 ) cat(paste("iteration", iter, "complete\n"))
} 






Trying Darwine:

http://mac.softpedia.com/dyn-postdownload.php?p=14184&t=0&i=2


re-installed the OpenBUGS package from there



Then went:

X11 terminal
cd /Users/nickm/Library/Application\ Support/Wine/prefixes/OpenBugs_WINE_v1/drive_c/Program\ Files/OpenBUGS/OpenBUGS321/
wine OpenBUGS.exe





UNIX

renaming files as a batch
http://staff.washington.edu/dittrich/misc/faqs/unix.rename.wildcard

The command sequence is then:

foreach file (*.sub22.t)
mv $file `basename $file .sub22.t`.sub2.t
end



R key bindings
To simply re-assign the key binding for the Edit:Execute command, you can just use the Keyboard System Preferences, which supports application-specific user bindings. If you don't mind interposing another layer of software (and any associated instability), you can also chain together a sequence of such keyboard shortcuts using a number of utilities. The shareware Butler (http://www.manytricks.com/butler/ ) allows this to be done in an application-specific way, and is easily configured to bind to the sequence Cmd-Return -> Option-DownArrow -> DownArrow, though you may need to insert delays between the individual steps. 





Testing for Phylogenetic Signal in R
(Pagel's Lambda)
http://bodegaphylo.wikispot.org/IV._Testing_Phylogenetic_Signal_in_R



# Multiple apply (mapply) allows counters
# apply with increments
counter=1:length(mbtrees_fixed2)
mbtrees_fixed2 = mapply(midpoint2, mbtrees_fixed, counter)



# apply gsub to every cell in a data frame, with mapply (multiple apply):

nexus_fn = "/Users/nickm/Desktop/__projects/_bird_phylo_Jessie/_data/AvFam_Prune-Merge_Data_v3simp.nex"
nexd = read_nexus_data2(nexus_fn, check_ambig_chars=TRUE)
# success

nexd_df_abcd = as.data.frame(nexd)
nexd_df_abcd = 
#printall(nex1_df)

# This file is in format a b c d e, convert to standard
uniq_chars = unique(as.character(unlist(nexd_df)))

# exclude all multiple character states; remove "?"
uniq_chars_1char = uniq_chars[nchar(uniq_chars) == 1]
uniq_chars_1char = uniq_chars_1char[uniq_chars_1char != "?"]
uniq_chars_1char = uniq_chars_1char[order(uniq_chars_1char)]

# Recode to "standard" format:
nexd_df = nexd_df_abcd
tmp_charstate_codes = charstate_codes()
for (i in 1:length(uniq_chars_1char))
	{
	cat("Converting '", uniq_chars_1char[i], "' --> '", tmp_charstate_codes[i], "'\n", sep="")
	
	# YES
	nexd_df = mapply(gsub, uniq_chars_1char[i], tmp_charstate_codes[i], nexd_df)

	# NO
	#nexd_df = gsub(uniq_chars_1char[i], tmp_charstate_codes[i], nexd_df)
	}


set seed:

set.seed(4)


###################################
#floating point strings
###################################

# make e.g. 00001
i =1 

for (i in 1:100)
	{
	tmpi = sprintf("%05.0f", i)
	tmpi
	
	new_fn1 = paste0("file_num_", i, ".pdf")
	new_fn2 = paste0("file_num_", tmpi, ".pdf")
	cat("My filenames: ", new_fn1, ",", new_fn2, "\n", sep="")
	}



significant digits
format a string to include e.g. 0004 instead of 4

weights_vals = c("-", paste("=", numbered_weights, sep=""), "[1 0")

txt_numbered_weights = c()
for (i in 1:length(numbered_weights))
	{
	tmp_wt = sprintf("%05.0f", 100*numbered_weights[i])
	txt_numbered_weights = c(txt_numbered_weights, tmp_wt)
	}


Number of digits:
number = 500
numdigits = nchar(as.character(number))

(sigdtf = signif_digits_df(newdtf))

# You can go down in significant digits...
(sigdtf = signif_digits_df(newdtf, numsig=2))

# But you can't go up...
(sigdtf = signif_digits_df(sigdtf, numsig=4))





# Set default time zone to UTC ALWAYS
Sys.setenv(TZ='UTC')



# Take all the tnttree files and convert to newick
tnttree_fns = list.files(pattern=".tnttree")

file_info = file.info(tnttree_fns)

# Do only files close in time (in minutes)
keepTF = (Sys.time() - file_info$mtime) < 30
tnttree_fns = tnttree_fns[keepTF]








# Ordered factors
theo<-c("cons", "mod", "cons", "cons", "lib", "mod")
table(theo)

theo<-ordered(theo, levels=c("lib", "mod", "cons"))
table(theo)

# GLM with ordered factors, i.e. ordinal/nominal data
glmobj = glm(formula=wetness~cloudiness+daynight+cloudiness*daynight, family=binomial(link="logit"), data=indata, weights=num_obs)





#######################################################
# Moving averages, smoothing with time series -- by timestamp, not just sequential
#######################################################

# time series analysis
times = as.POSIXct(modsite1$mod_POSIXct)
calc_mov_average_pt <- function(centerpoint, width_in_secs, times, vals)
	{
	defaults = '
	centerpoint = as.POSIXct("2011-03-12 10:25:00 UTC")
	width_in_secs = 48 * 60 * 60
	times = as.POSIXct(modsite1$mod_POSIXct)
	vals = modsite1$intvals
	'
	
	# produces difftime class
	timedifs = as.numeric(times) - as.numeric(centerpoint)
	valsTF = abs(timedifs) <= (width_in_secs/2)
	
	avg = mean(vals[valsTF], na.rm=TRUE)
	return(avg)
	}

avgs = sapply(X=times, FUN=calc_mov_average_pt, width_in_secs=96*60*60, times=times, vals=modsite1$intvals)

plot(times, avgs)





#######################################################
# wiggle points so they don't overlap
#######################################################
# to wiggle points:
jitter







vignette()

install.packages("phangorn")

install.packages(c("diversitree"), lib="/Library/Frameworks/R.framework/Resources/library/", contriburl=contrib.url(getOption("repos"), 'source'), type='source', dependencies=TRUE, configure.args=list(R_KEEP_PKG_SOURCE=TRUE))

library(phangorn)
?simSeq

vignette(Trees, package=phangorn)


install.packages("PHYLOGR")
library(PHYLOGR)
?PHYLOGR


vignette of the package
phylobase

vignette()
vignette("phylobase", package="phylobase")






# execute a string in python:
# evaluate a string in python -- eval perhaps
exec

#R:
.Last.value

.Options for options
?options to see options options

R, list all objects
ls()

# New commands
subset


# returns things matching the condition
which(is.na)


# approximate equality of floating point numbers
all.equal(x,y) is a utility to compare R objects x and y testing ‘near equality’. 


# get a specific thing out of an object
slot

# plot on top, not new plot
plot(new=FALSE)

# try
args(biogeomancer)
b = try( biogeomancer('Peru', locality=lonzero$locality[3], progress='') )
b

# vignette

# Search just the R website

colnames

# bring up a table
fix



# Data frame to list or matrix
c(as.matrix(rfd_m))
unlist(rfd_m)



# compare upper and lower triangles
# AAGH!  Library matlab screws up "sum" and maybe other standard functions...
install.packages("matlab")
library(matlab)
fliplr
flipud
sum
std
strcmp
rot90
reshape
ones
zeros
filesep
repmat
padarray
isempty
fix
eye
colorbar
ceil
ndims



# upper case
# lower case
tolower(x)
toupper(x)

#
chartr translates each character in x that is specified in old to the corresponding character specified in new. 


Help with functions in S4

#####################################################
# Phylobase intro
#####################################################
# getting help
library(phylobase)
rand_p4_tree <- as(rand_tree, "phylo4")
plot(rand_p4_tree)

# All fine and good, but how to we nd out about all the great features of the phylobase
# plotting function? R has two nifty ways to nd it, the rst is to simply put a question mark in
# front of the whole call:
`?`(plot(rand_p4_tree))

# R looks at the class of the rand p4 tree object and takes us to the correct help le (note:
# this only works with S4 objects). The second ways is handy if you already know the class of
# your object, or want to compare to generics for dierent classes:
`?`(method, plot("phylo4"))

# More information about how S4 documentation works can be found in the methods package,
# by running the following command.
help("Documentation", package = "methods")


# Get a list of datums
A=make_EPSG() 
A


================

Scripts dealing with factors & associated crapola!!!



dtf_classes = cls.df(pvals_results, printout=TRUE)
cls.df(dtf_classes, printout=TRUE)

# Convert numbers to numeric
newdtf = dfnums_to_numeric(dtf)

(sigdtf = signif_digits_df(newdtf))

# You can go down in significant digits...
(sigdtf = signif_digits_df(newdtf, numsig=2))

# But you can't go up...
(sigdtf = signif_digits_df(sigdtf, numsig=4))


(df_factors_to_char(sigdtf))

(x=df_everything_to_char(sigdtf))



=================

terminal: locate <filename>

/Library/Receipts/R-Framework.pkg/Contents/Info.plist
/Library/Receipts/R-Framework.pkg/Contents/Resources/BundleVersions.plist
/Library/Receipts/R-Framework.pkg/Contents/Resources/English.lproj/Description.plist
/Library/Receipts/R-GUI.pkg/Contents/Info.plist
/Library/Receipts/R-GUI.pkg/Contents/Resources/BundleVersions.plist
/Library/Receipts/R-GUI.pkg/Contents/Resources/English.lproj/Description.plist





# remove/clear all objects
rm(list=ls(all=TRUE))




# generic functions
flipdiag



# R, close all devices/windows/ quartz windows
graphics.off()

# Keep the source in loaded functions (so you can see comments)
options(keep.source=TRUE)

# Turn off all BS stringsAsFactors silliness
options(stringsAsFactors = FALSE)

options(width=200)


# Close all open graphics windows
graphics.off()

sourcedir = '/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))

sourcedir = '/Dropbox/_njm/'
source4 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source4, sep=""))



R, plot outside of subplots:

# outer margins
plot_NMMDS_list(NMDS_results_list, points_to_color)
title("Multiple figures in a plot", outer=TRUE, line=-2, cex=2)



R, print all attributes of an object:
attributes(z)
attributes(z)$data
head(attributes(z)$data)
# display the new table
str(meuse.pb)
rownames(meuse.pb)


# Read a text file into a list of strings
tmplines = scan(sptree_fn, what="character", sep="\n")


# split on non-whitespace (or words?)
unord_raw = scan(unordered_fn, what="character", sep="\n")
unord_words = strsplit(unord_raw, "\\w")[[1]]

# split on whitespace (or words?)
unord_raw = scan(unordered_fn, what="character", sep="\n")
unord_words = strsplit(unord_raw, "\\W")[[1]]

# number of characters
nchar(taxon_names[i], type="chars")





# XML comment node

# Comments list
tmp_comments_list = NULL
tmp_comments_list[[1]] = xmlCommentNode(" ")
tmp_comments_list[[2]] = xmlCommentNode(" ")
tmp_comments_list[[3]] = xmlCommentNode("Monophyly statistics")










# Example repeat loop
# be careful of infinite loops!
# http://stackoverflow.com/questions/4357827/do-while-loop-in-r
i = 1
repeat
	{
	i=i+1
	print(i)
	if (i > 10)
		{
		break
		}
	}










# Strip whitespace:
trim:
install.packages("gregmisc")
library(gregmisc)


# String splitting
strsplit
stringsplit

# Set operations:
# Performs set union, intersection, (asymmetric!) difference, equality and membership on two vectors.
union(x, y)
intersect(x, y)
setdiff(x, y)
setequal(x, y)
is.element(el, set)

# matching finding which elements match
match
%in%

see also:
which
search (not)
find



# check if a variable exists:
exists("t$node.label")




# close all graphics
graphics.off()

# means for each level in a factor

data(meuse.grid)
mg = meuse.grid
rowSums(table(mg$x, mg$y))

# or
tapply(mg$x, mg$y, length)




lib.loc <- .libPaths()
package = "gstat"
.find.package(package, lib.loc, verbose = verbose)


   
     
# PDF rotated to match page:
pdf(file="diagnostics_cophenetic_v1.pdf", width=10, height=7, paper="USr")

pdffn = "diagnostics_cophenetic_v1.pdf"
pdf(file=pdffn, width=10, height=7, paper="USr")

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


	#dev.off()
	graphics.off()

--> turn off all graphicsf


# paste lists: use collapse instead of sep
spaces_to_add = paste(rep(" ", numspaces_to_add), collapse="")




# versions of all packages
sessionInfo()

sessionInfo()






setwd("/Dropbox/_njm/__packages/")
Rcpp.package.skeleton( "RcppSkeleton" , force=TRUE, example_code=TRUE, module=TRUE)

#writeLines( system( "tree", intern = TRUE )
install.packages("/Dropbox/_njm/__packages/RcppSkeleton", lib="/Library/Frameworks/R.framework/Resources/library/", NULL, type='source', INSTALL_opts=c("--debug"))



# system time
unclass(Sys.time())


# PDF AWESOMENESS


# Get results
#pdf(file = "conus_cors_all_vs_all_v3.pdf", width=15, height=15)
doPDFs = TRUE
if (doPDFs == TRUE)
	{
	pdffn = "conus_cors_all_vs_all_v3.pdf"
	pdf(file=pdffn, width=15, height=15, paper="USr")
	}


# Turn off the PDF and open
if (doPDFs == TRUE)
	{
	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	}




# PDF awesomeness

# Setup the PDF filename
pdffn = "conus_cors_all_vs_all_v3.pdf"

# Open the PDF to write to
pdf(file=pdffn, width=15, height=15)

# If you just want standard quartz window on-screen
#quartz(title="Hey look at this title", width=15, height=15)

# Turn off the PDF / close quartz
dev.off()

# Open the PDF automatically from R
# make the system command-line command
cmdstr = paste("open ", pdffn, sep="")

# Execute the command in the system
# (i.e. "open filename.pdf"
system(cmdstr)





# Transparent colors via:
rgb(1,0,0, alpha=0.5)




# gdata: upperTriangle
# Extract or replace the upper/lower triangular portion of a matrix
tmpphyd = upperTriangle(phyd[1:numsp_to_get_G, 1:numsp_to_get_G], diag=FALSE)



x=c(0,0)
y=c(0,1)
d=c(1, 1)
junkdata = as.data.frame(cbind(x, y, d))
coordinates(junkdata) = ~x+y
(junkdistskm = spDists(junkdata, longlat=TRUE))

     
summ <- function(x)
	{
	cat("\n")
	cat("SUMM(): PROVIDING OVERALL SUMMARY OF OBJECT...\n")
	cat("\nCLASS OF OBJECT:\n")
	print(class(x))

	cat("\nDIMENSIONS OF OBJECT:\n")
	print(dim(x))

	cat("\nLENGTH OF OBJECT:\n")
	print(length(x))

	cat("\nATTRIBUTES OF OBJECT:\n")
	print(attributes(x))

	cat("\nSUMMARY() OF OBJECT:\n")
	print(summary(x))
	}




titles on outside of subplots
multiplots
multiple subplots
outer margins
oma

# logEMSY by bin
abS$logEMSY = log(abS$EMSY)

plot.new()
thing_tohist_str = "logEMSY"
plot_hists_of_rates_by_binsizes(abS, list_binsizes, thing_tohist_str)
titletxt = "logEMSY by binsize (fraction out of mean diversity extinct in bin,\ndivided by binsize in my)"
title(titletxt, outer=TRUE)


use outer=TRUE



R: 
size(array)


cumulative sums, etc
cumsum(1:10)
cumprod(1:10)
cummin(c(3:1, 2:0, 4:2))
cummax(c(3:1, 2:0, 4:2))


BEST R PAGE EVER
http://stackoverflow.com/questions/1189759/expert-r-users-whats-in-your-rprofile



wait: 
Sys.sleep(1)



prod -- multiply values in a list

The prod function is similar to sum, but does products rather than sums. Prod-
ucts can often over
ow or under
ow (a suburb of Circle 1)|taking logs and
doing sums is generally a more stable computation.




disappearing attributes
Most coercion functions strip the attributes from the object. For example, the
result of:
as.numeric(xmat)
will not be a matrix. A command that does the coercion but keeps the attributes
is:
storage.mode(xmat) <- 'numeric'


# Install from source rather than binary
# (rgl doesn't install at the moment, though...)
install.packages(c("rgl"), lib="/Library/Frameworks/R.framework/Resources/library/", contriburl=contrib.url(getOption("repos"), 'source'), type='source', dependencies=TRUE)

install.packages(c("Rcpp"), lib="/Library/Frameworks/R.framework/Resources/library/", contriburl=contrib.url(getOption("repos"), 'binary'), type='binary', dependencies=TRUE)

install.packages(c("Rcpp"), lib="/Library/Frameworks/R.framework/Resources/library/", contriburl=contrib.url("http://R-Forge.R-project.org", 'source'), type='source', dependencies=FALSE)

install.packages(c("Rcpp"), lib="/Library/Frameworks/R.framework/Resources/library/", contriburl=contrib.url("http://R-Forge.R-project.org", 'source'), type='source', dependencies=FALSE)




# Library, then un-library, a package:
library(sp)
detach("package:sp")

# Find where a package is installed:
find.package(package="cladoRcpp")

# Check the version:
packageVersion("cladoRcpp")

# Remove packages installed in different places:
remove.packages(pkgs="Rcpp", lib="/Library/Frameworks/R.framework/Resources/library/")
remove.packages(pkgs="Rcpp", lib="/Users/nickm/Library/R/2.10/library/")

# See a list of all installed packages (allegedly):
installed.packages()






find.package(package="cladoRcpp")
remove.packages(pkgs="cladoRcpp")
remove.packages(pkgs="cladoRcpp", lib="/Library/Frameworks/R.framework/Versions/3.4/Resources/library/")
remove.packages(pkgs="cladoRcpp", lib="/Library/Frameworks/R.framework/Resources/library/")

detach("package:cladoRcpp")

install.packages(pkgs="/drives/Dropbox/_njm/__packages/cladoRcpp_0.15.tar.gz", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')
library(cladoRcpp)




uninstall
reload
detach


detach("package:sp")


detach("package:sp")
library(zoo)

detach("package:zoo")
library(zoo)

# Install from gzip, or directory

# install.packages("bigmemory")  # old version, 3.0
# to remove it:
# 

# detach(unload = TRUE)

# remove.packages("bigmemory", lib="/Library/Frameworks/R.framework/Resources/library/")
# remove.packages("bigmemory", lib="/Users/nickm/Library/R/2.10/library/")
# remove.packages("bigmemory")
# remove.packages("bigalgebra")
# remove.packages("biganalytics")
# remove.packages("bigtabulate")
installed.packages()


These are some tools that help you find/remove packages:

# Library, then un-library, a package:
library(sp)
detach("package:sp")

# Find where a package is installed:
find.package(package="cladoRcpp")

# Check the version:
packageVersion("cladoRcpp")

# Remove packages installed in different places:
remove.packages(pkgs="Rcpp", lib="/Library/Frameworks/R.framework/Resources/library/")
remove.packages(pkgs="Rcpp", lib="/Users/nickm/Library/R/2.10/library/")

# See a list of all installed packages (allegedly at least):
installed.packages()



########################################
https://gist.github.com/wiltonkiss/54aa277a71206ed1a6b91f4c00919ee9
########################################
#' Remove everything except functions
#' 
#' Remove everything from the global environment except functions
#' @export 

rm.except.fn <- function(){
        rm(list = setdiff(ls(envir = .GlobalEnv), lsf.str(envir = .GlobalEnv)), envir = .GlobalEnv)
}

#' Remove only functions
#' 
#' Remove functions from the global environment
#' @export 

rm.fn <- function(){
        rm(list=lsf.str(envir = .GlobalEnv), envir = .GlobalEnv)
}


#' Remove everything except...
#' 
#' Remove everything from the global environment except named element. No default.
#' @param x Character vector of elements to be kept. No default
#' @export 

rm.all.except <- function(x){
        x <- c("rm.all.except",x)
        rm(list = ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv) %in% x], envir = .GlobalEnv)
}
########################################


detach(unload = TRUE)







http://r.789695.n4.nabble.com/Modify-base-R-functions-in-Rprofile-site-td905121.html


# Find out where the function lives.
(  env <- as.environment( 'package:phylobase' ))

  # Crack the binding.
  unlockBinding( '.phylobase', env )

  # Replace the function.
  assignInNamespace( 'parse', function(...){
 
   # Your function here, or an object that contains it.

   }, ns = 'base' )

  # Relock the binding.
  lockBinding( 'parse', env ) 











# These compiled without a problem from downloaded zipfiles
#install.packages("/bioinformatics/R/bigmemory_4.2.7.tar.gz", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')

#install.packages("/bioinformatics/R/biganalytics_1.0.14.tar.gz", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')

#install.packages("/bioinformatics/R/bigtabulate_1.0.14.tar.gz", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')

###################
# Installing bigalgebra is tougher:
# install.packages("bigalgebra")
#
# 1. Had to copy the "bigmemory" and "boost" libraries from the unzipped bigmemory package
# "include" directory, into the bigalgebra "src" directory, to get compile to work
#
# 2. compiled from local modified directory
#install.packages("/bioinformatics/R/bigalgebra", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')
#
#################################################################

library(bigmemory)
library(bigalgebra)
library(biganalytics)
library(bigtabulate)




#### We first load the twalk package:
rm(list=ls(all=TRUE)) # remove previous definitions

library(Rtwalk)






R:
rep() (replicate)
seq() (sequence of items)
# series
# sequence
# repeat



# Plot PDFs
x = seq(-2, 2, 0.1)
y = dnorm(x, mean=0, sd=2)
plot(x, y, pch="")
lines(x, y)
lines(x, 0.9*y)
legend()

R: evaluate a string without printing
eval(parse(text = {a string}))



eval() - evaluate a string

# Why are you using eval?  The following is equivalent:
for(name in names(x)) {
  y <- x[[name]][1]
}


repr(a) - represent a as a string



in gdata:
strsplit
substr -- subset of a string
trim -- remove outside whitespace

R: 
jitter -- add a little noise to data


The presence of UseMethod indicates this is a generic function. To see what methods are
available we can use methods()



Converting data.frame character to numeric
http://stackoverflow.com/questions/2288485/dataframe-coverting-column-type-to-numeric

as.data.frame()
NOT
data.frame()
!!!!













# Get a tree structure in the same order as the prt structure:


ordered_nodenames = get_nodenums(tree_to_chainsaw)
parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, tree_to_chainsaw$edge[,2])



odd
even
in: gtools



get slots

getSlots(class(mymap))



# Geography



library(laser)
library(ape)
library(dismo)
library(maptools)
library(sp)

wd = "/Users/nickm/Desktop/__projects/_2011-05-01_Lagrange_fossil_horsies_SVP/_data/"
setwd(wd)

fn = "NAbioboundaries.shp"

shapefile(fn)

sf = readShapePoly(fn)
spplot(sf)


data(wrld_simpl)


coordinates(wrld_simpl)
coordinates(sf)

proj4string(wrld_simpl)
proj4string(sf)


mymap<-readOGR(fn, "NAbioboundaries")
#"Mth030607.shp", layer="Mth030607")


plot(wrld_simpl)
plot(mymap)

mymap2 = spTransform(mymap, CRS(proj4string(wrld_simpl)))
plot(wrld_simpl)

attr(mymap2, "bbox") = attr(wrld_simpl, "bbox")
plot(mymap2, new=FALSE)

# convert to unprojected
unprojected_CRS = CRS("+proj=longlat +datum=WGS84")
mymap3 = spTransform(mymap2, unprojected_CRS)

# write to KML
writeOGR(mymap3, "NAbioboundaries.kml", "NAbioboundaries", driver="KML")


# display what's in the slots
getSlots(class(mymap))
attr(mymap, "data")
attr(mymap, "polygons")
attr(mymap, "plotOrder")
attr(mymap, "bbox")
attr(mymap, "proj4string")


# get the polygons
slot(mymap2, "polygons")







sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

sourcedir = '/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))

# Dear Nick,
# Log-likelihoods are typically not normalized. If you look at the density of the normal 
# distribution, and use this to calculate the likelihood of a really-great fitting model, 
# you'll see that the positive log-likelihoods are not uncommon:
#  
# Here's a tight-fitting model example

y = rnorm(10, 1:10, 0.000001)
x = rnorm(10, 1:10, 0.000001)



# Remove the space between axis and border with par like this:
# set the inside box ("i") ahead of time!!  -- removes the 4% extension
# this friggin' solution took an hour to find,
# solution here: http://tolstoy.newcastle.edu.au/R/help/06/08/32529.html
par(xaxs = "i")
par(yaxs = "i")

plot(grd_polygon)




par(mfrow=c(1, 2))
hist(x)
hist(y)
plot(x, y)

################################
# manual log-likelihood
################################
# fit linear model
lm1 = lm(y ~ x)
summ(lm1)

# this model fits really, really well!

# dnorm gets the density (the PDF value) when the PDF is a normal distribution
# (normally distributed errors)

# sapply gets the density for each residual from the model prediction
# (this assumes the errors are normally distributed)
(PDF_density_vals_for_data = sapply(lm1$resid, dnorm, mean = 0, sd = sd(lm(y ~ x)$resid), log = F))

(log_PDF_density_vals_for_data = log(PDF_density_vals_for_data))

# total log-likelihood for the model =
sum(log_PDF_density_vals_for_data)

# which is a very high log-likelihood!

# proper log-likelihood
(lm1)
(tmp_loglik = logLik(lm1))
summ(tmp_loglik)

# extract the values
attr(tmp_loglik, "df")
attr(tmp_loglik, "nall")
as.numeric(tmp_loglik)


# and then here's a model closer to what you might find in real life:
y = rnorm(10, 1:10, 1)
x = rnorm(10, 1:10, 1)
plot(x, y)
lm2 = lm(y ~ x)
(lm2)

# extracting values
# R-squared: simple way
summary(lm2)$r.squared

# R-squared: Here's the difficult way where "model" is the lm output:
(cor(predict(lm2), y))^2
 

logLik(lm2) # proper log-likelihood













Flatten a list of lists into 1 list of the atomic elements of the combined lists

fast:
unlist

slow:
stack(t2$clades)
(but has unstack)




unlist, then put in square matrix

rotations_fn = "/GIS_paleo/ODSN/ODSN_rotations.txt"
rotations = read_table_good(fn=rotations_fn, tmpskip=12)
# printall(rotations)
# head(rotations)

# Edits made to rotations file to fix minor errors (see notes)

# Parse the rotations file
platecodes = rotations$plate_wrt
platecodes2 = strsplit(platecodes, split="-")
platecodes3 = matrix(unlist(platecodes2), ncol=2, byrow=TRUE)
numcodes = rotations$numeric_code
numcodes2 = strsplit(numcodes, split="-")
numcodes3 = matrix(unlist(numcodes2), ncol=2, byrow=TRUE)
newcodes = cbind(platecodes3[,1], numcodes3[,1], platecodes3[,2], numcodes3[,2])

newcodes = paste(platecodes3[,1], numcodes3[,1], platecodes3[,2], numcodes3[,2], sep="-")
newcodes2 = strsplit(newcodes, split="-")
newcodes3 = matrix(unlist(newcodes2), nrow=4, byrow=TRUE)














F Duan wrote:
> Dear R people,
> 
> I am using par(mfrow=c()) to plot multi-figures in the same window. And I like 
> to put a common title (and xlab, ylab) for all of plots. I have already left 
> some margin by resetting omi values in par() and hided all (xlab, ylab) for 
> each sub-plot. Could anyone tell me how to do that?
> 
> Thanks a lot,
> 
> Frank

see ?title

argument outer is used to place titles in the outer margin.

HTH






> clades.from.polytomies
function (tree) 
{
    from <- tree$edge[, 1]
    to <- tree$edge[, 2]
    n.taxa <- length(tree$tip.label)
    is.node <- seq_len(max(tree$edge)) %in% from
    edge.counts <- tapply(to, from, length)
    poly.nodes <- as.integer(names(edge.counts[edge.counts > 
        2]))
    ans <- lapply(poly.nodes, ancestors2, tree)
    ans1 <- mapply(setdiff, ans, poly.nodes, SIMPLIFY = FALSE)
    clades <- lapply(ans1[!(poly.nodes %in% unlist(ans1))], function(x) x[x <= 
        n.taxa])
    clades.repr <- tree$tip.label[sapply(clades, "[", 1)]
    names(clades) <- clades.repr
    clades.spp <- lapply(clades, function(x) tree$tip.label[x])
    clades.drop <- sort(unlist(lapply(clades, "[", -1)))
    tree2 <- drop.tip.fixed(tree, clades.drop)
    make.clade.tree(tree2, clades.spp)
}
<environment: namespace:diversitree>
> make.clade.tree
function (tree, clades) 
{
    if (!identical(class(tree), "phylo")) 
        stop("tree must be a plain 'phylo' tree")
    if (!all(names(clades) %in% tree$tip.label)) 
        stop("Unknown clade representatives")
    if (!all(sapply(clades, is.character))) 
        stop("'clades' must be a list of character vectors")
    if (any(duplicated(unlist(clades)))) 
        stop("Duplicated species names")
    tree$clades <- clades
    class(tree) <- c("clade.tree", class(tree))
    tree
}

















############################################################################################
############################################################################################
############################################################################################
############################################################################################
#################################################################
# Compare all of the methods
#################################################################
sourcedir = "/Dropbox/_njm/"
source8 = '_matrix_utils_v1.R'
source(paste(sourcedir, source8, sep=""))


# By hand
#################################################################
# Make sparse matrices with nrow nonzero values (default)
matrix_dim = 100	# fast at 1000, malloc at 10000
num_to_make = 1000
tmp_newvals_str = "newvals=0.5"
tmp_numvals_to_change = matrix_dim

t_list = rep(1, length(A_list))

A_list = make_lotsa_sparse_matrices(num_to_make, matrix_dim, newvals_str=tmp_newvals_str,  numvals_to_change=matrix_dim)

Qmat = A_list[[1]]
t = t_list[[1]]


# works
z = try( (Pmat = expm::expm(Qmat*t, method='Pade')) )

# won't work
z = try( (Pmat = expm::expm(Qmat*t, method='R_Eigen')) )
# Error in eigen(x, sym = isSym) : infinite or missing values in 'x'

z
# [1] "Error in eigen(x, sym = isSym) : infinite or missing values in 'x'\n"
# attr(,"class")
# [1] "try-error"

class(z)
# [1] "try-error"

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################




array subsetting:
x = array([[1,2,3], [4,5,6], [7,8,9]], dtype=float)
compress([1,2], x, axis=1)

result:
array([[ 1.,  2.],
       [ 4.,  5.],
       [ 7.,  8.]])



# alter diagonals
from LR_run_functions_v2 import *
from numpy import diagonal, NaN
x = array([[1,2,3], [4,5,6], [7,8,9]], dtype=float)

#Set diagonals:
i=range(0,3)
i
x[i,i]
x[i,i]=0

x[i,i+1]

nick@mws2[phylocom]|38> x
                   <38> 
array([[ 0.,  2.,  3.],
       [ 4.,  0.,  6.],
       [ 7.,  8.,  0.]])



DIAGONALS WITH INDEXING 

# “Fancy” indexing also works.

>>> i = [0,1,2]

>>> a[i,i]

array([11, 22, 33]) 

# Indexing can also be used

# to set diagonal values

>>> a[i,i] = 2

>>> i = array([0,1])

# upper diagonal

>>> a[i,i+1] = 1

# lower diagonal

>>> a[i+1,i]= = -1

>>> a

array([[ 2,  1, 13],

       [-1,  2,  1],

       [31, -1,  2]])

triu(a)	Triangular, upper
tril(a)	Triangular, lower
http://www.cfa.harvard.edu/~jbattat/computer/python/science/idl-numpy.html





x = array([[1,2,3], [4,5,6], [7,8,9]], dtype=float)
square_array =x
dims = shape(square_array)

# Make a mask of the upper-right triangle, including diagonal
mask1 = tri(dims[0], dims[1], k=-1, dtype=bool) == False
#mask1 = reshape(mask1, dims[0]*dims[1])
mask1

c = square_array
c[mask1] = NaN






Color maps, heat maps

# Molecular tree
# heat map
# scale to 0-1, times 130, round down, +1, = scale from 20 to 120
tmp_colors = rev(heat.colors(130, alpha=1))
rescaled_colors = 20+floor( 99*((rates_for_edges-min(rates_for_edges)) / (max(rates_for_edges)-min(rates_for_edges))) )
br.col = tmp_colors[rescaled_colors]
plot(tr, edge.color=br.col, edge.width=8)

# blue-red
# scale to 0-1, times 100, round down, +1, = scale from 1 to 100
tmp_colors = rev(rainbow(100, start=0, end=4.5/6, alpha=1))
rescaled_colors = 1+floor( 99*((rates_for_edges-min(rates_for_edges)) / (max(rates_for_edges)-min(rates_for_edges))) )
br.col = tmp_colors[rescaled_colors]
plot(tr, edge.color=br.col, edge.width=8)






##########################################################
# Install bigmemory and related packages in R!!
# downloaded from:
# https://r-forge.r-project.org/R/?group_id=556
##########################################################
# install.packages("bigmemory")  # old version, 3.0
# to remove it:
# remove.packages("bigmemory", lib="/Library/Frameworks/R.framework/Resources/library/")
# uninstall package
# detach(unload = TRUE)


###################
# Installing bigalgebra is tougher:
# install.packages("bigalgebra")
#
# 1. Had to copy the "bigmemory" and "boost" libraries from the unzipped bigmemory package
# "include" directory, into the bigalgebra "src" directory, to get compile to work
#
# 2. compiled from local directory
install.packages("/bioinformatics/R/bigalgebra", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')
#
# These compiled without a problem 
install.packages("/bioinformatics/R/bigmemory_4.2.7.tar.gz", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')

install.packages("/bioinformatics/R/biganalytics_1.0.14.tar.gz", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')

install.packages("/bioinformatics/R/bigtabulate_1.0.14.tar.gz", "/Library/Frameworks/R.framework/Resources/library/", NULL, type='source')

library(bigmemory)
library(bigalgebra)
library(biganalytics)
library(bigtabulate)







> save("f.mean.Rdata",file="D:/Users/Ays/Documents/Results")
> it saves but when I open the file in notepad it is just some characters
> meaningless.

The purpose of save is to save R objects to a file so you can retrieve
them later or on another system via the load command.
pickle

On 4/02/2010, at 12:09 PM, mkna005 mkna005 wrote: 

> Hello all! 
> I was wondering if there is a way to pickle an R object into a file 
> like it is possible in python? Such as you have an complicated R 
> object(not a dataframe) , you use a function to write it to a file and 
> than you have a function where you can retrieve the object from that 
> file later on. 

?dput 
?dget 



string split
strsplit
on whitespace:
l = "Node	   Fix [Mod]	  Min"
strsplit(l, "[ \t]+")[[1]]


# on individual whitespace
numlines_str = "   75091 combined_75000_trees.trees"
strsplit(numlines_str, "[[:blank:]]")[[1]]


# This should be a faster list2 string
list2str_fast <- function(list1, spacer="")
	{
	# convert to character
	# split into list of lists of characters
	# merge lists with unlist
	# paste with collapse argument
	tmpstr = paste(unlist(strsplit(as.character(list1), split="")), collapse="")
	return(tmpstr)
	}


# Split a column into 2 columns
splits_left1 = strsplit(splits$split, split="\\|")
splits_left2 = unlist(splits_left1)
splits_left3 = matrix(data=splits_left2, ncol=2, byrow=TRUE)
head(splits_left3)




# replace whitespace or at least space
gsub("[: :]", ",", tntstr2)

gsub("[ \t]+", ",", tntstr2)


# Convert any spaces to underscores
dates_df$OTUname = gsub("\\+", "_", dates_df$OTUname)



#reverse:
rev


# Remove/delete files

# Wipe these files
system(paste("rm ", lgfn, sep=""))
system(paste("rm ", outfn, sep=""))
system(paste("rm ", keyfn, sep=""))
system(paste("rm ", splitsfn, sep=""))
system(paste("rm ", statesfn, sep=""))



R control statements if else for
break
nextncol(results_tables2$taxon_diversity_counts)




R:
break -- break out of innermost loop
next -- go to next iteration of loop (like python continue)
#(like pass in python)

R: 
rev: reverse

find rows with parts of items matching text:
GREPL
grepl
sp_hits = grepl("sp.", df$species)




methods(print)
methods(plot)




sprintf("%02.6f", 8.448334e-01)

grepl("\\d", tmpstr)
grep("\\d", tmpstr)
regexpr("\\d", tmpstr)
gregexpr("\\d", tmpstr)


great little online Regex builder at http://txt2re.com/
regular expression



EVEN AWESOMER


merge_words_nonwords <- function(words, nonwords)
	{
	if (length(nonwords) == ((length(words) + 1)))
		{
		words = c(words, "")
		}
	
	wordsmat = cbind(nonwords, words)
	paste1 = apply(X=wordsmat, MARGIN=1, FUN=paste, sep="", collapse="")
	paste1
	
	sentence = paste(paste1, sep="", collapse="")
	return(sentence)
	}



str = '1-$s-$y'
mstr = paste(rownames(BioGeoBEARS_model_object@params_table), sep="", collapse="|")
mstr

regmatches(x=str, m=gregexpr(mstr, str))

str = '1-s-y-mx01/mx01j-j'
tmpwords = c("d", "e", "a", "b", "x", "u", "j", "ys", "y", "s", "v", "mx01", "mx01j", "mx01y", "mx01s", "mx01v", "mx01r", "dp")
mstr = paste(tmpwords, sep="", collapse="|")

mstr

words=regmatches(x=str, m=gregexpr(mstr, str), invert=FALSE)[[1]]
words
nonwords=regmatches(x=str, m=gregexpr(mstr, str), invert=TRUE)[[1]]
nonwords

merge_words_nonwords(words, nonwords)


str = 's-y-mx01/mx01j-j'
tmpwords = c("d", "e", "a", "b", "x", "u", "j", "ys", "y", "s", "v", "mx01", "mx01j", "mx01y", "mx01s", "mx01v", "mx01r", "dp")
mstr = paste(tmpwords, sep="", collapse="|")

mstr

words=regmatches(x=str, m=gregexpr(mstr, str), invert=FALSE)[[1]]
words
nonwords=regmatches(x=str, m=gregexpr(mstr, str), invert=TRUE)[[1]]
nonwords
merge_words_nonwords(words, nonwords)


str = 'ys'
tmpwords = c("d", "e", "a", "b", "x", "u", "j", "ys", "y", "s", "v", "mx01", "mx01j", "mx01y", "mx01s", "mx01v", "mx01r", "dp")
mstr = paste(tmpwords, sep="", collapse="|")

mstr

words=regmatches(x=str, m=gregexpr(mstr, str), invert=FALSE)[[1]]
words
nonwords=regmatches(x=str, m=gregexpr(mstr, str), invert=TRUE)[[1]]
nonwords

merge_words_nonwords(words, nonwords)



str = 'ysys'
tmpwords = c("d", "e", "a", "b", "x", "u", "j", "ys", "y", "s", "v", "mx01", "mx01j", "mx01y", "mx01s", "mx01v", "mx01r", "dp")
mstr = paste(tmpwords, sep="", collapse="|")

mstr

words=regmatches(x=str, m=gregexpr(mstr, str), invert=FALSE)[[1]]
words
nonwords=regmatches(x=str, m=gregexpr(mstr, str), invert=TRUE)[[1]]
nonwords

merge_words_nonwords(words, nonwords)



# Remove everything up to the first =
# (e.g. tree TREE1 =  
#  in NEXUS files...)
tree <- gsub("^.*= *", "", tree)


# remove everything between square brackets
gsub("\\[[^]]*\\]", "", X[s])


'(tree)'	# Word 1
re2='( )'	# White Space 1
re3='(STATE)'	# Word 2
re4='(_)'	# Any Single Character 1
re5='(\\d+)'

'(tree)( )(STATE)(_)(\\d+)'

'\\(tree\\)\\( \\)\\(STATE\\)\\(_\\)\\(\\d+\\)'

$txt=',74a';

$re1='(,)';	# Any Single Character 1
$re2='(\\d+)';	# Integer Number 1
$re3='((?:[a-z][a-z0-9_]*))';	# Variable Name 1

$re3='((?:[a-z][a-z0-9_]*))';	# Variable Name 1




# http://txt2re.com/index-python.php3?s=tree%20STATE_2409000%20[%26lnP=-1684.3463984799218,lnP=-69316.05058319442]%20=%20[%26R]%20&-15&-27&-11&-155&7&-28&1&-29&-142&-30&-12&-31

txt='tree STATE_2409000 [&lnP=-1684.3463984799218,lnP=-69316.05058319442] = [&R] '

tmpre = NULL
ri = 0
tmpre[[(ri=ri+1)]] = '(tree)'	# Word 1
tmpre[[(ri=ri+1)]] = '( )'	# White Space 1
tmpre[[(ri=ri+1)]] = '(STATE)'	# Word 2
tmpre[[(ri=ri+1)]] = '(_)'	# Any Single Character 1
tmpre[[(ri=ri+1)]] = '(\\d+)'	# Integer Number 1
tmpre[[(ri=ri+1)]] = '( )'	# White Space 2
tmpre[[(ri=ri+1)]] = '(\\[.*?\\])'	# Square Braces 1
tmpre[[(ri=ri+1)]] = '( )'	# White Space 3
tmpre[[(ri=ri+1)]] = '(=)'	# Any Single Character 2
tmpre[[(ri=ri+1)]] = '( )'	# White Space 4
tmpre[[(ri=ri+1)]] = '(\\[&R\\])'	# Square Braces 2
tmpre[[(ri=ri+1)]] = '( )'	# White Space 5



tmpre_string = list2str(tmpre, spacer="")
(tmpre_string)



# Extract just the numbers, with (characters) (numbers) pattern
gregexpr("(?:\\w+(\\d+))+", tmpstr)
grep("(?:\\w+(\\d+))+", tmpstr)
grepl("(?:\\w+(\\d+))+", tmpstr)
regexpr("(?:\\w+(\\d+))+", tmpstr)

matches = gregexpr("(?:\\w+(\\d+))+", tmpstr)[[1]]
matches_end = matches-1+attr(matches,"match.length")
x = mapply(substr, tmpstr, matches, matches_end)

# pull out the numbers / extract the numbers / extract numbers
# just the numbers
tmpstr = "190.1Ma - 65Ma"
matches = gregexpr("(?:(\\d+))+", tmpstr)[[1]]
matches_end = matches-1+attr(matches,"match.length")
x = mapply(substr, tmpstr, matches, matches_end)
x





# pull out the numbers / extract the numbers / extract numbers
# just the numbers INCLUDING THOSE CONNECTED BY DECIMAL POINTS!!!!
# (but not negative symbols)
tmpstr = "190.1Ma - 65Ma"
matches = gregexpr("(?:([0-9\\.]+))+", tmpstr)[[1]]
matches_end = matches-1+attr(matches,"match.length")
x = mapply(substr, tmpstr, matches, matches_end)
x


#numbers,including decimals
tmpstr = "190.1Ma - 65Ma"
matches = gregexpr("(?:(\\d*\\.\\d+))+", tmpstr)[[1]]
matches_end = matches-1+attr(matches,"match.length")
x = mapply(substr, tmpstr, matches, matches_end)
x


#numbers,including decimals
tmpstr = "190.1Ma - 65Ma"
matches = gregexpr("([+-]?\\d\\*\\.\\d+)(?![-+0-9\\.])", tmpstr)[[1]]
matches_end = matches-1+attr(matches,"match.length")
x = mapply(substr, tmpstr, matches, matches_end)
x


([+-]?\\d*\\.\\d+)(?![-+0-9\\.])

\\d*\\.\\d+


# greedy, match everythign between 1st and last parens
z2 = "2(0 1)002???91011101001002100511?200?1??01?11?0?30030100000100?1??1???10201?10?1?10?1030101201000?10?1???10?1?1?1?????0???1??0000201?"
gregexpr('(\\(.*\\))', z2)


# http://r.789695.n4.nabble.com/Regular-expression-to-define-contents-between-parentheses-td848782.html
# The problem is that regular expressions are greedy,
# so you were matching everything between the first and
# last parens, as you noticed.  Putting the question
# mark there makes it a "minimal" matching operation.
# Apparently this is only implemented in perl regex's,
# or at least in that syntax.  Hence the 'perl=TRUE'. 

# find e.g. all (0 1)
gregexpr('(\\(.*?\\))', paste(z2, z2, sep=""))
gregexpr('(\\[.*?\\])', paste(z2, z2, sep=""), perl=TRUE)[[1]]

# find e.g. all [0 1]
z3 = gsub("\\(", "[", z2)
z3 = gsub("\\)", "]", z3)
gregexpr('(\\[(.*?)\\])', paste(z3, z3, sep=""), perl=TRUE)[[1]]









# WORKS
    
    tmpre = NULL
	ri = 0
	tmpre[[(ri=ri+1)]] = '(tree)'	# Word 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 1
	tmpre[[(ri=ri+1)]] = '(STATE)'	# Word 2
	tmpre[[(ri=ri+1)]] = '(_)'	# Any Single Character 1
	tmpre[[(ri=ri+1)]] = '(\\d+)'	# Integer Number 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 2
	tmpre[[(ri=ri+1)]] = '(\\[.*?\\])'	# Square Braces 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 3
	tmpre[[(ri=ri+1)]] = '(=)'	# Any Single Character 2
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 4
	tmpre[[(ri=ri+1)]] = '(\\[&R\\])'	# Square Braces 2
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 5
	
	tmpre_string = list2str(tmpre, spacer="")
	(tmpre_string)

    # Remove the header   
    X3 <- gsub(tmpre_string, "", X2)
    











# WORKS

    edges <- gsub("\\[(&[[:alnum:]_=%!.,{}-]+)\\]", "", X)
    edges <- unlist(strsplit(gsub("[()]*", "", edges), ":"))
    tips <- gsub("^[0-9]+.[0-9]+,*", "", edges)


	re1='(,)'										# Any Single Character 1
	#re2='(?:[a-z][a-z]*[0-9]+[a-z0-9]*)'		# Alphanum 1
	re2='(?:[a-zA-Z0-9]*)'		# Alphanum 1
    regstr = paste(re1, re2, sep="")
    
    gsub(regstr, "", "asd,234g")
    gsub(regstr, "", edges)











rev -- reverse











# SPAM: Sparse matrices
install.packages("spam")
library(spam)
vignette(package="spam")
demo(package="spam")

Use ‘demo(package = .packages(all.available = TRUE))’
to list the demos in all *available* packages.





###########################################

loose stuff:



paste list2str

paste(x, collapse=”,”) is very useful. Would you have a solution for the following



xaxt
nice tick marks
nice histogram bins
pretty



Axis for a generic interface.
axTicks returns the axis tick locations corresponding to at=NULL; pretty is more flexible for computing pretty tick coordinates and does not depend on (nor adapt to) the coordinate system in use.



# Blank plot
# Plot, no borders (bty="n"), no labels (xlab, ylab), no tick marks (xaxt, yaxt)
plot(1:10, 1:10, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

# Defaults
xaxt="s", yaxt="s"




# get tips
not_tips = mytree$edge[,2] %in% mytree$edge[,1]
tips = list_of_not_tips == FALSE

nodenum_edge1 = mytree$edge[not_tips,][,1] == edge1
nodenum_edge1_name = mytree$node.label[nodenum_edge1]

e
nodelabels("root", 116)
nodelabels("root", 115)



apply:
# apply a mini-function to columns of a data matrix
uniq_states = apply(x, 2, function(x) length(unique(x)))



#convert dataframe columns to factor:
apply(char_dtf, 2, as.factor)


# convert data to phyDat
d = apply(char_dtf, 2, paste)
d2 = as.data.frame(d)
tmplevels = as.character(sort(unique(unlist(d2))))
as.phyDat(d2, type="USER", levels=tmplevels)


# convert columns to strings
colnames_to_use = c(clad_chars_names)
char_dtf = all_params[, colnames_to_use]

# apply over overs
apply(char_dtf, 1, paste, collapse="")


margins:

par(oma=c(1,1,3,1), mfrow=c(2,3))

# set outer margins & number of subplots
par(oma=c(1,1,3,1), mfrow=c(2,3))
# then do outer title
mtext("observed and null distances between DNA tree and params tree", outer=TRUE)




SUBPLOTS WITH NO MARGINS:

# Subplots	
par(mfrow=c(3,3))
# Subplot margins: c(bottom, left, top, right) 
par(mar=c(0, 0, 0, 0))
# outer margins
par(oma=c(4, 4, 4, 2))	


Master global Title for graph with multiple plots
mtext("Densities", outer = TRUE, cex = 1.5)



oma
A vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text.


# Get just the files that end in .R
codefiles = list.files(path=".", pattern=".R$")


# Graphics fun...
source("http://research.stowers-institute.org/efg/R/Graphics/Basics/mar-oma/mar-oma.R")

plot(0:10, 0:10, type="n", xlab="X", ylab="Y")

text(5,5, ID, col="red", cex=size1)

box("plot", col="red")
mtext("Figure", SOUTH<-1, line=3, adj=1.0, cex=size2, col="blue")
box("figure", col="blue")



nf <- layout(matrix(c(2,0,1,3), 2, 2, byrow=TRUE), widths=c(3,1), heights=c(1,3), TRUE)
layout.show(nf)

par(mar=c(3,3,1,1))
plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
par(mar=c(0,3,1,1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
par(mar=c(3,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)

# find the edges existing at i:
TF_edge_exists_at_i = (TF_edge_top_above_i + TF_edge_bottom_below_i == 2)

print(paste("#edges existing = ", sum(TF_edge_exists_at_i), "", sep=""))

# indices of the edges (same diff)
indices_of_edges_existing = which(TF_edge_exists_at_i == TRUE)
print(paste("#edges existing = ", length(indices_of_edges_existing), "", sep=""))





# Lotsa histograms

#A vector of the form c(bottom, left, top, right)
par(mfrow=c(4,6), mar=c(2,1,3,1), oma=c(7,3,6,2))





# padj = 0 means bottom-justified
# padj = 1 means top-justified

####################################
# Histograms of each main continuous parameter
####################################
pdffn = "FigureCD_histograms_params_lnparams.pdf"
pdf(file=pdffn, width=10, height=7, paper="USr")

# Good settings for having a lot of histograms
par(mfrow=c(4,6), mar=c(2,1,3,1), oma=c(7,3,6,2))

###################
# Draw histograms
for (colnum in 1:ncol(params))
	{
	hist(params[,colnum], main=names(params)[colnum], xlab="", ylab="", breaks=10)
	}

# Add title
titletxt = paste("Histograms: params", sep="")
mtext(titletxt, side=3, outer=TRUE, line=1, padj=0, cex=1.5)

captiontxt = paste("Figure B: Histograms of the variable continuous shell-growth parameters.", sep="")
mtext(captiontxt, side=1, outer=TRUE, adj=0, padj=1, line=4, cex=1.0)	
# End histograms
###################

###################
# Draw histograms
for (colnum in 1:ncol(ln_params))
	{
	hist(ln_params[,colnum], main=names(ln_params)[colnum], xlab="", ylab="", breaks=10)
	}

# Add title
titletxt = paste("Histograms: ln_params", sep="")
mtext(titletxt, side=3, outer=TRUE, line=1, padj=0, cex=1.5)

captiontxt = paste("Figure B: Histograms of the natural log transformation of the parameters. If a parameter\ncontained some zero values, then before transformation, those values were reassigned to be 10% of the\nvalue of the minimum nonzero value for that parameter.", sep="")
mtext(captiontxt, side=1, outer=TRUE, adj=0, padj=1, line=4, cex=1.0)	
# End histograms
###################















# Plot leftwards-facing, diagonal tree
#plot(mytree, type = "c", use.edge.length = FALSE, direction = "leftwards")


# Plot unrooted tree
#plot(mytree, type = "u", font = 1, x.lim=c(0, 0.246), y.lim=c(1, 14),  no.margin = TRUE, lab4ut 
#= "horizontal")

#plot(mytree, type = "u", font = 1, x.lim=c(1, 0.123), y.lim=c(1, 0.246),  no.margin = TRUE, 
#lab4ut = "horizontal")






par(xpd = NA, mar = c(5, 4, 4, 5))   #, omd = c(0,0.9,0,1))


# start the plot
xdim = c(-10,130)
ydim = c(-10,130)
plot(xdim, ydim, ty="n", asp=1, xlab = "x", ylab = "y")

# and use colors; heat is one set of color options
color.range <- rainbow(length(wo$data))
# pch=21 means solid circle
# col = colors for corresponding points
points(cbind(wo$coords, wo$data), pch=21, col=color.range, bg=color.range)



zr<- range(data)
tempcolor = rainbow(128)
firstmatch = match(color.range[1], tempcolor)
lastmatch = match(color.range[length(color.range)], tempcolor)
tempcolor = tempcolor[firstmatch:lastmatch]

image.plot(legend.only=TRUE, zlim=zr, nlevel=length(tempcolor), col=tempcolor)





###################################
# Objects in R.
###################################
# Most things are just lists:
# Introduction to data types and objects in R 
# http://www.nealgroothuis.name/?page_id=67
# 
# R Programming/Introduction
# From Wikibooks, the open-content textbooks collection
# http://en.wikibooks.org/wiki/R_Programming/Introduction
#
# 

# Build a linear model
y<-c(1,2,4)
x<-c(1,2,3)
foo<-lm(y~x)

# See all of the ways of displaying what is in there
class(foo)
names(foo)

summary(foo)

# Names and class:
attributes(foo)
attr(foo, "name")
attr(foo, "class")
# set class:
attr(foo, "class") <- "gregion"

# give the structure of your data
str(foo)
str(head(toccs))

# Nah:
# slotNames(foo)
# getSlots(foo)
# showMethods(foo)
# methods(foo)


dimnames(toccs)
# works only for data.frames
# this doesn't really work at least in R.app
View(head(toccs))



Usage

setClass(Class, representation, prototype, contains=character(),
         validity, access, where, version, sealed, package)

removeClass(Class, where)

isClass(Class, formal=TRUE, where)

getClasses(where, inherits = missing(where))

findClass(Class, where, unique = "")

resetClass(Class, classDef, where)

sealClass(Class, where)





R classes

x=setClass("ds", contains = "data.frame",
representation(
	clade_name = "character", 
	taxon_level = "character",
	sp_included = "logical",
	region = "character",
	countries_list = "character",
	oldest_desired_age = "numeric",
	youngest_desired_age = "numeric",
	make_genus_table = "logical",
	run_minmax_hist = "numeric",
	occsfn = "character",  				# raw occurrences filename
	rangesfn = "character", 			# raw ranges filename
	all_mammal_names_fn = "character" 	# filename containing the list of all mammal species
	))


y=new("ds", subset_ranges)






LearnBayes
ML phylogeny R
network R


GARP etc.



# install from source:
install.packages(repos=NULL, pkgs='/Users/nickm/Desktop/downloads/phangorn_1.0-0.tgz', type='source')


# R 2.10.1 downloaded April 16 2010

Many binaries installed:
/private/tmp/RtmpNrNMfs/downloaded_packages


install.packages("LearnBayes")
library(LearnBayes)

library(ape)
library(phangorn)

# search for packages
search()

# detach a package
detach("package:multicore")


pr(data) doesn't cover:

attributes(data)

# specific attribute:
attr(data)



get all available demos:
demo()
example()


e.g.:
demo(cokriging, package="gstat", ask=FALSE)


see R Package tutorial (dowloaded PDF)


The Slots in an Object from a Formal Class
Description
These functions return or set information about the individual slots in an object.
Usage
object@name
object@name <- value

slot(object, name)
slot(object, name, check = TRUE) <- value

slotNames(x)
getSlots(x)


Get or set specific attributes of an object.
Usage
attr(x, which, exact = FALSE)

Description
These functions access an object's attributes. The first form below returns the object's attribute list. The replacement forms uses the list on the right-hand side of the assignment as the object's attributes (if appropriate).
Usage
attributes(obj)
attributes(obj) <- value
mostattributes(obj) <- value
Arguments
obj
an object
value
an appropriate named list of attributes, or NULL.



grep, grepl, regexpr and gregexpr search for matches to argument pattern within a character vector: they differ in the format of and amount of detail in the results.

sub and gsub perform replacement of the first and all matches respectively.


Matrix exponentiation:

https://stat.ethz.ch/pipermail/r-sig-phylo/2009-June/000393.html

expm <- function(...) expm:::expm(..., method="Ward77")




# Write strings to a file
make a 1-column matrix of strings, then write.table to file



# list free datasets
data()
# load simple world map
data(wrld_simpl)
# 
plot(wrld_simpl, xlim=c(-130,10), ylim=c(-60,60))

Data sets in package ‘maptools’:

h1pl (gpcholes)                            
h2pl (gpcholes)                            
wrld_simpl                                 Simplified world country polygons


# example issue with list of objects:

library(ape)
library(phangorn)
example(NJ)

# Jukes-Cantor (starting tree from NJ)  
fitJC1 <- pml(tree, Laurasiatherian)  

# optimize edge length parameter     
fitJC2 <- optim.pml(fitJC1)
fitJC2 
  
# search for a better tree using NNI rearrangements     
fitJC3 <- optim.pml(fitJC2, optNni=TRUE)
fitJC3

# Now, the function SH.test can allegedly take "objects of class 'pml' separated by commas, [or] a list containing such objects".  Since I'm going to have hundreds of these fits, I'd like to make submit a list of them to SH.test, something like this:

#list_of_fits = c(fitJC1, fitJC2, fitJC3)
list_of_fits = list(fitJC1, fitJC2, fitJC3)
SH.test(list_of_fits, B=100)

# ...but "list of fits" is something weird, e.g. 
(list_of_fits[1])

#...does not return the same thing as:
(fitJC1)

SOLUTION:
nextComp = length(list_of_lists)+1
list_of_lists[[nextComp]] = fitJC1


# New commands
subset


# returns things matching the condition
which(is.na)


# get a specific thing out of an object
slot

# plot on top, not new plot
plot(new=FALSE)

# try
args(biogeomancer)
b = try( biogeomancer('Peru', locality=lonzero$locality[3], progress='') )
b

# vignette

# Search just the R website

colnames

# bring up a table
fix

biogeomancer
# We recommend using a tool like BioGeomancer: http://bg.berkeley.edu/
# latest (Guralnick et al., 2006) to georeference textual locality descriptions. An
# important feature of BioGeomancer is that it attempts to capture the uncer-
# tainty associated with each georeference (Wieczorek et al., 2004). The dismo
# package has a function biogeomancer that you can use for this, and that we
# demonstrate below, but its use is generally not recommended because you really
# need a detailed map interface for accurate georeferencing.
# Here is an example for one of the records with longitude = 0. We put the
# biogeomacer function into a 'try' function, to assure elegant error handling if
# the computer is not connected to the Internet.
# Guralnick, R.P., J. Wieczorek, R. Beaman, R.J. Hijmans and the BioGeo-
# mancer Working Group, 2006. BioGeomancer: Automated georeferenc-
# ing to map the world's biodiversity data. PLoS Biology 4: 1908-1909.
# http://dx.doi.org/10.1371/journal.pbio.0040381



Google Maps:
gmap <- function(x, maptype='terrain', pch='x', cex=1, col='red') {
      require(RgoogleMaps)
      xr <- range(x$lon)
      yr <- range(x$lat)
      mykey <- "ABQIAAAAx_Zq0CG7Dz9YNSzDR0PYtxT2yXp_ZAY8_ufC3CFXhHIE1NvwkxQ9M3z-hbUeB-0ItTVP2WPiFXA8PA"
      gm <- GetMap.bbox(key=mykey, lonR=xr, latR=yr, maptype=maptype)
      tmp <- PlotOnStaticMap(gm, lon=gb$lon, lat=gb$lat, pch=pch, cex=cex, col=col, verbose=0)
}

gmap(gb, maptype='terrain')











data(singer, package = "lattice")

## using traditional graphics

singer.split <- with(singer, split(height, voice.part))
par(mfrow = c(2, 4))
for (i in names(singer.split))
    hist(singer.split[[i]], main = i, xlab = "height")

## using lattice

library(lattice)
histogram(~height | voice.part, singer)

## using ggplot2

library(ggplot2)
qplot(height, data = singer, geom = "histogram", facets = voice.part ~ .)

qplot(height, data = singer, geom = "points", facets = voice.part ~ .)



# Get all coordinates
# number of polygons
numpolys = length(wrld_simpl@polygons)
allcoords = NULL
for (i in 1:numpolys)
	{
	x = wrld_simpl@polygons[[i]]@Polygons
	for (j in 1:length(x))
		{
		cat("Polygon ", wrld_simpl@polygons[[i]]@ID, ": ", i, ".", j, " dims:\n", sep="")
		print(dim(x[[j]]@coords))
		allcoords = rbind(allcoords, x[[j]]@coords)
		}
	}
(head(allcoords))
(dim(allcoords))





	# T-test to see if empirical (hypothesized) value is different from the population mean
	# http://www.r-tutor.com/elementary-statistics/hypothesis-testing/two-tailed-test-population-mean-unknown-variance
	xbar = meanval				# sample mean 
	mu0 = empirical_val			# hypothesized value 
	s = sdval					# sample standard deviation 
	n = length(allvals)			# sample size 
	t = (xbar-mu0) / (s/sqrt(n)) 
	cat("t-statistic = ", t, "\n", sep="")		# test statistic 
	
	# get the test statistic for half-p-value
	alpha = .05
	t.half.alpha = qt(1-alpha/2, df=n-1)	
	cat("half-alphas: ", -t.half.alpha, ", ", t.half.alpha, "\n", sep="")
	
	pval = 2 * pt(t, df=n-1)			# lower tail 
	cat("pval = ", pval, sep="", "\n")	# two-tailed p-value 

	pval = 2 * pt(-t.half.alpha, df=n-1)			# lower tail 
	cat("pval = ", pval, sep="", "\n")	# two-tailed p-value 

	pval = 2 * pt(t.half.alpha, df=n-1)			# lower tail 
	cat("pval = ", pval, sep="", "\n")	# two-tailed p-value 



# Pretty numbering for plot axes...
pretty(allvals)



# ========================================================
# Fix axes on histograms:
# ========================================================

		# Hist of PDs
		f1 = paste(dir, prefix, '', geog_regions[i], '_emp_distmat.txt', sep="")
		emp_dists = read.table(f1)
		barbreaks <- seq(0, 800, by=50)
		
		# Remove space between line and actual 0 in histogram
		# set the inside box ("i") ahead of time!!  -- removes the 4% extension
		# this friggin' solution took an hour to find,
		# solution here: http://tolstoy.newcastle.edu.au/R/help/06/08/32529.html
		par(xaxs = "i")
		par(yaxs = "i") 
		h1 <- hist(emp_dists[,1], breaks=barbreaks, plot=FALSE)

		# label plots individually
		if (i==1)
			{
			plot(h1, xlim=c(0,800), ylim=c(0, 1.2*max(h1$counts)), col="darkgreen", main=geog_region_names[i], xlab="", ylab="count", bty="o", xaxt="n")

			# x-axis with tick-marks but no numbers
			axis(1, at=seq(0, 800, by=200), label=rep("", length(seq(0, 800, by=200))), cex.axis=3, tick=TRUE, tcl=ticklength, pos=0)			
			}
		else
			{
			plot(h1, xlim=c(0,800), ylim=c(0, 1.2*max(h1$counts)), col="darkgreen", main=geog_region_names[i], xlab="", ylab="", xaxt="n", bty="o")			
			axis(1, at=seq(0, 800, by=200), label=rep("", length(seq(0, 800, by=200))), cex.axis=3, tick=TRUE, tcl=ticklength, pos=0)
			}
			
		
		# Lines of average counts & 95% confidence intervals
		lines(h1$mids, simmeans[i,], type="l", lty=1, xpd=FALSE)
		lines(h1$mids, simmeans[i,] + 2*simsds[i,], type="l", lty=3)
		lines(h1$mids, simmeans[i,] - 2*simsds[i,], type="l", lty=3)

	
		# Vertical lines
		# Empirical
		abline(v=means_of_means[1,i], lty=2, lwd=3, col="darkgreen") 	# lty=2 (Line TYpe = dashed)
		# Simulated
		abline(v=means_of_means[3,i], lty=2, lwd=2, col="black") 	# lty=2 (Line TYpe = dashed)
		abline(v=means_of_means[3,i]+2*means_of_means[4,i], lty=3, lwd=1, col="black") 	# lty=2 (Line TYpe = dashed)
		abline(v=means_of_means[3,i]-2*means_of_means[4,i], lty=3, lwd=1, col="black") 	# lty=2 (Line TYpe = dashed)




str()
R Documentation
Compactly Display the Structure of an Arbitrary R Object




print plots after making/storing them:

pl1 <- spplot(z["zn.se"], main="log-zinc std.err.")
pl2 <- spplot(z["cu.se"], main="log-copper std.err.")
pl3 <- spplot(z["cd.se"], main="log-cadmium std.err.")
pl4 <- spplot(z["pb.se"], main="log-lead st.err.")

print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))


quantiles
percentiles, CDF, PDF
percentile_025 = aggregate(percent_extinct, by=list(timebin_sizes), FUN=quantile, probs=0.025, na.rm=TRUE)
	
percentile_975 = quantile(randK_sims_params, probs=0.975, na.rm=TRUE)

percentile_025 = apply(randK_sims_params, 2, quantile, probs=0.025, na.rm=TRUE)
percentile_025 = apply(randK_sims_params, 2, quantile, probs=0.975, na.rm=TRUE)



max.col
colmax
maxcol

colSums (x, na.rm = FALSE, dims = 1)
rowSums (x, na.rm = FALSE, dims = 1)
colMeans(x, na.rm = FALSE, dims = 1)
rowMeans(x, na.rm = FALSE, dims = 1)



title for overall plot with subplots

mtext allows this if outer = TRUE.
main text = mtext
mtext("observed and null distances between DNA tree and params tree", outer=TRUE)





adj
The value of adj determines the way in which text strings are justified in text, mtext and title. A value of 0 produces left-justified text, 0.5 (the default) centered text and 1 






http://bridgewater.wordpress.com/2010/12/21/my-favorite-r-packages-installed-with-one-command/

My favorite R packages (installed with one command)

Posted on December 21, 2010. Filed under: Uncategorized | Tags: rstats |

I just started a new job (working on social search awesomeness at Bing) and so I had to set up my “dev” environment with all of my usual tools (R, python,vim,etc). One thing that made this a bit easier is my habit of keeping an R script around that installs all of my common packages for me in one shot.

This is also a nice way to share your list of favorite packages with your friends. Please feel free to share your list of cool packages.

# R packages I use commonly: 12/21/2010 twitter: drbridgewater

#list all packages that I commonly use
p<-c()

#literate programming
p<-c(p,"R2HTML") #plus utils::sweave
p<-c(p,"Rpad")

#plotting
p<-c(p,"ggplot2")
p<-c(p,"YaleToolkit")

#data mining/machine learning/NLP
p<-c(p,"ElemStatLearn")
p<-c(p,"gbm")
p<-c(p,"bayesm")
p<-c(p,"RWeka")
p<-c(p,"lsa")
p<-c(p,"tm")

#graphs and networks
p<-c(p,"igraph")

#statistics
p<-c(p,"survival")
p<-c(p,"Hmisc")
p<-c(p,"ICSNP")
p<-c(p,"zipfR")

repositories<-c("http://cran.cnr.Berkeley.edu","http://www.stats.ox.ac.uk/pub/RWin")
install_package<-function(pack,repositories)
{
if(!(pack %in% row.names(installed.packages())))
{
update.packages(repos=repositories, ask=F)
install.packages(pack, repos=repositories, dependencies=T)
}
require(pack,character.only=TRUE)
}

for( pack in p)
{
install_package(pack,repositories)
}



kickass bar plots:


par(oma=c(1,1,1,1), mar=c(15.6, 4.1, 4.1, 7.1), mfrow = c(1,1))
par(xpd=TRUE)
tmp_colors = c("white", "cyan", "yellow", "green2", "orange", "red", "brown")
tmp_xnames = rep("", length(colnames(aic_weights)))
legend_text = rev(c("Brownian", "O-U", "lambda", "kappa", "delta", "early burst", "white noise"))
tmp_legend = list(x=36, y=0.9, pt.cex=4)
bar_x_positions = barplot(aic_weights, names.arg=tmp_xnames, col=tmp_colors, ylab="AIC weight", xlab="shell growth parameter", legend.text=legend_text, args.legend=tmp_legend)


## Create plot and get bar midpoints in 'mp'
#bar_x_positions = barplot(1:length(colnames(aic_weights)))

## Set up x axis with tick marks alone
#axis(1, at = bar_x_positions, labels = FALSE)

# Create some text labels
tmp_x_labels = colnames(aic_weights)

# Plot x axis labels at mp
text(x=bar_x_positions-0.1, y=-0.015, adj=0, srt=315, labels = tmp_x_labels, xpd=TRUE)

# Plot the best models across the top
text(x=bar_x_positions-0.1, y=1.015, adj=0, srt=45, labels = best_models, xpd=TRUE)

toptxt = "Best model under AIC:"
text(x=1, y=1.12, labels=toptxt)

caption = "Figure SX: Akaike weights of 7 models for the evolution of shell parameters.  The models compared are: \nBrownian motion (BM, brown); the Ornstein-Uhlenbeck (OU, red) stabilizing selection model; \nthe lambda model (orange) which rescales internal branch lengths by a linear fraction; the \nkappa model (green) which rescales each branch length by a power equal to the kappa parameter, \nand which becomes a speciational model as kappa approaches 0; delta (yellow) which focuses \nchange towards the base or tips; early burst (EB, cyan) which has an initial high rate of change \nthat then declines; and white noise (white), where observations are produced by a normal distribution \nwith no tree structure, which represents the situation of no phylogenetic signal. Brownian motion (BM) \nhas the highest AIC weight for 63% of the shell parameters 15/24.  White noise (no phylogenetic signal) \nis superior for 25% (6/24) shell parameters."

text(x=-2, y=-0.45, labels=caption, adj=0)


















#######################################################
# master_graphic_fromscratch_v02.R
# Hominin chart idea
#######################################################
# 1. Top: phylogeny and fossils
# 2. Middle: previous estimates by year (y-axis) and histogram
# 3. Bottom: phylogeny estimate(s) for that node
#######################################################

library(gdata)
library(ape)

library(phyloch) 	# for 95% HPDbars -- http://www.christophheibl.de/Rpackages.html,
# https://stat.ethz.ch/pipermail/r-sig-phylo/2010-November/000843.html

sourcedir = '/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))

sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))


# Set working directory
wd = "/Users/nickm/Desktop/__projects/_hominin_phylo_dating/_graphics"
setwd(wd)



#######################################################
# 1. Get the (consensus) phylogeny
#######################################################

# Load the consensus tree
connexus_FigTree_fn = "/Users/nickm/Desktop/__projects/_hominin_phylo_dating/BEAST_link_gene_morph_trees/v06_BEAST_fix_morph_monophyly_operators_longer/v06_morph_timetrees.nexus_FigTree_big_v02_noperiods.con"

# contr = read.nexus(connexus_FigTree_fn)


# Add a nice long root for display purposes
out_trfn = "/Users/nickm/Desktop/__projects/_hominin_phylo_dating/BEAST_link_gene_morph_trees/v06_BEAST_fix_morph_monophyly_operators_longer/v06_morph_timetrees.nexus_FigTree_big_v02_noperiods_wroot.con"
out_trfn = add_root_to_nexus_trees(trfn=connexus_FigTree_fn, out_trfn=out_trfn, root_brlen=1.5)
moref(out_trfn)


contr = read.beast_original(out_trfn)
condtf = read_beast_prt(connexus_FigTree_fn)


# Fix tip labels
new_tiplabels = contr$tip.label
new_tiplabels = gsub(pattern="\\.", replacement="_", x=new_tiplabels)
new_tiplabels = gsub(pattern="'", replacement="", x=new_tiplabels)
new_tiplabels
contr$tip.label = new_tiplabels



# Plot 3 subplots

# Subplot formatting
# http://www.statmethods.net/advgraphs/layout.html
# The layout( ) function has the form layout(mat) where
# mat is a matrix object specifying the location of the N figures to plot.
# One figure in row 1 and two figures in row 2
attach(mtcars)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
hist(wt)
hist(mpg)
hist(disp)

# Subplot margins: c(bottom, left, top, right) 
par(mar=c(1, 4, 1, 1))
# outer margins
#par(oma=c(5, 4, 4, 2))	
par(oma=c(3,0,3,0))	

# set the inside box ("i") ahead of time!!  -- removes the 4% extension
# this friggin' solution took an hour to find,
# solution here: http://tolstoy.newcastle.edu.au/R/help/06/08/32529.html
par(xaxs = "i")
par(yaxs = "i") 




# Plot the default tree
#plot(contr)
#nodelabels()

# Rotate the branches
tmptr = rotate(phy=contr, node=38)
tmptr = rotate(phy=tmptr, node=40)
tmptr = rotate(phy=tmptr, node=43)
tmptr = rotate(phy=tmptr, node=44)
tmptr = rotate(phy=tmptr, node=46)
tmptr = rotate(phy=tmptr, node=47)
tmptr = rotate(phy=tmptr, node=48)
tmptr = rotate(phy=tmptr, node=49)
tmptr = rotate(phy=tmptr, node=59)
tmptr = rotate(phy=tmptr, node=65)
tmptr = rotate(phy=tmptr, node=67)


# Get the node coordinates
coords = plot_phylo3_nodecoords(tmptr, root.edge=TRUE)
ylims = c(min(pretty(coords$yy)), max(pretty(coords$yy)))
ylims = ylims + c(-1, 1)

plot.new()
layout(mat=matrix(data=c(1,1,1,1,2,3,4), nrow=7, ncol=1, byrow=TRUE))


# plot(tmptr, root.edge=TRUE, y.lim=ylims) # with ylims
plot(tmptr, root.edge=TRUE)

# How to do colors in text, numeric, hex, and transparent colors
# last 2 digits of hex indicate transparency
# http://research.stowers-institute.org/efg/R/Color/Chart/
# https://stat.ethz.ch/pipermail/r-help/2007-October/142934.html
colorname = "blue"
colnums = col2rgb(colorname)
colhex = rgb(t(colnums/255))
colhex_transparent = paste(colhex, "50", sep="")

HPDbars(phy=tmptr, col=colhex_transparent, lwd = 5)
#nodelabels()
axisPhylo()
























trstr = "((((1:61.9979,2:61.9979):10.0994,3:72.0973):30.5957,(4:97.17,5:97.17):5.52294):165.913,(((((6:45.8837,7:42.7899):90.6427,8:11.748):9.81895,((9:71.9983,10:30.1075):59.5189,(((11:48.3447,12:48.3447):77.8739,(13:86.8531,(14:53.665,15:53.665):33.1881):39.3656):26.0816,((16:94.7537,17:94.7537):44.1509,(18:82.3045,(19:72.9737,20:72.9737):9.33074):56.6002):13.3956):27.9983):10.886):22.824,21:214.009):34.3554,(22:186.877,23:186.877):61.4866):20.2417);"


###################################
# Original node ordering (APE "cladewise")
###################################
library(ape)
tr = read.tree(file="", text=trstr)
plot(tr)
nodelabels()


###################################
# Postorder version 1 (APE "pruningwise")
###################################
tr2 = reorder(tr, order="pruningwise")
plot(tr2)
nodelabels()

internal_nodenums_in_postorder = unique(tr2$edge[,1])
internal_nodenums_in_postorder = internal_nodenums_in_postorder - length(tr2$tip.label)

tr3 = tr2
tr3$node.label[internal_nodenums_in_postorder] = 1:tr3$Nnode
plot(tr3)
nodelabels(tr3$node.label)


###################################
# Postorder version 2 (phylobase "postorder")
###################################
library(phylobase)
tr4 = as(tr2, "phylo4")
tr5 = reorder(tr4, "postorder")

tmp_edge = attr(tr5, "edge")
tmp_edge2 = tmp_edge[tmp_edge[,2] > length(tr2$tip.label), ]
tmp_edge2[,2]

nodenums_to_change = tmp_edge2[,2]-length(tr2$tip.label)
tr2$node.label[nodenums_to_change] = 1:tr2$Nnode
plot(tr2)
nodelabels(tr2$node.label)









Successful reprojection of a grid to match another grid:

#===================================================
# Starter source material
#===================================================
sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

sourcedir = '/Dropbox/_njm/'
source3 = '_genericR_v1.R'
source(paste(sourcedir, source3, sep=""))

source2 = 'pbdbtools_v2.R'
source(paste(sourcedir, source2, sep=""))

source1 = 'pdb_utils_v6.R'
source(paste(sourcedir, source1, sep=""))

source8 = '_biogeog_sim_utils_v1.R'
source(paste(sourcedir, source8, sep=""))

# For latlong
source5 = '_phylogeostats_v1.R'
source(paste(sourcedir, source5, sep=""))
#===================================================


#===================================================
# Run with these commands
#===================================================
library(ape)
library(phylobase)
library(gdata) # needed for trim (whitespace trimming) function
library(lattice) # for histogram

# for mapping
library(dismo)
library(maptools)
library(sp)
library(rgdal) 	# for readOGR


# for loading ASCII grid
library(adehabitat)
library(maptools)


##########################################
# Memory improvement (requires some of the above things for some reason)
##########################################
# based on comment in:
# https://r-forge.r-project.org/forum/message.php?msg_id=4025&group_id=294
#setOptions(chunksize = 1e+04, maxmemory = 1e+06)
setOptions(chunksize = 1e+05, maxmemory = 30e+06)

# http://r-sig-geo.2731867.n2.nabble.com/change-projection-of-large-raster-file-td6499121.html
beginCluster(type="SOCK")



#wd = "/Dropbox/_njm/BioGeoR"
wd = "/Users/nickm/Desktop/__2011-05-25_Hawaii/_code"
setwd(wd)



# Reproject the full-resolution grid to Behrmann Cylindrical Equal-Area projection

##################################
# Load DEM data -- as basemap
# source: http://www.prism.oregonstate.edu/pub/prism/pacisl/
##################################

# Coarsening a grid:
# https://stat.ethz.ch/pipermail/r-sig-geo/2011-May/011798.html
# output.dim
# The number of rows and columns to return in the created object 
# using GDAL's method to take care of image decimation / replication;
# presently ordered (y,x) - this may change
# for dem, default:
xdim = 10800
ydim = 21600

# Crashes at /1 or /2...up to /5, it looks like...
xdim_new = floor(xdim/6)
ydim_new = floor(ydim/6)

xdim_new_dem = xdim_new
ydim_new_dem = ydim_new

#fn1 = "hi_dem_15s.asc"
# ETOPO1: 1-minute elevation dataset, top of ice sheet
# http://www.ngdc.noaa.gov/mgg/global/global.html
fn1 = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/ETOPO1/ETOPO1_Ice_c_geotiff.tif"

# Projection
metadata = '
Geographic:

	Latitude Resolution:
		0.01666666667 
	Longitude Resolution:
		0.01666666667 
	Geographic Coordinate Units:
		Decimal degrees 

Geodetic Model:

	Horizontal Datum Name:
		WGS 84 
	Ellipsoid Name:
		WGS 84 
	Semi-major Axis:
		6378137.000000 
	Denominator of Flattening Ratio:
		298.257224 
'


# full resolution
# dem_orig = readGDAL(fn1)

# subsampling -- here, load the full dataset
dem = readGDAL(fn1, output.dim=c(ydim_new, xdim_new))

# Apply the Geographic (latlong) projection
proj4string(dem)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
demproj = CRS(proj4string(dem))
demproj
dembbox = attr(dem, "bbox")
dembbox

	
	

###########################################
# Load the environmental dataset
###########################################

# Bio-ORACLE Environmental data
# figured out projection here:
# http://groups.google.com/group/oracle_ugent/browse_thread/thread/d9652df4799cac75

#######################
tmptxt='> headf(fn2)
Read 1574 items
ncols         3709
nrows         1568
xllcorner     -17367529.33951
yllcorner     -7342230.1364987
cellsize      9365.622736
NODATA_value  -9999
-9999 -9999 -9999 -9999 -9999'
#######################

xdim = 1568
ydim = 3709
xdim_new = floor(xdim/2)
ydim_new = floor(ydim/2)

xdim_new_env = xdim_new
ydim_new_env = ydim_new
#
fn2 = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/Bio-ORACLE/BioOracle_EAgrid_9090RV/calcite.asc"

# full resolution
# dem_orig = readGDAL(fn1)

# subsampling (no subsampling here)
env = readGDAL(fn2, output.dim=c(ydim_new, xdim_new))

# Projection from:
# http://sites.google.com/site/spatialr/crsprojections
proj4string(env)=CRS("+proj=cea +lat_ts=30 +lat_ts=0 +ellps=WGS84 +datum=WGS84")
attr(env, "bbox")
envproj = CRS(proj4string(env))
envproj
envbbox = attr(env, "bbox")
envbbox




###########################################
# Projecting from Behrmann to other projections is apparently weird
# (Lame!!)
#
# So we will project the elevation grid to Berhmann
###########################################
demr = raster(dem)
envr = raster(env)

#ex <- projectExtent(demr, envproj) # reprojects the extent
demr_reproj <- projectRaster(from=demr, to=envr) # to project the actual GLC-data
dim(demr_reproj)

attr(demr, "extent")
attr(env, "bbox")
attr(demr_reproj, "extent")
image(demr_reproj)

wd = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/ETOPO1/"
setwd(wd)

save(demr_reproj, file="demr_reproj.Rnw")
save(env, file="env.Rnw")



# Get the latlong coords
dim(env)

coords_xy = adf(coordinates(env))
dim(coords_xy)
x = coords_xy$x
y = coords_xy$y
xy = cbind(x,y)
dim(xy)

wd = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/ETOPO1/"
setwd(wd)
#save(demr_reproj, file="demr_reproj.Rnw")
save(xy, file="env_xy.Rnw")


# Read in the resampled DEM file
load("demr_reproj.Rnw")
dim(demr_reproj)
attr(demr_reproj, "extent")
attr(env, "bbox")

# Extract the values from the elevation grid
elev_vals = extract(demr_reproj, xy)



#demr_reproj_fn = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/ETOPO1/ETOPO1_Ice_c_geotiff_reproj_to_Behrmann.asc"

#demr_reproj = readGDAL(demr_reproj_fn)
#demr_reproj_SpatialGrid = readGDAL(demr_reproj_fn)

#demr_reproj = raster(demr_reproj_SpatialGrid)
wd = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/ETOPO1/"
setwd(wd)
#save(demr_reproj, file="demr_reproj.Rnw")
save(xy, file="env_xy.Rnw")



wd = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/ETOPO1/"
setwd(wd)

load("demr_reproj.Rnw")
load("env_xy.Rnw")



# Extract the values from the elevation grid
elev_vals = extract(demr_reproj, xy)

dim(elev_vals)
save(elev_vals, file="elev_vals.Rnw")

length(elev_vals)
dim(elev_vals)


xdim = 1568
ydim = 3709
xdim_new = floor(xdim/2)
ydim_new = floor(ydim/2)

elev_vals2 = matrix(data=elev_vals, nrow=ydim_new, ncol=xdim_new, byrow=TRUE)
#dim(elev_vals) = c(ydim_new, xdim_new)
dim(elev_vals2)

attr(env, "grid")
env_bbox = adf(attr(env, "bbox"))
#demr_new = raster(x=elev_vals2)
demr_new = raster(x=elev_vals2, xmn=env_bbox$min[1], xmx=env_bbox$max[1], ymn=env_bbox$min[2], ymx=env_bbox$max[2], crs=envproj)
dim(demr_new)
image(demr_new)

save(demr_new, file="demr_new.Rnw")

demr_reproj_fn = "/Users/nickm/Desktop/__projects/___phylokriging_example/_data/ETOPO1/ETOPO1_Ice_c_geotiff_reproj_to_Behrmann_sample_to_env.asc"
demr_new_SpatialGrid = as(demr_reproj, "SpatialGridDataFrame")
writeAsciiGrid(demr_new_SpatialGrid, demr_reproj_fn) 










# Label below the x-axis
mtext(text="Time (Ma)", side=1, line=2)

blankraster_world = raster(nrows=numcells_lat, ncols=numcells_long, xmn=tmpext2@xmin, xmx=tmpext2@xmax, ymn=tmpext2@ymin, ymx=tmpext2@ymax, crs=envproj)










# A typical tree string

((((Aotearoa_magna_hw0050:57.55665505,Zearchaea_sp_hw0051:57.55665505):17.1867611,Mesarchaea_bell_hw0040_20A:74.74341615):27.71428065,(Chilarchaea_quellon_hw0029:51.72488185,Mecysmauchenius_segmentatus_hw0098:51.72488185):50.73281496):167.2431139,((Huttonia_sp_hwHutt41:177.3053354,Palpimanus_sp_hw0082:177.3053354):88.1331157,(Stenochilus_sp_hw0081:239.5191569,(Patarchaea_muralis:46.25785221,((((Myrmecarchaea_sp:52.52267389,Baltarchaea_conica:52.53323845):12.82754307,Archaea_paradoxa:65.31673034):32.81723581,Afrarchaea_grimaldii:53.14529283):38.20828846,(((Austrarchaea_nodosa_hw0074:49.40542201,Austrarchaea_davisae_hwAu066:49.40542201):26.12254178,Austrarchaea_mainae_hw0075:75.5279638):90.82062321,((Eriauchenius_workmani_hw0006:71.64058191,Eriauchenius_bourgini_hw0061:71.64058191):69.96785536,((Afrarchaea_pilgrimsrest_hw59Af3:49.13614886,Afrarchaea_woodae_hw57Af1:49.13614886):73.44416036,(NewGenus_legendrei_hw0014:82.4140833,(NewGenus_lavatenda_hw0003AB:57.30071564,NewGenus_jeanneli_hw0057:57.30071564):25.11336767):40.1662259):19.02812805):24.74014973):16.51404697):30.99101197):25.66551091):25.91929425):4.262359573);


# Add a root edge
library(ape)

trstr = "((((Aotearoa_magna_hw0050:57.55665505,Zearchaea_sp_hw0051:57.55665505):17.1867611,Mesarchaea_bell_hw0040_20A:74.74341615):27.71428065,(Chilarchaea_quellon_hw0029:51.72488185,Mecysmauchenius_segmentatus_hw0098:51.72488185):50.73281496):167.2431139,((Huttonia_sp_hwHutt41:177.3053354,Palpimanus_sp_hw0082:177.3053354):88.1331157,(Stenochilus_sp_hw0081:239.5191569,(Patarchaea_muralis:46.25785221,((((Myrmecarchaea_sp:52.52267389,Baltarchaea_conica:52.53323845):12.82754307,Archaea_paradoxa:65.31673034):32.81723581,Afrarchaea_grimaldii:53.14529283):38.20828846,(((Austrarchaea_nodosa_hw0074:49.40542201,Austrarchaea_davisae_hwAu066:49.40542201):26.12254178,Austrarchaea_mainae_hw0075:75.5279638):90.82062321,((Eriauchenius_workmani_hw0006:71.64058191,Eriauchenius_bourgini_hw0061:71.64058191):69.96785536,((Afrarchaea_pilgrimsrest_hw59Af3:49.13614886,Afrarchaea_woodae_hw57Af1:49.13614886):73.44416036,(NewGenus_legendrei_hw0014:82.4140833,(NewGenus_lavatenda_hw0003AB:57.30071564,NewGenus_jeanneli_hw0057:57.30071564):25.11336767):40.1662259):19.02812805):24.74014973):16.51404697):30.99101197):25.66551091):25.91929425):4.262359573)ingroup:200;"

t = read.tree(file="", text=trstr)

plot(t, root.edge=TRUE)


# Add an outgroup (doesn't work)

library(ape)

t = read.tree(file="", text="((((Aotearoa_magna_hw0050:57.55665505,Zearchaea_sp_hw0051:57.55665505):17.1867611,Mesarchaea_bell_hw0040_20A:74.74341615):27.71428065,(Chilarchaea_quellon_hw0029:51.72488185,Mecysmauchenius_segmentatus_hw0098:51.72488185):50.73281496):167.2431139,((Huttonia_sp_hwHutt41:177.3053354,Palpimanus_sp_hw0082:177.3053354):88.1331157,(Stenochilus_sp_hw0081:239.5191569,(Patarchaea_muralis:46.25785221,((((Myrmecarchaea_sp:52.52267389,Baltarchaea_conica:52.53323845):12.82754307,Archaea_paradoxa:65.31673034):32.81723581,Afrarchaea_grimaldii:53.14529283):38.20828846,(((Austrarchaea_nodosa_hw0074:49.40542201,Austrarchaea_davisae_hwAu066:49.40542201):26.12254178,Austrarchaea_mainae_hw0075:75.5279638):90.82062321,((Eriauchenius_workmani_hw0006:71.64058191,Eriauchenius_bourgini_hw0061:71.64058191):69.96785536,((Afrarchaea_pilgrimsrest_hw59Af3:49.13614886,Afrarchaea_woodae_hw57Af1:49.13614886):73.44416036,(NewGenus_legendrei_hw0014:82.4140833,(NewGenus_lavatenda_hw0003AB:57.30071564,NewGenus_jeanneli_hw0057:57.30071564):25.11336767):40.1662259):19.02812805):24.74014973):16.51404697):30.99101197):25.66551091):25.91929425):4.262359573)ingroup:100,outgroup:100;")

plot(t, root.edge=TRUE)








# Generic handy R functions
sourcedir = '/Dropbox/_njm/'
source2 = '_genericR_v1.R'
source(paste(sourcedir, source2, sep=""))

sourcedir = '/Dropbox/_njm/'
source3 = '_R_tree_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

sourcedir = '/Dropbox/_njm/'
source3 = '_R_XMLtools_v1.R'
source(paste(sourcedir, source3, sep=""))


# R functions for dealing with XML
library(XML)

# R functions for dealing with trees
library(gtools)		# These are Greg's miscellaneous functions; includes "trim"
					# Formerly package ('gregmisc')
# 

# Get the directories, file names etc.

# Rdata files
XMLdata_dir = "/Users/nickm/Desktop/__projects/2012-03-14_Markos_Alexandrou_BEAST/05_BEAST_v1/"

setwd(XMLdata_dir)





################################################
# Extract character matrix from XML file
################################################
sourcedir = '/Dropbox/_njm/'
source3 = '_R_xml_BEAST_functions_v1.R'
source(paste(sourcedir, source3, sep=""))

fn = "BDincFossilOTUsr8s1.xml"
dirfn = paste(XMLdata_dir, fn, sep="")


charmat_df2 = extract_character_matrix_from_BEAST_XML(dirfn)
#head(charmat_df2)
#names(charmat_df2)

partitions_df = charmat_df_to_alignment_df(charmat_df2)
#head(partitions_df)
#names(partitions_df)

nexus_df2 = partitions_df_to_nexus_df(partitions_df)
row.names(nexus_df2)
dim(nexus_df2)


charmatrix = nexus_df_to_charmatrix(nexus_df2)
row.names(charmatrix)
dim(charmatrix)

# Write to NEXUS file
tmp = row.names(charmatrix)
maxlengths = unlist(lapply(tmp, nchar))
maxlength = max(maxlengths) +2

charmatrix_to_nexus_file(charmatrix, nexus_outfn="nexus_fromBEAST_v2.nexus", namewidth=maxlength, ambiguity_check=FALSE)







# results of print to variable

# Print the first numlines part of the xml node list to screen
# hp = headprint
hp <- function(item, numlines=15)
	{
	txt = capture.output(print(xmlnodelist))
	
	cat("\n")
	cat("headxml(): \n\n")
	for (i in 1:numlines)
		{
		cat(txt[i], "\n", sep="")
		}
	return()
	}




#######################################################
# Scatter plot with axes drawn on the same scale
#######################################################
# http://onertipaday.blogspot.com/2007/05/scatter-plot-with-axes-drawn-on-same.html
# 
# I'd like to produce some scatter plots where N units on the X axis are > equal to N 
# units on the Y axis (as measured with a ruler, on screen or paper).
# x <- sample(10:200,40)
# y <- sample(20:100,40)
# windows(width = max(x),height = max(y))
# plot(x,y)

# try:
plot(x, y, asp = 1)

# or, better:
library(MASS)
eqscplot(x,y)

#or
library(lattice)
xyplot(y ~ x, aspect = "iso")





#######################################################
# STUPID BUGS
#######################################################

> source3 = '_R_tree_functions_v1.R'
> source(paste(sourcedir, source3, sep=""))
Warning messages:
1: In grepl("\n", lines, fixed = TRUE) :
  input string 27 is invalid in this locale
2: In grepl("\n", lines, fixed = TRUE) :
  input string 28 is invalid in this locale
3: In grepl("\n", lines, fixed = TRUE) :
  input string 32 is invalid in this locale
4: In grepl("\n", lines, fixed = TRUE) :
  input string 35 is invalid in this locale
5: In grepl("\n", lines, fixed = TRUE) :
  input string 37 is invalid in this locale
6: In grepl("\n", lines, fixed = TRUE) :
  input string 564 is invalid in this locale
7: In grepl("\n", lines, fixed = TRUE) :
  input string 620 is invalid in this locale
  
--> convert all code to ASCII to fix this...




















#######################################################
# Compiler switcheroo
#######################################################
junk='
switch_compiler_to_new()
switch_compiler_to_default()
junk='


# Change the compiler to e.g. g++46
switch_compiler_to_new <- function(old_compiler_symlink="/usr/bin/g++", new_compiler_fn="/my_gcc/bin/g++46", new_R_Makevars="/Users/nickm/.R/Makevars46")
	{
	
	# Replace the symlink to the compiler
	cmdstr = paste("rm ", old_compiler_symlink, "; ln -s ", new_compiler_fn, " ", old_compiler_symlink, sep="")
	system(cmdstr)
	
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	# Check the compiler
	cmdstr2 = "g++ -v"
	system(cmdstr2)
	
	return(cmdstr)
	}


# Change the compiler back to e.g. default g++42
switch_compiler_to_default <- function(old_compiler_symlink="/usr/bin/g++", new_compiler_fn="/usr/bin/llvm-g++-4.2", new_R_Makevars="/Users/nickm/.R/Makevars_default")
	{
	
	# Replace the symlink to the compiler
	cmdstr = paste("rm ", old_compiler_symlink, "; ln -s ", new_compiler_fn, " ", old_compiler_symlink, sep="")
	system(cmdstr)
	
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	# Check the compiler
	cmdstr2 = "g++ -v"
	system(cmdstr2)

	return(cmdstr)
	}



#######################################################
# SET/CHANGE ENVIRONMENT VARIABLES
#######################################################

# Running TNT in R requires something like this:
Sys.setenv("TERM"="dumb")
tntcmd = paste("./tnt.command proc ", x, ";", sep="")
system(tntcmd)




# OR THIS


final_script_fn = "master_scr.sh"

# Write the list of commands to a file
write.table(tntcmd, file=final_script_fn, quote=FALSE, append=FALSE, sep="", row.names = FALSE, col.names=FALSE)

cat("Run by pasting this:\n")
cat("\n")
#		chmod a+x *.sh
#		./master_scr.sh
cat("chmod a+x *.sh\n")
cat("./", final_script_fn, "", sep="")

Sys.setenv("TERM"="dumb")
system("chmod a+x *.sh\n")
system("./master_scr.sh")





# This requires modifying the environment variables for CPPFLAGS (-I flags in the compile statement to include headers)
# ...and perhaps -L flags (LIBS)

Sys.getenv("CPPFLAGS")
Sys.getenv("LIBS")

...this doesn't seem to get to install.packages/R CMD INSTALL, it is easier to modify the Makevars file at:

/Users/nickm/.R/Makevars_phyRmcmc

...e.g.:

new_Makevars()
old_Makevars()


# Change default Makevars for install
new_Makevars <- function(new_R_Makevars="/Users/nickm/.R/Makevars_phyRmcmc")
	{
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	return(cmdstr)
	}

# Change default Makevars for install
old_Makevars <- function(new_R_Makevars="/Users/nickm/.R/Makevars_default")
	{
	
	# Replace the $HOME/.R/Makevars file
	cmdstr = paste("rm ~/.R/Makevars; cp ", new_R_Makevars, " ~/.R/Makevars", sep="")
	system(cmdstr)
	
	return(cmdstr)
	}










#######################################################
# TNT commands
#######################################################


Hi! As the error message says, you are running out of memory.

Ideas:

1. Try different settings for the maximum number of trees, e.g.:
http://tnt.insectmuseum.org/index.php/Commands/hold

hold 178;
hold 400;
hold 10000;

(The stats.run script on tntwiki seems to crash for me about 400, but for other purposes I have used up to hold  is about the maximum I can consistently get to work in all my scripts, but sometimes one can go much higher; my maximum has been about:

hold 500000;

)



2. Try increasing the amount of RAM you allow TNT to use.  The default is only 16 MB:

http://tnt.insectmuseum.org/index.php/Commands/mxram

The highest I can go is:

mxram 4096;

...for a 4-gigabyte machine.  You should start with smaller values first, you can crash your computer if you go too big. e.g.:

mxram 32;
mxram 64;
mxram 128;
mxram 256;
mxram 512;
mxram 1024;
mxram 2048;
mxram 4096;



3. Finally, try running on both 32-bit and 64-bit TNT.  Some scripts work on one but not the other.

Cheers!
Nick





#######################################################
# RNCL -- Access to NEXUS Class Library
#######################################################
library(rncl)

wd = "/drives/Dropbox/_njm/"
setwd(wd)

trfn = "dinotree.newick"

tr = rncl(file=trfn, file.format="newick", spacesAsUnderscores=TRUE)
tr










#######################################################
# Nice tables with R package pixiedust
#######################################################

Tables so Beautifully Fine-Tuned You Will Believe It's Magic [R logo]

[Up] [Top]
Documentation for package ‘pixiedust’ version 0.7.0

DESCRIPTION file.
User guides, package vignettes and other documentation.
Package NEWS.
Help Pages

pixiedust-package	Tables So Beautifully Fine-Tuned You Will Believe It's Magic.
%<>%	Chain together multiple operations.
%>%	Chain together multiple operations.
as.data.frame.dust	Convert 'dust' Object to Data Frame
as.data.frame.dust_list	Convert 'dust' Object to Data Frame
assert_match_arg	Checkmate Compatible Version of 'match.arg'
dust	Dust Table Construction
dust.default	Dust Table Construction
dust.grouped_df	Dust Table Construction
dust.list	Dust Table Construction
get_pixie_count	Access and manipulate table numbers counters
glance_foot	Prepare Glance Statistics for 'pixiedust' Table Footer
increment_pixie_count	Access and manipulate table numbers counters
medley	Sprinkle Medleys
medley_all_borders	Apply Cell Borders to All Cells in a Region
medley_bw	Sprinkle Medleys
medley_model	Sprinkle Medleys
pixiedust	Tables So Beautifully Fine-Tuned You Will Believe It's Magic.
pixieply	Apply Functions Over 'dust_list' Objects
pixie_count	Access and manipulate table numbers counters
print.dust	Print A 'dust' Table
print.dust_list	Print A 'dust' Table
pvalString	Format P-values for Reports
redust	Dust Table Construction
redust.default	Dust Table Construction
redust.dust_list	Dust Table Construction
set_pixie_count	Access and manipulate table numbers counters
sprinkle	Define Customizations to a Table
sprinkle.default	Define Customizations to a Table
sprinkle.dust_list	Define Customizations to a Table
sprinkle_colnames	Column Names for 'dust' Tables
sprinkle_colnames.default	Column Names for 'dust' Tables
sprinkle_colnames.dust_list	Column Names for 'dust' Tables
sprinkle_print_method	Define Customizations to a Table
sprinkle_print_method.default	Define Customizations to a Table
sprinkle_print_method.dust_list	Define Customizations to a Table
sprinkle_table	Define Customizations to a Table
sprinkle_table.default	Define Customizations to a Table
sprinkle_table.dust_list	Define Customizations to a Table
tidy_levels_labels	Term and Level Descriptions for 'pixiedust' Tables



https://cran.r-project.org/web/packages/pixiedust/vignettes/pixiedust.html

x <- dust(lm(mpg ~ qsec + factor(am), data = mtcars))
x


fit <- lm(mpg ~ qsec + factor(am) + wt + factor(gear), 
          data = mtcars)
          
summary(fit)
broom::tidy(fit)

library(pixiedust)
dust(fit)


dust(fit) %>% 
  sprinkle(cols = c("estimate", "std.error", "statistic"), round = 2)
  
  
dust(fit) %>% 
  sprinkle(cols = c("estimate", "std.error", "statistic"),
           round = 3) %>% 
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) 



dust(fit) %>% 
  sprinkle(cols = c("estimate", "std.error", "statistic"),
           round = 3) %>% 
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  sprinkle_colnames(term = "Term", p.value = "P-value")
  
  

dust(fit) %>% 
  sprinkle(cols = c("estimate", "std.error", "statistic"),
           round = 3) %>% 
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  sprinkle_colnames(term = "Term", p.value = "P-value", 
                    std.error = "SE", statistic = "T-statistic",
                    estimate = "Coefficient")


dust(fit) %>% 
  sprinkle(cols = c("estimate", "std.error", "statistic"),
           round = 3) %>% 
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  sprinkle_colnames("Term", "Coefficient", "SE", "T-statistic", "P-value")



dust(fit) %>% 
  sprinkle(cols = "term", 
           replace = c("Intercept", "Quarter Mile Time", "Automatic vs. Manual",
                       "Weight", "Gears: 4 vs. 3", "Gears: 5 vs 3")) %>%
  sprinkle(cols = c("estimate", "std.error", "statistic"),
           round = 3) %>% 
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  sprinkle_colnames("Term", "Coefficient", "SE", "T-statistic", "P-value")
  
  

basetable <- dust(fit) %>% 
  sprinkle(cols = c("estimate", "std.error", "statistic"),
           round = 3) %>% 
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  sprinkle_colnames(term = "Term", estimate = "Coefficient", 
                    std.error = "SE", statistic = "T-statistic", 
                    p.value = "P-value") %>% 
  sprinkle_print_method("html")



basetable %>% 
  sprinkle(rows = c(2, 4), bold = TRUE, italic=TRUE)

