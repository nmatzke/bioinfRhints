#######################################################
# Scripts: manipulating alignments, sequence names / tip names, etc.
#######################################################
library(ape)
library(phytools)
library(seqinr)				# for read.fasta
library(BioGeoBEARS)	# for cls.df
library(gdata)				# for trim
library(stringr) 			# for regmatches

sourceall("/GitHub/bioinfRhints/R/BEASTmasteR/") 
# for remove_equals_from_tips()
# for read.beast.table_original


#wd = "/GitHub/bioinfRhints/flag/530_MotA_rename/"
wd = "/GitHub/bioinfRhints/flag/589_MotA_rename/"
setwd(wd)

#alnfn = "motA_alignment_refined.fasta"
alnfn = "motA_hmmalign589_full.fasta"
aln = read.fasta(alnfn)

#trfn = "530_sequences_Alignment_contree_reRootLadder_gIDs.newick"
#nexfn = "530_sequences_Alignment_contree_reRootLadder_gIDs.nexus"
trfn = "motA589_IQtree_midpoint_FigTree.newick"
nexfn = "motA589_IQtree_midpoint_FigTree.nexus"



tr = read.tree(trfn)
tr2 = read.nexus(nexfn)



# "aln" is an R list, each element is a sequence, with a name and some other attributes:
attributes(aln[[1]])


# Get the full names of all the sequences
# (lapply = "list apply" = apply the function "attr" to each element in the list "aln")
# (unlist turns the list of names into a vector)
fullnames = unlist(lapply(X=aln, FUN=attr, which="Annot"))
gid = unlist(lapply(X=aln, FUN=attr, which="name")) # gid = Genbank IDs

# Extract the species (between brackets)
# Examples:
tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium]"
regmatches(tmptxt, gregexpr("\\[.+?\\]", tmptxt))

tmptxt = ">QQS07318.1 MAG: MotA/TolQ/ExbB proton channel family protein [Fibrobacteres bacterium] [hello]"
regmatches(tmptxt, gregexpr("\\[.+?\\]", tmptxt))

# Function to extract last item in [square brackets like this] in a string
# ...repeat for each item in a list
list_of_strings = fullnames; replace_spaces=TRUE
extract_last_brackets <- function(list_of_strings, replace_spaces=TRUE)
	{
	species_names = rep("", length(list_of_strings))
	txt = paste0("\nextract_last_brackets() is processing ", length(list_of_strings), " strings. String #")
	cat(txt)
	for (i in 1:length(list_of_strings))
		{
		cat(i, ",", sep="")
		tmptxt = list_of_strings[i]
		tmpstrs = unlist(regmatches(tmptxt, gregexpr("\\[.+?\\]", tmptxt)))
		tmpstr = tmpstrs[length(tmpstrs)] # take the last bracketed text, if more than 1
		
		# Remove "[", "]"
		tmpstr = gsub(pattern="\\[", replacement="", x=tmpstr)
		tmpstr = gsub(pattern="\\]", replacement="", x=tmpstr)
		
		# Replace spaces with "_"
		if (replace_spaces == TRUE)
			{
			tmpstr = gsub(pattern=" ", replacement="_", x=tmpstr)
			}
		species_names[i] = tmpstr
		}
	return(species_names)
	}

# Run the function
species_names = extract_last_brackets(list_of_strings=fullnames)
species_names


#######################################################
# Get the sequence lengths (ignoring indels "-")
#######################################################
get_seqlengths <- function(aln)
	{
	seqlengths = rep(0, times=length(aln))
	for (i in 1:length(aln))
		{
		# Length of this row of the alignment, minus indels ("-")
		seqlengths[i] = length(aln[[i]]) - sum(aln[[i]] == "-")
		}
	return(seqlengths)
	}


seqlengths = get_seqlengths(aln)
hist(seqlengths, breaks=50)




#######################################################
# Parse the names for protein name info
#######################################################
protein_name = rep("", times=length(aln))
species_names_wSpaces = extract_last_brackets(list_of_strings=fullnames, replace_spaces=FALSE)

for (i in 1:length(aln))
	{
	orig_name = fullnames[i]

	# Remove ">", ], [
	subname = gsub(pattern=">", replacement="", x=orig_name)
	subname = gsub(pattern="\\[", replacement="", x=subname)
	subname = gsub(pattern="\\]", replacement="", x=subname)

	# Remove GenBank ID:
	subname = gsub(pattern=gid[i], replacement="", x=subname)
	
	# Remove species name:
	subname = gsub(pattern=species_names_wSpaces[i], replacement="", x=subname)

	# REMOVE SEMICOLONS, GD IT
	subname = gsub(pattern="\\;", replacement="", x=subname, ignore.case=TRUE)
	# REMOVE EQUALS, GD IT
	subname = gsub(pattern="\\=", replacement="EQ", x=subname, ignore.case=TRUE)

	
	# fix spaces
	subname = gsub(pattern="  ", replacement=" ", x=subname)
	subname = gsub(pattern="  ", replacement=" ", x=subname)
	subname = gsub(pattern="  ", replacement=" ", x=subname)
	subname = gsub(pattern="  ", replacement=" ", x=subname)
	subname = gdata::trim(subname)
		
	protein_name[i] = subname
	}

sort(table(protein_name))


list_of_strings = protein_name
sum(grepl(pattern="flagellar stator protein MotA", x=list_of_strings, ignore.case=TRUE))

classify_MotAfam_labels <- function(list_of_strings)
	{
	short_protname = rep("", length(list_of_strings))
	for (i in 1:length(list_of_strings))
		{
		tmpstr = list_of_strings[i]
		
		# Add semicolon on the end
		
		# Many are just MotA/TolQ/ExbB
		if (grepl(pattern="MotA/TolQ/ExbB", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "AQB"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="gliding", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "gliding"
			short_protname[i] = tmpstr
			next()
			}
		
		# Remove obvious words - flagellum
		tmpstr = gsub(pattern="transport membrane proton channel, TolQ-related protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport membrane proton channel, TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)


		tmpstr = gsub(pattern="flagellar motor rotation protein ", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar motor protein ", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MAG: flagellar motor protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="chemotaxis protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="chemotaxis", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="proton channel family", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar motor component", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MotA protein", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="motility protein A", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar basal body stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="endoflagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="flagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="MAG\\: flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="MotA\\; MotA component of the H\\+\\-coupled stator flagellum complex", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="Mot family proton \\(H\\+\\) or sodium \\(Na\\+\\) transporter MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="sodium channel stator-force generator subunit of flagellar rotation", replacement="sodium_MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="signal recognition particle-docking protein FtsY", replacement="SRP_FtsY", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="sodium-driven polar flag MotA", replacement="sodiumPolarMotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="motA\\, \\(MotA\\)", replacement="MotA", x=tmpstr, ignore.case=TRUE)
 		tmpstr = gsub(pattern="MAG: flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="pomA protein", replacement="PomA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="protein of unassigned function", replacement="unk", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="multi-sensor hybrid histidine kinase", replacement="HistKin", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MotA\\; MotA component of the H\\+\\-coupled stator flagellum complex", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Mot family proton \\(H+\\) or sodium \\(Na+\\) transporter MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="membrane spanning protein in TonB-ExbB-ExbD complex", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		#tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)
		#tmpstr = gsub(pattern="", replacement="", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="flagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="putative", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="motor", replacement="", x=tmpstr, ignore.case=TRUE)



		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)
		
		# Any blank here is just a motor protein
		if ((tmpstr == "") || (tmpstr == "protein"))
			{
			tmpstr = "motor"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="tonB system transport protein ExbB/TolQ", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "TolQ_ExbB_wTonB"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="tonB system transport protein ExbB/TolQ", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "TolQ_ExbB_wTonB"
			short_protname[i] = tmpstr
			next()
			}

		if (grepl(pattern="ExbB/TolQ", x=tmpstr, ignore.case=TRUE) == TRUE)
			{
			tmpstr = "ExbB_TolQ"
			short_protname[i] = tmpstr
			next()
			}
		tmpstr = gsub(pattern="system transport protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="accessory protein", replacement="", x=tmpstr, ignore.case=TRUE)


		if (tolower(tmpstr) == "probable biopolymer transport protein")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}

		if (tolower(tmpstr) == "biopolymer transporter")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}
		if (tolower(tmpstr) == "biopolymer transport proteins")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}
		if (tolower(tmpstr) == "biopolymer transport protein")
			{
			tmpstr = "biopoly_transp"
			short_protname[i] = tmpstr
			next()
			}


		tmpstr = gsub(pattern="hypothetical protein ", replacement="HYP_", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transport protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="outer membrane transport energization protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Biopolymer", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="domain-containing protein", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="colicin uptake protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbB-related protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ferric siderophore transport system,", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable tolQ-type", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable tolQ protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Cell division and transport-associated protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolQ protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Tol-Pal system protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Protein TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="(ExbB-like)", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Ton complex subunit ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQ-like protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQ transporter", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Tol-Pal system subunit TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="\\(ExbB\\)", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbB protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tonB\\-system energizer ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TonB exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)


		tmpstr = gsub(pattern="Cell division and transportassociated TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tonBsystem energizer ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolPal system TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQtype", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolPal system subunit TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQlike protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolQrelated protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="sodium channel statorforce generator subunit of flagellar rotation", replacement="sodium MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="related to TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="related to \\(TolQ\\)", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolerance to group A colicins, singlestranded DNA filamentous phage, required for OM integrity", replacement="TolA_OM", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="\\(ExbBlike\\)", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBrelated protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBlike", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBrelated protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="ExbBlike protein", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Cell division and transportassociated TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tonBsystem energizer ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TolPal system TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolpal systemassociated acylCoA thioesterase", replacement="tolpal thioest", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="sodiumdriven polar flag MotA", replacement="sodium MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MotA MotA component of the H+coupled stator flagellum complex", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical proteintransmembrane prediction", replacement="HYP", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical protein", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="probable TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="membrane spanning protein in TonBExbBExbD complex", replacement="ExbB", x=tmpstr, ignore.case=TRUE)




		tmpstr = gsub(pattern="MAG: TolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical membrane protein", replacement="HYP", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="hypothetical protein-transmembrane prediction", replacement="", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern=", tolQ", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tolerance to group A colicins, single-stranded DNA filamentous phage, required for OM integrity", replacement="colicin-tol", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="tol-pal system-associated acyl-CoA thioesterase ", replacement="TolPal_thioester", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="two-component sensor histidine kinase YdfI", replacement="2compSens_hisKin_YdfI", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="protein of unassigned function", replacement="unk", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="conserved hypothetical protein", replacement="HYP", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="motA\\, \\(MotA\\)", replacement="MotA", x=tmpstr, ignore.case=TRUE)


		tmpstr = gsub(pattern="flagellar stator protein MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MAG: MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar protein", replacement="flag", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="flagellar proton channel", replacement="flag", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="MAG\\: MotA", replacement="MotA", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)

		tmpstr = gsub(pattern="tol\\-pal system\\-associated acyl\\-CoA thioesterase", replacement="tolPal_thioest", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="probable ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="TonB ExbB2", replacement="ExbB2", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TonB exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="TonB ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="LafT protein", replacement="LafT", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="Tolq", replacement="TolQ", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="exbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="transporter ExbB", replacement="ExbB", x=tmpstr, ignore.case=TRUE)
		tmpstr = gsub(pattern="LafT protein", replacement="LafT", x=tmpstr, ignore.case=TRUE)

		tmpstr = gsub(pattern="uncharacterized protein", replacement="unk", x=tmpstr)


		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gsub(pattern="  ", replacement=" ", x=tmpstr)
		tmpstr = gdata::trim(tmpstr)

		if ((tmpstr == "") || (tmpstr == "transporter"))
			{
			tmpstr = "transport"
			short_protname[i] = tmpstr
			next()
			}

		
		short_protname[i] = tmpstr
		}
	return(short_protname)
	}

short_protname = classify_MotAfam_labels(list_of_strings=protein_name)
short_protname

sort(table(short_protname))

# Remove "=" from tipnames before renaming
species_names_wSpaces = gsub(pattern="=", replacement="EQ", x=species_names_wSpaces) # remove "="
species_names = gsub(pattern="=", replacement="EQ", x=species_names) # remove "="


newnames1_df = cbind(species_names, short_protname, gid, seqlengths)
newnames2_df = cbind(short_protname, species_names, gid, seqlengths)
newnames3_df = cbind(seqlengths, short_protname, species_names, gid)
newnames1 = apply(X=newnames1_df, MARGIN=1, paste0, collapse="|")
newnames2 = apply(X=newnames2_df, MARGIN=1, paste0, collapse="|")
newnames3 = apply(X=newnames3_df, MARGIN=1, paste0, collapse="|")


nexfn2 = gsub(pattern=".nexus", replacement="_noEQ.nexus", x=nexfn)
remove_equals_from_tips(nexfn=nexfn, outfn=nexfn2, format="raxml")


# Write the new names into the NEXUS file; this keeps the bootstraps etc.
outfn = gsub(pattern="_noEQ.nexus", replacement="_newNames1.nexus", x=nexfn2)
tmpstrs = readLines(nexfn)

for (i in 1:length(tmpstrs))
	{
	tmpstr = tmpstrs[i]
#	tmpstr = gsub(pattern="'", replacement="", x=tmpstr)
	for (j in 1:length(gid))
		{
		tmpstr = gsub(pattern=gid[j], replacement=newnames1[j], x=tmpstr)
		}
	tmpstrs[i] = tmpstr
	}
writeLines(text=tmpstrs, con=outfn)
outfn1 = gsub(pattern=".nexus", replacement="_nameFirst.nexus", x=nexfn)
writeLines(text=tmpstrs, con=outfn1)


# Write the new names into the NEXUS file; this keeps the bootstraps etc.
outfn = gsub(pattern=".nexus", replacement="_newNames2.nexus", x=nexfn)
tmpstrs = readLines(nexfn)

for (i in 1:length(tmpstrs))
	{
	tmpstr = tmpstrs[i]
	tmpstr = gsub(pattern="'", replacement="", x=tmpstr)
	
	for (j in 1:length(gid))
		{
		tmpstr = gsub(pattern=gid[j], replacement=newnames2[j], x=tmpstr)
		}
	tmpstrs[i] = tmpstr
	}
writeLines(text=tmpstrs, con=outfn)
outfn2 = gsub(pattern=".nexus", replacement="_protFirst.nexus", x=nexfn)
writeLines(text=tmpstrs, con=outfn2)



# Write the new names into the NEXUS file; this keeps the bootstraps etc.
outfn = gsub(pattern=".nexus", replacement="_newNames3.nexus", x=nexfn)
tmpstrs = readLines(nexfn)

for (i in 1:length(tmpstrs))
	{
	tmpstr = tmpstrs[i]
	tmpstr = gsub(pattern="'", replacement="", x=tmpstr)
	for (j in 1:length(gid))
		{
		tmpstr = gsub(pattern=gid[j], replacement=newnames3[j], x=tmpstr)
		}
	tmpstrs[i] = tmpstr
	}
writeLines(text=tmpstrs, con=outfn)
outfn3 = gsub(pattern=".nexus", replacement="_lengthFirst.nexus", x=nexfn)
writeLines(text=tmpstrs, con=outfn3)


sort(table(short_protname))







#######################################################
# Coloring in the branches by size
#######################################################

# Classify the sizes
size_classes = rep(0, times=length(seqlengths))
categories = c(0, 100, 200, 300, 400, 500, 600, 700, 800, (max(seqlengths)+1))
charval = c(0, 1, 2, 3, 4, 5, 6, 7, 8)

for (i in 1:(length(categories)-1))
	{
	minval = categories[i]
	maxval = categories[i+1]
	
	TF1 = seqlengths >= minval
	TF2 = seqlengths < maxval
	TF = (TF1 + TF2) == 2
	size_classes[TF] = charval[i]
	}

hist(size_classes)

names(size_classes) = newnames3

library(ape)
library(phytools)

# Ladderize will flip up/down, so it looks like FigTree

tr2 = phytools::readNexus(outfn3, format="standard")

if (is.binary(tr2) == FALSE)
	{
	tr2 = multi2di(tr2)
	
	brlens = tr2$edge.length
	TF1 = is.na(brlens)
	TF2 = is.nan(brlens)
	brlens[TF1] = 0.0
	brlens[TF2] = 0.0
	TF3 = brlens <= 0.0
	brlens[TF3] = min(brlens > 0.0)
	tr2$edge.length = brlens
	}

TF = sort(tr2$tip.label) == sort(newnames3)

namevals = cbind(sort(tr2$tip.label), sort(newnames3))
namevals[TF==FALSE,]
sum(TF == FALSE)



#stats_table = read.beast.table_original(file=outfn3, digits=4) 


tr2 = ladderize(tr2, right=FALSE)
tr2$tip.label = gsub(pattern="'", replacement="", x=tr2$tip.label)
#plot(tr2, show.tip.label=FALSE)


tip_rownum_in_sorted = rep(0, times=length(size_classes))
char_rownum_in_charnames = rep(0, times=length(size_classes))
for (i in 1:nrow(namevals))
	{
	tip_rownum_in_sorted[i] = match(namevals[i,1], table=tr2$tip.label)
	char_rownum_in_charnames[i] = match(namevals[i,2], table=newnames3)
	}
orders = cbind(tip_rownum_in_sorted, char_rownum_in_charnames)
tiporder = order(orders[,1])
char_into_tip_order = orders[tiporder,2]

size_classes = size_classes[char_into_tip_order]
names(size_classes)[1:5]
tr2$tip.label[1:5]


# Create an ordered character
params_matrix = matrix(data=0, nrow=length(unique(size_classes)), ncol=length(unique(size_classes)))
for (i in 1:nrow(params_matrix))
	{
	for (j in 1:ncol(params_matrix))
		{
		if (j == i+1)
			{
			params_matrix[i,j] = 1
			}
		if (i == j+1)
			{
			params_matrix[i,j] = 1
			}
		}
	}
params_matrix	

res = ace(x=size_classes, phy=tr2, type="discrete", model=params_matrix)


#pdffn = "530_sequences_Alignment_contree_reRootLadder_gIDs_protFirst_sizeColors.pdf"
pdffn = gsub(pattern=".nexus", replacement="_sizeColors.pdf", x=outfn2)
pdf(file=pdffn, width=12, height=48)

titletxt = "Colored plot of sequence lengths\n(red=low, yellow=high"
plot.phylo(tr2, type="phylogram", label.offset=0.0, cex=0.45, font=1, align.tip.label=TRUE)#, tip.color="white")
title(titletxt)

colors <- heat.colors(n=length(unique(size_classes)))
tiplabels(pch=22, bg=colors[as.numeric(size_classes)], cex=2, adj=0.5)
res$lik.anc[res$lik.anc < 0] = 0.0
nodelabels(pie=res$lik.anc, piecol=colors, cex=0.4)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


