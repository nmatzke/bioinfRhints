library(seqinr)  # Required for reading FASTA files

wd = "/Users/cpue388/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Flagellum/Alignments/motA/hmmer_alignments/relabel"
setwd(wd)
wd


args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("motA589_subset_to_EcoliK12.fasta, accessions.tsv")
}





label_file = args[1]

# Read the CSV file
csv_file = args[2]
metadata.df = read.table(csv_file, sep="\t", header=T)




# Rename descriptions
tmpstr = metadata.df$description
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
metadata.df$description = tmpstr




# Return new name of an accession
rename = function(acc, accessions.df){




	# Find the row in the metadata
	match = which(metadata.df$oldname == acc)
	if (length(match) == 0){

		if (length(match) == 0){
			stop(paste("Error: could not find'", acc, "'in the metadata file under column 'oldname'"))
		}


		
	}

	if (length(match) > 1){
		stop(paste("Error: found", length(match), "matches'", acc, "'in the metadata file under column 'oldname'"))
	}

	domain = substr(metadata.df[match,"domain"], 1, 4)
	phylum = metadata.df[match,"phylum"]
	species = metadata.df[match,"species"] # only take the first two tokens of species, to avoid '=' or ',' characters
	species = gsub(" ", "_", species)
	species = sapply(strsplit(species, "_"), function(ele) paste(ele[1], ele[2], sep="_"))
	desc = metadata.df[match,"description"]
	desc = gsub(" ", "_", desc)
	desc = gsub("/", "_", desc)
	accession = metadata.df[match,"protein"]
	if (accession == "") accession = metadata.df[match,"prokka"]

	newname = paste(domain, species, desc, accession, sep="_")

	newname 

}





if (length(grep(".fasta", label_file)) == 1){

	# Read the FASTA file
	sequences = read.fasta(label_file)
	outfile = gsub(".fasta", ".tidy.fasta", label_file)

	for (acc in names(sequences)){


		newname = rename(acc, metadata.df)

		names(sequences)[names(sequences) == acc] = newname

		cat(paste("Renaming", acc, "to", newname, "\n"))

	}


	cat(paste("Saving fasta to", outfile, "\n"))


	write.fasta(sequences, names(sequences), outfile)



}else{
	cat(paste("Please ensure the input file is .fasta format\n"))
}

