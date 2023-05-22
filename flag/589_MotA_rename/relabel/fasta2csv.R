library(seqinr)
library(taxizedb)
library(ape)
library(rentrez)
library(Biostrings)




fasta = read.fasta("motA589_subset_to_EcoliK12.fasta")
n = names(fasta)



accessions = sapply(strsplit(n, "[.]1"), function(ele) ele[1])



#accessions.df = data.frame(domain = character(0),  phylum = character(0), class = character(0), order = character(0),  family = character(0),  genus=character(0),
						#species = character(0), genbank = character(0), protein = character(0), description = character(0), oldname = character(0))

accessions.df = read.table("accessions.tsv", sep="\t", header=T)
accessions.df = accessions.df[accessions.df$species != "",]
for (i in 1:length(n)){

	acc = accessions[i]
	accessions.df2 = data.frame(domain = "", phylum = "", class = "", order="", family="", genus="", species = "", genbank = "" , protein = acc, prokka="", description = "", oldname = n[i])


	if (any(accessions.df$protein == acc)){
		cat(paste("Skipping", acc, "because it's already in outfile\n"))
		next
	}


	cat(paste0(n[i], " -> ", acc, "\n"))




	out <- tryCatch(
	        {

		result = entrez_fetch(db="protein", id=acc, rettype="gp")

		# Get genbank
		result = strsplit(result, "\n")[[1]]
		gb = grep("^DBSOURCE    accession", result)

		


		if (length(gb) == 0){
			cat(paste("Warning: cannot find genbank for", acc, "\n"))

			gb = ""
			


		}else{
			gb = strsplit(result[gb], " +")[[1]][3]


			# Get taxonomy
			#seq_DNAbin = read.GenBank(gb, as.character = TRUE)


			# Get the species and its phylum
			#species = attr(seq_DNAbin, "species")

			

		}

		species = gsub("^ +", "", gsub(".+GANISM ", "", result[grep("ORGANISM ", result)][1]))
			

		# Description
		desc = result[grep("^DEFINITION", result)]
		desc = gsub("^DEFINITION +", "", desc)
		desc = gsub(" [[].+", "", desc)
		desc = gsub(":", "", desc)
		desc = gsub(";", "", desc)
		desc = gsub(",", "", desc)


		species = gsub("_", " ", species)
		cat(paste0("\t\tSpecies: " , species, "\n"))
		species.result = name2taxid(species, out_type="summary")
		class.df = classification(species.result$id)[[1]]
		genus = class.df[class.df$rank == "genus", "name"] 
		phylum = class.df[class.df$rank == "phylum", "name"] 
		family = class.df[class.df$rank == "family", "name"] 
		class = class.df[class.df$rank == "class", "name"] 
		order = class.df[class.df$rank == "order", "name"] 
		domain = class.df[class.df$rank == "superkingdom", "name"] 
		cat(paste0("\t\tGenus: ", genus, "\n"))
		cat(paste0("\t\tFamily: ", family, "\n"))
		cat(paste0("\t\tOrder: ", order, "\n"))
		cat(paste0("\t\tClass: ", class, "\n"))
		cat(paste0("\t\tPhylum: ", phylum, "\n"))
		cat(paste0("\t\tDomain: ", domain, "\n"))

		if (length(genus) == 0) {
			genus = ""
		}
		if (length(phylum) == 0) {
			phylum = ""
		}
		if (length(family) == 0) {
			family = ""
		}
		if (length(class) == 0) {
			class = ""
		}
		if (length(order) == 0) {
			order = ""
		}
		if (length(domain) == 0) {
			domain = ""
		}




		# Append to dataframe
		accessions.df2 = data.frame(domain = domain,  phylum = phylum, class = class,  order = order, family = family, genus = genus,  species = species,  genbank = gb , protein = acc, prokka="", description = desc, oldname = n[i])

		#write.table(accessions.df, "accessions.tsv", sep="\t", row.names=F, quote=F)

		rm(result)
		#rm(seq_DNAbin)


		


	}, error = function(cond){
		cat(paste("Warning: encountered error for", acc, "\n"))
		

		
	})


	accessions.df = rbind(accessions.df, accessions.df2)
	write.table(accessions.df, "accessions.tsv", sep="\t", row.names=F, quote=F)


}
