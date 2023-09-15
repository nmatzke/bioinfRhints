library(ape)
library(msa)
library(seqinr)
library(BioGeoBEARS)
library(openxlsx)

sourceall("/GitHub/bioinfRhints/Rsrc/")

wd = "/GitHub/bioinfRhints/flag/get_MotBs/"
setwd(wd)



#######################################################
# Get potential MotBs by gene order
#######################################################

# Excel table of MotAs, with genome listed
xlsfn = "/GitHub/bioinfRhints/flag/AQB_classification/groupTax_1282_mafftConstr_2023-08-07_edit.xlsx"
genomes_master_dir = "~/Downloads/Full_genomes/genomes"
#setwd(genomes_dir)

xls = read.xlsx(xlsfn)
head(xls)
names(xls)


head(xls$genome_id)
tail(sort(xls$genome_id), n=20)

genome_dirs = list.files(genomes_master_dir)

unzipTF = TRUE  # If true, unzip the zipfiles if needed
i=1
hitscount = rep(0, times=nrow(xls))
genome_idnums = rep("", times=nrow(xls))
for (i in 1:nrow(xls))
	{
	tmprow = xls[i,]
	# Remove GCA_, GCF, or VER
	words = strsplit(tmprow$genome_id, split="_")[[1]]
	genome_idnum = strsplit(words[2], split="\\.")[[1]]
	#genome_idnum = gsub(pattern="VER_", replacement="", x=genome_idnum, ignore.case=FALSE)
	#genome_idnum = gsub(pattern="GCA_", replacement="", x=genome_idnum, ignore.case=FALSE)
	#genome_idnum = gsub(pattern="GCF_", replacement="", x=genome_idnum, ignore.case=FALSE)
	genome_idnums[i] = genome_idnum
	
	TF = grepl(pattern=genome_idnum, x=genome_dirs)
	hitscount[i] = sum(TF)
	} # END for (i in 1:nrow(xls))
table(hitscount)

xls[hitscount==0,]

cat("\nSearching (gene order tables) from ", length(genome_dirs), " genome directories...\n")
for (i in 1:length(genome_dirs))
	{
	if (i == 1)
		{
		cat("\ngenome_dirs #", i, ",", sep="")
		} else {
		cat(i, ",", sep="")
		}
	genome_dir = genome_dirs[i]
	# Access that directory, look for protein in list
	gene_order_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_feature_table.txt")
	if ((file.exists(gene_order_table_fn) == FALSE) && (unzipTF == TRUE))
		{
		gene_order_table_zipfn = slashslash(paste0(wd, "/", "genomes/", genome_dir, "/", genome_dir, "_feature_table.txt.gz"))
		cmdstr = paste0("gunzip ", gene_order_table_zipfn)
		system(cmdstr)
		}

	prot_feature_table_df = read_gene_order_table(gene_order_table_fn=gene_order_table_fn)
	dim(prot_feature_table_df)

	# Also get the species name from the 
	# GCA_000005845.2_ASM584v2_assembly_report.txt
	assembly_report_table_fn = paste0("genomes/", genome_dir, "/", genome_dir, "_assembly_report.txt")
	spname = get_spname_from_assembly_report(assembly_report_table_fn)

	# Add some stuff to the front of the line
	id = genome_dir
	prot_feature_table_df = cbind(id, prot_feature_table_df)
	prot_feature_table_df = cbind(spname, prot_feature_table_df)
	
	#prot_feature_tables_all_df = rbind(prot_feature_tables_all_df, prot_feature_table_df)

	if (i == 1)
		{
		write.table(prot_feature_table_df, file=prot_feature_tables_all_fn, sep="\t", row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE)
		} else {
		write.table(prot_feature_table_df, file=prot_feature_tables_all_fn, sep="\t", row.names=FALSE, append=TRUE, quote=FALSE, col.names=FALSE)
		}
	} # END for (i in 1:length(genome_dirs))

