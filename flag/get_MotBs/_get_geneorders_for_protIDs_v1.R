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
genome_ids = rep("", times=nrow(xls))
genome_dir_hits = rep("", times=nrow(xls))
for (i in 1:nrow(xls))
	{
	tmprow = xls[i,]
	
	# Try original version
	TF = grepl(pattern=tmprow$genome_id, x=genome_dirs)
	if (sum(TF) == 1)
		{
		genome_ids[i] = tmprow$genome_id
		genome_dir_hits[i] = genome_dirs[TF]
		hitscount[i] = sum(TF)
		next()
		}
	
	# Try GCF instead
	gcf_idnum = gsub(pattern="GCA_", replacement="GCF_", x=tmprow$genome_id)
	TF = grepl(pattern=gcf_idnum, x=genome_dirs)
	if (sum(TF) == 1)
		{
		genome_ids[i] = gcf_idnum
		genome_dir_hits[i] = genome_dirs[TF]
		hitscount[i] = sum(TF)
		next()
		}
	
	hitscount[i] = sum(TF)
	} # END for (i in 1:nrow(xls))
table(hitscount)

head(genome_dir_hits)


#######################################################
# Get the gene orders near each input MotA homolog
#######################################################
i = 1
sourceall("/GitHub/bioinfRhints/Rsrc/")

gene_neighbors_df_all = NULL

for (i in 1:nrow(xls))
	{
	cat("\n")
	cat(i, "/", nrow(xls))
	cat("\n")
	
	genome_dir_hit = slashslash(paste0(genomes_master_dir, "/", genome_dir_hits[i]))
	gene_order_table_fn = slashslash(paste0(genome_dir_hit, "/", genome_dir_hits[i], "_feature_table.txt"))
	gene_order_table_fn
	
	# Unzip if needed
	if ((file.exists(gene_order_table_fn) == FALSE) && (unzipTF == TRUE))
		{
		gene_order_table_zipfn = slashslash(paste0(genome_dir_hit, "/", genome_dir_hits[i], "_feature_table.txt.gz"))
		if (file.exists(gene_order_table_zipfn) == FALSE)
			{
			txt = paste0("Stop on i=", i, ". Gene order / feature table '", gene_order_table_zipfn, "' not found. Fix and repeat.")
			cat("\n")
			cat(txt)
			cat("\n")
			stop(txt)
			}
		cmdstr = paste0("gunzip ", gene_order_table_zipfn)
		system(cmdstr)
		}
	prot_feature_table_df = read_gene_order_table(gene_order_table_fn=gene_order_table_fn)
	print(dim(prot_feature_table_df))
	
	list_of_protIDs = xls$tipnames3_uniq[i]
	prot_feature_tables_all_df = prot_feature_table_df
	gene_neighbors_df = get_adjacent_genes(list_of_protIDs, prot_feature_tables_all_df, genomes_to_spnames_df=NULL, printwarnings=FALSE, prot_feature_tables_all_fn=gene_order_table_fn)
	
	# Build up table
	if (i == 1)
		{
		gene_neighbors_df_all = gene_neighbors_df
		} else {
		gene_neighbors_df_all = rbind(gene_neighbors_df_all, gene_neighbors_df)
		}
	} # END for (i in 1:nrow(xls))

head(gene_neighbors_df_all)
tail(gene_neighbors_df_all)
dim(gene_neighbors_df_all)



# Merge and save

outfn = gsub(pattern=".xlsx", replacement="_wGeneOrder.xlsx", x=xlsfn)
xls2 = cbind(xls, gene_neighbors_df_all)
openxlsx::write.xlsx(xls2, file=outfn)

file.copy(from=outfn, to="/GitHub/bioinfRhints/flag/get_MotBs/")




#######################################################
# Terminal commands to run HMMER searches against MotB / ExbD / TolR precalculated Interpro profiles:
#######################################################
# /Users/nmat471/HD/GitHub/bioinfRhints/flag/get_MotBs/_cmds_hmmer_motBs_v1.txt

wd = "/Users/nmat471/HD/GitHub/bioinfRhints/flag/get_MotBs/"
getwd(wd)

xlsfn = "groupTax_1282_mafftConstr_2023-08-07_edit_wGeneOrder.xlsx"
xls = openxlsx::read.xlsx(xlsfn)

gids = c(xls$acc0, xls$acc1, xls$acc2)

gids2 = readLines("unique_MotB_ExbD_TolR_hits_IDs.txt")

# Overlap between the two lists
TF1 = gids %in% gids2
TF2 = gids2 %in% gids

sum(TF1)
sum(TF2)

gids1a = gids[TF1]
gids2a = gids2[TF2]

TF3 = gids2a %in% gids1a

length(gids1a)
length(gids2a)

gids2a[TF3]
length(gids2a)
