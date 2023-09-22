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
setwd(wd)

xlsfn = "groupTax_1282_mafftConstr_2023-08-07_edit_wGeneOrder.xlsx"
xls = openxlsx::read.xlsx(xlsfn)

gids = c(xls$accM3, xls$accM2, xls$accM1, xls$acc1, xls$acc2, xls$acc3)

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
length(TF3)

#######################################################
# Add *matching* MotBs from gene order, and from hmmsearch, to new columns in the Excel file
#######################################################
TF_M3 = xls$accM3 %in% gids2a
TFtmp = grepl(pattern="MotB", x=xls$nameM3, ignore.case=TRUE); TF_M3[TFtmp] = TRUE
TFtmp = grepl(pattern="ExbD", x=xls$nameM3, ignore.case=TRUE); TF_M3[TFtmp] = TRUE
TFtmp = grepl(pattern="TolR", x=xls$nameM3, ignore.case=TRUE); TF_M3[TFtmp] = TRUE
TFtmp = grepl(pattern="AglS", x=xls$nameM3, ignore.case=TRUE); TF_M3[TFtmp] = TRUE
length(TF_M3)
sum(TF_M3)
xls[TF_M3,]

TF_M2 = xls$accM2 %in% gids2a
TFtmp = grepl(pattern="MotB", x=xls$nameM2, ignore.case=TRUE); TF_M2[TFtmp] = TRUE
TFtmp = grepl(pattern="ExbD", x=xls$nameM2, ignore.case=TRUE); TF_M2[TFtmp] = TRUE
TFtmp = grepl(pattern="TolR", x=xls$nameM2, ignore.case=TRUE); TF_M2[TFtmp] = TRUE
TFtmp = grepl(pattern="AglS", x=xls$nameM2, ignore.case=TRUE); TF_M2[TFtmp] = TRUE
length(TF_M2)
sum(TF_M2)

TF_M1 = xls$accM1 %in% gids2a
TFtmp = grepl(pattern="MotB", x=xls$nameM1, ignore.case=TRUE); TF_M1[TFtmp] = TRUE
TFtmp = grepl(pattern="ExbD", x=xls$nameM1, ignore.case=TRUE); TF_M1[TFtmp] = TRUE
TFtmp = grepl(pattern="TolR", x=xls$nameM1, ignore.case=TRUE); TF_M1[TFtmp] = TRUE
TFtmp = grepl(pattern="AglS", x=xls$nameM1, ignore.case=TRUE); TF_M1[TFtmp] = TRUE
length(TF_M1)
sum(TF_M1)
# 1094


TF1 = xls$acc1 %in% gids2a
TFtmp = grepl(pattern="MotB", x=xls$name1, ignore.case=TRUE); TF1[TFtmp] = TRUE
TFtmp = grepl(pattern="ExbD", x=xls$name1, ignore.case=TRUE); TF1[TFtmp] = TRUE
TFtmp = grepl(pattern="TolR", x=xls$name1, ignore.case=TRUE); TF1[TFtmp] = TRUE
TFtmp = grepl(pattern="AglS", x=xls$name1, ignore.case=TRUE); TF1[TFtmp] = TRUE
length(TF1)
sum(TF1)
# 1094

TF2 = xls$acc2 %in% gids2a
TFtmp = grepl(pattern="MotB", x=xls$name2, ignore.case=TRUE); TF2[TFtmp] = TRUE
TFtmp = grepl(pattern="ExbD", x=xls$name2, ignore.case=TRUE); TF2[TFtmp] = TRUE
TFtmp = grepl(pattern="TolR", x=xls$name2, ignore.case=TRUE); TF2[TFtmp] = TRUE
TFtmp = grepl(pattern="AglS", x=xls$name2, ignore.case=TRUE); TF2[TFtmp] = TRUE
length(TF2)
sum(TF2)

TF3 = xls$acc3 %in% gids2a
TFtmp = grepl(pattern="MotB", x=xls$name3, ignore.case=TRUE); TF3[TFtmp] = TRUE
TFtmp = grepl(pattern="ExbD", x=xls$name3, ignore.case=TRUE); TF3[TFtmp] = TRUE
TFtmp = grepl(pattern="TolR", x=xls$name3, ignore.case=TRUE); TF3[TFtmp] = TRUE
TFtmp = grepl(pattern="AglS", x=xls$name3, ignore.case=TRUE); TF3[TFtmp] = TRUE
length(TF3)
sum(TF3)
xls[TF3,]

# Choose 1st, 2nd, or 3rd neighbor (choose closest)
TFs_table = cbind(TF_M3, TF_M2, TF_M1, TF1, TF2, TF3)
colnums = rep(NA, times=nrow(TFs_table))
MotB_bestMatch_gid = rep("", times=nrow(TFs_table))
MotB_bestMatch_name = rep("", times=nrow(TFs_table))
MotB_bestMatch_symbol = rep("", times=nrow(TFs_table))
nums = 1:6
score = c(1,2,3,6,5,4)
xls_acc = xls[,c("accM3","accM2","accM1","acc1","acc2","acc3")]
xls_name = xls[,c("nameM3","nameM2","nameM1","name1","name2","name3")]
xls_symbol = xls[,c("sym3","symM2","symM1","sym1","sym2","sym3")]
for (i in 1:nrow(TFs_table))
	{
	tmpTFs = TFs_table[i,]
	if (sum(tmpTFs) == 0)
		{
		colnums[i] = NA
		MotB_bestMatch_gid[i] = ""
		MotB_bestMatch_name[i] = ""
		MotB_bestMatch_symbol[i] = ""
		} else {
		max_score = max(score[tmpTFs])
		max_score_TF = score == max_score
		colnums[i] = nums[max_score_TF]
		MotB_bestMatch_gid[i] = xls_acc[i,colnums[i]]
		MotB_bestMatch_name[i] = xls_name[i,colnums[i]]
		MotB_bestMatch_symbol[i] = xls_symbol[i,colnums[i]]
		} # END if (is.infinite(colnums[i]) == TRUE)
	} # END for (i in 1:nrow(TFs_table))

MotB_bestmatches_df = as.data.frame(cbind(MotB_bestMatch_gid, MotB_bestMatch_name, MotB_bestMatch_symbol), stringsAsFactors=FALSE)
head(MotB_bestmatches_df)
dim(MotB_bestmatches_df)
sum(MotB_bestmatches_df$MotB_bestMatch_gid == "")

# Add back to the excel spreadsheet
xls2 = cbind(xls, MotB_bestmatches_df)

outfn2 = "groupTax_1282_mafftConstr_2023-08-07_edit_wGeneOrder_BestMotBs.xlsx"
openxlsx::write.xlsx(x=xls2, file=outfn2)
system(paste0("open ", outfn2))



