

# 4.2NCBI identical protein retrieval
# https://bioconductor.riken.jp/packages/3.21/bioc/vignettes/ginmappeR/inst/doc/ginmappeR.html
library(ape)
library(BioGeoBEARS)
library(ginmappeR)

sourceall("/GitHub/str2phy/Rsrc/")
source("/GitHub/str2phy/Rsrc/blast/blastR_setup/blastsequences_v4.R")

wd = "/GitHub/bioinfRhints/flag/2025-02-19_CPL_mBio/10_seqIDs_to_classifications/"
setwd(wd)

trfn = "A379tree.nexus"
tr = read.nexus(trfn)
seqnames = tr$tip.label
head(seqnames)
seqids = get_leading_seqids_from_name(strings=seqnames, split="_")
catn(seqids[1:10])


runslow = TRUE
base_fn = "379_AQBs" # includes directory
version = "v1"

out_table_fn = paste0(base_fn, "_wUniProtdata_", version, ".Rdata")


# Get identical protein seqids for each seqid, using getNCBIIdenticalProteins
motA_identical_protIDs_fn = paste0(base_fn, "_identical_proteins_", version, ".Rdata")
identical_protein_seqids_list = get_IDs_identical_proteins(seqids, runslow=runslow, motA_identical_protIDs_fn =motA_identical_protIDs_fn)
# Counts
counts_identical_proteins = sapply(FUN=length, X=identical_protein_seqids_list)
rev(sort(counts_identical_proteins))
sort(counts_identical_proteins)

# Take a list of identical proteins, return just the ones that match genbank_prefixes()

# Extract IDs with genbank prefixes
recodes = genbank_prefixes()
identical_protein_seqids_list_wGenBank_codes = sapply(X=identical_protein_seqids_list, FUN=reduce_identical_proteins_to_matching_codes, codes=recodes)
seqids_wGenBank_codes_counts = sapply(X=identical_protein_seqids_list_wGenBank_codes, FUN=length)

seqids_wGenBank_codes_txt = sapply(X=identical_protein_seqids_list_wGenBank_codes, FUN=paste0, split="", collapse=",")

# Use queryup to search for UniProt for each entry
library(queryup)
query_fields$field
uniprot_xref_list = get_uniprot_xrefs(seqids)

# fields you can return:
printall(return_fields)

query <- list("xref" = seqids)

runslow = runslow
#seqids_to_uniprot_accession_entries_fn = "seqids_to_uniprot_accession_entries_df_v1.Rdata"
seqids_to_uniprot_accession_entries_fn = paste0(base_fn, "_seqids_to_uniprot_accession_entries_df_", version, ".Rdata")

if (runslow)
	{
	seqids_to_uniprot_accession_entries_df = get_uniprot_accession_by_seqid(seqids)
	save(seqids_to_uniprot_accession_entries_df, file=seqids_to_uniprot_accession_entries_fn)
	} else {
	# Loads to: seqids_to_uniprot_accession_entries_df
	load(file=seqids_to_uniprot_accession_entries_fn)
	} # END if (runslow)

# Get lots of fields, but not linked to seqid
bigdf = queryup::query_uniprot(query, columns=c("accession", "id", "gene_names", "protein_name", "organism_name", "organism_id", "reviewed", "xref_refseq", "xref_cdd", "xref_pfam", "organelle", "length", "mass"), show_progress=TRUE)
# Do these separately, they crash when linked to other fields
bigdf2_pdbsum = queryup::query_uniprot(query, columns=c("xref_pdbsum"), show_progress=TRUE)
bigdf2_pdb = queryup::query_uniprot(query, columns=c("xref_pdb"), show_progress=TRUE)
bigdf2_alphafold = queryup::query_uniprot(query, columns=c("xref_alphafolddb"), show_progress=TRUE)

dim(seqids_to_uniprot_accession_entries_df)
bigdf_outdf = link_query_uniprot_table_to_accession_table(seqids_to_uniprot_accession_entries_df, bigdf)
dim(bigdf_outdf)

bigdf_outdf = link_query_uniprot_table_to_accession_table(bigdf_outdf, bigdf2_pdbsum)
dim(bigdf_outdf)
bigdf_outdf = link_query_uniprot_table_to_accession_table(bigdf_outdf, bigdf2_pdb)
dim(bigdf_outdf)
bigdf_outdf = link_query_uniprot_table_to_accession_table(bigdf_outdf, bigdf2_alphafold)
dim(bigdf_outdf)


bigdf_outdf = cbind(bigdf_outdf, counts_identical_proteins, seqids_wGenBank_codes_txt, seqids_wGenBank_codes_counts)
bigdf_outdf







bigdf_outdf = get_uniprot_data_on_seqids(seqids, runslow=TRUE, base_fn="379_AQBs", version="v1")
head(bigdf_outdf)
dim(bigdf_outdf)

