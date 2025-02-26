
# Search UniProt with queryup
# https://cran.r-project.org/web/packages/queryup/readme/README.html

# Example: P0ABV0 - TOLQ_ECO57
# https://www.uniprot.org/uniprotkb/P0ABV0/entry
# https://www.uniprot.org/help/query-fields
# query <- list("gene_exact" = "Pik3r1")

library(queryup)
query_fields$field

seqid = "NP_308799.1"
query <- list("xref" = seqid)
df <- query_uniprot(query, show_progress = FALSE)
head(df)

query <- list("xref" = seqid, database="pfam")
df <- query_uniprot(query, show_progress = FALSE)
head(df)



seqid = "WP_000131314.1"
query <- list("xref" = seqid)
df <- query_uniprot(query, show_progress = FALSE)
head(df)

seqid = "WP_000131314.1"
query <- list("xref" = seqid)
df <- query_uniprot(query, columns=c("xref_pdbsum", "xref_pdb"), show_progress = FALSE)
head(df)


seqid = "AAC73831"
query <- list("xref" = seqid)
df <- query_uniprot(query, show_progress = FALSE)
head(df)


seqid = "ABC76751"
query <- list("xref" = seqid)
df <- query_uniprot(query, show_progress = FALSE)
head(df)


seqid = "USF25003"
query <- list("xref" = seqid)
df <- query_uniprot(query, show_progress = FALSE)
head(df)


query <- list("xrefcount" = seqid)
df <- query_uniprot(query, show_progress = FALSE)
head(df)

# Others you might use:
# accession_id
# protein_name
# organism_name, organism_id
# reviewed


