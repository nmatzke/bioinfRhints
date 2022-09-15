#######################################################
# Converting a list of FASTA files to NEXUS
#######################################################

library(ape)
library(BioGeoBEARS)
library(gdata)	# for trim
library(seqinr)
library(stringr)
library(pegas)

# Get the files from here:
alns_dir = "/drives/GDrive/__GDrive_projects/2016-07-31_divide_and_conquer_starBEAST/nana5_BEAST_50alignments/"

fns = list.files(path=alns_dir)

# These are FASTA files
fns = fns[endsWith(x=fns, suffix="fas")]
fns


# Make a for-loop, convert each file to NEXUS
cat("\nProcessing ", length(fns), " files...\n")
for (i in 1:length(fns))
	{
	cat(i, " ", sep="")

	# Read in FASTA file
	fasta_fn = fns[i]
	out_nexus_aln_fn = paste0(fasta_fn, ".nexus")
	fasta_fn = slashslash(paste0(addslash(alns_dir), fns[i]))
	out_nexus_aln_fn = slashslash(paste0(addslash(alns_dir), out_nexus_aln_fn))
		
	# Write out to NEXUS
	seqs = seqinr::read.fasta(fasta_fn)
	ape::write.nexus.data(x=seqs, file=out_nexus_aln_fn)
	
	}



# These are NEXUS files
fns = list.files(path=alns_dir)
fns = fns[endsWith(x=fns, suffix="nexus")]
fns

cat(fns, sep="\n")

convert_NEXUS_to_FASTA <- function(fasta_fn, out_nexus_aln_fn=NULL)
	{
	require(ape)
	require(seqinr)
	
	if (is.null(out_nexus_aln_fn))
		{
		out_nexus_aln_fn = paste0(fasta_fn, ".nex")
		}
		
	# Write out to NEXUS
	seqs = seqinr::read.fasta(fasta_fn)
	ape::write.nexus.data(x=seqs, file=out_nexus_aln_fn)
	txt = paste0("\nconvert_NEXUS_to_FASTA() converted a FASTA file named '", fn, "' to a NEXUS data file named '", out_nexus_aln_fn, "'.\n")
	cat(txt)
	
	return(out_nexus_aln_fn)
	} # END convert_NEXUS_to_FASTA <- function(fasta_fn, out_nexus_aln_fn=NULL)
