#######################################################
# TNT commands in R
#######################################################

# Attempt to change:
#cd /drives/GDrive/__GDrive_projects/2014-11-21_Randy_Irmis_autapomorphies/_02_tnt/
#tnt bground proc dataSEA_df.nex, log dataSEA_df.tnt, echo =, echo ], xread *, quote .,, quote proc /.,, quote comments 0, quote .,, echo -, quit, &


# Ordering, inactivate uninformative characters
#tnt proc supermatrix_v8_simp.tnt, outgroup Morotopithecus, ccodes, xinact, auto, zzz, ; open auto_logfile.txt;
#Rscript _02_read_autorun_v1.R


# Do basic TNT analysis
# WORKED on 2015-08-27
cd /drives/GDrive/__GDrive_projects/2016-06-16_venerid3/02_TNT/
tnt proc v130613mNNs_nmo_morph_wdata.tnt, xinact, auto, zzz, ; open auto_logfile.txt;
Rscript _02_read_autorun_v1.R







# 2023-05-29_
# MotA refined alignment


cd /Users/nmat471/HD/GitHub/bioinfRhints/flag/589_MotA_rename/CD-search/mafft_constrained_3/02_tnt/
xinact


# taxname +200 ; /* allow taxon names-tags of at 200 */


tnt mxram 128, sectsch:slack 50, taxname +200, proc data.tnt, auto, zzz, ; open auto_logfile.txt;
Rscript _02_read_autorun_v1.R





tnt
mxram 128;
sectsc:slack 28;
proc data.tnt;
auto;
zzz;


tnt mxram 128, sectsc:slack 28, proc data.tnt, auto, ;
zzz;