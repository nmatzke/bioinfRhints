
library(ape)

newickstr2tntstr <- function(newickstr)
	{
	# 1. Change any ; to *
	tntstr1 = gsub(pattern=";", replacement="*", x=newickstr)
	
	# 2. Change any commas to whitespace
	tntstr2 = gsub(pattern=",", replacement=" ", x=tntstr1)

	# 3. # replace "),(" with ")("
	tntstr3 = gsub(pattern="\\),\\(", replacement="\\)\\(", x=tntstr2)

	# 3. # replace ") (" with ")("
	tntstr3a = gsub(pattern="\\) \\(", replacement="\\)\\(", x=tntstr3)

	# 4. replace ")" with " )"
	tntstr4 = gsub("\\)", " \\)", tntstr3a)
		
	# 4. replace ") )" with ")"
	tntstr5 = gsub("\\) \\)", "\\)\\)", tntstr4)
	tntstr6 = gsub("\\) \\)", "\\)\\)", tntstr5)
	
	tntstr = tntstr6
	
	return(tntstr6)
	}




newickstr = "('Wuhan SARs-CoV-2':0.0002654909,('USA SARs-CoV-2':0.0016976448,'2022 SARs-CoV-2':0.0030846578):0.0034879598,(Bat_cov_RaTG13:0.0117870014,((Pangolin_cov_1:0.0012653324,Pangolin_cov_2:0.0010022307):0.0705829809,('Bat sars-like cov':0.037907956,((('SARs-CoV BJ01':0.0002396983,'SARs-CoV TWS':0.0005233065):0.0009289262,'SARs-CoV civet':0.0037344069):0.0880014726,(Rousettus_bat_cov:0.8486810264,(('MERs-CoV 1':0.002137137,'MERs-CoV 2':0.006704995):0.8413245682,((((((Bovine_cov_2014:0.0075475269,Bovine_cov_2020:0.0078516395):0.0037348025,Canine_respiratory_cov:0.0153630643):0.0074039671,Human_cov_OC43:0.0418889172):0.0619887321,Equine_cov:0.0750373621):0.1287045342,(Murine_cov:0.1766622784,Human_cov_HKU1:0.2433501012):0.0660268207):0.9476123939,((((((Porcine_deltacov:0.0874721006,Quail_deltacov:0.0918300683):0.2692631739,Common_moorhen_cov:0.302041813):0.4002992604,Wigeon_cov:0.6822117753):0.9319133582,(Breda_virus:9.7852800017,Atlantic_salmon_bafinivirus:6.6040375559):9.9999987954):0.4524232245,((Turkey_cov:0.0466033214,Guinea_fowl_cov:0.0683218043):0.8737894405,Bottlenose_dolphin_cov:1.1480907414):0.450777981):0.3883492555,((Canine_cov:0.1330088547,Feline_cov:0.1502421848):0.4902183526,(Human_Cov_NL63:0.4722117347,Porcine_epidemic_diarrhea:0.4190487373):0.2553774461):1.0002527327):0.6027792331):0.2462945656):0.1788795367):0.690465841):0.0658034033):0.0283779573):0.0161770895):0.0085110003);"

tr = read.tree(file="", text=newickstr)

# Remove branch lengths
tr$edge.length = NULL

tr$tip.label



# Replace spaces with "_"
newnames1 = gsub(pattern=" ", replacement="_", x=tr$tip.label)

# Replace "-" with ""
newnames2 = gsub(pattern="-", replacement="", x=newnames1)

# Replace "'" with ""
newnames3 = gsub(pattern="'", replacement="", x=newnames2)


tr$tip.label = newnames3

newickstr2 = write.tree(tr, file="")
tntstr = newickstr2tntstr(newickstr2)
tntstr


# THen place that tntstr into a .tnttree file

# Check that no tip names start with a number

# Run the TNT commands

# Process the output logfile;
wd = "~/Downloads/Brynn2/"
setwd(wd)
#fn = "apomorphies_v1.txt"
wd = "/Users/nickm/Downloads/retnthelp/"
setwd(wd)
sars_cov2_apomorphies_from_tnt_fn = "apomorphies_v3.txt"
coronaviridae_apomorphies_from_tnt_fn = "apomorphies_v4.txt"


#######################################################
# SARS-CoV-2 apomorphies
#######################################################
fn = sars_cov2_apomorphies_from_tnt_fn
tmplines = readLines(fn)

# First, let's eliminate the junk

tmplines2 = gsub(pattern="      Char. ", replacement="", x=tmplines)
tmplines2

# Make a blank dataframe
tdf = data.frame(matrix(data="", nrow=length(tmplines2), ncol=4))

tdf_row = 0
for (i in 1:length(tmplines2))
	{
	cat(i, ",", sep="")
	tmpline = tmplines2[i]
	# Get first character of tmpline
	first_char = substr(tmpline, 1,1)
	
	is_numberTF = !is.na(as.numeric(first_char))
	if (is_numberTF == FALSE)
		{
		tmpline_nospaces = gsub(pattern=" ", replacement="", x=tmpline)
		tmpline_nospaces = gsub(pattern="\\:", replacement="", x=tmpline_nospaces)
		if (tmpline_nospaces != "")
			{
			node_ID = tmpline_nospaces
			}
		}

	if (is_numberTF == TRUE)
		{
		# Split the string
		# tmpline "929: K --> R "
		tmpline = gsub(pattern="\\:", replacement=" ", x=tmpline)
		words = strsplit(x=tmpline, split=" ")[[1]]
		site_position = as.numeric(words[1])
		starting_AA = words[3]
		ending_AA = words[5]
		
		# Store the data from this line
		tdf_row = tdf_row + 1
		tdf[tdf_row, 1] = node_ID
		tdf[tdf_row, 2] = site_position
		tdf[tdf_row, 3] = starting_AA
		tdf[tdf_row, 4] = ending_AA
		}
	}

head(tdf)
names(tdf) = c("node_ID", "site_position", "starting_AA", "ending_AA")
tdf$site_position = as.numeric(tdf$site_position)


# Sort by position
new_order = order(tdf$site_position)
tdf2 = tdf[new_order,]
head(tdf2)

outfn = gsub(pattern="\\.txt", replacement="_subsdf.txt", x=fn)
write.table(tdf2, file=outfn, quote=FALSE, sep="\t", row.names=FALSE)


# Convert site positions to start with #1, go to max site position
tdf$site_position = tdf$site_position+1
max_site_position = max(tdf$site_position, na.rm=TRUE)

# Actually, let's MANUALLY SET a maximum site position, based on TNT data file
max_site_position = 3975

# Count # of changes
counts = rep(0, times=max_site_position)

for (i in 1:max_site_position)
	{
	TF = tdf$site_position == i
	counts[i] = sum(TF, na.rm=TRUE)
	}
pos = 1:max_site_position
num_changes_df = cbind(pos, counts)
num_changes_df = as.data.frame(num_changes_df, stringsAsFactors=FALSE)

outfn = gsub(pattern="\\.txt", replacement="_num_changes_sdf.txt", x=fn)
write.table(num_changes_df, file=outfn, quote=FALSE, sep="\t", row.names=FALSE)

sars_cov2_substitutions = tdf2
sars_cov2_num_changes_df = num_changes_df





#######################################################
# Coronaviridae apomorphies
#######################################################
fn = coronaviridae_apomorphies_from_tnt_fn
tmplines = readLines(fn)

# First, let's eliminate the junk

tmplines2 = gsub(pattern="      Char. ", replacement="", x=tmplines)
tmplines2

# Make a blank dataframe
tdf = data.frame(matrix(data="", nrow=length(tmplines2), ncol=4))

tdf_row = 0
for (i in 1:length(tmplines2))
	{
	cat(i, ",", sep="")
	tmpline = tmplines2[i]
	# Get first character of tmpline
	first_char = substr(tmpline, 1,1)
	
	is_numberTF = !is.na(as.numeric(first_char))
	if (is_numberTF == FALSE)
		{
		tmpline_nospaces = gsub(pattern=" ", replacement="", x=tmpline)
		tmpline_nospaces = gsub(pattern="\\:", replacement="", x=tmpline_nospaces)
		if (tmpline_nospaces != "")
			{
			node_ID = tmpline_nospaces
			}
		}

	if (is_numberTF == TRUE)
		{
		# Split the string
		# tmpline "929: K --> R "
		tmpline = gsub(pattern="\\:", replacement=" ", x=tmpline)
		words = strsplit(x=tmpline, split=" ")[[1]]
		site_position = as.numeric(words[1])
		starting_AA = words[3]
		ending_AA = words[5]
		
		# Store the data from this line
		tdf_row = tdf_row + 1
		tdf[tdf_row, 1] = node_ID
		tdf[tdf_row, 2] = site_position
		tdf[tdf_row, 3] = starting_AA
		tdf[tdf_row, 4] = ending_AA
		}
	}

head(tdf)
names(tdf) = c("node_ID", "site_position", "starting_AA", "ending_AA")
tdf$site_position = as.numeric(tdf$site_position)


# Sort by position
new_order = order(tdf$site_position)
tdf2 = tdf[new_order,]
head(tdf2)

outfn = gsub(pattern="\\.txt", replacement="_subsdf.txt", x=fn)
write.table(tdf2, file=outfn, quote=FALSE, sep="\t", row.names=FALSE)


# Convert site positions to start with #1, go to max site position
tdf$site_position = tdf$site_position+1
max_site_position = max(tdf$site_position, na.rm=TRUE)

# Actually, let's MANUALLY SET a maximum site position, based on TNT data file
max_site_position = 3975

# Count # of changes
counts = rep(0, times=max_site_position)

for (i in 1:max_site_position)
	{
	TF = tdf$site_position == i
	counts[i] = sum(TF, na.rm=TRUE)
	}
pos = 1:max_site_position
num_changes_df = cbind(pos, counts)
num_changes_df = as.data.frame(num_changes_df, stringsAsFactors=FALSE)

outfn = gsub(pattern="\\.txt", replacement="_num_changes_sdf.txt", x=fn)
write.table(num_changes_df, file=outfn, quote=FALSE, sep="\t", row.names=FALSE)

coronaviridae_substitutions = tdf2
coronaviridae_num_changes_df = num_changes_df



#######################################################
# Plot comparing the two 
#######################################################
xvals = jitter(coronaviridae_num_changes_df$counts)
yvals = jitter(sars_cov2_num_changes_df$counts)
plot(xvals, yvals)

# Pretty messy, let's try...
combined_df = as.data.frame(cbind(sars_cov2_num_changes_df$counts, coronaviridae_num_changes_df$counts), stringsAsFactors=FALSE)
names(combined_df) = c("sars_cov2_changes", "coronaviridae_changes")
head(combined_df)


# This suggests something!
means_by_num_sars2_substitutions = aggregate(x=combined_df$coronaviridae_changes, by=list(combined_df$sars_cov2_changes), FUN=mean)
means_by_num_sars2_substitutions

# Lots of noise though:
lower025_by_num_sars2_substitutions = aggregate(x=combined_df$coronaviridae_changes, by=list(combined_df$sars_cov2_changes), FUN=quantile, probs=0.025)
lower025_by_num_sars2_substitutions

upper025_by_num_sars2_substitutions = aggregate(x=combined_df$coronaviridae_changes, by=list(combined_df$sars_cov2_changes), FUN=quantile, probs=0.975)
upper025_by_num_sars2_substitutions


# Generalized linear model for count data (responses are poisson
# Model: the rate of SARS2 substitutions (0-3 observed) is predicted by the 
#        number of Coronaviridae substitutions (0-12 subs. observed)
res = glm(formula=sars_cov2_changes~coronaviridae_changes, data=combined_df, family=poisson(link = "log"))
summary(res) # not significant

res = glm(formula=coronaviridae_changes~sars_cov2_changes, data=combined_df, family=poisson(link = "log"))
summary(res) # not significant


# Does having substitutions predict having substitutions?
cutoff_sars = 0
cutoff_coronoviridae = 0

sars_mutTF = combined_df$sars_cov2_changes > cutoff_sars
coro_mutTF = combined_df$coronaviridae_changes > cutoff_coronoviridae
TFdf = as.data.frame(cbind(sars_mutTF, coro_mutTF), row.names=NULL)
res = aov(formula=coro_mutTF~sars_mutTF, data=TFdf)
summary(res) # Significant!

# Logistic regression (0/1 data)
res = glm(formula=coro_mutTF~sars_mutTF, data=TFdf, family=binomial(link=logit))
summary(res)

# "The logistic regression coefficients give the change in the 
# log odds of the outcome for a one unit increase in the predictor variable."
# https://stats.oarc.ucla.edu/r/dae/logit-regression/

# Go from probability to log-odds
logit <- function(p)
	{
	log(p/(1-p))	# It is the inverse of the expit function.
	}

# Go from log-odds to probability
expit <- function(x)
	{
	1/(1+exp(-x)) # It is the inverse of the logit function.
	}

logit(0.5)
expit(0.0)

# Coefficients
coef(res)
confint(res)

# Base probability of having at least one substitution in Coronoviridae
# (when sars_mutTF==0, i.e., when SARS-CoV-2 has 0 substutitions)
expit(coef(res)[1])
# 95% CI:
expit(confint(res)[1,])

# Modifier to probability of having at least one substitution in Coronoviridae
# (when sars_mutTF==1, i.e., when SARS-CoV-2 has 1 or more substutitions)
expit(coef(res)[2])
# 95% CI:
expit(confint(res)[2,])


odds_ratio_intercept = exp(coef(res)[1])
odds_ratio_coro1 = exp(coef(res)[2])
odds_ratio_ttl = odds_ratio_intercept * odds_ratio_coro1
prob_total = expit(log(odds_ratio_ttl))
prob_total

# Same as
expit(coef(res)[1] + coef(res)[2])

# So, observing at least 1 substitution in Coronoviridae at a location
# changes our probability of seeing at least 1 substitution in SARS-CoV2 from
expit(coef(res)[1])
# ...to...
expit(coef(res)[1] + coef(res)[2])


#######################################################
# You can try variations on all this, ie:
#######################################################
# * Different cutoffs
# * predicting the reverse
# * more detailed categories of predictors (0 vs. 1. vs more)







