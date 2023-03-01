

# Example data
sars_cov2_changes = c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1)
coronaviridae_changes = c(0,1,0,0,0,0,0,0,0,1,0,1,2,1,1,1,1,0,2,5)

# Make a data.frame
df = as.data.frame(cbind(sars_cov2_changes, coronaviridae_changes))
df


# Subset and histograms

# 0 changes in coronaviridae
subset_df_when_coronaviridae_EQ_0 = subset(x=df, subset=(coronaviridae_changes==0))

# 1 or more changes in coronaviridae (GT_0 = # changes greater than 0)
subset_df_when_coronaviridae_GT_0 = subset(x=df, subset=(coronaviridae_changes>0))

# Histograms of each, on same plot

# Blank plot with 2 rows and 1 column
par(mfrow=c(2,1))

xmax = max(df$sars_cov2_changes)
breakvals = 0:(xmax+1)

# Histogram when Coronaviridae changes == 0 
titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae=0"
hist(x=subset_df_when_coronaviridae_EQ_0$sars_cov2_changes, breaks=breakvals, right=FALSE, freq=TRUE, xlab="# of SARS-CoV-2 changes", main=titletxt)


# Histogram when Coronaviridae changes > 0 
titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae>0"
hist(x=subset_df_when_coronaviridae_GT_0$sars_cov2_changes, breaks=breakvals, right=FALSE, freq=TRUE, xlab="# of SARS-CoV-2 changes", main=titletxt)

# Another way, with barplot, looks better when the values are just counts of 0 and 1
titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae=0"
counts = table(subset_df_when_coronaviridae_EQ_0$sars_cov2_changes)
barplot(counts, main=titletxt)

titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae>0"
counts = table(subset_df_when_coronaviridae_GT_0$sars_cov2_changes)
barplot(counts, main=titletxt)



# Do it the reverse way -- subset based on the # of SARS-CoV-2 changes being 0 or >0

# 0 changes in SARS-CoV-2
subset_df_when_sarscov2_EQ_0 = subset(x=df, subset=(sars_cov2_changes==0))

# >0 changes in SARS-CoV-2
subset_df_when_sarscov2_GT_0 = subset(x=df, subset=(sars_cov2_changes>0))

# Histograms of each, on same plot

# Blank plot with 2 rows and 1 column
par(mfrow=c(2,1))

# Histogram when Coronaviridae changes == 0 
xmax = max(df$coronaviridae_changes)
breakvals = 0:xmax
titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae=0"
hist(x=subset_df_when_sarscov2_EQ_0$coronaviridae_changes, breaks=breakvals, right=FALSE, freq=TRUE, xlab="# of Coronaviridae changes", xlim=c(0,xmax), main=titletxt)


# Histogram when Coronaviridae changes > 0 
titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae>0"
hist(x=subset_df_when_sarscov2_GT_0$coronaviridae_changes, breaks=breakvals, right=FALSE, freq=TRUE, xlab="# of Coronaviridae changes", xlim=c(0,xmax), main=titletxt)




# Histograms of each, on same plot - DENSITY VERSION

# Blank plot with 2 rows and 1 column
par(mfrow=c(2,1))

# Histogram when Coronaviridae changes == 0 
titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae=0"
hist(x=subset_df_when_sarscov2_EQ_0$coronaviridae_changes, breaks=breakvals, right=FALSE, freq=FALSE, xlab="# of Coronaviridae changes", xlim=c(0,xmax), main=titletxt)


# Histogram when Coronaviridae changes > 0 
# freq=FALSE makes the y axis the relative probability density rather than the count

titletxt = "# of SARS-CoV-2 changes when # of changes in Coronaviridae>0"
hist(x=subset_df_when_sarscov2_GT_0$coronaviridae_changes, breaks=breakvals, right=FALSE, freq=FALSE, xlab="# of Coronaviridae changes", xlim=c(0,xmax), main=titletxt)
