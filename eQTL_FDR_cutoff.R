## eQTl_FDR_cutoff.R
## by Rachel M Gittelman

## In this script I'll take all of the eQTL data, including true p values and the null p values
## based on permutations, and determine what P value cut off corresponds to the desired FDR

##########################################################################################
## read in data
##########################################################################################
args <- commandArgs()
data_file <- args[4]

data <- read.table(data_file, sep="\t", header=F)

GTEx_sig_permute <- function(info, dataframe, fdr=0.05)
{
	## This function will calculate a p value cut off for the given fdr based on 
	## permutation data. It will return just the significant results.
	
	actual_fdr=0.5
	cutoffs <- c(0.05,sort(info[,4][info[,4] <= 0.05],dec=TRUE))                               ## The only possible P value cutoffs are the collection of real P values below 0.05
	cutoff_ind <- 1
	cutoff <- cutoffs[cutoff_ind]
	while(actual_fdr > fdr & sum(info[,4] <= cutoff) > 0)
	{
		sig <- sum(info[,4]<=cutoff)
		prop_sig_permute <- sum(unlist(dataframe[,2:ncol(dataframe)]) <= cutoff) /
						    ncol(dataframe)
		actual_fdr <- prop_sig_permute/sig
		cutoff_ind <- cutoff_ind + 1
		cutoff <- cutoffs[cutoff_ind]
	}
	cutoff_ind <- max(1, cutoff_ind-1)
	if(actual_fdr <= 0.05)
	{
		return(info[info[,4] <= cutoffs[cutoff_ind],1:4, drop=F])
	}
}

sig_associations <- GTEx_sig_permute(data[,1:4],data[,5:ncol(data)])

if(length(sig_associations) != 0)
{
	write.table(sig_associations, file=paste(data_file,".sig_FDR_0.05",sep=""), append=T,
			    sep="\t", quote=F, col.names=F, row.names=F)
}