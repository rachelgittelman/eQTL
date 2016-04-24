## eQTL.R
## by Rachel M Gittelman

## This script will take an expression matrix and a genotype matrix and perform eQTL analysis.
## It assumes just one eQTL per gene, so it retains the lowest P value obtained per gene as
## the test statistic. Then it will shuffle the genotype matrix (currently set to 1000 times)
## and obtain null P values for each gene (again by calculating associations with all SNPs and
## retaining just the smallest one). 

## The expression matrix is in a standard format: m genes x n samples, tab delimited

## The genotype matrix is assumed to be converted from the VCFTOOLs --012 option, so it 
## actually requires three files:

## .012: m samples x n loci (SNPs). rownames, but no header
## .012.indv: single column where each entry corresponds to the sample name
## .012.pos: two columns that correspond to the chromosome and 1-based position of each locus

## Finally, the script requires a "gene file", a list of genes from the expression matrix to 
## include in the analysis.

## Importantly, this script does not assume the two matrices have the same set of samples, or
## that they're in the correct order

## Data from multiple eQTL analyses (different tissues and haplotypes) can be combined in to
## a single dataset for FDR analysis.

##########################################################################################
## read in data
##########################################################################################
args <- commandArgs()

## genotype data
genotype_file <- args[4]
pos_file <- args[5]
indv_file <- args[6]

## expression data
exp_file <- args[7]
gene_file <- args[8]

## name of the file to write to, and the particular cell type for this set of analyses
outfile <- args[9]
cell_type <- args[10]

nper=1000

write(paste("beginning associations for ",cell_type,sep=""),stderr())

data <- read.table(genotype_file, sep="\t", header=F, row.names=1)
pos <- read.table(pos_file, sep="\t", header=F)
indv <- read.table(indv_file, sep="\t", header=F)
exp <- read.table(exp_file, header=T, row.names=1)
genes <- scan(gene_file, what=character())


##########################################################################################
## prep data (get rid of missing genotypes), subset just the right genes to test, eliminate
## samples that aren't in the expression data, match correct samples and order
##########################################################################################

rownames(data) <- indv[,1]
colnames(data) <- paste(pos[,1],pos[,2],sep=".")

data <- data[,!apply(data,2,function(x)any(-1%in%x)), drop=F] 			  ## get rid of SNPs with any missing data (-1) 
exp <- exp[which(substr(rownames(exp),1,15) %in% genes),]                                 ## subset correct genes

data <- data[which(substr(rownames(data),1,9) %in% substr(colnames(exp),1,9)),]           ## only keep genotypes for samples that are in the expression data

order_of_exp <- match(substr(rownames(data), 1,9), substr(names(exp), 1,9))
missing_in_data <- which(is.na(order_of_exp) == TRUE)
exp <- exp[,order_of_exp[!(is.na(order_of_exp))]]                                         ## now reorder the expression data so that columns correspond to rows in the genotype data

print(dim(data))
print(dim(exp))

if(is.null(dim(exp)) || is.null(dim(data)) || nrow(exp)==0 || ncol(data)==0)              ## If, after filtering genes and missing data theres nothing left, quit.
{
	print("didn't even get started")
	write(c(cell_type, "NA", "NA", "NA"), file=outfile, append=T, ncolumns=4, sep="\t")
	quit()
}


##########################################################################################
## functions for running associations
##########################################################################################

single_association <- function(geno, exp)
{
	## Build a standard linear model between matched genotype and expression vectors
	
	a <- lm(unlist(exp) ~ unlist(geno))
	return(summary(a)$coefficients[8])
}	

gene_associations <- function(genotypes, expressions, real=T)
{
	## This function will perform all of the associations for a single gene, returning the
	## most significant SNP/Pvalue. Thus it takes the entire genotype matrix, and a single
	## vector of expression levels. real indicates whether this is permutated data or not
	
	pvals <- as.numeric(unlist(apply(genotypes, 2, single_association, exp=expressions)))
	min_index <- which(pvals == min(pvals, na.rm=T))[1]
	pval <- pvals[min_index]
	snp <- colnames(genotypes)[min_index]

 	if (real) {return(list("snp"=snp,"pval"=pval))}                                       ## Only care about the SNP in the real data
 	else return(pval)
}

##########################################################################################
## Now actually run the real associations
##########################################################################################

associations <- do.call(rbind, apply(exp, 1, gene_associations, genotypes=data))

##########################################################################################
## now permute genotype vector nper times and run associations again
##########################################################################################

set.seed(1)
empirical_null <- replicate(nper, apply(exp, 1, gene_associations, genotypes=data[sample(1:nrow(data)),], real=F), simplify=T)

##########################################################################################
## make sure data is in correct format to write out
##########################################################################################

final <- data.frame(rep(celltype,nrow(associations)),associations, empirical_null)
final <- data.frame(lapply(final, as.character), stringsAsFactors=FALSE)
rownames(final) <- rownames(associations)

write.table(final, file=outfile, append=T, sep="\t", row.names=T, col.names=F, quote=F)
