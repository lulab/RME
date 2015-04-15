# Function: calculation of the posterior pairing probability based on fitted parameters
# Author: Yang Wu (wuyang.bnu@gmail.com)
# Version: 1.0 (Jan-14, 2015)

library('argparse')
library('evir')

# Get arguments
parser <- ArgumentParser(description="Calculation of the posterior pairing probability")
parser$add_argument('-i', '--infile', type="character", help="input file for normalized reactivity, e.g. SHAPE/SHAPE.train.23SrRNA.data.normed")
parser$add_argument('-o', '--outfile', type="character", help="outfile for posterior pairing probability, e.g. SHAPE.txt")
parser$add_argument('-O', '--outdir', type="character", help="output directory for probability files, e.g. SHAPE")
parser$add_argument('-r', '--pairRatio', type="character", help="file containing the pair_ratio information, either a constant or different values for each RNA")
parser$add_argument('-p', '--pairParam', type="character", help="parameter file for paired bases, e.g. SHAPE/pairParameter.txt")
parser$add_argument('-l', '--loopParam', type="character", help="parameter file for single-stranded bases, e.g. SHAPE/loopParameter.txt")
parser$add_argument('-m', '--modeOutlier', type="character", help="mode for dealing with outliers, DEFAULT is use.  'use': use outliers as normal data; 'cap': cap outliers to specified boundary (e.g. small outliers are reset to bottomOutlier and large outliers are reset to topOutlier; filter': filter outliers defined by -t and -m, OPTIONAL")
parser$add_argument('-b', '--bottomOutlier', type="double", help="data smaller than bottomOutlier is defined as outliers, OPTIONAL")
parser$add_argument('-t', '--topOutlier', type="double", help="data larger than topOutlier is defined as outliers, OPTIONAL")

args <- parser$parse_args()
if (is.null(args$infile) || is.null(args$outfile) || is.null(args$outdir) || is.null(args$pairRatio) || is.null(args$pairParam) || is.null(args$loopParam) || is.null(args$modeOutlier) || is.null(args$bottomOutlier) || is.null(args$topOutlier)) {  # all parameters required
	parser$print_help()
	q(status=1)
}
if (!file.exists(args$infile)) {  # infile should exist
	cat("The input file does not exists, -infile\n")
	q(status=1)
}
if (!file.exists(args$pairRatio)) {  # pairRatio should exists
	cat("The pair_ratio file does not exists, -pairRatio\n")
	q(status=1)
}
if (!file.exists(args$pairParam)) {  # pairParam should exists
	cat("The param file for paired bases does not exists, -pairParam\n")
	q(status=1)
}
if (!file.exists(args$loopParam)) {  # loopParam should exists
	cat("The param file for single-stranded bases does not exists, -loopParam\n")
	q(status=1)		
}
if (args$bottomOutlier>args$topOutlier) {  # check relative size for define outliers
	cat("bottomOutlier should be smaller than topOutlier\n")
	q(status=1)
}

# Function, get constant pair_ratio or variable pair_ratio for each RNA
getPairRatio <- function(pairRatioFile, RNA) {
	c <- read.table(pairRatioFile)
	if (nrow(c)==1 && ncol(c)==1) {  # P(S=1) as a constant
		return(c[1,1])  # one value for all RNAs
	} else if (ncol(c==2)) {  # variable P(S=1) for each RNA
		return(as.numeric(c[c[,1]==RNA,2]))  # two-column file specify a value for each RNA
	} else {
		cat(paste("Incorrect format for file", pairRatioFile, "\n", sep=""))
		q()
	}
}

# Function, deal with outliers
capData <- function (data) {
	data[data<args$bottomOutlier] <- args$bottomOutlier
	data[data>args$topOutlier] <- args$topOutlier
	return(data)
}
filterData <- function(indata) {
	filter <- indata$Reactivity>=args$bottomOutlier & indata$Reactivity<=args$topOutlier
	return(filter)
}


# Section 1, read input data
indata <- read.table(args$infile, header=T)
ncol <- ncol(indata)
if (args$modeOutlier=='cap') {  # cap outliers to specified outliers
	indata$Reactivity <- capData(indata$Reactivity)
} else if (args$modeOutlier=='filter') {  # filter outliers
	filter <- filterData(indata)
	indata <- indata[filter,]
} 


# Section 2, calculate likelihoods for paired bases
source(args$pairParam)
if (distri=='gev') {  # calculate likelihoods based on generalized extreme value distribution
	indata$pair_likeli <- dgev(indata$Reactivity, xi, mu, sigma)
} else if (distri=='gaussian') {  # calculate likelihoods based on gaussian distribution
	indata$pair_likeli <- dnorm(indata$Reactivity, mu, sigma)
} else {
	cat(paste("The distribution has not been implemented yet, ", distri, "\n", sep=""))
	q(status=1)
}

# Section 3, calculate likelihoods for single-stranded bases
source(args$loopParam)
if (distri=='gamma') {  # calculate likelihoods based on gamma distribution
	indata$loop_likeli <- dgamma(indata$Reactivity, shape, rate)
} else if (distri=='gaussian') {  # calculate likelihoods based on gaussian distribution
	indata$loop_likeli <- dnorm(indata$Reactivity, mu, sigma)
} else if (distri=='2gaussian') { # calculate likelihoods based on gaussian mixture of 2 components 
	loop_likeli1 <- dnorm(indata$Reactivity, mu1, sigma1)
	loop_likeli2 <- dnorm(indata$Reactivity, mu2, sigma2)
	indata$loop_likeli <- loop_likeli1 * lambda1 + loop_likeli2 * lambda2
} else {                                                          
	cat(paste("The distribution has not been implemented yet, ", distri, "\n", sep=""))
	q(status=1)                                                   
}

# Section 4, calculate posterior pairing probability
out <- indata[0,]
for (i in 1:length(unique(indata$RNA))) {
    RNA=as.character(unique(indata$RNA)[i])
	data <- indata[indata$RNA==RNA,]
	pair_ratio <- getPairRatio(args$pairRatio, RNA)
	
	data$prob <- data$pair_likeli * pair_ratio / (data$pair_likeli * pair_ratio + data$loop_likeli * (1-pair_ratio))	
	write.table(data.frame(data$Index, data$prob), file=paste(args$outdir, '/', unique(indata$RNA)[i], '.prob', sep=""), row.names=F, col.names=F, quote=F, sep="\t")  # seperate .prob files
	out <- rbind(out, data)
}
write.table(out, file=paste(args$infile, '.probDetail', sep=""), row.names=F, col.names=T, quote=F, sep="\t")  # detail file

if (colnames(out)[3] == 'Reactivity') {
	colnames(out)[3] = 'Probability'
}
out$Probability <- out$prob
write.table(out[,c(1:ncol)], file=args$outfile, row.names=F, col.names=T, quote=F, sep="\t")  # merged file 
