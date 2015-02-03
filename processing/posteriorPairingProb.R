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
parser$add_argument('-r', '--pairRatio', type="double", help="pair_ratio, the prior probability of paired bases, e.g. 0.534558833156579")
parser$add_argument('-p', '--pairParam', type="character", help="parameter file for paired bases, e.g. SHAPE/pairParameter.txt")
parser$add_argument('-l', '--loopParam', type="character", help="parameter file for single-stranded bases, e.g. SHAPE/loopParameter.txt")

args <- parser$parse_args()
if (is.null(args$infile) || is.null(args$outfile) || is.null(args$outdir) || is.null(args$pairRatio) || is.null(args$pairParam) || is.null(args$loopParam)) {  # all parameters required
	parser$print_help()
	q(status=1)
}
if (!file.exists(args$infile)) {  # infile should exist
	cat("The input file does not exists, -infile\n")
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


# Section 1, read input data
indata <- read.table(args$infile, header=T)
ncol <- ncol(indata)

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
indata$prob <- indata$pair_likeli * args$pairRatio / (indata$pair_likeli * args$pairRatio + indata$loop_likeli * (1-args$pairRatio))
write.table(indata, file=paste(args$infile, '.probDetail', sep=""), row.names=F, col.names=T, quote=F, sep="\t")  # detail file

for (i in 1:length(unique(indata$RNA))) {
	data <- subset(indata, indata$RNA==unique(indata$RNA)[i])
	write.table(data.frame(data$Index, data$prob), file=paste(args$outdir, '/', unique(indata$RNA)[i], '.prob', sep=""), row.names=F, col.names=F, quote=F, sep="\t")  # seperate .prob files
}

if (colnames(indata)[3] == 'Reactivity') {
	colnames(indata)[3] = 'Probability'
}
indata$Probability <- indata$prob
write.table(indata[,c(1:ncol)], file=args$outfile, row.names=F, col.names=T, quote=F, sep="\t")  # merged file 
