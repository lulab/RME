# Function: distribution fitting for normalized reactivity on training RNAs
# Author: Yang Wu (wuyang.bnu@gmail.com)
# Version: 1.0 (Jan-14, 2015)

library('argparse')
library('mixtools')
library('MASS')
library('evir')

# Get arguments
parser <- ArgumentParser(description="Distribution fitting on normalized reactivity")
parser$add_argument('-i', '--infile', type="character", help="input file for normalized reactivity, e.g. SHAPE/SHAPE.train.23SrRNA.data.normed")
parser$add_argument('-o', '--outfile', type="character", help="outfile for distribution parameters, e.g. SHAPE/pairParameter.txt")
parser$add_argument('-d', '--distribution', type="character", help="chosen distribution, e.g. gev, gamma, gaussian, 2gaussian")
parser$add_argument('-c', '--class', type="character", help="structure class, e.g. pair or loop")

args <- parser$parse_args()
if (is.null(args$infile) || is.null(args$outfile) ||is.null(args$distribution) ||is.null(args$class)) {  # all parameters required
	parser$print_help()
	q(status=1)
}
if (!file.exists(args$infile)) {  # infile should exist
	cat("The input file does not exists, -infile\n")
	q(status=1)
}


# Section 1, read input data
indata <- read.table(args$infile, header=T)
if (args$class=='pair') {
	data <- subset(indata$Reactivity, indata$Structure==1)  # distribution fitting on paired set
} else if (args$class=='loop') {
	data <- subset(indata$Reactivity, indata$Structure==0)  # distribution fitting on single-stranded set
} else {
	cat("Illegal class, pair or loop for -class\n")
	q(status=1)
}

# Section 2, fit distribution
distri <- paste("'", args$distribution, "'", sep="")
if (args$distribution=='gev') {  # fit generalized extreme value distribution
	fit <- gev(data)
	test <- ks.test(data, "pgev", fit$par.ests[1], fit$par.ests[3], fit$par.ests[2])
	param <- data.frame(c('distri', 'xi', 'mu', 'sigma', 'pvalue'), c(distri, fit$par.ests[1], fit$par.ests[3], fit$par.ests[2], test$p.value))

} else if (args$distribution=='gamma') {  # fit gamma distribuion
	fit <- fitdistr(data, "gamma")
	test <- ks.test(data, "pgamma", fit$estimate[[1]], fit$estimate[[2]])
	param <- data.frame(c('distri', 'shape', 'rate', 'pvalue'), c(distri, fit$estimate[[1]], fit$estimate[[2]], test$p.value))

} else if (args$distribution=='gaussian') {  # fit gaussain distribuion
	mu <- mean(data)
	sigma <- sd(data)
	test <- ks.test(data, "pnorm", mu, sigma)
	param <- data.frame(c('distri', 'mu', 'sigma', 'pvalue'), c(distri, mu, sigma, test$p.value))

} else if (args$distribution=='2gaussian') {  # fit gaussain mixture of 2 components
	gaussain_mix <- normalmixEM(data, k=2)
	param <- data.frame(c('distri', 'mu1', 'mu2', 'sigma1', 'sigma2', 'lambda1', 'lambda2'), c(distri, gaussain_mix$mu, gaussain_mix$sigma, gaussain_mix$lambda))

} else {
	cat("The distribution has not been implemented yet, -distribution\n")
	q(status=1)
}

# Section 3: output parameters
write.table(param, file=args$outfile, row.names=F, col.names=F, quote=F, sep="=")
