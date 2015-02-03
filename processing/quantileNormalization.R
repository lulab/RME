# Function: normalize SHAPE reactivity for different RNAs by Quantile normalization
# Author: Yang Wu (wuyang.bnu@gmail.com)
# Version: 1.0 (Jan-14, 2015)

library('preprocessCore')
library('argparse')
library('reshape')

# Get arguments
parser <- ArgumentParser(description="Quantile normalization on SHAPE data")
parser$add_argument('-i', '--infile', type="character", help="input file for raw reactivity, e.g. ../00.dat/data/SHAPE.train.23SrRNA.data")
parser$add_argument('-o', '--outfile', type="character", help="outfile for normalized reactivity, e.g. SHAPE/SHAPE.train.23SrRNA.data.normed")
parser$add_argument('-n', '--normfile', type="character", help="file storing the normalization target, e.g. SHAPE/quantile-target.txt")
parser$add_argument('-m', '--mode', type="character", help="run mode, train or test")

args <- parser$parse_args()
if (is.null(args$infile) || is.null(args$outfile) ||is.null(args$normfile) ||is.null(mode)) {  # all parameters required
	parser$print_help()
	q(status=1)
}
if (!file.exists(args$infile)) {  # infile should exist
	cat("The input file does not exists, -infile\n")
	q(status=1)
}
if (args$mode != 'train' & !file.exists(args$normfile)) {  # normfile should exist unless using train mode
	cat("Please prepare the target file, -normfile, unless in train mode\n")
	q(status=1)
}


# Section 1, read input data
indata <- read.table(args$infile, header=T)

# Section 2, quantile normalization, get target based on training RNAs
if (args$mode == 'train') {
	tmp <- data.frame(Index=indata$Index, Name=indata$RNA, Value=indata$Reactivity)
	m <- cast(tmp, Index ~ Name)
	target <- normalize.quantiles.determine.target(as.matrix(m[,-1]))  # get normalization target using training RNAs
	write.table(target, file=args$normfile, quote=F, row.names=F, col.names=F)
}
target <- read.table(args$normfile)

# Section 3, quantile normalization, apply to individual RNA
out <-indata[0,]  
for (i in 1:length(unique(indata$RNA))) {
	data <- subset(indata, indata$RNA==unique(indata$RNA)[i])
	tmp <- normalize.quantiles.use.target(as.matrix(data$Reactivity), target$V1)  # quantile normalization on each RNA
	data$Reactivity <- tmp
	out <- rbind(out, data)
}
write.table(out, file=args$outfile, quote=F, row.names=F, col.names=T, sep="\t")
