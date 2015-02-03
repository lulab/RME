# Function: normalize raw read counts for different RNA by log-abundance of the RNA it belongs to
# Author: Yang Wu (wuyang.bnu@gmail.com)
# Version: 1.0 (Jan-14, 2015)

library('argparse')

# Get arguments
parser <- ArgumentParser(description="Read count normalization on high-throughput data")
parser$add_argument('-i', '--infile', type="character", help="input file for raw read counts, e.g. ../00.dat/data/PARS.train.25SrRNA.data")
parser$add_argument('-o', '--outfile', type="character", help="outfile for normalized reactivity, e.g. PARS/PARS.train.25SrRNA.data.normed")

args <- parser$parse_args()
if (is.null(args$infile) || is.null(args$outfile)) {  # all parameters required
	parser$print_help()
	q(status=1)
}
if (!file.exists(args$infile)) {  # infile should exist
	cat("The input file does not exists, -infile\n")
	q(status=1)
}

# Define function
normAbunLen <- function(data) {
	data <- log(data + 1)
	sum <- sum(data)
	scale <- sum/length(data)
	data <- data/scale
	return(data)
}


# Section 1, read input data
indata <- read.table(args$infile, header=T)

# Section 2, normalization for each RNA, normed count 1 and 2, then subtraction
normed <- data.frame(indata[0,], norm1=double(), norm2=double(), Reactivity=double())
for (i in 1:length(unique(indata$RNA))) {
	data <- subset(indata, indata$RNA==unique(indata$RNA)[i])
	data$norm1 <- normAbunLen(data[,3])
	data$norm2 <- normAbunLen(data[,4])
	data$Reactivity <- data$norm2-data$norm1
	normed <- rbind(normed, data)
}

# Section 3, output normed reactivity
if (ncol(indata)==5) {
	out <- data.frame(RNA=normed$RNA, Index=normed$Index, Reactivity=normed$Reactivity, Base=normed$Base)
} else if (ncol(indata)==6) {
	out <- data.frame(RNA=normed$RNA, Index=normed$Index, Reactivity=normed$Reactivity, Base=normed$Base, Structure=normed$Structure)
} else {
	cat("Illegal input format, at least 5 columns (RNA, Index, Count1, Count2, Base) and 1 optional columns (Structure)\n")
	q(status=1)                                                     
}
write.table(out, file=args$outfile, row.names=F, col.names=T, sep="\t", quote=F)  
write.table(normed, file=paste(args$outfile, '.scaleDetail', sep=""), row.names=F, col.names=T, sep="\t", quote=F)        
