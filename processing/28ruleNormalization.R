# Function: normalize SHAPE reactivity for different RNAs by 2-8% rule
# Author: Yang Wu (wuyang.bnu@gmail.com)
# Version: 1.0 (Mar-2, 2015)

library('argparse')
library('reshape')

# Get arguments
parser <- ArgumentParser(description="2-8% normalization on SHAPE data")
parser$add_argument('-i', '--infile', type="character", help="input file for raw reactivity, e.g. ../00.dat/data/SHAPE.train.23SrRNA.data")
parser$add_argument('-o', '--outfile', type="character", help="outfile for normalized reactivity, e.g. SHAPE/SHAPE.train.23SrRNA.data.normed")

args <- parser$parse_args()
if (is.null(args$infile) || is.null(args$outfile)) {  # all parameters required
    parser$print_help()
    q(status=1)
}
if (!file.exists(args$infile)) {  # infile should exist
    cat("The input file does not exists, -infile\n")
    q(status=1)
}

scaleData28 <- function (data) {
    order <- data[order(data, decreasing=T)]
	upper2 <- floor(length(order)*0.02)+1
	upper10 <- floor(length(order)*0.1)
    scale <- sum(order[upper2:upper10])/(upper10-upper2+1)
#    print(scale)
    data <- round(data/scale, 4)
    return(data)
}



# Section 1, read input data
indata <- read.table(args$infile, header=T)

# Section 2, 2-8% normalization, apply to individual RNA
out <-indata[0,]  
for (i in 1:length(unique(indata$RNA))) {
    data <- subset(indata, indata$RNA==unique(indata$RNA)[i])
#    print(as.character(unique(indata$RNA)[i]))
    data$Reactivity <- scaleData28(data$Reactivity)
    out <- rbind(out, data)
}
write.table(out, file=args$outfile, quote=F, row.names=F, col.names=T, sep="\t")
