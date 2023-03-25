#!/usr/bin/env Rscript

# Author: Ulrich Lautenschlager
# Basic usage: Rscript 4_index3entropy_index4silhouette.r -k [desired k-mer length, default=8] -l [desired linkage criterion, default=single] -o [output_prefix] input_file
# For recommendations on choosing k-mer lengths and linkage criteria, please refer to the supplementary extended methods document of the paper
# NOTE: The individuals' IDs in the input sequence names ('BC02', 'BC04', 'BC05', 'BC06' in the present study) should conform to those used in script 3_index2_KBS.r
# Further information is provided by the help option: Rscript 4_index3entropy_index4silhouette.r -h

library(ape)
library(argparse)
library(cluster)  # silhouette()
library(kmer)

plogp <- function(p, base=2) {
    ifelse(p==0, 0, p * log(p, base=base))
}

# Shannon entropy
entropy <- function(x, base=2) {
    if (length(x) == 1) {
        return(NA)
    }
    p <- x/sum(x)
    H <- -sum(plogp(p, base=base))
    return(H)
}

# computation of k-mer distances considering sequence orientation of input consensuses
kdist_min <- function(alignment, ...) {
    n_seq <- length(alignment)
    m <- matrix(NA, nrow=n_seq, ncol=n_seq, dimnames=list(names(alignment), names(alignment)))
    if (n_seq < 2) {
        return(NA)
    }
    for (i in 1:(n_seq-1)) {
        for (j in (i+1):n_seq) {
            kd1 <- kmer::kdistance(alignment[c(i,j)], ...)
            kd2 <- kmer::kdistance(c(alignment[i], ape::complement(alignment[j])), ...)
            m[j,i] <- min(kd1, kd2) # lower triangle
        }
    }
    return(as.dist(m))
}

# create argument parser
parser <- ArgumentParser(description="This script allows to hierarchically cluster sequences based on kmer distances, estimating the optimal number of clusters based on the average silhouette width.")
parser$add_argument("-k", "--kmer_length", help="kmer length", type="integer", default=8)
parser$add_argument("-l", "--linkage", help="linkage criterion", type="character", default="single")
parser$add_argument("-o", "--output_prefix", help="prefix for output files (allows shorter file names)", type="character", default="")
parser$add_argument("-r", "--round", help="round average silhouette widths to a specified number of decimal places - this may affect the proposed values of k (disabled by default)", type="integer", default=-1)
parser$add_argument("infile", help="input file (FASTA)", type="character")

# parse command line arguments
argv <- parser$parse_args()

if (argv$output_prefix == "") {
    argv$output_prefix <- sub("\\.fa$|\\.fna$|\\.fsa$|\\.fasta$", "", argv$infile)
}

# read input
alignment <- read.FASTA(argv$infile, type="DNA")
if (is.null(alignment) || length(alignment)<=2) {
    cat("skip \"", infile, "\" (less than 3 input sequences)\n", sep="")
    quit()
}

# distance matrix
d <- kdist_min(alignment, k=argv$kmer_length)
m <- as.matrix(d)

# hierarchical clustering:
h <- hclust(d, method = argv$linkage)

# calculate silhouette coefficients for different cuts through the dendrogram
sil_coeffs <- NULL
for (k in 2:(length(alignment)-1)) {
    si <- silhouette(cutree(h, k), d)
    sil_avg <- mean(si[,3])
    if (argv$round >= 0) {
        sil_avg <- round(sil_avg, argv$round)
    }
    sil_coeffs <- c(sil_coeffs, sil_avg)
}

# plot graphical representation of silhouette coefficients as PDF file
pdf(paste(argv$output_prefix, "_silh_", argv$linkage, "_", argv$kmer_length, "mer", ".pdf", sep=""))
plot(2:(length(h$labels)-1), sil_coeffs, type="l", xlab="number of clusters", ylab="silhouette coefficient")
invisible(dev.off())

# write values of k and silhouette coefficients as CSV file 
outfile <- paste(argv$output_prefix, "_silh_", argv$linkage, "_", argv$kmer_length, "mer", ".csv", sep="")
outtable <- cbind(2:(length(h$labels)-1), sil_coeffs)
colnames(outtable)[1] <- "k"
write.table(outtable, file=outfile, quote=F, sep="\t", row.names=F, col.names=T)

print(outtable)

best_ks <- which(sil_coeffs==max(sil_coeffs)) + 1
cat("\nbest values of k: ")
cat(best_ks, "\n")

cat("\ncluster frequencies for k=", best_ks[1], ":\n", sep="")

# calculate entropy (mean cluster uncertainty of individuals)
# sequence names are assumed to begin with the individual id, everything after the first "_" is ignored
partition <- cutree(h, best_ks[1])
individuals <- sub("_.*$", "", attr(partition, "names"))
ent_weighted_sum <- 0
for (ind in unique(individuals)) {
    cluster_counts <- as.vector(table(factor(partition[individuals == ind], levels=1:best_ks[1])))
    ent_weighted_sum <- ent_weighted_sum + sum(cluster_counts) * entropy(cluster_counts)
    cat(ind, ": ", cluster_counts, "\n")
}
cat("entropy: ", ent_weighted_sum/length(alignment), " (weighted mean across individuals)\n", sep="")
