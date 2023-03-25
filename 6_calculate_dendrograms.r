# Author: Tankred Ott
# Usage: Replace "/path/to/candidate_loci" at the beginning with your own input path
# In lines 36 and 37, "k=8" or "method = 'single'" can be accordingly modified in case a different k-mer length or linkage criterion for distance/dendrogram computation is desired, but refer to the supplementary extended methods document of the paper for recommendations; be sure to rename the output suffix in line 32 accordingly

library(kmer)
library(ape)
library(tools)

in_dir <- '/path/to/candidate_loci'
in_files <- list.files(in_dir, '*', full.names = T)

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

for (i in seq_along(in_files)) {
  in_file <- in_files[i]
  out_file <- paste0(
    file_path_sans_ext(in_file),
    '_modkmerdist_k8_hclust_singlelinkage.tre'
  )

  alignment <- read.dna(in_file, format = 'fasta')
  d <- kdist_min(alignment, k=8)
  h <- hclust(d, method = 'single')
  p <- ape::as.phylo(h)
  ape::write.tree(p, file=out_file)
}
