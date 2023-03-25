# Author: Tankred Ott
# Usage: Replace "/path/to/" in lines 10 and 83 with your own input and output paths
# In line 11 ("ind_names"), the individuals' ids  ('BC02', 'BC04', 'BC05', 'BC06') are to be complemented and modified as needed, however they must be identical to the start of the sequence names within the locus FASTA files in the input folder!
# In line 50, "k=8" can be modified in case a different k-mer length for distance computation is desired, but refer to the supplementary extended methods document of the paper for recommendations

library(ape)
library(kmer)


in_dir <- '/path/to/prechoiceloci_clusterconsensuses_locuswise'
ind_names <- c('BC02', 'BC04', 'BC05', 'BC06')


ind_names_sorted <- sort(ind_names)
in_files <- list.files(in_dir, '*', full.names = T)

locus_names <- sapply(strsplit(basename(in_files), '_'), function(x) x[1])

res_mat <- matrix(NA, nrow = length(in_files), ncol = length(ind_names) * 2 + 3)
row.names(res_mat) <- locus_names
colnames(res_mat) <- c(
  'd_mean',
  paste(ind_names, 'd_mean', sep = '_'),
  paste(ind_names, 'KBS', sep = '_'),
  'KBS_mean',
  'KBS_sd'
)

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

for (i in seq_along(in_files)) {
  in_file <- in_files[i]

  alignment <- read.dna(in_file, format = 'fasta')
  d <- as.matrix(kdist_min(alignment, k=8))

  d_mean <- mean(d[upper.tri(d, diag = F)])

  ind_prefixes <- sapply(strsplit(row.names(d), '_'), function(x) x[1])
  print(ind_prefixes)
  ord <- order(ind_prefixes, decreasing = F)

  ind_prefixes <- ind_prefixes[ord]


  ind_means <- sapply(unique(ind_prefixes), function(cur_prefix) {
    cur_prefix_idcs <- ind_prefixes == cur_prefix
    d_subset <- d[cur_prefix_idcs, cur_prefix_idcs, drop=F]

    if (ncol(d_subset) == 1) {return(0.0)}

    return(mean(d_subset[upper.tri(d_subset, diag = F)]))
  })

  kbs_scores <- (ind_means / d_mean - 1) * 100
  kbs_scores_mean <- mean(kbs_scores)
  kbs_scores_sd <- sd(kbs_scores)

  idcs <- match(names(kbs_scores), ind_names_sorted)

  res_mat[i, 1] <- d_mean
  res_mat[i, 1 + idcs] <- ind_means
  res_mat[i, 1 + length(ind_names) + idcs] <- kbs_scores
  res_mat[i, 1 + length(ind_names) * 2 + 1] <- kbs_scores_mean
  res_mat[i, 1 + length(ind_names) * 2 + 2] <- kbs_scores_sd  
}

write.csv(res_mat, file = '/path/to/KBS_prechoiceloci.csv', row.names=T)
