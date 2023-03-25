# Author: Tankred Ott
# Usage: Replace "/path/to/5_input_example.xlsx" at the beginning with your own input path and file, preserving the single quotation marks; also replace "sheet_name" with the name of the sheet that contains the data in your own Excel input file

library(readxl)

in_file <- '/path/to/5_input_example.xlsx'

x <- read_xlsx(in_file, sheet = 'sheet_name')

plot(seq(0.01, 0.99, 0.01), sapply(
  seq(0.01, 0.99, 0.01),
  function(q) {
    nm_q <- quantile(x$norm.mean, probs = 1 - q)
    nm_idx <- x$norm.mean > nm_q
    
    kbs_q <- quantile(x$KBS_mean, probs = q)
    kbs_idx <- x$KBS_mean < kbs_q
    
    kbs_sd_q <- quantile(x$KBS_sd, probs = q)
    kbs_sd_idx <- x$KBS_sd < kbs_sd_q
    
    ent_q <- quantile(x$entropy, probs = q)
    ent_idx <- x$entropy < ent_q
    
    return(
      sum(nm_idx & kbs_idx & kbs_sd_idx & ent_idx)
    )
  }
),
    ylab = 'number of loci',
    xlab = 'proportion best values',
)

q <- 0.5
nm_q <- quantile(x$norm.mean, probs = 1 - q)
nm_idx <- x$norm.mean > nm_q

kbs_q <- quantile(x$KBS_mean, probs = q)
kbs_idx <- x$KBS_mean < kbs_q

kbs_sd_q <- quantile(x$KBS_sd, probs = q)
kbs_sd_idx <- x$KBS_sd < kbs_sd_q

ent_q <- quantile(x$entropy, probs = q)
ent_idx <- x$entropy < ent_q

sum(nm_idx & kbs_idx & kbs_sd_idx & ent_idx)

qs <- seq(0.01, 0.99, 0.01)
r_2 <- vector(mode = 'list', length = length(qs))
for (i in seq_along(qs)) {
  q <- qs[i]
  nm_q <- quantile(x$norm.mean, probs = 1 - q)
  nm_idx <- x$norm.mean > nm_q
  
  kbs_q <- quantile(x$KBS_mean, probs = q)
  kbs_idx <- x$KBS_mean < kbs_q
  
  kbs_sd_q <- quantile(x$KBS_sd, probs = q)
  kbs_sd_idx <- x$KBS_sd < kbs_sd_q
  
  ent_q <- quantile(x$entropy, probs = q)
  ent_idx <- x$entropy < ent_q
  
  idx <- nm_idx & kbs_idx & kbs_sd_idx & ent_idx
  
  if (i == 1) {
    r_2[[i]] <- x[idx,]$locus
  } else {
    new_r <- x[idx,]$locus
    new_r <- new_r[!(new_r %in% r_2[[i-1]])]
    r_2[[i]] <- c(r_2[[i-1]], new_r)
  }
}

nrValues <- max(sapply(r_2, length))
m <- data.frame(matrix('',nrow = nrValues, ncol=length(r_2)))
colnames(m) <- qs
for (i in 1:length(r_2)) {
  v <- c(r_2[[i]], rep('', nrValues - length(r_2[[i]])))
  m[,i] <- v
}
m

write.csv(m, 'best_loci.csv')
