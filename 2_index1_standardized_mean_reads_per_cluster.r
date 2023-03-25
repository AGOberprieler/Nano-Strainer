# Author: Ulrich Lautenschlager
# Usage: Replace "2_input_example.csv" at the beginning with your own input, adapt the cbind command to the columns/column names of the CSV file and run: Rscript 2_index1_standardized_mean_reads_per_cluster.r

df <- read.csv("2_input_example.csv", header=T, sep=";")

attach(df)

zscores <- cbind(
    scale(readsperlocus_2 / clusterperlocus_2),
    scale(readsperlocus_4 / clusterperlocus_4),
    scale(readsperlocus_5 / clusterperlocus_5),
    scale(readsperlocus_6 / clusterperlocus_6)
)

df$mean_standardized_readsperclus <- apply(zscores, 1, mean, na.rm=T)
df$mean_standardized_readsperclus <- round(df$mean_standardized_readsperclus, 4)

write.csv(df, file="scores_standardized_mean_reads_per_cluster.csv", quote=F, row.names=F)
