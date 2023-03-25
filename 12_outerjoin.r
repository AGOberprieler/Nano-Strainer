# Author: Tankred Ott
# Usage: Replace [path/to/statistics_sample_x] in lines 6, 7 and 10 with your own input and output paths and files; all input files must contain headers in their first line
# In line 9, "by='locus'" states the name/header of the column used as basis for joining (i.e., the column containing the locus names in the two files). Change to your own column header if necessary. Entries in the "locus" column of the files to be joined must not contain sample IDs!
# NOTE: Only two files can be joined at a time, i.e., the result file has to be used as input for another round of joining the next sample's statistics file, until all samples have been joined

a <- read.csv('/[path/to/statistics_sample_1].csv', sep = ';')
b <- read.csv('/[path/to/statistics_sample_2].csv', sep = ';')

r <- merge(a, b, by='locus', all=TRUE)
write.table(r, '/[path/to/statistics_sample_1and2].csv', sep = ';', row.names=F)
