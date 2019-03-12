# Rscript for adjusted p-value calculation per transcript

arguments <- commandArgs(trailingOnly = TRUE)
filename <- arguments[1]
output <- arguments[2]

df <- read.table(filename, header = FALSE)

# fdr1 = p.adjust(df$V1, 'fdr')
fdr1 = p.adjust(df$V1, 'fdr', n = 2200000)

write.table(format(-log10(fdr1), digits = 6, scientific = F), file = output, quote = F, sep = "\t", col.names = F, row.names = F)