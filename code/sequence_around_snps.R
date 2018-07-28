library(GenomicRanges) # essential package
library(BSgenome.Mmusculus.UCSC.mm10) # this will specify the reference genome that we are interested in
#library(BSgenome.Hsapiens.UCSC.hg19) # if you want to do human, consider using this one

# Make a data frame of the SNPs of interest... convention is chr1, chr2... chrX
snps_df <- data.frame(
 chr = c("chrX", "chrX", "chrX"),
 position = c(13291199, 13291286, 94209650),
 ID = c("rs33663591", "rs33663592", "rs29071854")
)

# Determine the number of basepairs to pad on each side
pad <- 20

# Make a GenomicRanges object of the regions of interest
snps_df$start <- snps_df$position - pad
snps_df$end <- snps_df$position + pad
snps_gr <- makeGRangesFromDataFrame(snps_df, keep.extra.columns = TRUE)

# Extract sequences from a reference genome
sequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10, snps_gr)

# Put it all together
final_df <- data.frame(snps_gr, sequence)
final_df$reference_allele <- substr(as.character(final_df$sequence), pad+1, pad + 1) # will extract the 1 bp for the SNP of interest

# Write to a TSV file 
write.table(final_df, file = "SNPs-mm10-sequences-output.tsv", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
