library(GenomicRanges) # essential package
library(BSgenome.Mmusculus.UCSC.mm10) # this will specify the reference genome that we are interested in
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # package for mm10 reference genes
#library(BSgenome.Hsapiens.UCSC.hg19) # if you want to do human, consider using this one

# Make a data frame of the SNPs of interest... convention is chr1, chr2... chrX
snps_df <- data.frame(
  chr = c("chrX", "chrX", "chrX"),
  position = c(13291199, 94209650, 94209648),
  ID = c("rs33663591", "PROBLEM1", "PROBLEM2")
)

# Determine the number of basepairs to pad on each side
pad <- 20

# Make a GenomicRanges object of the regions of interest
snps_df$start <- snps_df$position - pad
snps_df$end <- snps_df$position + pad
snps_gr <- makeGRangesFromDataFrame(snps_df, keep.extra.columns = TRUE)

# Extract sequences from a reference genome-- define the naive sequence as the genomic DNA
dna_sequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10, snps_gr)

# Put it all together
final_df <- data.frame(snps_gr, dna_sequence)
final_df$reference_allele <- substr(as.character(final_df$dna_sequence), pad+1, pad + 1) # will extract the 1 bp for the SNP of interest

# Overlap with exons to see if the nominated sequences are fully contained
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
exons_snps_overlap <- exonsByOverlaps(txdb, snps_gr, minoverlap = width(snps_gr)[1]) 

# See which SNPs were contained in one exon
overlaps <- findOverlaps(snps_gr, exons_snps_overlap)
final_df$in_one_exon <- 1:dim(final_df)[1] %in% unique(queryHits(overlaps))

# Filter for problem snps -- do these again for RNA sequence
snps_gr_problem <- snps_gr[!final_df$in_one_exon]

# Have to loop over these manually...
rna_sequence_combined_exons <- sapply(1:length(snps_gr_problem), function(idx){
  gr_one <- snps_gr_problem[idx]
  # Grab gene associated with variant
  ov1 <- transcriptsByOverlaps(txdb, gr_one)
  tx_id <- head(mcols(ov1)[,"tx_name"], 1)
  
  # Get exons associated with gene and then the snp
  exons_gr <- exonsBy(txdb, by = "tx", use.names = TRUE)[[tx_id]]
  ov_exon <- findOverlaps(exons_gr, gr_one)
  
  # Pull exon number; grab this exon
  exon_snp <- queryHits(ov_exon)
  
  # Determine the overlapping parts of the query sequence
  in_exon <- pintersect(gr_one, exons_gr[exon_snp])
  missing_from_exon <- psetdiff(gr_one, exons_gr[exon_snp])
  missing_bp_n <- width(missing_from_exon)
  
  # If what is missing is bigger than what we have, then we need to look one exon further and vice-versa
  idx_shift <- ifelse(start(missing_from_exon) > start(in_exon), 1, -1)
  new_exon <- exons_gr[idx_shift + exon_snp]
  
  # get missing bases from the left side of the exon if it comes after the sequence
  if(idx_shift == 1) {
    gr_new <- GRanges(seqnames = c(seqnames(in_exon), seqnames(new_exon)),
                      ranges = c(
                        IRanges(start = start(in_exon),
                                end = end(in_exon)),
                        IRanges(start = start(new_exon),
                                end = start(new_exon) + missing_bp_n - 1)
                      ))
    totalSeq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr_new)
    
  } else{
    # Grab from the right side of the preceding exon
    gr_new <- GRanges(seqnames = c(seqnames(new_exon), seqnames(in_exon)),
                      ranges = c(IRanges(start = end(new_exon) - missing_bp_n + 1,
                                         end = end(new_exon)),
                                 IRanges(start = start(in_exon),
                                         end = end(in_exon))
                      ))
    totalSeq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr_new)
  }
  
  # Return new sequence
  paste0(totalSeq, collapse = "")
})

# Copy the DNA to RNA sequence and fix problems
final_df$rna_sequence <- final_df$dna_sequence
final_df$rna_sequence[!final_df$in_one_exon] <- rna_sequence_combined_exons

# Write to a TSV file 
write.table(final_df, file = "SNPs-mm10-sequences-RNA-DNA.tsv", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
