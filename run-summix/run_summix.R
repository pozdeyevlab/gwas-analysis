library(optparse)
library(data.table)
library(glue)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Summix")

library(Summix)

# Define command-line arguments
option_list <- list(
  make_option(c("--input"), type="character", help="Aligned Result"),
  make_option(c("--reference"), type="character", help="Gnomad Reference"),
  make_option(c("--output"), type="character", help="Output File")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

print(opt$input)
print(opt$reference)
print(opt$output)

# reading summary statistics
# subsetting for chr21
header = glue('AF_afr\tAF_amr\tAF_eas\tAF_nfe\tAF_mid\tAF_sas')
write.table(header, file=opt$output, quote=FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
ss <- fread(opt$input)
ss$MarkerName <- ss$STUDY_ID
ss <- ss[grepl("^21", ss$MarkerName), ]
ss$MarkerName <- paste0("chr", ss$MarkerName)
  
# loading gnomad 4.0.0 reference for chr21
gnomad <- fread(opt$reference)

gnomad <- gnomad[ , c("# [1]CHROM", "[2]POS", "[3]REF", "[4]ALT", "[11]AF_afr", "[17]AF_amr", "[23]AF_eas","[32]AF_nfe", "[29]AF_mid", "[39]AF_sas")]

colnames(gnomad) <-  c("CHROM", "POS", "REF", "ALT", "AF_afr", "AF_amr", "AF_eas","AF_nfe", "AF_mid", "AF_sas")

gnomad$MarkerName <- paste0(gnomad$CHROM, ":", gnomad$POS, ":", gnomad$REF, ":", gnomad$ALT)

merged <- merge(gnomad, ss, by = "MarkerName")

# filtering for variantsd that matched
merged <- merged[!is.na(merged$Aligned_AF), ]

# downsampling to random 10000 SNPs
for (i in 1:5) {
  merged_ds <- merged[sample.int(nrow(merged), 10000, replace=TRUE)]

  # inverting freq when needed
  merged_ds$true_freq <- merged_ds$Aligned_AF
  
  merged_ds <- merged_ds[!is.na(merged_ds$true_freq), ] # all frequencies are calculated
  ##write.table(header, file=glue('results_3/{ancestry}_summix_{phenotype}.tsv'), append=TRUE)
  c <-summix(data = merged_ds,
         reference=c("AF_afr",
                     "AF_amr",
                     "AF_eas",
                     "AF_nfe",
                     "AF_mid",
                     "AF_sas"),
         observed="true_freq") #,pi.start = c(.167, .167, .167, .167, .167, .167))
  line = glue('{c$AF_afr}\t{c$AF_amr}\t{c$AF_eas}\t{c$AF_nfe}\t{c$AF_mid}\t{c$AF_sas}')
  write.table(line, file=opt$output, append=TRUE, quote=FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
}
