rm(list=ls())
library(Summix)
library(data.table)
# reading summary statistics
# subsetting for chr21
ss <- fread("C:/Users/Nikitos/Downloads/thyroid_cancer_1.tbl")
ss <- ss[grepl("^21", ss$MarkerName), ]
ss$MarkerName <- paste0("chr", ss$MarkerName)
write.csv(ss, "~/Downloads/tmp.csv", row.names = F)
# loading gnomad 4.0.0 reference for chr21
gnomad <- fread("C:/Users/Nikitos/Downloads/gnomad_ref_chr21.tsv")
gnomad <- gnomad[ , c("# [1]CHROM", "[2]POS", "[3]REF", "[4]ALT", "[11]AF_afr", "[17]AF_amr", "[23]AF_eas",
                      "[32]AF_nfe", "[29]AF_mid", "[39]AF_sas")]
colnames(gnomad) <-  c("CHROM", "POS", "REF", "ALT", "AF_afr", "AF_amr", "AF_eas",
                       "AF_nfe", "AF_mid", "AF_sas")
gnomad$MarkerName <- paste0(gnomad$CHROM, ":", gnomad$POS, ":", gnomad$REF, ":", gnomad$ALT)
write.csv(gnomad, "~/Downloads/gnomad_tmp.csv", row.names = F)
#ss <- fread("~/Downloads/tmp.csv")
#gnomad <- fread("~/Downloads/gnomad_tmp.csv")
merged <- merge(gnomad, ss, by = "MarkerName")
# filtering for variantsd that matched
merged <- merged[!is.na(merged$Freq1), ]
# downsampling to random 10000 SNPs
merged_ds <- merged[sample.int(nrow(merged), 10000)]
# inverting freq when needed
merged_ds$true_freq <- NA
for (i in 1:nrow(merged_ds)){
  if ((toupper(merged_ds$REF[i]) == toupper(merged_ds$Allele1[i])) & (toupper(merged_ds$ALT[i]) == toupper(merged_ds$Allele2[i])))
    merged_ds$true_freq[i] <- merged_ds$Freq1[i]
  if ((toupper(merged_ds$REF[i]) == toupper(merged_ds$Allele2[i])) & (toupper(merged_ds$ALT[i]) == toupper(merged_ds$Allele1[i])))
    merged_ds$true_freq[i] <- 1 - merged_ds$Freq1[i]
}
merged_ds <- merged_ds[!is.na(merged_ds$true_freq), ] # all frequencies are calculated
summix(data = merged_ds,
       reference=c("AF_afr",
                   "AF_amr",
                   "AF_eas",
                   "AF_nfe",
                   "AF_mid",
                   "AF_sas"),
       observed="true_freq",
#       pi.start = c(.2, .2, .2, .2, .2),
       goodness.of.fit=TRUE)
