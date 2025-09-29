## performing TAD analysis 
suppressPackageStartupMessages({
  library(SpectralTAD)
})
DIR ="HIC/"
setwd(DIR)

cool_mat = read.table(paste0("GSE240410_data/MB288_HL.10kb.txt"))
#Convert to sparse 3-column matrix using cooler2sparse from HiCcompare
sparse_mats = HiCcompare::cooler2sparse(cool_mat)
spec_hier = SpectralTAD(sparse_mats[["chr11"]], chr = "chr11", 
                        resolution = 25000, qual_filter = TRUE, levels = 3)
head(spec_hier$Level_1)
df <- data.frame(spec_hier$Level_1)
head(df)
write.table(df, file = paste0("beds/",sam,"_chr11_TAD_lv1_10kb.bed"),
            row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
q()
df <- data.frame(spec_hier$Level_2)
df <- df[!is.na(df$Sil_Score),]
head(df)
write.table(df, file = "MB288_chr11_TAD_lv2_25kb.bed",
            row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
df <- data.frame(spec_hier$Level_3)
df <- df[!is.na(df$Sil_Score),]
head(df)
write.table(df, file = "MB288_chr11_TAD_lv3_25kb.bed",
            row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

q()