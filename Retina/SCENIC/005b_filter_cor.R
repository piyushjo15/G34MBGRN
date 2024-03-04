args = commandArgs(trailingOnly=TRUE)
DIR="Retina/SCENIC/"
df <- read.delim(paste0(DIR,args,"/cor",args,".tsv"))
keep <- is.na(df$rho)
df <- df[!keep,]
write.table(df, file =paste0(DIR,args,"/cor",args,"_cleaned.tsv"), 
            sep="\t",quote = FALSE, row.names = FALSE )
