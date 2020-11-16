library(rtracklayer)
library(dplyr)

options(stringsAsFactors = FALSE)

pos_in_region <- function(pos, start, end){
  greater <- pos >= start
  smaller <- pos <= end
  row_n <- which(smaller & greater)
  return(row_n)
}

#read input files
target_list <- read.csv("Input/20201015 gene list_314.csv")
target_trackers <- c()

#set input for UCSC track
position <- "chr19:29,925,114-29,960,718"
hg <- "hg19"

position <- gsub(",", "", position)
chr <- gsub("(\\w+)?:(\\d+)-(\\d+)", "\\1", position)
pos_start <- as.numeric(gsub("(\\w+)?:(\\d+)-(\\d+)", "\\2", position))
pos_end <- as.numeric(gsub("(\\w+)?:(\\d+)-(\\d+)", "\\3", position))

#extract TFBS
ucsc <- browserSession()
genome(ucsc) <- hg
query <- ucscTableQuery(ucsc, track="knownGene", GRangesForUCSCGenome(hg, chr, IRanges(pos_start, pos_end)))
tf <- data.frame(track(query))

target <- read.table("D:/Han/Biotools/Input/SNP_info.input", header = TRUE)
target$flag <- FALSE
for (i in 1:nrow(target)) {
  for (j in 1:nrow(tf)) {
    if(target[i, 'POS'] >= tf[j, 'start'] & target[i, 'POS'] <= tf[j, 'end']){
      target[i, 'flag'] <- TRUE 
      break()
    }
  }
}

#extract CREs
CREs <- read.table("D:/Han/Annotation/mm10-ccREs.bed", col.names = c("chr", "start", "end",  "cCRE_name1", "cCRE_name2", "Type"))
target <- read.csv("D:/Han/Biotools/Input/20201015 gene list_314.csv", header = TRUE)
target$cCRE <- "No"
for (i in 1:nrow(target)) {
  chr_CREs <- filter(CREs, chr==target[i,'chr'])
  row_n <- pos_in_region(target[i, 'pos'], chr_CREs$start, chr_CREs$end)
  if (length(row_n) > 0) {
    print(row_n)
    target[i, "cCRE"] <- chr_CREs[row_n, "Type"]
  }
}


#get SNP position
snp_mart = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
test <- getBM(attributes = c('refsnp_id','allele','consequence_type_tv'), 
              filters = c('snp_filter'), values = target$SNP, mart = snp_mart)

grch37_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", 
                 path="/biomart/martservice", dataset="hsapiens_snp")
snp_pos <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end'),
                 filters = c('snp_filter'), values = target, mart = grch37_snp)
