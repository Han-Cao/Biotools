library(rtracklayer)
library(biomaRt)
library(dplyr)

#set input for UCSC track
position <- "chr9:6,215,786-6,257,983"
hg <- "hg19"

position <- gsub(",", "", position)
chr <- gsub("(\\w+)?:(\\d+)-(\\d+)", "\\1", position)
pos_start <- as.numeric(gsub("(\\w+)?:(\\d+)-(\\d+)", "\\2", position))
pos_end <- as.numeric(gsub("(\\w+)?:(\\d+)-(\\d+)", "\\3", position))

#extract TFBS
ucsc <- browserSession()
genome(ucsc) <- hg
query <- ucscTableQuery(ucsc, track="wgEncodeRegTfbsClusteredV3", GRangesForUCSCGenome(hg, chr, IRanges(pos_start, pos_end)))
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
CREs <- read.table("Input/IL33_CRE.bed", col.names = c("chr", "start", "end",  "CRE", "z"))
target <- read.table("D:/Han/Biotools/Input/SNP_info.input", header = TRUE)
target$flag <- FALSE
for (i in 1:nrow(target)) {
  for (j in 1:nrow(CREs)) {
    if(target[i, 'POS'] >= CREs[j, 'start'] & target[i, 'POS'] <= CREs[j, 'end']){
      target[i, 'flag'] <- TRUE 
      break()
    }
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
