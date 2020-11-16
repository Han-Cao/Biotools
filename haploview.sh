#!/bin/bash
if [[ $# -eq 0 ]] ; then
	echo "Generate Haploview input file from vcf"
	echo "Similar with 'plink --recode HV', but phased genotype is used"
	echo ""
	echo "Usage:"
	echo "-v input vcf file"
	echo "-o output name"
	echo "-l SNP list file"
	echo "-S sample list file"
    exit 0
fi

while getopts ":v:o:l:S:" opt; do
  case $opt in
    v) vcf="$OPTARG"
    ;;
    o) output="$OPTARG"
	;;
   	l) list="$OPTARG"
	;;
	S) sample="$OPTARG"
	;;
  esac
done

if [[ -n $list ]]; then
	if [[ -n $sample ]]; then
		bcftools convert --hapsample tmp --vcf-ids -i "ID=@./$list" -S $sample $vcf
	else
		bcftools convert --hapsample tmp --vcf-ids -i "ID=@./$list" $vcf
	fi
elif [[ -n $sample ]]; then
	bcftools convert --hapsample tmp --vcf-ids -S $sample $vcf
else
	bcftools convert --hapsample tmp --vcf-ids -S $vcf
fi
bgzip -d tmp.hap.gz

echo '
convert_code <- function(data){
  ref <- data[4]
  alt <- data[5]
  if (nchar(ref)>1 | nchar(alt)>1){
  	ref <- "A"
  	alt <- "T"
  }
  data[data==0] <- ref
  data[data==1] <- alt
  return(data)
}
dedup_snp <- function(hap){
  dup_list <- hap[duplicated(hap[,2]),2]
  if (length(dup_list) == 0){
    return(hap)
  }
  remove_snp <- c()
  for (dup_snp in dup_list) {
    dup_snp_index <- which(hap[,2]==dup_snp)
    keep_snp <- which.max(apply(hap[dup_snp_index, 6:ncol(hap)], 1, sum))
    dup_snp_index <- dup_snp_index[-keep_snp]
    remove_snp <- c(remove_snp, dup_snp_index)
  }
  return(hap[-remove_snp,])
}
args <- commandArgs(trailingOnly=TRUE)
hap_file <- args[1]
sample_file <- args[2]
out_ped <- args[3]
out_info <- args[4]
hap <- read.table(hap_file)
hap <- dedup_snp(hap)
sample <- read.table(sample_file)
sample <- sample[3:nrow(sample),]
write.table(hap[,2:3], out_info, col.names = FALSE, row.names= FALSE, quote = FALSE)
hap <- apply(hap, 1, convert_code)
hap_1 <- hap[seq(6,nrow(hap)-1,2), ]
hap_2 <- hap[seq(7,nrow(hap),2), ]
output <- cbind(sample, 0, 0, 0)
for (i in 1:ncol(hap_1)){
	output <- cbind(output, hap_1[,i], hap_2[,i])
}
write.table(output, out_ped, col.names = FALSE, row.names= FALSE, quote = FALSE)
' >tmp.R

Rscript tmp.R tmp.hap tmp.sample $output.ped $output.info
rm tmp.R tmp.hap tmp.sample