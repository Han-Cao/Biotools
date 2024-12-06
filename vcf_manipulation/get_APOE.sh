#!/bin/bash
genome=hg19
prefix=""

vcf=$1
genome=$2
prefix=$3

if [[ $genome == "hg19" ]]; then
	bcftools convert --hapsample APOE -r ${prefix}19:45411941,${prefix}19:45412079 --vcf-ids $vcf
elif [[ $genome == "hg38" ]]; then
	bcftools convert --hapsample APOE -r ${prefix}19:44908684,${prefix}19:44908822 --vcf-ids $vcf
fi

bgzip -df APOE.hap.gz

echo '
convert_code <- function(data){
  ref <- data[4]
  alt <- data[5]
  data[data==0] <- ref
  data[data==1] <- alt
  return(data)
}

CountAllele <- function(target, AlleleTable){
  counts <- ifelse(AlleleTable[,1]==target,1,0) + ifelse(AlleleTable[,2]==target,1,0)
  return(counts)
}

hap=read.table("APOE.hap")
sample=read.table("APOE.samples", header=TRUE)[-1,-3]
hap=apply(hap, 1,convert_code)
hap1_index <- seq(6, nrow(hap)-1, 2)
hap2_index <- seq(7, nrow(hap), 2)
sample$rs429358_1=hap[hap1_index,1]
sample$rs429358_2=hap[hap2_index,1]
sample$rs7412_1=hap[hap1_index,2]
sample$rs7412_2=hap[hap2_index,2]
sample$hap1=paste0(sample$rs429358_1, sample$rs7412_1)
sample$hap2=paste0(sample$rs429358_2, sample$rs7412_2)

sample$E2=CountAllele("TT",sample[,c("hap1","hap2")])
sample$E3=CountAllele("TC",sample[,c("hap1","hap2")])
sample$E4=CountAllele("CC",sample[,c("hap1","hap2")])
write.table(sample,"APOE.txt",row.names=FALSE, quote=FALSE, sep="\t")
' >tmp_get_APOE.R
Rscript --vanilla tmp_get_APOE.R

rm tmp_get_APOE.R APOE.hap APOE.samples
