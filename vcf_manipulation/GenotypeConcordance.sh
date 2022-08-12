#!/bin/bash

if [[ $# -eq 0 ]] ; then
	echo "Compare the genotype of 2 samples in one vcf file"
	echo "usage: bash GenotypeConcordance.sh vcf sample1 sample2"
    exit 0
fi

vcf=$1
sample1=$2
sample2=$3

bcftools query -f '[%GT ]\n' -s ${sample1},${sample2} -r 1 $vcf > tmp.txt

echo '
df = read.table("tmp.txt")
df$ref1 = df[,1] == "0/0"
df$ref2 = df[,2] == "0/0"
df$het1 = df[,1] %in% c("1/0", "0/1")
df$het2 = df[,2] %in% c("1/0", "0/1")
df$alt1 = df[,1] == "1/1"
df$alt2 = df[,2] == "1/1"
output = matrix(0,3,3, dimnames=list(c("ref/ref","ref/alt","alt/alt"),c("ref/ref","ref/alt","alt/alt")))

compare_geno <- function(df, col1, col2){
	sum(df[df[,col1], col2])
}

output[1,1] = compare_geno(df, "ref1" ,"ref2")
output[1,2] = compare_geno(df, "ref1" ,"het2")
output[1,3] = compare_geno(df, "ref1" ,"alt2")
output[2,1] = compare_geno(df, "het1" ,"ref2")
output[2,2] = compare_geno(df, "het1" ,"het2")
output[2,3] = compare_geno(df, "het1" ,"alt2")
output[3,1] = compare_geno(df, "alt1" ,"ref2")
output[3,2] = compare_geno(df, "alt1" ,"het2")
output[3,3] = compare_geno(df, "alt1" ,"alt2")
print(output)
' > tmp.R 
Rscript --vanilla tmp.R
rm tmp.R tmp.txt