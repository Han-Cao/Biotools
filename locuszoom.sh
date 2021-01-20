#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Usage: bash locuszoom.sh -p plink_file -v vcf_file -o outpath -L region"
    printf "\n"

	echo "Required arguments:"
	echo "-p plink result (first line contain SNP BP and P column names)"
    echo "-v vcf for calculate LD (or -N for plot without ld, -V for ld file)"
    echo "-V pre-calculated ld file instead of vcf file "
    echo "-o outpath (absolute path only)"
    echo "-L plot interval (e.g. 1:1000-2000)"
    printf "\n"

    echo "Optional arguments:"
    echo "-r refsnp (default: auto detect by P value)"
    echo "-R additional refsnp"
    echo "-e exclude snps, seperate by comma"
    echo "-y max y axis (default=10)"
    echo "-n remove reference SNP name from plots"
    echo "-N no ld plot"
    echo "-l add significant line (-log10 value)"
    exit 0
fi

#set default value
y_max=10
no_SNP_name=FALSE
abline=NULL
rplot_drawMarkerNames=TRUE
rplot_ldColors="gray50,gray50,#f18c8d,#ea5354,#e41a1c,#e41a1c"
ref2=NULL

while getopts ":p:v:V:r:R:L:e:o:l:y:nN" opt; do
  case $opt in
    p) plink="$OPTARG"
    ;;
    v) vcf="$OPTARG"
    ;;
    V) ld="$OPTARG"
    ;;
    o) outpath="$OPTARG"
    ;;
    r) ref="$OPTARG"
    ;;
    R) ref2="$OPTARG"
    ;;
    L) region="$OPTARG"
    ;;
    e) exclude="$OPTARG"
    ;;
    y) y_max="$OPTARG"
    ;;
    n) rplot_drawMarkerNames="FALSE"
    ;;
    l) abline="$OPTARG"
    ;;
    N) no_LD="TRUE"
    ;;
  esac
done

[[ ! -d $outpath ]] && mkdir -p $outpath
chr=${region%:*}
region=${region#*:}
start=${region%-*}
end=${region#*-}

awk '
BEGIN {OFS="\t"}
NR==1 {
    for (i=1; i<=NF; i++) {
        ix[$i] = i
    }
    print "SNP", "BP", "P"
}
NR>1 {
    print $ix["SNP"], $ix["BP"], $ix["P"]
}' $plink >tmp.plink

if [[ -n $exclude ]] ; then
    echo -e "${exclude//,/\n}" >tmp.exclude
    cat tmp.exclude | while read snp
    do
        sed -i "/${snp}\s/d" tmp.plink
    done
fi

if [[ -z "$ref" ]]; then
    ref=$(awk -v start="$start" -v end="$end" '{if ($2>=start && $2<=end) print $0}' tmp.plink | \
          sort -g -k 3 | awk '{if ($3!="NA") print $1}' | head -n 1)
fi

if [[ ! -z "$vcf" ]]; then
    plink --vcf $vcf \
    --ld-snp $ref \
    --maf 0.05 \
    --r2 dprime \
    --ld-window-r2 0 \
    --ld-window 10000 \
    --out tmp
    printf '%s\t%s\t%s\t%s\n' "snp1" 'snp2' 'dprime' 'rsquare' > tmp.ld2
    awk 'FNR>1 {print $6 "\t" $3 "\t" $8 "\t" $7}' tmp.ld >>tmp.ld2
    ld=tmp.ld2
fi

name=$(basename $plink)
name=${name%%.*}

#set additional plot parameters
plot_par="ymax=${y_max} drawMarkerNames=${rplot_drawMarkerNames}"
if [[ ${abline} != "NULL" ]]; then
    plot_par="${plot_par} signifLine=${abline}"
fi
if [[ ${ref2} != "NULL" ]]; then
    plot_par="${plot_par} condLdLow='gray50' condPch='16,16,16,16,16,16,16,16,16' ldCuts = '0,.5,.7,.9,1'"
fi


if [[ $no_LD == "TRUE" ]]; then
    cat tmp.plink | \
    locuszoom --metal - \
    --markercol SNP \
    --pvalcol P \
    --refsnp $ref \
    --chr $chr \
    --start $start \
    --end $end \
    --no-ld \
    --build hg19 \
    --plotonly \
    -p ${outpath}${name} \
    ${plot_par}

elif [[ ${ref2} == "NULL" ]]; then
    locuszoom --metal tmp.plink \
    --markercol SNP \
    --pvalcol P \
    --refsnp $ref \
    --chr $chr \
    --start $start \
    --end $end \
    --ld $ld \
    --build hg19 \
    --plotonly \
    -p ${outpath}${name} \
    ${plot_par}
else
    locuszoom --metal tmp.plink \
    --markercol SNP \
    --pvalcol P \
    --refsnp $ref \
    --chr $chr \
    --start $start \
    --end $end \
    --ld $ld \
    --build hg19 \
    --plotonly \
    --add-refsnps $ref2 \
    -p ${outpath}${name} \
    ${plot_par} 
fi

rm tmp.plink 
[[ -n $exclude ]] && rm tmp.exclude
[[ -e tmp.ld ]] && rm tmp.ld
[[ -e tmp.ld2 ]] && rm tmp.ld2
[[ -e tmp.log ]] && rm tmp.log
[[ -e tmp.nosex ]] && rm tmp.nosex