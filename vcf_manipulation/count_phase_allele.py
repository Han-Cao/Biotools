import pysam
import argparse

# parse arguments
parser = argparse.ArgumentParser(description="Annotate vcf with phased GT count")
parser.add_argument("invcf", metavar="vcf", type=str, nargs=1,
                    help="input vcf")
args = parser.parse_args()

# read the input file
invcf = args.invcf[0]
outvcf = invcf.replace(".vcf.gz", ".phase_count.vcf.gz")
myvcf = pysam.VariantFile(invcf, "r")

# Add AN_phased to header
myvcf.header.info.add("AN_phased", "1", "Integer", "Phased allele number")
myvcf.header.info.add("AC_phased", "1", "Integer", "Phased allele count")
# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(outvcf, "w", header=myvcf.header)

for variant in myvcf:
    AN_phased = 0
    AC_phased = 0
    for sample in variant.samples:
        if variant.samples[sample].phased:
            AN_phased += 2
            if None not in variant.samples[sample]["GT"]:
                AC_phased += sum(variant.samples[sample]["GT"])
            
    variant.info["AN_phased"] = AN_phased
    variant.info["AC_phased"] = AC_phased
    vcf_out.write(variant)