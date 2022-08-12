#!/usr/bin/env python3

# Extract INFO fields from vcf
# Created: 12/8/2022
# Author: Han Cao

import pysam
import argparse

# parse arguments
parser = argparse.ArgumentParser(description="Extract, rename, and merge INFO fields from vcf")
parser.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                    help="input vcf")
parser.add_argument("-o", "--outvcf", metavar="vcf", type=str, required=True,
                    help="input vcf")
parser.add_argument("--info", metavar="TAG", type=str, required=True,
                    help="Comma seperated INFO tags to extract, rename by NEW=OLD. Can merge multiple INFO into 1 new INFO, priority from high to low: AF,NEW=High_priority_old,NEW=Low_priority_old")
parser.add_argument("--keep-old", action="store_true", default=False,
                    help="Keep raw INFO when rename, useful when merging (default: False)")
parser.add_argument("--all", action="store_true", default=False,
                    help="Whether keep all INFO fields. If specified, only rename work in --info (default: False)")
args = parser.parse_args()

# parse args.info pattern string
class INFO_parser:
    def __init__(self, pattern_str) -> None:
        self.pattern_str = pattern_str
        self.pattern_dict = self.parse_pattern_str(pattern_str)
        self.new_tags = self.pattern_dict.keys()
        self.old_tags = [tag for old_tag_list in self.pattern_dict.values() for tag in old_tag_list]

    def parse_pattern_str(self, pattern_str) -> dict:
        """ Parse pattern string """
        pattern_dict = {}
        for x in pattern_str.split(","):
            if "=" in x:
                new_tag, old_tag = x.split("=")
            else:
                new_tag = x
                old_tag = x

            # add or rename
            if new_tag not in pattern_dict:
                pattern_dict[new_tag] = [old_tag]
            # merge multiple old tags into one new tag
            else:
                pattern_dict[new_tag].append(old_tag)

        return pattern_dict
    
    def map_tag(self, new_tag) -> list:
        """ Map new tag to old tag list """
        return self.pattern_dict[new_tag]

# read input
invcf = pysam.VariantFile(args.invcf, "r")
invcf_parser = INFO_parser(args.info)

# create new vcf header
new_header = invcf.header
for new_tag in invcf_parser.new_tags:
    if new_tag in new_header.info:
        continue
    else:
        old_tag = invcf_parser.map_tag(new_tag)[0]
        old_header = new_header.info[old_tag]
        new_header.info.add(new_tag, old_header.number, old_header.type, old_header.description)

# write new vcf
outvcf = pysam.VariantFile(args.outvcf, "w", header=new_header)
for variant in invcf:
    new_variant = variant.copy()
    
    # clear all old INFO
    if not args.all:
        new_variant.info.clear()

        # keep old tags if specified
        if args.keep_old:
            for old_tag in invcf_parser.old_tags:
                if old_tag in variant.info:
                    new_variant.info[old_tag] = variant.info[old_tag]

    # add new tags
    for new_tag in invcf_parser.new_tags:
        old_tag_list = invcf_parser.map_tag(new_tag)

        # set new_tag as the first existing old_tag
        for old_tag in old_tag_list:
            if old_tag in variant.info:
                new_variant.info[new_tag] = variant.info[old_tag]
                break
    
    outvcf.write(new_variant)

outvcf.close()
