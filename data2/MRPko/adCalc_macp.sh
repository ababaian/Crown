#!/bin/bash
# ADcalc.sh
# Allelic Depth Calculator
# for a range of coordinates

# Controls -----------------
REGION='chr13:1004900-1005000'
OUTPUT='macp_MRPko.gvcf'
DEPTH='100000'
BAMLIST='bam.list'

# Iterate through every bam file in directory
# look-up position and return VCF

bcftools mpileup -f hgr1.fa \
 --max-depth $DEPTH --min-BQ 30 \
 -a FORMAT/DP,AD \
 -r $REGION \
 --ignore-RG \
 -b $BAMLIST |
 bcftools annotate -x INFO,FORMAT/PL - |
 bcftools view -O v -H - >> $OUTPUT

