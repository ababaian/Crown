#!/bin/bash
# ADcalc.sh
# Allelic Depth Calculator
# for a range of coordinates

# Controls -----------------
REGION1='chr13:1003661-1005529'
REGION2='chr13:1007948-1013560'
# chr13:1003661-1005529 18S
# chr13:1007948-1013560 28S

OUTPUT1='18S_siMyc.gvcf'
OUTPUT2='28S_siMyc.gvcf'

DEPTH='100000'
BAMLIST='bam.list'

# Iterate through every bam file in directory
# look-up position and return VCF

## Region 1: 18S
    bcftools mpileup -f hgr1.fa \
  --max-depth $DEPTH --min-BQ 30 \
  -a FORMAT/DP,AD \
  -r $REGION1 \
  --ignore-RG \
  -b $BAMLIST |
  bcftools annotate -x INFO,FORMAT/PL - |
  bcftools view -O v -H - >> $OUTPUT1

## Region 2: 28S
    bcftools mpileup -f hgr1.fa \
  --max-depth $DEPTH --min-BQ 30 \
  -a FORMAT/DP,AD \
  -r $REGION2 \
  --ignore-RG \
  -b $BAMLIST |
  bcftools annotate -x INFO,FORMAT/PL - |
  bcftools view -O v -H - >> $OUTPUT2

