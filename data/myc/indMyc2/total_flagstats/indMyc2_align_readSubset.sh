#!/bin/bash
#
# readSubset.sh
# moved output bam files from 
# siMyc_align_v0 to ~/rawbam

mkdir -p bam

for BAM in $(ls *.bam)
do
   # Read Subset ------------------------------
    # Extract mapped reads, and their unmapped pairs

      # Extract Header
      samtools view -H $BAM > align.header.tmp

      # Unmapped reads with mapped pairs
      # Extract Mapped Reads
      # and their unmapped pairs
      samtools view -b -F 4 $BAM > align.F4.bam #mapped

        samtools sort align.F4.bam -o align.hgr1.bam
        samtools index align.hgr1.bam
        samtools flagstat align.hgr1.bam > align.hgr1.flagstat

      # Clean up 
      rm *tmp* align.F4.bam align.f4F8.bam

    # Rename the hgr Bam files
      mv align.hgr1.bam bam/$BAM.hgr1.bam
      mv align.hgr1.bam.bai bam/$BAM.hgr1.bam.bai
      mv align.hgr1.flagstat bam/$BAM.hgr1.flagstat

done

