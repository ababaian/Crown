#!/bin/bash
#
# readSubset.sh
# moved output bam files from 
# siMyc_align_v0 to ~/rawbam

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
      samtools view -b -f 4 -F 8 $BAM > align.f4F8.bam #unmapped pairs

      # What are the mapped readnames
      samtools view align.F4.bam | cut -f1 - > read.names.tmp

      # Extract mapped reads
      samtools view align.F4.bam | grep -Ff read.names.tmp - > align.F4.tmp.sam


      # Extract cases of read pairs mapped on edge of region of interest
      # -------|======= R O I ======| ----------
      # read:                  ====---====
      samtools view align.F4.bam | grep -Ff read.names.tmp - > align.F4.tmp.sam

      # Complete mapped reads list
      #cut -f1 align.F4.tmp.sam > read.names.45s.long.tmp

      # Extract unmapped reads with a mapped pair
      samtools view align.f4F8.bam | grep -Ff read.names.tmp - > align.f4F8.tmp.sam

      # Re-compile bam file
      cat align.header.tmp align.F4.tmp.sam align.f4F8.tmp.sam | samtools view -bS - > align.hgr1.tmp.bam

        samtools sort align.hgr1.tmp.bam -o align.hgr1.bam
        samtools index align.hgr1.bam
        samtools flagstat align.hgr1.bam > align.hgr1.flagstat

      # Clean up 
      rm *tmp* align.F4.bam align.f4F8.bam

    # Rename the hgr Bam files
      mv align.hgr1.bam $BAM.hgr1.bam
      mv align.hgr1.bam.bai $BAM.hgr1.bam.bai
      mv align.hgr1.flagstat $BAM.hgr1.flagstat

done

