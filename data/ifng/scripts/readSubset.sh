#!/bin/bash
#
# readSubset.sh
#

for LIBRARY in $(ls *.bam | cut -f1 -d'.' -)
do

  # Read Subset ------------------------------
    # Extract mapped reads, and their unmapped pairs

      # Extract Header
      samtools view -H $LIBRARY.bam > align.header.tmp

      # Extract Mapped Reads (single-end)
      samtools view -b -F 4 $LIBRARY.bam > align.F4.bam #mapped

        samtools sort align.F4.bam -o align.hgr1.bam
        samtools index align.hgr1.bam
        samtools flagstat align.hgr1.bam > align.hgr1.flagstat

      # Clean up 
      rm *tmp* align.F4.bam align.f4F8.bam

    # Rename/move the hgr1 Bam files
      mv align.hgr1.bam ../bam/$LIBRARY.hgr1.bam
      mv align.hgr1.bam.bai ../bam/$LIBRARY.hgr1.bam.bai
      mv align.hgr1.flagstat ../bam/$LIBRARY.hgr1.flagstat

done
