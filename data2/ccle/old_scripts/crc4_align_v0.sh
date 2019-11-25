#!/bin/bash
# crc_align_hgr1.fa
# rDNA alignment pipeline
# for CRC Hep data on GSC to hgr1
# 170614 -- 2156 build
# xhost10 

# Control Panel -------------------------------

# Project Dir
  BASE='/home/ababaian/projects/crc4'
  cd $BASE

# Sequencing Data
  CRC_DIR='/home/ababaian/projects/crc4'
  LIB_LIST='crc4_data.txt' # list of crc data fastq files
  
# CPU
  THREADS='2'
  
# Initialize start-up sequence ----------------
# Make working directory
  mkdir -p align

#Resources
  #aws s3 cp s3://crownproject/resources/hgr1.fa ./
  samtools faidx hgr1.fa
  bowtie2-build hgr1.fa hgr1


# ---------------------------------------------
# SCRIPT LOOP ---------------------------------
# ---------------------------------------------
# For each line in input LIB_LIST; run the pipeline

cat $LIB_LIST | while read LINE
do
    #Initialize Run
    echo "Start Iteration:"
    echo "  $LINE"
    echo ''
    
    LIBRARY=$(echo $LINE | cut -f1 -d' ' -) # Library Name
    RGSM=$(echo $LINE | cut -f2 -d' ' -)    # Sample / Patient Identifer
    RGID=$(echo $LINE | cut -f3 -d' ' -)    # Read Group ID
    RGLB=$(echo $LINE | cut -f3 -d' ' -)    # Library Name. Accession Number
    RGPL='ILLUMINA'                   # Sequencing Platform.
    RGPO=$(echo $LINE | cut -f4 -d' '  -)    # Patient Population

    FASTQ1=$(echo $LINE | cut -f5 -d' ' -)  # Filename Read 1
    FASTQ2=$(echo $LINE | cut -f6 -d' ' -)  # Filename Read 2
    
    FQ1="$CRC_DIR/$FASTQ1"            # Fastq1 Filepath
    FQ2="$CRC_DIR/$FASTQ2"            # Fastq2 Filepath
    
    # Extract Sequencing Run Info
    RGPU=$RGID
    
    # Bowtie2: align to genome
    bowtie2 --very-sensitive-local -p $THREADS --rg-id $RGID \
      --rg LB:$RGLB --rg SM:$RGSM \
      --rg PL:$RGPL --rg PU:$RGPU \
      -x hgr1 -1 $FQ1 -2 $FQ2 |\
      samtools view -bS - > aligned_unsorted.bam
      
    # Calcualte library flagstats
    samtools flagstat aligned_unsorted.bam > aligned_unsorted.flagstat

    # Read Subset ------------------------------
    # Extract mapped reads, and their unmapped pairs

      # Extract Header
      samtools view -H aligned_unsorted.bam > align.header.tmp

      # Unmapped reads with mapped pairs
      # Extract Mapped Reads
      # and their unmapped pairs
      samtools view -b -F 4 aligned_unsorted.bam > align.F4.bam #mapped
      samtools view -b -f 4 -F 8 aligned_unsorted.bam > align.f4F8.bam #unmapped pairs

      # Extract just the 45S unit
      #aws s3 cp s3://crownproject/resources/rDNA_45s.bed ./
      #samtools view -b -L rDNA_45s.bed align.F4.bam > align.F4.45s.bam

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

    # Rename the total Bam Files
      rm aligned_unsorted.bam
      rm aligned_unsorted.bam.bai $LIBRARY.bam.bai
      rm aligned_unsorted.flagstat $LIBRARY.flagstat

    # Rename the hgr Bam files
      mv align.hgr1.bam $LIBRARY.hgr1.bam
      mv align.hgr1.bam.bai $LIBRARY.hgr1.bam.bai
      mv align.hgr1.flagstat $LIBRARY.hgr1.flagstat

done

# Primary VCF ----------------------------
# N/A

# Script complete
