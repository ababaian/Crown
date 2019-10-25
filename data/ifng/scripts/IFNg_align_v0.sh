#!/bin/bash
# IFNg_align_v0.sh
# rDNA alignment pipeline
# for IFNg RNAseq data on AWS

# Control Panel -------------------------------

# Project Dir
  mkdir -p /home/ubuntu/align
  BASE='/home/ubuntu/align'
  cd $BASE

# Sequencing Data
  CRC_DIR='/home/ubuntu/align'
  LIB_LIST='/home/ubuntu/IFNg_data_1.txt' # list of IFNg data for the run
  
# CPU
  THREADS='3'
  
# Initialize start-up sequence ----------------
# Make working directory
mkdir -p bam
mkdir -p rawbam
mkdir -p fq
mkdir -p flagstat

#Resources
  aws s3 cp s3://crownproject/resources/hgr1.fa ./
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
    RGPO=$(echo $LINE | cut -f4 -d' ' -)    # Patient Population

    FASTQ1=$(echo $LINE | cut -f5 -d' ' -)  # Filename Read 1   
    FQ1="$CRC_DIR/$FASTQ1"            # Fastq1 Filepath

    echo "Lib: $LIBRARY"
    echo "RGSM: $RGSM"
    echo "RGID: $RGID"
    echo "RGLB: $RGLB"
    echo "RGPL: $RGPL"
    echo "RGPO: $RGPO"
    echo "FQ: $FQ1"
    echo ''
    
    # Extract Sequencing Run Info
    RGPU=$(gzip -dc $FQ1 | head -n1 - | cut -f1 -d'.' | cut -f2 -d'@')
    
    # Bowtie2: align to genome
    bowtie2 --very-sensitive-local -p $THREADS --rg-id $RGID \
      --rg LB:$RGLB --rg SM:$RGSM \
      --rg PL:$RGPL --rg PU:$RGPU \
      -x hgr1 -U $FQ1 |\
      samtools view -bS - > aligned_unsorted.bam
      
    # Calcualte library flagstats
    samtools flagstat aligned_unsorted.bam > aligned_unsorted.flagstat

    # Rename the total Bam Files
      mv aligned_unsorted.bam $LIBRARY.bam
      mv aligned_unsorted.bam.bai $LIBRARY.bam.bai
      mv aligned_unsorted.flagstat $LIBRARY.flagstat


  # Read Subset ------------------------------
    # Extract mapped reads, and their unmapped pairs

      # Extract Header
      samtools view -H $LIBRARY > align.header.tmp

      # Extract Mapped Reads (single-end)
      samtools view -b -F 4 $LIBRARY > align.F4.bam #mapped

        samtools sort align.F4.bam -o align.hgr1.bam
        samtools index align.hgr1.bam
        samtools flagstat align.hgr1.bam > align.hgr1.flagstat

      # Clean up 
      rm *tmp* align.F4.bam align.f4F8.bam

    # Rename/move the hgr1 Bam files
      mv align.hgr1.bam bam/$LIBRARY.hgr1.bam
      mv align.hgr1.bam.bai bam/$LIBRARY.hgr1.bam.bai
      mv align.hgr1.flagstat bam/$LIBRARY.hgr1.flagstat

   # Rename/move the raw total bam file
      mv $LIBRARY.bam rawbam/
      mv $LIBRARY.bam.bai rawbam/
      mv $LIBRARY.flagstat flagstat/

done

# Script complete
