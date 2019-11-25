#!/bin/bash
# mrpko_align_hgr1.fa
# rDNA alignment pipeline
# for MRP KO data on GSC to hgr1
# 170614 -- 2156 build
# xhost10 

# Control Panel -------------------------------

# Project Dir
  BASE='/home/ababaian/projects/MRPko'
  cd $BASE

# Sequencing Data
  CRC_DIR='/home/ababaian/projects/MRPko'
  LIB_LIST='mrpko_data.txt' # list of crc data fastq files
  
# CPU
  THREADS='2'
  
# Initialize start-up sequence ----------------
# Make working directory
  mkdir -p align

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
      

    # Sort the library
    samtools sort aligned_unsorted.bam $LIBRARY.hgr1.bam

    # Index the library
    samtools index $LIBRARY.hgr1.bam

    # Calculate flagstats for library
    samtools flagstat $LIBRARY.hgr1.bam > $LIBRARY.flagstat
    
    
    # Rename the total Bam Files
      rm aligned_unsorted.bam


done

# Primary VCF ----------------------------
# N/A

# Script complete
