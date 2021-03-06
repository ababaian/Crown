#!/bin/bash
# siMyc_align_hgr1.fa
# rDNA alignment pipeline
# for CRC data on local machine

# Control Panel -------------------------------

# Project Dir
  BASE='/home/artem/Desktop/siMyc'
  cd $BASE

# Sequencing Data
  CRC_DIR='/home/artem/Desktop/siMyc'
  LIB_LIST='siMyc_data.txt' # list of crc data fastq files
  
# CPU
  THREADS='1'
  
# Initialize start-up sequence ----------------
# Make working directory
  mkdir -p align

#Resources
  cp ~/Crown/resources/hgr1/hgr1.fa ./
  samtools faidx hgr1.fa
  bowtie2-build hgr1.fa hgr1

# GATK variant calling resources
  cp ~/Crown/resources/hgr1/hgr1.gatk.fa ./
  cp ~/Crown/resources/hgr1/hgr1.gatk.fa.fai ./
  cp ~/Crown/resources/hgr1/hgr1.gatk.dict ./


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


done

# Primary VCF ----------------------------
# N/A

# Script complete
