#!/bin/bash
# hgr1_align_v4.ccle.sh
# rDNA alignment pipeline - SRA version
PIPE_VERSION='190531 build -- CCLE'
AMI_VERSION='crown-190601 - ami-0b375c9c58cb4a7a2'
# EC2: c4.2xlarge (8cpu / 15 gb)
# EC2: c4.xlarge  (4cpu / 8  gb)
# Storage: 200 Gb
#

# Input Requirements --------------------------

# $1 : Library name + Output name(unique)
# $2 : Seq-read type (wgs|rna)
# $3 : BioSample ID
# $4 : Library SRA Accession

# Control Panel -------------------------------
# Amazon AWS S3 Home URL
  S3URL='s3://crownproject/ccle'

# CPU
	THREADS='3'

# Terminate instances upon completion (for debuggin)
  TERMINATE='TRUE'
    
# Read Group Data
  LIBRARY=$1    # Library Name / File prefix
  TYPE=$2       # wgs OR rna data-type
	RGPO='ccle'   # Patient Population - CCLE
	RGSM=$3       # Sample / Patient Identifer
	RGID=$4       # Read Group ID. SRA Accession Number
  RGLB=$LIBRARY # Library Name. Accession Number
  RGPL='ILLUMINA' # Seq Platform
  RGPU=$5       # Read Group. Platform Unit (SRA Experiment)


  SRA=$4
  OUTPUT="$LIBRARY.$TYPE" # Output filename

  # Fastq files from fastq-dump
  # FQ1=<SRR>_1.fastq.gz
  # FQ2=<SRR>_2.fastq.gz

echo " -- hgr1 Alignment Pipeline -- "
echo " version: $PIPE_VERSION "
echo " ami:     $AMI_VERSION  "
echo " s3:      $S3URL  "
echo " library: $LIBRARY -- $TYPE"
echo " date:    $(date)"
echo ''

# Initialize wordir ---------------------------
echo 'Initializing ...'

# Make working directory
  mkdir -p align; cd align

# Copy hgrX genome and create bowtie2 index
  cp ~/resources/hgr1/* ./

# Load updated GDC token to instance
#aws s3 cp $S3URL/scripts/gdc.token gdc.token
#chmod 600 gdc.token

echo "Download SRA file: $SRA"
echo "  cmd: prefetch -X 200G --ascp-path <PATH> $SRA"

$HOME/bin/prefetch -X 100G --ascp-path \
  "$HOME/.aspera/connect/bin/ascp|$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh"\
  $SRA

# SRA INPUT ====================================================

#if [ -s $FQ1 ]
#then
    echo "SRA Input Pipe"
    echo ''
    echo "Starting hgr1 alignment"
    # Paired-End Extracted Reads Alignment

    # Bowtie2: align to genome
    ~/bin/bowtie2 --very-sensitive-local -p $THREADS \
      --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \
      --rg PL:$RGPL --rg PU:$RGPU \
      -x hgr1 --sra-acc $SRA | \
      ~/bin/samtools view -bS - > aligned_unsorted.bam
      
    echo "Alignment complete."
    echo "Calculate flagstats."

    # Calcualte library flagstats
      ~/bin/samtools flagstat aligned_unsorted.bam > aligned_unsorted.flagstat
      
    # Read Subset ------------------------------
    echo "Subset reads (retain mapped & their pairs, remove unmapped)."

    # Extract mapped reads, and their unmapped pairs

      # Extract Header
      ~/bin/samtools view -H aligned_unsorted.bam > align.header.tmp

      # Unmapped reads with mapped pairs
      # Extract Mapped Reads
      # and their unmapped pairs
      ~/bin/samtools view -b -F 4 aligned_unsorted.bam > align.F4.bam #mapped
      ~/bin/samtools view -b -f 4 -F 8 aligned_unsorted.bam > align.f4F8.bam #unmapped pairs

      echo "Recompiling mapped bam file."

      # Re-compile bam file
      ~/bin/samtools cat -h align.header.tmp -o align.hgr1.tmp.bam align.F4.bam align.f4F8.bam
        ~/bin/samtools sort -@ $THREADS -O BAM align.hgr1.tmp.bam > align.hgr1.bam
        ~/bin/samtools index align.hgr1.bam
        ~/bin/samtools flagstat align.hgr1.bam > align.hgr1.flagstat

      # Clean up 
      rm *tmp* align.F4.bam align.f4F8.bam

    echo "Processing complete. Copy files to S3"
    # Rename the total Bam Files
      mv aligned_unsorted.bam $OUTPUT.bam
      #mv aligned_unsorted.bam.bai $OUTPUT.bam.bai
      mv aligned_unsorted.flagstat $OUTPUT.flagstat

    # Rename the hgr-aligned Bam files
      mv align.hgr1.bam $OUTPUT.hgr1.bam
      mv align.hgr1.bam.bai $OUTPUT.hgr1.bam.bai
      mv align.hgr1.flagstat $OUTPUT.hgr1.flagstat
      
    # Copy output to AWS S3 
    # Alignments (Full)
      aws s3 cp $OUTPUT.flagstat $S3URL/hgr1/

    # Alignments (Aligned)
      aws s3 cp $OUTPUT.hgr1.bam $S3URL/hgr1/
      aws s3 cp $OUTPUT.hgr1.bam.bai $S3URL/hgr1/
      aws s3 cp $OUTPUT.hgr1.flagstat $S3URL/hgr1/

#fi
 
echo "Create log files and copy to S3"
# Copy screen log file to AWS S3
cp ~/screenlog.0 ./$OUTPUT.screenlog
aws s3 cp $OUTPUT.screenlog $S3URL/logs/

# Shutdown and Terminate instance
EC2ID=$(ec2metadata --instance-id)
sleep 20s # to catch errors

if [ "$TERMINATE" = TRUE ]
then
  echo "Run Complete -- Shutting down instance."
  aws ec2 terminate-instances --instance-ids $EC2ID
else
  echo "Run Complete -- Instance is online."
fi

# Script complete
