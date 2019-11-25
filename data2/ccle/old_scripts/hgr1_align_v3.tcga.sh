#!/bin/bash
# hgr1_align_v3.tcga.sh
# rDNA alignment pipeline
PIPE_VERSION='190506 build -- TCGA'
AMI_VERSION='crown-180813 - ami-0031fd61f932bdef9'
# EC2: c4.2xlarge (8cpu / 15 gb)
# EC2: c4.xlarge  (4cpu / 8  gb)
# Storage: 200 Gb
#

# Input Requirements --------------------------

# $1 : Library name and file-output name (unique)
# $2 : Library population/analysis set
# $3 : Library UUID

# Control Panel -------------------------------
# Amazon AWS S3 Home URL
  S3URL='s3://crownproject/tcga2'

# CPU
	THREADS='3'

# Terminate instances upon completion (for debuggin)
  TERMINATE='TRUE'

# Sequencing Data
	LIBRARY=$1 # Library/ File name

# TCGA FILE UUID
  UUID=$3

 # FastQ File-names
    FQ0="$LIBRARY.tmp.sort.0.fq"
    FQ1="$LIBRARY.tmp.sort.1.fq"
    FQ2="$LIBRARY.tmp.sort.2.fq"
    
# Read Group Data
# Extract from downloaded BAM file / input
	RGPO=$2 # Patient Population

	#RGSM= # Sample. Patient Identifer
	#RGID= # Read Group ID. Accession Number
    
	RGLB=$LIBRARY # Library Name. Accession Number
	RGPL='ILLUMINA'  # Sequencing Platform.
    
	# Extract Sequencing Run Info
	#  RGPU=$(gzip -dc $FQ1 | head -n1 - | cut -f1 -d':' | cut -f2 -d' ')

echo " -- hgr1 Alignment Pipeline -- "
echo " version: $PIPE_VERSION "
echo " ami:     $AMI_VERSION  "
echo " s3:      $S3URL  "
echo " library: $LIBRARY  "
echo " date:    $(date)"
echo ''

# Initialize wordir ---------------------------
echo 'Initializing ...'

# Make working directory
  mkdir -p align; cd align

# Copy hgrX genome and create bowtie2 index
  cp ~/resources/hgr1/* ./

# Load updated GDC token to instance
aws s3 cp $S3URL/scripts/gdc.token gdc.token
chmod 600 gdc.token

echo "Download BAM file: $UUID"
  
# Download RNAseq BAM file
# with a UUID as input
  ~/bin/gdc-client download -t gdc.token -d ./ \
  -n $THREADS $UUID
  
# Link the RNA-seq bamfile which is called by its UID to workdir
  ln -s */*.bam input.bam
  
# Extract ReadGroup Sample Name (SM)
  RGSM=$(~/bin/samtools view -H input.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq )

# Extract ReadGroup identifier (ID)
  #RGID=$(~/bin/samtools view -H input.bam | grep '^@RG' | sed "s/.*ID:\([^\t]*\).*/\1/g" | uniq )
  RGID=$(~/bin/samtools view -H input.bam | grep '^@RG' | sed "s/.*ID:\([^\t]*\).*/\1/g" | head -n1 - )

echo " Convert BAM file to fastq file"

# Convert input bam file to fastq files for re-alignment
~/bin/samtools sort -@ $THREADS-n input.bam | \
    ~/bin/samtools fastq -@ $THREADS -O \
    -0 $FQ0 \
    -1 $FQ1 \
    -2 $FQ2 -

# SINGLE END READS ====================================================

if [ -s $FQ0 ]
then
    echo "Single-end Reads Detected."
    echo ''
    echo "Starting hgr1 alignment"
    # Single-End Extracted Reads Alignment

    # Extract Sequencing Run Info
    #RGPU=$(gzip -dc $FQ0| head -n1 - | cut -f1 -d':' | cut -f2 -d' ')
    RGPU=$(head -n1 $FQ0 | cut -f1 -d':' | cut -f2 -d' ')

    # Bowtie2: align to genome
    bowtie2 --very-sensitive-local -p $THREADS \
      --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \
      --rg PL:$RGPL --rg PU:$RGPU \
      -x hgr1 -U $FQ0 | \
      ~/bin/samtools view -bS - > aligned_unsorted.bam
     
    echo "Alignment complete."
    echo "Calculate flagstats."

    # Calculate library flagstats
      ~/bin/samtools flagstat aligned_unsorted.bam > aligned_unsorted.flagstat
      rm $FQ0 # Remove fastq files to save space

    echo "Subset reads (retain mapped, remove unmapped)."

    # Read Subset ------------------------------
    # Extract mapped reads, and their unmapped pairs

      # Extract Header
      ~/bin/samtools view -H aligned_unsorted.bam > align.header.tmp

      # Extract Mapped Reads
      ~/bin/samtools view -b -F 4 aligned_unsorted.bam | \
      ~/bin/samtools sort -@ $THREADS -O BAM - > align.hgr1.bam #mapped
      
    # Calcualte library flagstats
    ~/bin/samtools index align.hgr1.bam
    ~/bin/samtools flagstat align.hgr1.bam > align.hgr1.flagstat

    echo "Processing complete. Copy files to S3"

    # Rename the total Bam Files
      #mv aligned_unsorted.bam $LIBRARY.se.bam
      #mv aligned_unsorted.bam.bai $LIBRARY.se.bam.bai
      mv aligned_unsorted.flagstat $LIBRARY.se.flagstat

    # Rename the hgr-aligned Bam files
      mv align.hgr1.bam $LIBRARY.hgr1.se.bam
      mv align.hgr1.bam.bai $LIBRARY.hgr1.se.bam.bai
      mv align.hgr1.flagstat $LIBRARY.hgr1.se.flagstat
      
    # Alignments (Full)
    aws s3 cp $LIBRARY.se.flagstat $S3URL/$RGPO/

    # Alignments (Aligned)
    aws s3 cp $LIBRARY.hgr1.se.bam $S3URL/$RGPO/
    aws s3 cp $LIBRARY.hgr1.se.bam.bai $S3URL/$RGPO/
    aws s3 cp $LIBRARY.hgr1.se.flagstat $S3URL/$RGPO/


fi

# PAIRED END READS ====================================================


if [ -s $FQ1 ]
then
    echo "Paired-end Reads Detected."
    echo ''
    echo "Starting hgr1 alignment"
    # Paired-End Extracted Reads Alignment

    # Extract Sequencing Run Info
    #RGPU=$(gzip -dc $FQ1| head -n1 $FQ1 | cut -f1 -d':' | cut -f2 -d' ')
    RGPU=$(head -n1 $FQ1 | cut -f1 -d':' | cut -f2 -d' ')
    
    # Bowtie2: align to genome
    bowtie2 --very-sensitive-local -p $THREADS \
      --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \
      --rg PL:$RGPL --rg PU:$RGPU \
      -x hgr1 -1 $FQ1 -2 $FQ2 | \
      ~/bin/samtools view -bS - > aligned_unsorted.bam
      
    echo "Alignment complete."
    echo "Calculate flagstats."

    # Calcualte library flagstats
      ~/bin/samtools flagstat aligned_unsorted.bam > aligned_unsorted.flagstat
      
      rm $FQ1 $FQ2 # Remove fastq files to save space

      
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
      mv aligned_unsorted.bam $LIBRARY.bam
      #mv aligned_unsorted.bam.bai $LIBRARY.bam.bai
      mv aligned_unsorted.flagstat $LIBRARY.flagstat

    # Rename the hgr-aligned Bam files
      mv align.hgr1.bam $LIBRARY.hgr1.bam
      mv align.hgr1.bam.bai $LIBRARY.hgr1.bam.bai
      mv align.hgr1.flagstat $LIBRARY.hgr1.flagstat
      
    # Copy output to AWS S3 
    # Alignments (Full)
      aws s3 cp $LIBRARY.flagstat $S3URL/$RGPO/

    # Alignments (Aligned)
      aws s3 cp $LIBRARY.hgr1.bam $S3URL/$RGPO/
      aws s3 cp $LIBRARY.hgr1.bam.bai $S3URL/$RGPO/
      aws s3 cp $LIBRARY.hgr1.flagstat $S3URL/$RGPO/

fi
 
echo "Create log files and copy to S3"
# Copy screen log file to AWS S3
cp ~/screenlog.0 ./$LIBRARY.screenlog
aws s3 cp $LIBRARY.screenlog $S3URL/logs/

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
