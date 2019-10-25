#!/bin/bash
# 1kgpilot_2
# DNS: ec2-52-34-12-139.us-west-2.compute.amazonaws.com
# AMI: crown-170220 - ami-66129306
# EC2: c4.2xlarge (8cpu / 15 gb)
# Storage: 1000 Gb
# Start: 
# Alignment done: 
# Align.subset done: 
# End:
#
# CMD:
# ec2-52-34-12-139.us-west-2.compute.amazonaws.com
# sh 1kgpilot_2.sh NA12878_pp NA12878 ERR194147 CEU ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

# ec2-52-27-70-31.us-west-2.compute.amazonaws.com
# sh 1kgpilot_2.sh NA12891_pp NA12891 ERR194160 CEU ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194160/ERR194160_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194160/ERR194160_2.fastq.gz



# Control Panel -------------------------------
# CPU
	THREADS='7'

# Sequencing Data
	LIBRARY=$1 # Library/ File name
	FASTQ1=$5
	FASTQ2=$6

    # File-names
    FQ1=$(basename $FASTQ1)
    FQ2=$(basename $FASTQ2)

# Read Group Data
	RGSM=$2   # Sample. Patient Identifer
	RGID=$3 # Read Group ID. Accession Number
	RGLB=$LIBRARY # Library Name. Accession Number
	RGPL='ILLUMINA'  # Sequencing Platform.
	RGPO=$4 # Patient Population
	# Extract Sequencing Run Info
	#  RGPU=$(gzip -dc $FQ1 | head -n1 - | cut -f1 -d':' | cut -f2 -d' ')


# Initialize wordir ---------------------------

# Make working directory
  mkdir -p align; cd align

# Copy hgr genome and create bowtie2 index
  cp ~/resources/hgr_45s.fa ./
  samtools faidx hgr_45s.fa
  
  bowtie2-build hgr_45s.fa hgr
  
# Download Genome Sequencing Data
  wget $FASTQ1
  wget $FASTQ2

    # Extract Sequencing Run Info
    RGPU=$(gzip -dc $FQ1| head -n1 - | cut -f1 -d':' | cut -f2 -d' ')

# Primary Alignment -------------------------

# Bowtie2: align to genome

bowtie2 --very-sensitive-local -p $THREADS --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \
--rg PL:$RGPL --rg PU:$RGPU -x hgr -1 $FQ1 -2 $FQ2 | samtools view -bS - > aligned_unsorted.bam

#rm $FQ1 $FQ2 # Remove fastq files to save space

# Sort alignment file
#  samtools sort -@ $THREADS aligned_unsorted.bam aligned
#  samtools index aligned.bam
#  rm aligned_unsorted.bam

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
  cat align.header.tmp align.F4.tmp.sam align.f4F8.tmp.sam | samtools view -bS - > align.hgr.tmp.bam
    samtools sort align.hgr.tmp.bam align.hgr
    samtools index align.hgr.bam
    samtools flagstat align.hgr.bam > align.hgr.flagstat
    
  # Read Counts: align.hgr0.bam (NA19240_pcr)
    # 651340 + 0 in total (QC-passed reads + QC-failed reads)
    # 0 + 0 duplicates
    # 614264 + 0 mapped (94.31%:-nan%)
    # 651340 + 0 paired in sequencing
    # 325670 + 0 read1
    # 325670 + 0 read2
    # 166576 + 0 properly paired (25.57%:-nan%)
    # 577188 + 0 with itself and mate mapped
    # 37076 + 0 singletons (5.69%:-nan%)
    # 0 + 0 with mate mapped to a different chr
    # 0 + 0 with mate mapped to a different chr (mapQ>=5)
  
  rm *tmp* align.F4.bam align.f4F8.bam # Clean-up

# Rename the total Bam Files
  mv aligned_unsorted.bam $LIBRARY.bam
  mv aligned_unsorted.bam.bai $LIBRARY.bam.bai
  mv aligned_unsorted.flagstat $LIBRARY.flagstat

# Rename the hgr Bam files
  mv align.hgr.bam $LIBRARY.hgr.bam
  mv align.hgr.bam.bai $LIBRARY.hgr.bam.bai
  mv align.hgr.flagstat $LIBRARY.hgr.flagstat
  
# Primary VCF ----------------------------

# GATK variant calling over 45S region
#  aws s3 cp s3://crownproject/resources/hgr.gatk.fa ./
#  aws s3 cp s3://crownproject/resources/hgr.gatk.fa.fai ./
#  aws s3 cp s3://crownproject/resources/hgr.gatk.dict ./
  
#  java -Xmx12G -jar /home/ubuntu/software/GenomeAnalysisTK.jar \
#  -R hgr.gatk.fa -T HaplotypeCaller \
#  -ploidy 2 --max_alternate_alleles 6 \
#  -I $LIBRARY.bam -o $LIBRARY.hgr.vcf
   # Memory issues, restrict to 45S region only
     # -ploidy 100, 50, 20 failed... do 2 and analyze 45S further
     
# Upload final output files to S3
 
# Alignments (Full)
 #aws s3 cp $LIBRARY.bam s3://crownproject/1kg_hgr0/
 #aws s3 cp $LIBRARY.bam.bai s3://crownproject/1kg_hgr0/
 aws s3 cp $LIBRARY.flagstat s3://crownproject/1kg_pilot/

# Alignments (Aligned)
  aws s3 cp $LIBRARY.hgr.bam s3://crownproject/1kg_pilot/
  aws s3 cp $LIBRARY.hgr.bam.bai s3://crownproject/1kg_pilot/
  aws s3 cp $LIBRARY.hgr.flagstat s3://crownproject/1kg_pilot/

# VCF
 aws s3 cp $LIBRARY.hgr.vcf s3://crownproject/1kg_pilot/
 aws s3 cp $LIBRARY.hgr.vcf.idx s3://crownproject/1kg_pilot/
 
# Shutdown and Terminate instance
EC2ID=$(ec2metadata --instance-id)
#aws ec2 terminate-instances --instance-ids $EC2ID

# Script complete
