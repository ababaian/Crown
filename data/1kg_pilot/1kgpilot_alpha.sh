#!/bin/bash
# 1kgpilot.sh script

# Control Panel -------------------------------
# CPU
	THREADS='2'
	AWSID='#########'

# Sequencing Data
	LIBRARY='HG03118_lowcov' # Library/ File name
	FASTQ1='ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/ERR187387/ERR187387_1.fastq.gz'
	FASTQ2='ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/ERR187387/ERR187387_2.fastq.gz'

# Read Group Data
	RGID='ERR187387' # Read Group ID. Accession Number
	RGSM='HG03118'   # Sample. Patient Identifer
	RGLB=$LIBRARY # Library Name. Accession Number
	RGPL='ILLUMINA'  # Sequencing Platform.
	RGPO='ESN' # Patient Population
	# Extract Sequencing Run Info
	#  RGPU=$(gzip -dc SRR794330_1.fastq.gz | head -n1 - | cut -f1 -d':' | cut -f2 -d' ')


# Initialize wordir ---------------------------

# Make working directory
  mkdir align; cd align

# Copy hgr genome and create bowtie2 index
  cp ~/resources/hgr.fa ./
  cp ~/resources/hgr.fa.fai ./
  
  bowtie2-build hgr.fa hgr
  
# Download Genome Sequencing Data
  wget $FASTQ1
  wget $FASTQ2

    # Extract Sequencing Run Info
    RGPU=$(gzip -dc SRR794330_1.fastq.gz | head -n1 - | cut -f1 -d':' | cut -f2 -d' ')


# Treehash calculators / Glacier Upload script
aws s3 cp s3://crownproject/scripts/treehash.py ./
 chmod 755 treehash.py
aws s3 cp s3://crownproject/scripts/glacierupload.sh ./
 chmod 755 glacierupload.sh

# Primary Alignment -------------------------

# Bowtie2: align to genome

bowtie2 --very-sensitive-local -p $THREADS --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM --rg PL:$RGPL --rg PU:$RGPU --rg PO:$RGPO -x hgr -1 SRR794330_1.fastq.gz -2 SRR794330_2.fastq.gz | samtools view -bS - > aligned_unsorted.bam

# Sort alignment file
  samtools sort aligned_unsorted.bam aligned
  samtools index aligned.bam
  rm aligned_unsorted.bam


# Read Subset ------------------------------
# Extract mapped reads, and their unmapped pairs

  # Extract Header
  samtools view -H aligned.bam > align.header.tmp

  # Unmapped reads with mapped pairs
  # Extract Mapped Reads
  # and their unmapped pairs
  samtools view -b -F 4 aligned.bam > align.F4.bam #mapped
  samtools view -b -f 4 -F 8 aligned.bam > align.f4F8.bam #unmapped pairs
  
  # Extract just the 45S unit
  echo -e "chr13\t999000\t1013500" > rDNA_45s.bed
  samtools view -b -L rDNA-45s.bed align.F4.bam > align.F4.45s.bam
  
  # What are the 45S mapped readnames
  samtools view align.F4.45s.bam | cut -f1 - > read.names.45s.tmp
  
  # Extract cases of read pairs mapped on edge of region of interest
  # -------|======= R O I ======| ----------
  # read:                  ====---====

  samtools view align.F4.bam | grep -Ff read.names.45s.txt - > align.F4.tmp.sam

  # Complete mapped reads list
  cut -f1 align.F4.tmp.sam > read.names.45s.long.tmp

  # Extract unmapped reads with a mapped pair in 45s.long
  samtools view align.f4F8.bam | grep -Ff read.names.45s.long.tmp - > align.f4F8.tmp.sam

  # Re-compile bam file
  cat align.header.tmp align.F4.tmp.sam align.f4F8.tmp.sam | samtools view -bS - > align.45s.tmp.bam
    samtools sort align.45s.tmp.bam align.45s
    samtools index align.45s.bam
    rm align.45s.tmp.bam
    
  # Read Counts: align.45s
  # 66736 Reads
  # 60877 mapped (91.22%)
  # 33368 read1; 33368 read2
  
  # Recompile the total rDNA mapped read pairs
  samtools cat -h align.header.tmp -o align.rDNA.bam align.F4.bam align.f4F8.bam
  samtools index align.rDNA.bam
  
# Calcualte library flagstats
  samtools flagstat aligned.bam > aligned.flagstat
  samtools flagstat align.45s.bam > align.45s.flagstat
  samtools flagstat align.rDNA.bam > align.rDNA.flagstat
  
  rm *tmp* # Clean-up

# Rename the final Bam Files
  mv aligned.bam $LIBRARY.bam
  mv aligned.bam.bai $LIBRARY.bam.bai
  mv aligned.flagstat $LIBRARY.flagstat

  mv align.rDNA.bam $LIBRARY.rDNA.bam
  mv align.rDNA.bam.bai $LIBRARY.rDNA.bam.bai
  mv align.rDNA.flagstat $LIBRARY.rDNA.flagstat

  mv align.45s.bam $LIBRARY.45s.bam
  mv align.45s.bam.bai $LIBRARY.45s.bam.bai
  mv align.45s.flagstat $LIBRARY.45s.flagstat


# GATK variant calling over 45S region
  aws s3 cp s3://crownproject/resources/hgr.gatk.fa ./
  aws s3 cp s3://crownproject/resources/hgr.gatk.fa.fai ./
  aws s3 cp s3://crownproject/resources/hgr.gatk.dict ./
  
  java -Xmx12G -jar /home/ubuntu/software/GenomeAnalysisTK.jar \
  -R hgr.gatk.fa -T HaplotypeCaller \
  -ploidy 2 --max_alternate_alleles 6 \
  -I $LIBRARY.bam -o NA19240.rDNA_p2.vcf
   # Memory issues, restrict to 45S region only
     # -ploidy 100, 50, 20 failed... do 2 and analyze 45S further
     
  java -Xmx12G -jar /home/ubuntu/software/GenomeAnalysisTK.jar \
  -R hgr.gatk.fa -T HaplotypeCaller \
  -ploidy 50 --max_alternate_alleles 10 \
  -I $LIBRARY.45s.bam -o NA19240.45s_p50.vcf


# Upload final output files to S3
 
# 45S Only Alignent
 aws s3 cp $LIBRARY.45s.bam s3://crownproject/1kg_pilot/
 aws s3 cp $LIBRARY.45s.bam.bai s3://crownproject/1kg_pilot/
 
# rDNA Only Alignent
 aws s3 cp $LIBRARY.rDNA.bam s3://crownproject/1kg_pilot/
 aws s3 cp $LIBRARY.rDNA.bam.bai s3://crownproject/1kg_pilot/

#VCF over 45S
 aws s3 cp NA19240.45s_p50.vcf s3://crownproject/1kg_pilot/
 aws s3 cp NA19240.45s_p50.vcf.idx s3://crownproject/1kg_pilot/
 
#VCF over rDNA
 aws s3 cp NA19240.rDNA_p2.vcf s3://crownproject/1kg_pilot/
 aws s3 cp NA19240.rDNA_p2.vcf.idx s3://crownproject/1kg_pilot/

  
# Upload full bam file to amazon Glacier vault
  ./glacierupload.sh $LIBRARY.bam $AWSID
    aws s3 cp $LIBRARY.bam.glacer s3://crownproject/1kg_pilot/


# Script complete
