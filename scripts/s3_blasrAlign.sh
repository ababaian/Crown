#!/bin/sh
# s3_blasrAlign.sh
#
# Usage: s3_blasrAlign.sh <s3_dir> <file list in s3> <genome.fa>
# For each file in file_list in the s3_dir
# download the fastq file;
# align it to the genome.fa
# compress output sam into a bam file
# remove fastq file and begin with next one
#

# INPUT =====================

# S3 directory containing fastq files
S3_DIR=$1

# List of FASTQ filenames in S3_DIR to download
# iteratively
FILE_LIST=$2

# Genome to align to
GENOME=$3

# Count Variable
COUNT='1'

# SCRIPT ====================
for FILE in $(cat $FILE_LIST)
do
	FILEPATH=$(echo $S3_DIR/$FILE)

	echo "Starting download of $FILEPATH"
	echo ""
	
	# aws configure must be run beforehand
	aws s3 cp $FILEPATH ./

	# Use blasr to align $FILE to $GENOME
	# Can set different number of processors
	blasr $FILE $GENOME --nproc 2 --sam --out alignedTMP.sam

	# Remove fastq file
	rm $FILE

	# Sort, index and bam-file the output
	samtools view -bS alignedTMP.sam | samtools sort - aligned_$COUNT

	samtools index aligned_$COUNT.bam

	# Remove sam file
	rm alignedTMP.sam

	# Make a subset bam file of the
	# transcribed rDNA locus (18-28S)
	samtools view -b -L rDNA_main.bed aligned_$COUNT.bam | samtools sort - hgMain_$COUNT

	# Count
	COUNT=$((COUNT+1))
done
	 
# Concatenate all the individual bam files for hgrMain
samtools cat $(ls hgMain_*.bam) | samtools sort - na12878.pb.hgr_N

#~~~ End of Script ~~~ 
