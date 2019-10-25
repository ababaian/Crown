#!/bin/bash
#
# Initialize data from SRA
# for IFNg experiment
# Donor 23 PBMC only for now

# Prefetch
prefetch SRR3623724
prefetch SRR3623727
prefetch SRR3623730
prefetch SRR3623722
prefetch SRR3623725
prefetch SRR3623728
prefetch SRR3623723
prefetch SRR3623726
prefetch SRR3623729
prefetch SRR3623733
prefetch SRR3623736
prefetch SRR3623739
prefetch SRR3623731
prefetch SRR3623734
prefetch SRR3623737
prefetch SRR3623732
prefetch SRR3623735
prefetch SRR3623738


# Extract to Fastq
fastq-dump --gzip SRR3623724
fastq-dump --gzip SRR3623727
fastq-dump --gzip SRR3623730
fastq-dump --gzip SRR3623722
fastq-dump --gzip SRR3623725
fastq-dump --gzip SRR3623728
fastq-dump --gzip SRR3623723
fastq-dump --gzip SRR3623726
fastq-dump --gzip SRR3623729
fastq-dump --gzip SRR3623733
fastq-dump --gzip SRR3623736
fastq-dump --gzip SRR3623739
fastq-dump --gzip SRR3623731
fastq-dump --gzip SRR3623734
fastq-dump --gzip SRR3623737
fastq-dump --gzip SRR3623732
fastq-dump --gzip SRR3623735
fastq-dump --gzip SRR3623738


# Run IFNg_align_v0.sh

