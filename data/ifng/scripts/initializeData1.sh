#!/bin/bash
#
# Initialize data from SRA
# for IFNg experiment
# Donor 1 PBMC only for now
mkdir -p align
cd align


# Prefetch
prefetch SRR3623715
prefetch SRR3623718
prefetch SRR3623721
prefetch SRR3623713
prefetch SRR3623716
prefetch SRR3623719
prefetch SRR3623714
prefetch SRR3623717
prefetch SRR3623720

# Extract to Fastq
fastq-dump --gzip SRR3623715
fastq-dump --gzip SRR3623718
fastq-dump --gzip SRR3623721
fastq-dump --gzip SRR3623713
fastq-dump --gzip SRR3623716
fastq-dump --gzip SRR3623719
fastq-dump --gzip SRR3623714
fastq-dump --gzip SRR3623717
fastq-dump --gzip SRR3623720

# Run IFNg_align_v0.sh

