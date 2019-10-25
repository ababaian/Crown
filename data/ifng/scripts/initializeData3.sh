#!/bin/bash
#
# Initialize data from SRA
# for IFNg experiment
# Donor 4-5 PBMC only for now
cd align

# Prefetch
prefetch SRR3623742
prefetch SRR3623745
prefetch SRR3623748
prefetch SRR3623740
prefetch SRR3623743
prefetch SRR3623746
prefetch SRR3623741
prefetch SRR3623744
prefetch SRR3623747
prefetch SRR3623751
prefetch SRR3623754
prefetch SRR3623757
prefetch SRR3623749
prefetch SRR3623752
prefetch SRR3623755
prefetch SRR3623750
prefetch SRR3623753
prefetch SRR3623756

# Extract to Fastq
fastq-dump --gzip SRR3623742
fastq-dump --gzip SRR3623745
fastq-dump --gzip SRR3623748
fastq-dump --gzip SRR3623740
fastq-dump --gzip SRR3623743
fastq-dump --gzip SRR3623746
fastq-dump --gzip SRR3623741
fastq-dump --gzip SRR3623744
fastq-dump --gzip SRR3623747
fastq-dump --gzip SRR3623751
fastq-dump --gzip SRR3623754
fastq-dump --gzip SRR3623757
fastq-dump --gzip SRR3623749
fastq-dump --gzip SRR3623752
fastq-dump --gzip SRR3623755
fastq-dump --gzip SRR3623750
fastq-dump --gzip SRR3623753
fastq-dump --gzip SRR3623756

# Run IFNg_align_v0.sh

