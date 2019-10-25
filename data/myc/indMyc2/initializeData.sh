#!/bin/bash
# initialize Data for indMyc2 experiment

# Prefetch SRA data
prefetch SRR2939547
prefetch SRR2939541
prefetch SRR2939537
prefetch SRR2939538
prefetch SRR2939548
prefetch SRR2939542
# Convert to Fastq
fastq-dump --gzip SRR2939547
fastq-dump --gzip SRR2939541
fastq-dump --gzip SRR2939537
fastq-dump --gzip SRR2939538
fastq-dump --gzip SRR2939548
fastq-dump --gzip SRR2939542

