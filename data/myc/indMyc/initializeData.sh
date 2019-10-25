#!/bin/bash

# Prefetch SRA data
prefetch SRR2069978
prefetch SRR2069979
prefetch SRR2069980
prefetch SRR2069983
prefetch SRR2069984
prefetch SRR2069985
# Convert to Fastq
fastq-dump --gzip SRR2069978
fastq-dump --gzip SRR2069979
fastq-dump --gzip SRR2069980
fastq-dump --gzip SRR2069983
fastq-dump --gzip SRR2069984
fastq-dump --gzip SRR2069985

