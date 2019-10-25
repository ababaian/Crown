#!/bin/bash

# Prefetch data from SRA
prefetch SRR2939550
prefetch SRR2939551
prefetch SRR2939552
prefetch SRR2939554
prefetch SRR2939555

# Dump to FQ files
fastq-dump --gzip  SRR2939550
fastq-dump --gzip  SRR2939551
fastq-dump --gzip  SRR2939552
fastq-dump --gzip  SRR2939554
fastq-dump --gzip  SRR2939555
