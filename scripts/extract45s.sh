#!/bin/sh

BAM=$1

# Read Subset ------------------------------
# Extract mapped reads, and their unmapped pairs

  # Extract Header
  samtools view -H $BAM > align.header.tmp

  # Unmapped reads with mapped pairs
  # Extract Mapped Reads
  # and their unmapped pairs
  samtools view -b -F 4 $BAM > align.F4.bam #mapped
  samtools view -b -f 4 -F 8 $BAM > align.f4F8.bam #unmapped pairs
  
  # Extract just the 45S unit
  aws s3 cp s3://crownproject/resources/rDNA_45s.bed ./
  samtools view -b -L rDNA_45s.bed align.F4.bam > align.F4.45s.bam
  
  # What are the 45S mapped readnames
  samtools view align.F4.45s.bam | cut -f1 - > read.names.45s.tmp
  
  # Extract cases of read pairs mapped on edge of region of interest
  # -------|======= R O I ======| ----------
  # read:                  ====---====

  samtools view align.F4.bam | grep -Ff read.names.45s.tmp - > align.F4.tmp.sam

  # Complete mapped reads list
  cut -f1 align.F4.tmp.sam > read.names.45s.long.tmp

  # Extract unmapped reads with a mapped pair in 45s.long
  samtools view align.f4F8.bam | grep -Ff read.names.45s.long.tmp - > align.f4F8.tmp.sam

  # Re-compile bam file
  cat align.header.tmp align.F4.tmp.sam align.f4F8.tmp.sam | samtools view -bS - > align.45s.tmp.bam
    samtools sort align.45s.tmp.bam FINAL.45s
    samtools index FINAL.45s.bam
    rm *.tmp*
    rm align.*.bam
    

#  sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/p11_1/alignment/p11_1.bam 
#  sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/p18_1/alignment/p18_1.bam 
#  sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/natSen_1/alignment/natSen_1.bam 


sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/p11_2/alignment/p11_2.bam 
mv FINAL.45s.bam p11_2.45s.bam
mv FINAL.45s.bam.bai p11_2.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/p11_3/alignment/p11_3.bam 
mv FINAL.45s.bam p11_3.45s.bam
mv FINAL.45s.bam.bai p11_3.45s.bam.bai


sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/p18_2/alignment/p18_2.bam 
mv FINAL.45s.bam p18_2.45s.bam
mv FINAL.45s.bam.bai p18_2.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/p18_3/alignment/p18_3.bam 
mv FINAL.45s.bam p18_3.45s.bam
mv FINAL.45s.bam.bai p18_3.45s.bam.bai


sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/natSen_2/alignment/natSen_2.bam 
mv FINAL.45s.bam natSen_2.45s.bam
mv FINAL.45s.bam.bai natSen_2.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/natSen_3/alignment/natSen_3.bam 
mv FINAL.45s.bam natSen_3.45s.bam
mv FINAL.45s.bam.bai natSen_3.45s.bam.bai

# done: Fri Jan 13 15:42:24 PST 2017

# Look at p53 -/- immortalized cells for variation too
# Mon Jan 16 11:44:28 PST 2017

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/imm_1/alignment/imm_1.bam 
mv FINAL.45s.bam imm_1.45s.bam
mv FINAL.45s.bam.bai imm_1.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/imm_2/alignment/imm_2.bam 
mv FINAL.45s.bam imm_2.45s.bam
mv FINAL.45s.bam.bai imm_2.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/imm_3/alignment/imm_3.bam 
mv FINAL.45s.bam imm_3.45s.bam
mv FINAL.45s.bam.bai imm_3.45s.bam.bai


# H2O2 Treatment
sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/h202_1/alignment/h202_1.bam 
mv FINAL.45s.bam h202_1.45s.bam
mv FINAL.45s.bam.bai h202_1.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/h202_2/alignment/h202_2.bam 
mv FINAL.45s.bam h202_2.45s.bam
mv FINAL.45s.bam.bai h202_2.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/h202_3/alignment/h202_3.bam 
mv FINAL.45s.bam h202_3.45s.bam
mv FINAL.45s.bam.bai h202_3.45s.bam.bai

# aza Treatment

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/5aza_1/alignment/5aza_1.bam 
mv FINAL.45s.bam 5aza_1.45s.bam
mv FINAL.45s.bam.bai 5aza_1.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/5aza_2/alignment/5aza_2.bam 
mv FINAL.45s.bam 5aza_2.45s.bam
mv FINAL.45s.bam.bai 5aza_2.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/5aza_3/alignment/5aza_3.bam 
mv FINAL.45s.bam 5aza_3.45s.bam
mv FINAL.45s.bam.bai 5aza_3.45s.bam.bai

# quiescent cells

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/qui_1/alignment/qui_1.bam 
mv FINAL.45s.bam qui_1.45s.bam
mv FINAL.45s.bam.bai qui_1.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/qui_2/alignment/qui_2.bam 
mv FINAL.45s.bam qui_2.45s.bam
mv FINAL.45s.bam.bai qui_2.45s.bam.bai

sh extract45s.sh ~/Gprojects/LIONS/projects/azaSen/qui_3/alignment/qui_3.bam 
mv FINAL.45s.bam qui_3.45s.bam
mv FINAL.45s.bam.bai qui_3.45s.bam.bai

ls

