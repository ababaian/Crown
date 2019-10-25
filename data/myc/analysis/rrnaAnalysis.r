# rrnaAnalysis.R
#
# Analysis of adCalc.sh
# output gvcf files
#
library(ggplot2)

# Import
GVCF = read.table('28S_siMyc.gvcf')
  GVCF = data.frame(t(GVCF))
  colnames(GVCF) = seq(1, length(GVCF[1,]))

# Cut 28S to just 28S (5071 bases)
# chr13    1007948    1013018    28S
  # trim = which( (apply(GVCF[2,], 1, as.numeric) < 1007949) |
  #                 (apply(GVCF[2,], 1 , as.numeric) > 1013018) )
  # GVCF = GVCF[, -trim]

  
# Cut 18S to just 18S (1869 bases)
# chr13    1003661    1005529    18S
  # trim = which( (apply(GVCF[2,], 1, as.numeric) < 1003661) |
  #               (apply(GVCF[2,], 1 , as.numeric) > 1005529) )
  # GVCF = GVCF[, -trim]

# Cut 5.8S to just 5.8S (1869 bases)
# chr13    1006622    1006779    18S
  # trim = which( (apply(GVCF[2,], 1, as.numeric) < 1006623) |
  #                 (apply(GVCF[2,], 1 , as.numeric) > 1006779) )
  # GVCF = GVCF[, -trim]
  

refAllele = GVCF[4,]
altAllele = GVCF[5,]
genCoord  = GVCF[2,]
rnaCoord  = seq(1, length(GVCF[2,]))

sampleN = length(GVCF[,1]) - 9 # remove 9 header vcf rows
bpN     = length(genCoord) # 1869 for 18S; 5071 for 28S

# Functions =========================================================

# Convert DP:AD string to numeric DP (Total Depth)
dpCalc = function(inSTR){
# inSTR is from vcf
# in format DP:AD
# 2000:1500,400,50,50
# extract 2000
inSTR = as.character(inSTR)
as.numeric(unlist(strsplit(inSTR,split=':'))[1])

}

# Convert DP:AD string to numeric RD for the REFERENCE ALLELE DEPTH
# Thus Alternative_Allele_Depth = Total_Depth - Reference_Allele_Depth
# for all alternative alleles.
rdCalc = function(inSTR){
  # inSTR is from vcf
  # in format DP:AD
  # 2000:1500,400,50,50
  # extract 1500
  inSTR = as.character(inSTR)
  as.numeric(unlist(strsplit(unlist(strsplit(inSTR,split=":"))[2], split = ","))[1])
  
}


# Calculations ======================================================
# Calculate Depth of Coverage (baq > 30)
# for all positions

#Initialize DP vector
DP = vapply( GVCF[-c(1:9),1], dpCalc, 1)

#Extend the DP vector for all positions
for (i in 2:bpN){
DP = cbind(DP,
           vapply( GVCF[-c(1:9),i], dpCalc, 1) )
}

# Calculate Reference Depth of Coverage (baq > 30)
# for all positions
#
#Initialize
RD = vapply( GVCF[-c(1:9),1], rdCalc, 1)

#The rest
for (i in 2:bpN){
  RD = cbind(RD,
             vapply( GVCF[-c(1:9),i], rdCalc, 1) )
}


# Reference Allele Frequency
# Intra-Library
# RD / DP
RAF = RD / DP

# NOTE: division by zero is possible here and will introduce NAs


# Deconvolute Sample A from Sample B
# siCTRL = Row 1 - 3
idxA = 1:3
# siMyc  = Row 4 - 5
idxB = 4:5

RAF_A = RAF[idxA,]
RAF_B = RAF[idxB,]
              
# Change in Reference Allele Frequency
 # (For Paired-Samples) ----------------------

# dRAF = RAF_A - RAF_B
# 
# # Calculate some descriptive statistics
# # about the change in Reference Allele Frequency
# # Remove NA from calculations (no sequencing depth in a library)
# mean_dRAF = apply(dRAF,2,mean, na.rm = TRUE)
# sd_dRAF   = apply(dRAF,2,sd, na.rm = TRUE)
# var_dRAF  = apply(dRAF,2,var, na.rm = TRUE)
# mean_DP = apply(DP,2,mean, na.rm = TRUE)

# Remove poorly 'covered' positions (i.e. less then 1000x coverage on average)
# the magnitude of bias is simply too high at such regions
mean_DP = apply(DP,2, mean, na.rm = TRUE)
dropPOS = (mean_DP < 100)

RAF_A[,dropPOS]  = 0
RAF_B[,dropPOS] = 0

 # (For unpaired samples) --------------------
 mean_dRAF = apply(RAF_A, 2, mean, na.rm = TRUE) - 
  apply(RAF_B, 2, mean, na.rm = TRUE)

dRAF = mean_dRAF

plot(mean_dRAF)
