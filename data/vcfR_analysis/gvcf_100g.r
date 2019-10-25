# vcfR analysis
# 100 genomes hgr1
# 170311
#
# 100g_gvcf  (g variant input only) analysis

# Install vcfR
#install.packages("vcfR")
library("vcfR")
library("ggplot2")
library("reshape2")

# File pointers

  vcf_file = '100g.hgr1.g.vcf'
  dna_file = 'rDNA.fa'
  gff_file = 'rDNA.gff'

# output Prefix (for plots)
  outPrefix='100g_gvcf'
  
# FUNCTIONS =========================================================
  
## Count number of NODATA per variant position
N_NODATA = function(inputROW){
  # For an input ROW of DP (or any data)
  # returns how many in the row are NA
  
  return(  length(which(is.na(inputROW))) )
  
}

## Count Variant positions (iVAF > 0)
# N_VAR = function(inputROW){
#   # For an input ROW of DP (or any data)
#   # returns how many in the row are NA
#   
#   return( length(which(inputROW > 0 )) )
#   
# }

## Count Variant positions (iVAF > 0)
N_VAR = function(inputROW){
  # For an input ROW of DP (or any data)
  # returns how many in the row are NA
  
  return( length(which(inputROW > 0.02 )) )
  
}

# SCRIPT ============================================================
#====================================================================
# Import VCF / DNA / GFF
  VCF = read.vcfR(vcf_file)
  DNA = ape::read.dna(dna_file, format = 'fasta')
  GFF = read.table(gff_file, sep="\t", quote = "")
  
# Create chromR Object
  chrom = create.chromR(name="rDNA", vcf=VCF, seq=DNA, ann=GFF)

#plot(chrom)

## Quality Metrics ==================================================
 plot(chrom) # Standard chromR plot for VCF
 chromoqc(chrom, dp.alpha = 22)


 # Quality of each variant Vector
 #system(paste0("cut -f 6 ",vcf_file," > tmp.qual"))
 #vcf.quality = read.csv('tmp.qual', sep = "\t", quote ="", col.names = T)

# # Describe Quality stats (log scale)
# ggplot(as.data.frame(vcf.quality), aes(vcf.quality)) + geom_histogram() +
#    theme_minimal() + scale_x_log10()
#    xlab('Alternative Genotype Quality, PHRED') +
#    ylab('Count')

  
# Variant Statistics ================================================

# Intra-individual variation (i*)
  
  # Depth of coverage at each variant (used reads)
  DP = extract.gt(chrom, element="DP", as.numeric=TRUE)

  # Extract Genotype Quality per sample
  GQ = extract.gt(chrom, element="GQ", as.numeric=TRUE)
  
  # Reference allele-only depth
  RAD = extract.gt(chrom, element="AD", as.numeric=TRUE)

  # Variant allele frequency (intra-individual)
  # Reference Allele Frequency (RAD / DP)
  iVAF = (DP-RAD)/DP
  
# Population variation (called variants)
  
  # Number of individual genomes in total population
  N_pop = length(DP[1,])
  
  # Column-length vector of how many NODATA there are in each
    # column of DP
  N_nd = apply(DP, 1, N_NODATA)
  
  # Number of individuals with measured data at each variant
  N_measured = N_pop - N_nd
  
  # Variant allele count (population)
  # number of people carrying measured variant
    # pV = rowSums(!is.na(iVAF))
  pV = apply(iVAF, 1, N_VAR)
  
  # Population Variant Allele Frequency
  pVAF = pV / N_measured
  
  # Remove 'zero' from the data
  # when measuring intra-genomic VAF, only consider
  # when the variant is present
  noZero_iVAF = iVAF
  noZero_iVAF[which(iVAF == 0)] = NA
  
  # Average intra-genomic variant allele frequency
  mean_iVAF = apply(iVAF,1,mean, na.rm = T)
  sd_iVAF = apply(iVAF,1,sd, na.rm = T)
  
  # Average intra-genomic variant allele frequency in only variant called
  mean_noZero_iVAF = apply(noZero_iVAF,1,mean, na.rm = T)
  sd_noZero_iVAF = apply(noZero_iVAF, 1, sd, na.rm = T)
  
# Population Variant Statistics
  POPVAR = data.frame(pV, N_measured, pVAF, mean_iVAF, sd_iVAF, mean_noZero_iVAF, sd_noZero_iVAF)
  
# Highlight variants at a mean_iVAF between two values ( 0.33 < iVAF < 0.66)
# Common, abundant variants
  CAV = which(0.33 < POPVAR$mean_noZero_iVAF & POPVAR$mean_noZero_iVAF < 0.66 & POPVAR$pVAF > 0.5)
  #highlight2 = which(POPVAR$mean_noZero_iVAF >= 0.66)
  
  POPVAR$colr = 'black'
  POPVAR$colr[CAV] = 'red'
  #POPVAR$colr[highlight2] = 'purple'

# Plot--------------------------------------------------------------

# Describe population stats
ggplot(as.data.frame(pV), aes(pV)) + geom_histogram() +
  theme_minimal() + scale_x_log10() +
  xlab('# Genomes variant was measured in') + 
  ylab('Count')

ggplot(as.data.frame(mean_noZero_iVAF), aes(mean_noZero_iVAF)) + geom_histogram() +
  theme_minimal() + 
  xlab('Intra-genomic VAF') + 
  ylab('Count')

ggplot(as.data.frame(pVAF), aes(pVAF)) + geom_histogram() +
  theme_minimal() +
  xlab('Population VAF') + 
  ylab('Count')

# Mean Intra-individual Variant Allele Frequency vs.
# Population-level Variant Allele Frequency (called vs not called)
# open PDF device
pdf(file = paste0(outPrefix,".MeaniVAF_pVAF_Nm.pdf"), width = 6.5, height = 5)

PLOT1 = ggplot(POPVAR, aes(pVAF, mean_noZero_iVAF ))
#PLOT1 = PLOT1 + geom_point(alpha = 0.5, stroke = 0, aes(size = pV/2, color = colr))
PLOT1 = PLOT1 + geom_point(alpha = 0.25, stroke = 0, aes(size = N_measured/2, color = colr))
PLOT1 = PLOT1 + scale_color_manual(values =  c('gray80','red') )
PLOT1 = PLOT1 + theme_bw()
#PLOT1 = PLOT1 + theme(legend.position="none")
#PLOT1 = PLOT1 + theme(legend.scale)
PLOT1 = PLOT1 + ylab('Mean Intra-individual Variant Allele Frequency')
PLOT1 = PLOT1 + xlab('Population Variant Allele Frequency')
PLOT1

dev.off()

# ------------------------------------------------------------------ 
# ------------------------------------------------------------------ 
# ------------------------------------------------------------------ 

# Stratify variants by their mean Intra-individual allele frequency
# the hypothesis was that there are two classes of variants;
# Directional Variants, which can range from 0 - 1
# Stabilized Variants, which are maintained in a narrower range

Stratify = function(iVAF, minVAF, maxVAF, noZero = TRUE){
  
  if (noZero){
    iVAF[which(iVAF == 0)] = NA #noZero
  }
  
  # calcualte the average intra-individual VAF for each variant
  mean_iVAF = apply(iVAF,1, mean, na.rm = T)

  # Subselect (index) mean iVAF between a range of values
  strat_lim = which(mean_iVAF > minVAF & mean_iVAF <= maxVAF)
  
  # Actually subselect by stratification
  iVAF_str = apply(iVAF[strat_lim,], 1, sort, decreasing = TRUE)
  
  N_strats = length(iVAF_str) # number of variants in this stratification
  N_samples = length(iVAF[1,]) # number of samples in the analysis (people)
  
  iVAF_str_matrix = matrix( rep(0, N_strats * N_samples), nrow = N_strats)
  
  for (X in 1:N_strats){
    LINE_VALUES = unlist(iVAF_str[X])
    length(LINE_VALUES) = N_samples
    LINE_VALUES[is.na(LINE_VALUES)] <- 0
    
    iVAF_str_matrix[X,] = iVAF_str_matrix[X,] + LINE_VALUES
  }
  
  return(iVAF_str_matrix)
  
}

Stratify_rownames = function(iVAF, minVAF, maxVAF, noZero = TRUE){
  
  if (noZero){
    iVAF[which(iVAF == 0)] = NA #noZero
  }
  
  # calcualte the average intra-individual VAF for each variant
  mean_iVAF = apply(iVAF,1, mean, na.rm = T)
  
  # Subselect (index) mean iVAF between a range of values
  strat_lim = which(mean_iVAF > minVAF & mean_iVAF <= maxVAF)
  
  return(names(strat_lim))
}

StratifyPlot = function(iVAF, minVAF, maxVAF){
  
  ## Stratify variants by intra-individual variant allele frequency
  iVAF_SUB_matrix = data.frame(Stratify(iVAF, minVAF, maxVAF), row.names = Stratify_rownames(iVAF, minVAF, maxVAF))
  
  # Shape data to be plotted
  iVAF_SUB_matrix = melt(iVAF_SUB_matrix)
  iVAF_SUB_matrix$rowid = Stratify_rownames(iVAF, minVAF, maxVAF)

    # Plot--------------------------------------------------------------
    # Mean Intra-individual Variant Allele Frequency vs.
    # Population-level Variant Allele Frequency (called vs not called)
    plotTitle = paste0("Mean VAF: ", minVAF," - ",maxVAF)
  
    PLOT2 = ggplot(iVAF_SUB_matrix, aes(variable, value, group = factor(rowid)))
    PLOT2 = PLOT2 + geom_line(aes(color = factor(rowid)))
    #PLOT2 = PLOT2 + geom_line(aes(color = )) COLOR BY STANDRD DEVIATION OF VARIANT ACROSS ALL MEAUSRED SAMPLES
    PLOT2 = PLOT2 + theme_bw()
    PLOT2 = PLOT2 + theme(legend.position="none")
    PLOT2 = PLOT2 + ylim(0,1)
    PLOT2 = PLOT2 + scale_x_discrete(breaks = c(1,N_pop))
    PLOT2 = PLOT2 + ggtitle(plotTitle)
    #PLOT2 = PLOT2 + ylab('Intra-individual Variant Allele Frequency')
    #PLOT2 = PLOT2 + xlab('Genomes')
    ## -----------------------------------------------------------------
    #PLOT2
    
    return(PLOT2)
}

StratifyBoxPlot = function(iVAF, minVAF, maxVAF){
  
  ## Stratify variants by intra-individual variant allele frequency
  iVAF_SUB_matrix = data.frame(Stratify(iVAF, minVAF, maxVAF), row.names = Stratify_rownames(iVAF, minVAF, maxVAF))
  
  # Shape data to be plotted
  iVAF_SUB_matrix = melt(iVAF_SUB_matrix)
  iVAF_SUB_matrix$rowid = Stratify_rownames(iVAF, minVAF, maxVAF)
  
  # Plot--------------------------------------------------------------
  # Mean Intra-individual Variant Allele Frequency vs.
  # Population-level Variant Allele Frequency (called vs not called)
  plotTitle = paste0("Mean VAF: ", minVAF," - ",maxVAF)
  
  PLOT2 = ggplot(iVAF_SUB_matrix, aes(variable, value, group = factor(rowid)))
  PLOT2 = PLOT2 + geom_boxplot(aes(color = factor(rowid)))
  #PLOT2 = PLOT2 + geom_line(aes(color = )) COLOR BY STANDRD DEVIATION OF VARIANT ACROSS ALL MEAUSRED SAMPLES
  PLOT2 = PLOT2 + theme_bw()
  PLOT2 = PLOT2 + theme(legend.position="none")
  PLOT2 = PLOT2 + ylim(0,1)
  PLOT2 = PLOT2 + scale_x_discrete(breaks = c(1,N_pop))
  PLOT2 = PLOT2 + ggtitle(plotTitle)
  #PLOT2 = PLOT2 + ylab('Intra-individual Variant Allele Frequency')
  #PLOT2 = PLOT2 + xlab('Genomes')
  ## -----------------------------------------------------------------
  #PLOT2
  
  return(PLOT2)
}


# Plot Stratify for multiple ranges
minVAF=0 # initialize with min VAF at 0
nStratification=3 # How many equal parts to divide the VAF from 0-1

# range width of each stratification
widthStratification=1/nStratification

for (nS in 1:nStratification){

  maxVAF = minVAF + widthStratification
  
  # Generate Plot
  PLOT_VAF = StratifyPlot(iVAF, minVAF, maxVAF)
  
  # Open PDF file to write to
  pdf(file = paste0(outPrefix,".",nS,".strat.pdf"),
     width = 5, height = 5)
  
    print(PLOT_VAF)
  
  dev.off()
  
  minVAF = minVAF + widthStratification
}

# Stratify plot for only CAV

StratifyPlot(iVAF[CAV,], 0.3, 0.9)
StratifyBoxPlot(iVAF[CAV,], 0.3, 0.9)



# Heatmap ================================================
# Log transform depth-data at each variant
logDP = log(DP)

# # chrom from tutotrial
heatmap.bp(iVAF)
heatmap.bp(iVAF, cbarplot = F, rbarplot = F, rlabels = F)

heatmap.bp(logDP, cbarplot = F, rbarplot = F, rlabels = F)

## Table of Variants ====================================

# Function to report rRNA variant stats from imported VCF

REGIONNAME= c("up_5S","5S","down_5S","Nmask_Alu","spacer1",
              "Nmask_SR","spacer2","RNA45S_promoter","RNA45S",
              "5ETS","18S","ITS1","5.8S","ITS2","28S","3ETS",
              "IGS_start")

REGIONSTART = c(10000,10219,10340,10376,10650,11590,11897,999000,
                1000000,1000000,1003660,1005529,1006622,1006779,
                1007947,1013018,1013408)

REGIONEND   = c(10219,10340,10376,10650,11590,11897,12241,1000000,
                1013408,1003660,1005529,1006622,1006779,1007947,
                1013018,1013408,1013558)
           

RDNA_LOCUS = data.frame(REGIONNAME, REGIONSTART, REGIONEND)





# Count per Genome ==================================================

# How much data is missing?------------------------------------------

NA_DP = apply(DP,1,is.na) # bin vector of present/missing

# How data is missing
NA_freq = length(which(NA_DP)) / length(NA_DP)

NA_perSite =  data.frame( colSums(NA_DP) / length(NA_DP[,1]) )
colnames(NA_perSite) = c('percentNA')


ggplot(as.data.frame(NA_perSite), aes(percentNA)) + geom_histogram() +
  theme_minimal() + 
  xlab('Intra-genomic VAF') + 
  ylab('Count')

# Number of variantions per individual ------------------------------
df_iVAF = data.frame(iVAF) # iVAF as data.frame

#length( which( ( colSums(NA_DP) / length(NA_DP[,1])) > .5 ) ) / 926


# Go through rDNA_LOCUS (apply) and return the count of how many variants are in a range

countVCF = function( RDNA_ROW, inVCF){
  # for RDNA_ROW from RDNA_LOCUS
  # how many variants of inVCF are within the coordinates
  # of this region
  variantID = rownames(inVCF)
  variantPOS = read.table(text = variantID, header = T, sep = "_")[,2] # position of variant
  
  CNT = which(variantPOS >= as.integer(RDNA_ROW[2]) &
              variantPOS <  as.integer(RDNA_ROW[3]) )
  
  return(length(CNT))
}

# Which variants are in each of the CEPH-1436 libraries
  var_78 = which(df_iVAF$NA12878 > 0.02 & df_iVAF$NA12878 < 0.998)
  var_91 = which(df_iVAF$NA12891 > 0.02 & df_iVAF$NA12891 < 0.998)
  var_92 = which(df_iVAF$NA12892 > 0.02 & df_iVAF$NA12892 < 0.998)
  
  var_trio = sort(unique(c(var_78, var_91, var_92)))

# Count per region in the Utah Trio
VariantsPerRegion = 
data.frame(as.character(RDNA_LOCUS$REGIONNAME),
  apply(RDNA_LOCUS, 1, countVCF, inVCF = df_iVAF[var_78,]),
  apply(RDNA_LOCUS, 1, countVCF, inVCF = df_iVAF[var_91,]),
  apply(RDNA_LOCUS, 1, countVCF, inVCF = df_iVAF[var_92,]))
  colnames(VariantsPerRegion) = c('Region','78','91','92')
  

VariantsPerRegion[1,2:4] = 
  VariantsPerRegion[1,2:4] +
  colSums(VariantsPerRegion[3:7, 2:4])

VariantsPerRegion = 
  VariantsPerRegion[-c(3,4,5,6,7,9,17),]

VariantsPerRegion$mean = apply(VariantsPerRegion[,2:4],1,mean)
VariantsPerRegion$sd = apply(VariantsPerRegion[,2:4],1,sd)

VariantsPerRegion$forceSORT = c('A','B','C','D','E','F','G','H','I','J')
VariantsPerRegion$Region = as.character(VariantsPerRegion$Region)

pdf('VariantsPerRegion.pdf',width = 5, height = 5)


PLOT3 = ggplot(VariantsPerRegion, aes(x = forceSORT, mean))
PLOT3 = PLOT3 + geom_bar(stat = 'identity')
PLOT3 = PLOT3 + geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), width = 0)
PLOT3 = PLOT3 + scale_x_discrete(labels=VariantsPerRegion$Region)
PLOT3 = PLOT3 + theme_bw()
PLOT3

dev.off()

PLOT1 = PLOT1 + geom_point(alpha = 0.5, stroke = 0, aes(size = pV/2, color = colr))
PLOT1 = PLOT1 + scale_color_manual(values =  c('gray80','gray80','red') )
PLOT1 = PLOT1 + theme_bw()
#PLOT1 = PLOT1 + theme(legend.position="none")
PLOT1 = PLOT1 + theme(legend.scale)
PLOT1 = PLOT1 + ylab('Mean Intra-individual Variant Allele Frequency')
PLOT1 = PLOT1 + xlab('Population Variant Allele Frequency')
PLOT1


# # Which variants are less than 2% minor allele frequency within the genome
# A = length(df_iVAF$NA12878[-which(df_iVAF$NA12878 < 0.02 | df_iVAF$NA12878 > 0.98)])
# B = length(df_iVAF$NA12891[-which(df_iVAF$NA12891 < 0.02 | df_iVAF$NA12878 > 0.98)])
# C = length(df_iVAF$NA12892[-which(df_iVAF$NA12892 < 0.02 | df_iVAF$NA12878 > 0.98)])
# 
# # Which variants are less than 10% minor allele frequency
# A1 = length(df_iVAF$NA12878[-which(df_iVAF$NA12878 < 0.02 | df_iVAF$NA12878 > 0.98)])
# B1 = length(df_iVAF$NA12891[-which(df_iVAF$NA12891 < 0.02 | df_iVAF$NA12878 > 0.98)])
# C1 = length(df_iVAF$NA12892[-which(df_iVAF$NA12892 < 0.02 | df_iVAF$NA12878 > 0.98)])
# 
# # Fraction of variants in sample which are between 2-10%
# mA = (A-A1)/A
# mB = (B-B1)/B
# mC = (C-C1)/C


## CDF Plots per Individual ----------------------------------------

# Extract iVAF from TRIO, correct for consensus (i.e. Minor Allele Frequency)
NA12878_iVAF = df_iVAF$NA12878
  #NA12878_iVAF[which(NA12878_iVAF > 0.5)] = 1 - NA12878_iVAF[which(NA12878_iVAF > 0.5)]

NA12891_iVAF = df_iVAF$NA12891
  #NA12891_iVAF[which(NA12891_iVAF > 0.5)] = 1 - NA12891_iVAF[which(NA12891_iVAF > 0.5)]
  
NA12892_iVAF = df_iVAF$NA12892
  #NA12892_iVAF[which(NA12892_iVAF > 0.5)] = 1 - NA12892_iVAF[which(NA12892_iVAF > 0.5)]
  
TRIO_iVAF = data.frame(NA12878_iVAF, NA12891_iVAF, NA12892_iVAF)
TRIO_iVAF = melt(TRIO_iVAF)

# PLOT4 = ggplot(TRIO_iVAF, aes(x = value))
# PLOT4 = PLOT4 + geom_histogram(data=subset(TRIO_iVAF,variable == 'NA12878_iVAF'),fill = "blue", alpha = 0.2, bins = 10)
# PLOT4 = PLOT4 + geom_histogram(data=subset(TRIO_iVAF,variable == 'NA12891_iVAF'),fill = "green", alpha = 0.2, bins = 10)
# PLOT4 = PLOT4 + geom_histogram(data=subset(TRIO_iVAF,variable == 'NA12892_iVAF'),fill = "red", alpha = 0.2, bins = 10)
# PLOT4 = PLOT4 + theme_bw() + theme(legend.position="none")
# #PLOT4 = PLOT4 + ylim(c(0,150))
# #PLOT4 = PLOT4 + geom_vline(xintercept = 0.02)
# PLOT4

pdf('CEPH_TRIO_iVAF_CDF_2.pdf', height = 5, width = 5)
PLOT4 = ggplot(TRIO_iVAF, aes(value, color = variable )) + stat_ecdf(geom = "step")
PLOT4 = PLOT4 + theme_bw() + theme(legend.position="none")
PLOT4 = PLOT4 + geom_vline(xintercept = 0.02)
PLOT4
dev.off()

NA12891_VCF = df_iVAF[which(df_iVAF$NA12891 > 0.02 & df_iVAF$NA12891 < 0.998),]
NA12892_VCF = df_iVAF[which(df_iVAF$NA12892 > 0.02 & df_iVAF$NA12892 < 0.998),]


## Evolutionary Conservation fo Variants ---------------------------
ENTROPY = read.table('hgr1_entropy.bed', header = F, sep = '\t')
  colnames(ENTROPY) = c('chr','start','end','S')
  ENTROPY$info = 2 - ENTROPY$S

df_trio = df_iVAF[var_trio, c('NA12878','NA12891','NA12892')]  
  
TrioVariants = rownames(df_trio)
TrioVariants = read.table(text = TrioVariants, header = F, sep = "_")[,2] # position of variant

TrioENTROPY = ENTROPY[which(ENTROPY$end %in% TrioVariants),]
trio_mean_iVAF = apply(df_trio[which(TrioVariants %in% ENTROPY$end),],1 , mean)
TrioENTROPY$mean_iVAF = trio_mean_iVAF

TrioENTROPY = TrioENTROPY[-which(TrioENTROPY$start == 1013018),]
# Removing the last position variants. These are insertions at the last position
# of rRNA, I can't define them as being inside 28S since I don't know if the mature
# 28S molecule retains this 'long tail' or not.
# worth looking into in the RNAseq datas.

pdf('CEPH_TRIO_Evolution.pdf', height = 5, width = 5)
  PLOT5 = ggplot(TrioENTROPY, aes(info, mean_iVAF))
  PLOT5 = PLOT5 + geom_point(alpha = 0.25, stroke = 0, size = 4)
  PLOT5 = PLOT5 + scale_color_manual(values =  c('gray80') )
  PLOT5 = PLOT5 + theme_bw() + theme(legend.position="none")
  PLOT5
dev.off()


pdf(file = paste0(outPrefix,".MeaniVAF_pVAF_Nm.pdf"), width = 6.5, height = 5)

PLOT1 = ggplot(POPVAR, aes(pVAF, mean_noZero_iVAF ))
#PLOT1 = PLOT1 + geom_point(alpha = 0.5, stroke = 0, aes(size = pV/2, color = colr))
PLOT1 = PLOT1 + geom_point(alpha = 0.25, stroke = 0, aes(size = N_measured/2, color = colr))
PLOT1 = PLOT1 + scale_color_manual(values =  c('gray80','red') )
PLOT1 = PLOT1 + theme_bw()
#PLOT1 = PLOT1 + theme(legend.position="none")
PLOT1 = PLOT1 + ylab('Mean Intra-individual Variant Allele Frequency')
PLOT1 = PLOT1 + xlab('Population Variant Allele Frequency')


# Table 1 Variants --------------------------------------------------

varSelection = which(0.33 < POPVAR$mean_noZero_iVAF & POPVAR$mean_noZero_iVAF < 0.66 & pVAF > 0.5)

TABLE1 = data.frame(names(varSelection))

  TABLE1 = cbind( TABLE1, POPVAR$mean_noZero_iVAF[varSelection],
                  POPVAR$sd_noZero_iVAF[varSelection],
                  POPVAR$pVAF[varSelection],
                  POPVAR$pV[varSelection])
  colnames(TABLE1) = c('Coordinate', 'iVAF','sd_iVAF','pVAF','pV')

# Supplementary Table 1 ---------------------------------------------
  
  write.table(POPVAR,file = 'supplementary_Table_1.csv', quote = F, sep = "\t")

  