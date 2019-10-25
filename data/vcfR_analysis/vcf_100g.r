# vcfR analysis
# 100 genomes hgr1
# 170308
#
# hgr1_v1 vcf (variant input only) analysis

# Install vcfR
#install.packages("vcfR")
library("vcfR")
library("ggplot2")
library("reshape2")

# File pointers

  vcf_file = '100g.hgr1_all.vcf'
  dna_file = 'rDNA.fa'
  gff_file = 'rDNA.gff'

# output Prefix (for plots)
  
  outPrefix='100g_gvcf'

# SCRIPT ============================================================
#====================================================================
# Import VCF / DNA / GFF
  VCF = read.vcfR(vcf_file)
  DNA = ape::read.dna(dna_file, format = 'fasta')
  GFF = read.table(gff_file, sep="\t", quote = "")
  
  
# Create chromR Object
  chrom = create.chromR(name="rDNA", vcf=VCF, seq=DNA, ann=GFF)

  #plot(chrom)
 #chromoqc(chrom, dp.alpha = 22)


# Variant Statistics ================================================

# Intra-individual variation (_i)
  
  # Depth of coverage at each variant
  DP = extract.gt(chrom, element="DP", as.numeric=TRUE)

  # Reference allele-only depth
  RAD = extract.gt(chrom, element="AD", as.numeric=TRUE)

  # Variant allele frequency (intra-individual)
  iVAF = (DP-RAD)/DP
  
  
# Population variation (called variants)
  
  # Number of individual genomes in total population
  N_pop = length(DP[1,])
  
  # Called Variant allele count (population)
  pV = rowSums(!is.na(iVAF))
  pVAF = pV / N_pop
  
  # Average intra-genomic variant allele frequency
  mean_iVAF = apply(iVAF,1,mean, na.rm = T)
  sd_iVAF = apply(iVAF,1,sd, na.rm = T)

# Population Variant Statistics
  
POPVAR = data.frame(pVAF, mean_iVAF, sd_iVAF)

# Plot--------------------------------------------------------------
# Mean Intra-individual Variant Allele Frequency vs.
# Population-level Variant Allele Frequency (called vs not called)

# open PDF device
pdf(file = paste0(outPrefix,".MeaniVAF_pVAF.pdf"), width = 5, height = 5)

PLOT1 = ggplot(POPVAR, aes(mean_iVAF, pVAF))
PLOT1 = PLOT1 + geom_point(alpha = 1.5/10, stroke = 0, aes(size = 0.5))
PLOT1 = PLOT1 + theme_bw()
PLOT1 = PLOT1 + theme(legend.position="none")
PLOT1 = PLOT1 + xlab('Mean Intra-individual Variant Allele Frequency')
PLOT1 = PLOT1 + ylab('Population Variant Allele Frequency')
PLOT1

dev.off()
# ------------------------------------------------------------------ 
  
# Stratify variants by their mean Intra-individual allele frequency
# the hypothesis was that there are two classes of variants;
# Directional Variants, which can range from 0 - 1
# Stabilized Variants, which are maintained in a narrower range

# strat0_2 =  which(mean_iVAF <= 0.2)
# strat2_4 =  which(mean_iVAF > 0.2 & mean_iVAF <= 0.4)
# strat4_6 =  which(mean_iVAF > 0.4 & mean_iVAF <= 0.6)
# strat6_8 =  which(mean_iVAF > 0.6 & mean_iVAF <= 0.8)
# strat8_1 =  which(mean_iVAF > 0.8)


Stratify = function(iVAF, minVAF, maxVAF){
  
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

Stratify_rownames = function(iVAF, minVAF, maxVAF){
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
    plotTitle = paste0("Mean VAF: ",minVAF," - ",maxVAF)
  
    PLOT2 = ggplot(iVAF_SUB_matrix, aes(variable, value, group = factor(rowid)))
    PLOT2 = PLOT2 + geom_line(aes(color = factor(rowid)))
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
nStratification=5 # How many equal parts to divide the VAF from 0-1

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


# Heatmap ================================================
# Log transform depth-data at each variant
logDP = log(DP)

# # chrom from tutotrial

heatmap.bp(iVAF, cbarplot = F, rbarplot = F, rlabels = F)

heatmap.bp(logDP, cbarplot = F, rbarplot = F, rlabels = F)

#   







  ## From tutorial