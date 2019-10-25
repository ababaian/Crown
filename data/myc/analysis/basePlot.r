# basePlot.r
#
# Analysis and plot the data
# for a single base (after running crcAnalysis.r)
#

library(ggplot2)
library(reshape2)

# Position of interest
# in hgr1 / chr13 coordinates 
# 28S.59G / hgr1:1008007
# 28S.470A / hgr1:1008418
# 28S.4532U / hgr1:1012480
# 18S.1248U / hgr1:1004908

pos ='1012480'
  pos = which(genCoord == pos)

# paired analysis ----------------------------------
    
# Change in allele frequency between sampleA-sampleB
#POS = data.frame(dRAF[,pos]) # paired
# # Raw reference allele frequency
# POS$sampleA_RAF = RAF_A[,pos]
# POS$sampleB_RAF = RAF_B[,pos]
# 
# # Depth of coverage at position of interest
# POS$sampleA_DP = DP[idxA, pos]
# POS$sampleB_DP = DP[idxB, pos]
# 
# # Plot RAF_A vs. RAF_B
# POSDATA = melt(POS[,2:3])
  
# unpaired analysis -------------------------------
# Raw reference allele frequency
POSDATA = data.frame( c(RAF_A[,pos], RAF_B[,pos] ))
  colnames(POSDATA) = 'value'
  
POSDATA$variable = c( rep('sampleA', length(RAF_A[,pos]) ),
                      rep('sampleB', length(RAF_B[,pos]) ))

# # Depth of coverage at position of interest
# POS$sampleA_DP = DP[idxA, pos]
# 
# POS$sampleB_DP = DP[idxB, pos]

PLOT = ggplot(POSDATA, aes(variable, value)) +
  geom_boxplot(stat = 'boxplot') + 
  geom_jitter( width = 0.2)
PLOT

# Test for signifiance
t.test(RAF_A[,pos], RAF_B[,pos], paired = FALSE)
var.test(RAF_A[,pos], RAF_B[,pos])


# Plot change in reference allele frequency
# comapred to a sampleB distribution with the same
# standard deviation (null hypothesis)

# POS$sim_dRAF = rnorm(sampleN/2, mean = 0,
#                      sd = sd(POS$dRAF))
# 
# POSDATA = melt(POS[,c(1,6)])

# Plot dRAF alone
# PLOT = ggplot(POS, aes('delta 28S.r.59G', dRAF)) +
#   geom_boxplot(stat = 'boxplot') +
#   geom_jitter(width = 0.2)
# PLOT

# Plot dRAF vs. Normal Distribution
# PLOT = ggplot(POSDATA, aes(variable, value)) +
#   geom_boxplot(stat = 'boxplot') + 
#   geom_jitter( width = 0.2)
# PLOT
# 
# t.test(POS$dRAF, POS$sim_dRAF, paired = FALSE)
# var.test(POS$dRAF, POS$sim_dRAF)
