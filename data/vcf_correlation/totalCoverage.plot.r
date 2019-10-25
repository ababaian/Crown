# readCoverage.plot.r
#
# Read Coverage over an area

library(ggplot2)
library(reshape2)


bedFile = 'data/na12892_cov.w1s1.bed'
OUTPUTNAME = 'NA12892'

BEDCOVERAGE = read.csv(bedFile, header = F, sep = '\t')
  colnames(BEDCOVERAGE) = c('chrom','start','end','coverage')


pdf( paste0(OUTPUTNAME,'.5S_coverage.pdf'), height = 3, width = 10 )
PLOT = ggplot(BEDCOVERAGE, aes(start,coverage))
PLOT = PLOT + geom_line()
PLOT = PLOT + theme_minimal() + scale_x_continuous( limits = c(10000,12500))
PLOT = PLOT + scale_y_log10(limits = c(1,100000))
print(PLOT)
dev.off()

pdf( paste0(OUTPUTNAME,'.45S_coverage.pdf'), height = 3, width = 10)
PLOT = ggplot(BEDCOVERAGE, aes(start,coverage))
PLOT = PLOT + geom_line()
PLOT = PLOT + theme_minimal() + scale_x_continuous( limits = c(998500,1015000))
PLOT = PLOT + scale_y_log10(limits = c(1,100000))
print(PLOT)
dev.off()

