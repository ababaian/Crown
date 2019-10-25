# readCoverage.plot.r
#
# Read Coverage over an area

library(ggplot2)
library(reshape2)

# ATGC Colors
ATGC=c('#64f73f','#eb413c','#ffb340','#3c88ee','magenta','gray50','black')

input = 'data/NA12878.28s.59.readCoverage'

DATA = read.csv(input, header = TRUE, sep = "\t")
  rownames(DATA) = DATA[,1]
  DATA = DATA[,-1]
  DATA = data.frame(DATA)
  DATA = t(DATA)
  
  posDATA = melt(DATA[ ,1:7])
  posDATA$value = as.numeric(posDATA$value)
  
  
  negDATA = melt(DATA[ ,8:12])
  
  plotTitle = "28S.59G>A Variant Coverage"
  

pdf('NA12878_pos.pdf',width = 10, height = 3)    
  PLOT = ggplot(posDATA, aes(Var1,value, group = factor(Var1)))
  PLOT = PLOT + geom_bar(stat = "identity", aes(fill = factor(Var2)))
  PLOT = PLOT + scale_fill_manual(values = ATGC)
  PLOT = PLOT + theme_minimal() + scale_y_continuous( limits = c(0,8000))
  PLOT = PLOT + theme(legend.position="none")
  print(PLOT)
dev.off()
  

pdf('NA12878_neg.pdf',width = 10, height = 3)  
  PLOT2 = ggplot(negDATA, aes(Var1,value, group = factor(Var1)))
  PLOT2 = PLOT2 + geom_bar(stat = "identity", aes(fill = factor(Var2)))
  PLOT2 = PLOT2 + scale_fill_manual(values = ATGC[1:5])
  PLOT2 = PLOT2 + theme_minimal() + scale_y_continuous( limits = c(0,8000))
  PLOT2 = PLOT2 + theme(legend.position="none")
  print(PLOT2) 
dev.off()