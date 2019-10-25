# screenPlot.r
#
# Screen for changes between allele frequency
# between two biological samples
# using adCalc.sh / rrnaAnalysis.r
#

library(ggplot2)
library(rgl)
library(grid)

# Calculate P-value for difference of means in RAF
# between cancer and normal
# use as a 'score' 

# t.test based
Pval = t.test(canRAF[,1], normRAF[,1], paired = TRUE)$p.value

for (i in 2:bpN){
  Pval = cbind(Pval,
               t.test(canRAF[,i], normRAF[,i], paired = TRUE)$p.value)
}

# Score
Pscore = -log(Pval)
Pscore[is.na(Pscore)] = 0
Pscore = as.numeric(Pscore)


# Bonferonni correction
# p < alpha / m
# p < signifiance_cutoff / numberTests
# p * numberTests < significance_cutoff

Pscore_bon = -log(Pval*bpN)
Pscore_bon[is.na(Pscore_bon)] = 0
Pscore_bon[Pscore_bon < 0] = 0
Pscore_bon = as.numeric(Pscore_bon)

# NOTE: this should be re-calculated as a Manhatten plot
# for the publication. T-test is a bit basic.


# PLOT ==================================================

DATA = data.frame(1:bpN) # for 18S
  colnames(DATA) = c('RNA')

DATA$mean_dRAF = mean_dRAF    
DATA$var_dRAF = var_dRAF
DATA$mean_DP = mean_DP
DATA$Pscore = Pscore
DATA$Pscore_bon = Pscore_bon

# plot(mean_dRAF)
#  plot(var_dRAF)
#  plot(mean_DP)
#  plot(log(mean_DP))
#  plot(Pscore)
#  plot(Pscore_bon)
# plot3d(mean_DP, Pscore)
 
# PLOTS for 18S
PLOT1 = ggplot(DATA, aes(RNA, mean_dRAF)) +
  geom_point() + theme_minimal() +
  ylim(c(-0.10,0.10))
#PLOT1

# PLOT1A = ggplot(DATA, aes(RNA, mean_dRAF)) +
#   geom_point() + theme_minimal() +
#   ylim(c(-0.12,0.12))

PLOT2 = ggplot(DATA, aes(RNA, var_dRAF)) +
  geom_point() + theme_minimal()
#PLOT2

PLOT3 = ggplot(DATA, aes(RNA, mean_DP)) +
  geom_point() + theme_minimal() +
  scale_y_log10()
#PLOT3

PLOT4 = ggplot(DATA, aes(RNA, Pscore)) +
  geom_point() + theme_minimal() +
  geom_hline(yintercept = -log(0.001/bpN))
#PLOT4

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(PLOT1, PLOT2, PLOT3, PLOT4, cols=1)
