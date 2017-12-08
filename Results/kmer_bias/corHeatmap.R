
#################################################
# Import
#################################################

library(ggplot2)
library(reshape)
library(gplots)
library(gtools)
library(RColorBrewer)

#################################################
# Parameters
#################################################

cellData = read.table('/home/egg/eg474423_Projects/trunk/Dream/exp/kmer_bias/bias/CellFull.txt', sep='\t', header=F)
cellTranslate = as.character(cellData[,2])
names(cellTranslate) = as.character(cellData[,1])
colorTranslate = as.character(cellData[,3])
names(colorTranslate) = as.character(cellData[,1])
myDist = function(p1) dist(p1, method="euclidean")
myHclust = function(p2) hclust(p2, method="ward.D")
source('/home/egg/eg474423_Projects/trunk/TfbsPrediction/Results/LabelTranslation/newHeatmap2S.R')

#################################################
# Function
#################################################

createCorHeatmapDendrogram <- function(inFileName,outFileName){

  # Fetching data
  hmbreaks = c(seq(0.7,0.85,length=10),seq(0.85,1.0,length=10))
  hmcol = colorRampPalette(c("black", "green"))(n = 19)
  data = read.table(inFileName, sep='\t', header=T, row.names = 1)
  data = as.matrix(data)
  rscol = c()
  cscol = c()

  for (i in 1:length(colnames(data))) {
    cscol = c(cscol,colorTranslate[colnames(data)[i]])
    colnames(data)[i] = cellTranslate[colnames(data)[i]]
  }
  for (i in 1:length(rownames(data))) {
    rscol = c(rscol,colorTranslate[rownames(data)[i]])
    rownames(data)[i] = cellTranslate[rownames(data)[i]]
  }

  #lmat=rbind(c(5,4,0), c(3,2,1))
  #lhei=c(0.5, 4.5)
  #lwid=c(1,4,0.25)
  #lwid=c(1.5,3.5,0.25)

  # Plot
  postscript(outFileName,width=8.0,height=7.0,horizontal=FALSE,paper='special')
  par(cex.axis=1.0, cex.main=1.0, mar=c(5,1,1,5))
  heatmap.2(data, col=hmcol, breaks=hmbreaks, main="", cexCol=1.0, cexRow=1.0, trace='none',
            sepwidth=c(2,2), sepcolor='black', Rowv=TRUE, Colv=TRUE, density.info = 'none',
            #lmat = lmat, lhei = lhei, lwid = lwid,
            distfun=myDist, hclustfun=myHclust, keysize = 1)
  dev.off()
  system(paste('epstopdf',outFileName,sep=' '))

}

#################################################
# Execution
#################################################

# Input
loc = "/home/egg/eg474423_Projects/trunk/Dream/exp/kmer_bias/correlation/"
inList = c("spearman_All_corr")

# Strand Loop
for (i in 1:length(inList)) {

  # Parameters
  inFileName = paste(loc,inList[i],".txt",sep="")
  outFileName = paste(loc,inList[i],".eps",sep="")

  # Creating heatmap
  createCorHeatmapDendrogram(inFileName,outFileName)

}


