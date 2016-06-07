coordPlot = function(betaData, betaDataSets, startCoord, stopCoord, chr, betaCut = NULL, plotType = "boxplot") {
  require(ggplot2)
  require(reshape2)

  #extract setnames and columnNames
  setNames = as.vector(betaDataSets$setNames)
  columnNames = as.vector(betaDataSets$columnNames)
  
  #extract only data from the selected probes
  subdata = subset(betaData, CHR == chr & MAPINFO > startCoord & MAPINFO < stopCoord)
  colnames(subdata) = columnNames
  
  #get the column names with the percentages
  percentageNames = rev(rev(columnNames)[1:11])
  
  if(plotType == "boxplot") {
    if(is.null(betaCut)) {
      print("betaCut not set for boxplot; picking for you")
      betaCut = 6
    } else if(betaCut < 0 | betaCut > 10) {
      print("betaCut outside of range for boxplot; setting it for you")
      betaCut = 6
    } else {
      betaCut = betaCut + 1
    }
    
    valuePicker = which(columnNames == percentageNames[betaCut])
    graphData = subdata[,valuePicker]
    colnames(graphData) = unique(setNames)[-1]
    
    graphData = melt(t(graphData))
    colnames(graphData) = c("Dataset", "ProbeID", "Value")
    graphData = cbind(graphData, x = subdata[as.character(graphData$ProbeID), "MAPINFO"])
    print(graphData)
    g = ggplot(graphData, aes(x = as.factor(x), y = Value)) + 
      geom_boxplot() + theme_bw() + geom_jitter(width = 0.1) + ggtitle(percentageNames[betaCut]) +
      xlab("Coordinates") + ylab("Percentage")
  }

  
  
  return(g)
}