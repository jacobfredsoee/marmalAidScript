probePlot = function(probes, betaData, betaDataSets, betaCut = NULL, plotType = "boxplot") {
  require(ggplot2)
  require(reshape2)
  
  probes = as.character(probes)

  setNames = as.vector(betaDataSets$setNames)
  columnNames = as.vector(betaDataSets$columnNames)
  
  betaCut = betaCut + 1
  
  subdata = betaData[probes,]
  colnames(subdata) = columnNames
  
  percentageNames = rev(rev(columnNames)[1:11])
  
  if(plotType == "boxplot") {
    valuePicker = which(columnNames == percentageNames[betaCut])
    
    graphData = subdata[,valuePicker]
    colnames(graphData) = unique(setNames)[-1]
    
    
    g = ggplot(melt(graphData), aes(x = as.factor(variable), y = value)) +
      geom_boxplot() + theme_bw() + geom_jitter(width = 0.1) + ggtitle(percentageNames[betaCut]) +
      xlab("Dataset") + ylab("Percentage")
  }
  if(plotType == "lineplot") {
    
    valuePicker = which(columnNames %in% percentageNames)
    
    graphData = subdata[,valuePicker]
    colnames(graphData) = as.vector(sapply(unique(setNames)[-1], function(setName) rep(setName, length(percentageNames))))
    
    graphData = melt(t(graphData))
    colnames(graphData) = c("Dataset", "ProbeID", "Value")
    graphData = cbind(graphData, x = rep(1:length(percentageNames), length(unique(graphData$ProbeID))))
    
    gplots = list()
    length(gplots) = length(unique(graphData$ProbeID))
    
    for(i in 1:length(unique(graphData$ProbeID))) {
      g = ggplot(subset(graphData, ProbeID == unique(graphData$ProbeID)[i]), aes(x = x, y = Value, color = Dataset, shape = Dataset)) + 
        geom_line() + 
        geom_point() +
        scale_x_discrete(limits=percentageNames) + xlab("") + ylab("Percentage") +
        ggtitle(as.character(unique(graphData$ProbeID)[i]))
      gplots[[i]] = g + theme_bw()
    }
    
    g = multiplot(plotlist = gplots, cols = sqrt(length(unique(graphData$ProbeID))))
  }
  return(g)
}
