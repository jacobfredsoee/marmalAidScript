coordPlot = function(betaData, betaDataSets, startCoord, stopCoord, chr, betaCut = NULL, plotType = "boxplot") {
  require(ggplot2)
  require(reshape2)

  startCoord = 1
  stopCoord = 100000
  chr = 1
  
  
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
    
    g = ggplot(graphData, aes(x = as.factor(x), y = Value)) + 
      geom_boxplot() + theme_bw() + geom_jitter(width = 0.1) + ggtitle(percentageNames[betaCut]) +
      xlab("Coordinates") + ylab("Percentage") + ylim(c(0,1))
  } else if(plotType == "lineplot") {
    
    #select the columns to use
    valuePicker = which(columnNames %in% percentageNames)
    graphData = subdata[,valuePicker]
    colnames(graphData) = as.vector(sapply(unique(setNames)[-1], function(setName) rep(setName, length(percentageNames))))
    
    #create data for graph
    graphData = melt(t(graphData))
    colnames(graphData) = c("Dataset", "ProbeID", "Value")
    graphData = cbind(graphData, x = subdata[as.character(graphData$ProbeID), "MAPINFO"])
    graphData = cbind(graphData, betaValue = as.factor(rep(1:length(percentageNames), length(unique(graphData$ProbeID)) * length(unique(graphData$Dataset)))))
    
    #create the list to hold the plots
    gplots = list()
    length(gplots) = length(unique(graphData$Dataset))
    
    for(i in 1:length(gplots)) {
      subGraphData = subset(graphData, Dataset == unique(graphData$Dataset)[i] & betaValue != 1 & betaValue != 11)
      
      g = ggplot(subGraphData, aes(x = x, y = Value, color = betaValue, shape = betaValue)) + 
        geom_line(size = 1) + scale_color_manual(values = colorRampPalette(c("lightgrey", "firebrick"))(nlevels(subGraphData$betaValue)),
                                                 breaks=c(2:10),
                                                 labels=percentageNames[2:10]) +
        geom_point(size = 3) + scale_shape_manual(values=1:nlevels(subGraphData$betaValue),
                                                  breaks=c(2:10),
                                                  labels=percentageNames[2:10]) +
        xlab("Coordinates") + ylab("Percentage") +
        ggtitle(unique(graphData$Dataset)[i]) + theme_bw()
      
      gplots[[i]] = g
    }
    
    g = multiplot(plotlist = gplots, cols = sqrt(length(unique(graphData$ProbeID))))
  }

  
  
  return(g)
}