coordPlot = function(betaData, betaDataSets, startCoord, stopCoord, chr, betaCut = NULL, plotType = "boxplot", showGeneNames = FALSE) {
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
    graphData = cbind(graphData, x = as.factor(subdata[as.character(graphData$ProbeID), "MAPINFO"]))
    graphData = cbind(graphData, gene = subdata[as.character(graphData$ProbeID), "GENE_BY_JBBR"])
    
    g = ggplot(graphData, aes(x = x, y = Value)) + 
      geom_boxplot() + theme_bw() + geom_jitter(width = 0.1) + ggtitle(percentageNames[betaCut]) +
      xlab("Coordinates") + ylab("Percentage") + ylim(c(0,1))
    
    if(showGeneNames) {
      geneInfo = unique(graphData[,c("x", "gene")])
      geneInfo = geneInfo[order(geneInfo$x),]
      geneInfo$x = as.numeric(geneInfo$x)

      for(geneName in unique(geneInfo$gene)) {
        if(geneName != "") {
          singleGeneInfo = subset(geneInfo, gene == geneName)
          g = g + annotate("rect", xmin = min(singleGeneInfo$x), xmax = max(singleGeneInfo$x), ymin = 0, ymax = 0.05, alpha = 0.2) +
            annotate(geom = "text", x = (min(singleGeneInfo$x) + max(singleGeneInfo$x)) / 2, y = 0.025, label = unique(singleGeneInfo$gene, hjust = 0.5, vjust = 0.5))
        }
      }
    }
  } else if(plotType == "lineplot") {
    
    #select the columns to use
    valuePicker = which(columnNames %in% percentageNames)
    graphData = subdata[,valuePicker]
    colnames(graphData) = as.vector(sapply(unique(setNames)[-1], function(setName) rep(setName, length(percentageNames))))
    
    #create data for graph
    graphData = melt(t(graphData))
    colnames(graphData) = c("Dataset", "ProbeID", "Value")
    graphData = cbind(graphData, x = subdata[as.character(graphData$ProbeID), "MAPINFO"])
    graphData = cbind(graphData, gene = subdata[as.character(graphData$ProbeID), "GENE_BY_JBBR"])
    graphData = cbind(graphData, betaValue = as.factor(rep(1:length(percentageNames), length(unique(graphData$ProbeID)) * length(unique(graphData$Dataset)))))
    
    if(!is.null(betaCut)) {
      betaCut = betaCut + 1
      graphData = subset(graphData, betaValue %in% betaCut)
    }
    
    if(is.null(betaCut) | length(betaCut) != 1) {
      #create the list to hold the plots
      gplots = list()
      length(gplots) = length(unique(graphData$Dataset))
      
      for(i in 1:length(gplots)) {
        subGraphData = subset(graphData, Dataset == unique(graphData$Dataset)[i] & betaValue != 1 & betaValue != 11)
        
        g = ggplot(subGraphData, aes(x = x, y = Value, color = betaValue, shape = betaValue)) + 
          geom_line(size = 1) + scale_color_manual(values = colorRampPalette(c("lightgrey", "firebrick"))(length(unique(subGraphData$betaValue))),
                                                   breaks=c(2:10),
                                                   labels=percentageNames[2:10]) +
          geom_point(size = 3) + scale_shape_manual(values=1:length(unique(subGraphData$betaValue)),
                                                    breaks=c(2:10),
                                                    labels=percentageNames[2:10]) +
          xlab("Coordinates") + ylab("Percentage") +
          ggtitle(unique(graphData$Dataset)[i]) + theme_bw() + ylim(c(0,1))
        
        if(showGeneNames) {
          geneInfo = unique(graphData[,c("x", "gene")])
          geneInfo = geneInfo[order(geneInfo$x),]
          
          for(geneName in unique(geneInfo$gene)) {
            if(geneName != "") {
              singleGeneInfo = subset(geneInfo, gene == geneName)
              g = g + annotate("rect", xmin = min(singleGeneInfo$x), xmax = max(singleGeneInfo$x), ymin = 0, ymax = 0.05, alpha = 0.2) +
                annotate(geom = "text", x = (min(singleGeneInfo$x) + max(singleGeneInfo$x)) / 2, y = 0.025, label = unique(singleGeneInfo$gene, hjust = 0.5, vjust = 0.5))
            }
          }
        }
        
        gplots[[i]] = g
      }
      
      g = multiplot(plotlist = gplots, cols = sqrt(length(unique(graphData$Dataset))))
      
    } else {
      g = ggplot(graphData, aes(x = x, y = Value, color = Dataset, shape = Dataset)) +
        geom_line(size = 1) + scale_color_manual(values = colorRampPalette(c("lightgrey", "darkblue"))(length(unique(graphData$Dataset)))) +
        geom_point(size = 3) + scale_shape_manual(values=1:length(unique(graphData$Dataset))) +
        xlab("Coordinates") + ylab("Percentage") +
        ggtitle(percentageNames[betaCut]) + theme_bw() + ylim(c(0,1))
      
      if(showGeneNames) {
        geneInfo = unique(graphData[,c("x", "gene")])
        geneInfo = geneInfo[order(geneInfo$x),]
        
        for(geneName in unique(geneInfo$gene)) {
          if(geneName != "") {
            singleGeneInfo = subset(geneInfo, gene == geneName)
            g = g + annotate("rect", xmin = min(singleGeneInfo$x), xmax = max(singleGeneInfo$x), ymin = 0, ymax = 0.05, alpha = 0.2) +
              annotate(geom = "text", x = (min(singleGeneInfo$x) + max(singleGeneInfo$x)) / 2, y = 0.025, label = unique(singleGeneInfo$gene, hjust = 0.5, vjust = 0.5))
          }
        }
      }
        
    }

  }
  return(g)
}