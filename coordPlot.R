## coordPlot
## function to create plots from a given set of coordinates and a corresponding beta summary set (created 
## from extractSelectedGroups.R).
## Two types of plots are available. One that shows all probes for each dataset at a give beta cutoff 
## value (boxplot) and one that shows all dataset for each probe at all beta cutoffs (lineplot)
## 
## betaData: a dataFrame with the beta summery info (created by extractSelectedGroups.R)
## betaDataSets: a list with the information to betaData (created by extractSelectedGroups.R)
## startCoord: numeric value of the start coordinate, minimum of 1
## stopCoord: numeric value of the stop coordinate, has to be greater than startCoord
## chr: chromosome, either numeric for chromosome number or character for X or Y chromosome
## betaCut: numeric value for the boxplot to determine the cufoff for betaValues to use (between 0 and 10)
## plotType: type of plot (either boxplot or lineplot)
## showGeneNames: boolean to display gene names on the plots
## setSelection: numeric vector with the number of sets to display graphs for

coordPlot = function(betaData, betaDataSets, startCoord, stopCoord, chr, betaCut = NULL, plotType = "boxplot", showGeneNames = FALSE, setSelection = NULL) {
  require(ggplot2)
  require(reshape2)

  #extract setnames and columnNames
  setNames = as.vector(betaDataSets$setNames)
  columnNames = as.vector(betaDataSets$columnNames)
  
  #Remove sets not needed
  if(!is.null(setSelection)) {
    
    #Adjust for probeinfo
    setSelection = c(0, setSelection) + 1
    selectedSets = unique(setNames)[setSelection]
    
    selectionVector = setNames %in% selectedSets
  } else {
    selectionVector = TRUE
  }
  
  setNames = setNames[selectionVector]
  columnNames = columnNames[selectionVector]
  
  #extract only data from the selected probes
  subdata = subset(betaData, CHR == chr & MAPINFO > startCoord & MAPINFO < stopCoord)
  subdata = subdata[,selectionVector]
  colnames(subdata) = columnNames
  
  #get the column names with the percentages
  percentageNames = rev(rev(columnNames)[1:11])

  if(plotType == "boxplot") {
    #boxplot start
    
    #betaCut checks
    if(is.null(betaCut)) {
      print("betaCut not set for boxplot; picking for you")
      betaCut = 6
    } else if(betaCut < 0 | betaCut > 10) {
      print("betaCut outside of range for boxplot; setting it for you")
      betaCut = 6
    } else {
      betaCut = betaCut + 1
    }
    
    #select the columns with the given cutoff
    valuePicker = which(columnNames == percentageNames[betaCut])
    graphData = subdata[,valuePicker]
    colnames(graphData) = unique(setNames)[-1]
    
    #Create the data for plot
    graphData = melt(t(graphData))
    colnames(graphData) = c("Dataset", "ProbeID", "Value")
    graphData = cbind(graphData, x = as.factor(subdata[as.character(graphData$ProbeID), "MAPINFO"]))
    graphData = cbind(graphData, gene = subdata[as.character(graphData$ProbeID), "GENE_BY_JBBR"])
    
    #Create the plot
    g = ggplot(graphData, aes(x = x, y = Value)) + 
      geom_boxplot() + theme_bw() + geom_jitter(width = 0.1) + ggtitle(percentageNames[betaCut]) +
      xlab("Coordinates") + ylab("Percentage") + ylim(c(0,1))
    
    #Paste on gene names
    if(showGeneNames) {
      #showGeneNames start
      
      #Get the gene info
      geneInfo = unique(graphData[,c("x", "gene")])
      geneInfo = geneInfo[order(geneInfo$x),]
      geneInfo$x = as.numeric(geneInfo$x)

      #Calculate the size of the rectangle for each gene
      for(geneName in unique(geneInfo$gene)) {
        if(geneName != "") {
          singleGeneInfo = subset(geneInfo, gene == geneName)
          g = g + annotate("rect", xmin = min(singleGeneInfo$x), xmax = max(singleGeneInfo$x), ymin = 0, ymax = 0.05, alpha = 0.2) +
            annotate(geom = "text", x = (min(singleGeneInfo$x) + max(singleGeneInfo$x)) / 2, y = 0.025, label = unique(singleGeneInfo$gene, hjust = 0.5, vjust = 0.5))
        }
      }
      #showGeneNames end
    }
    
    #boxplot end
  } else if(plotType == "lineplot") {
    #lineplot start
    
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
    
    #Adjust the betaValue and select the subset of data with the given betavalues
    if(!is.null(betaCut)) {
      betaCut = betaCut + 1
      graphData = subset(graphData, betaValue %in% betaCut)
    }
    
    #If no betavalue is set or more than one value is parsed
    if(is.null(betaCut) | length(betaCut) != 1) {
      #multiBeta start
      
      #create the list to hold the plots
      gplots = list()
      length(gplots) = length(unique(graphData$Dataset))
      
      #create the plots
      for(i in 1:length(gplots)) {
        #createPlots start
        
        #Not interesting to see 1 and 11
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
        
        #Paste on gene names
        if(showGeneNames) {
          #showGeneNames start
          
          #Get the gene info
          geneInfo = unique(graphData[,c("x", "gene")])
          geneInfo = geneInfo[order(geneInfo$x),]
          
          #Calculate the size of the rectangle for each gene
          for(geneName in unique(geneInfo$gene)) {
            if(geneName != "") {
              singleGeneInfo = subset(geneInfo, gene == geneName)
              g = g + annotate("rect", xmin = min(singleGeneInfo$x), xmax = max(singleGeneInfo$x), ymin = 0, ymax = 0.05, alpha = 0.2) +
                annotate(geom = "text", x = (min(singleGeneInfo$x) + max(singleGeneInfo$x)) / 2, y = 0.025, label = unique(singleGeneInfo$gene, hjust = 0.5, vjust = 0.5))
            }
          }
          #showGeneNames end
        }
        gplots[[i]] = g
        #createPlots end
      }
      
      g = multiplot(plotlist = gplots, cols = sqrt(length(unique(graphData$Dataset))))
      
      #multiBeta end
    } else {
      #singleBeta start
      
      #Create the plot
      g = ggplot(graphData, aes(x = x, y = Value, color = Dataset, shape = Dataset)) +
        geom_line(size = 1) + scale_color_manual(values = colorRampPalette(c("lightgrey", "darkblue"))(length(unique(graphData$Dataset)))) +
        geom_point(size = 3) + scale_shape_manual(values=1:length(unique(graphData$Dataset))) +
        xlab("Coordinates") + ylab("Percentage") +
        ggtitle(percentageNames[betaCut]) + theme_bw() + ylim(c(0,1))
      
      #Paste on gene names
      if(showGeneNames) {
        #showGeneNames start
        
        #Get the gene info
        geneInfo = unique(graphData[,c("x", "gene")])
        geneInfo = geneInfo[order(geneInfo$x),]
        
        #Calculate the size of the rectangle for each gene
        for(geneName in unique(geneInfo$gene)) {
          if(geneName != "") {
            singleGeneInfo = subset(geneInfo, gene == geneName)
            g = g + annotate("rect", xmin = min(singleGeneInfo$x), xmax = max(singleGeneInfo$x), ymin = 0, ymax = 0.05, alpha = 0.2) +
              annotate(geom = "text", x = (min(singleGeneInfo$x) + max(singleGeneInfo$x)) / 2, y = 0.025, label = unique(singleGeneInfo$gene, hjust = 0.5, vjust = 0.5))
          }
        }
        #showGeneNames end
      }
      #singleBeta end
    }
    #lineplot end
  }
  return(g)
}