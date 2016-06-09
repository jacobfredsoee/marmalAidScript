## probePlot
## function to create plots from a given set of probes and a corresponding beta summary set (created from
## extractSelectedGroups.R).
## Two types of plots are available. One that shows all probes for each dataset at a give beta cutoff 
## value (boxplot) and one that shows all dataset for each probe at all beta cutoffs (lineplot)
## 
## probes: a vector with the names of the probes
## betaData: a dataFrame with the beta summery info (created by extractSelectedGroups.R)
## betaDataSets: a list with the information to betaData (created by extractSelectedGroups.R)
## betaCut: numeric value for the boxplot to determine the cufoff for betaValues to use (between 0 and 10)
## plotType: type of plot (either boxplot or lineplot)
## force: boolean value. lineplot will display a maximum of 9 plots unless forced is set to true

probePlot = function(probes, betaData, betaDataSets, betaCut = NULL, plotType = "boxplot", force = FALSE, setSelection = NULL) {
  require(ggplot2)
  require(reshape2)
  
  #ensure it is not a factor
  probes = as.character(probes)

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
  subdata = betaData[probes, selectionVector]
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

    #Create the plot
    g = ggplot(melt(graphData), aes(x = as.factor(variable), y = value)) +
      geom_boxplot() + theme_bw() + geom_jitter(width = 0.1) + ggtitle(percentageNames[betaCut]) +
      xlab("Dataset") + ylab("Percentage") + ylim(c(0,1))
    
    #boxplot end
  } else if(plotType == "lineplot") {
    #lineplot start
    
    #check if too many probes were parsed
    if(length(probes) > 9 & !force) {
      stop("too many probes (max 9)! set force to TRUE in order to force it, but use with caution")
    }
    
    #select the columns to use
    valuePicker = which(columnNames %in% percentageNames)
    graphData = subdata[,valuePicker]
    colnames(graphData) = as.vector(sapply(unique(setNames)[-1], function(setName) rep(setName, length(percentageNames))))
    
    #create data for graph
    graphData = melt(t(graphData))
    colnames(graphData) = c("Dataset", "ProbeID", "Value")
    graphData = cbind(graphData, x = rep(1:length(percentageNames), length(unique(graphData$ProbeID)) * length(unique(graphData$Dataset))))
    
    #create the list to hold the plots
    gplots = list()
    length(gplots) = length(unique(graphData$ProbeID))
    
    #create the plots
    for(i in 1:length(unique(graphData$ProbeID))) {
      g = ggplot(subset(graphData, ProbeID == unique(graphData$ProbeID)[i]), aes(x = x, y = Value, color = Dataset, shape = Dataset)) + 
        geom_line() + 
        geom_point() +
        scale_x_discrete(limits=percentageNames) + xlab("") + ylab("Percentage") +
        ggtitle(as.character(unique(graphData$ProbeID)[i])) + ylim(c(0,1))
      gplots[[i]] = g + theme_bw()
    }
    
    #plot them all
    g = multiplot(plotlist = gplots, cols = sqrt(length(unique(graphData$ProbeID))))
    
    #lineplot end
  } else {
    stop("plotType not recognized")
  }
  
  #return the ggplot
  return(g)
}