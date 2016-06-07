
if(!exists("betaData")) {
  betaData = readRDS("O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/output.rds")
}
if(!exists("betaDataSets")) {
  betaDataSets = readRDS("O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/outputSets.rds")
}

#betaDataSets = list(setNames = betaDataSets[[1]], columnNames = betaDataSets[[2]])

probes = rownames(betaData)[grep("GABRE", betaData$GENE_BY_JBBR)]

columnSelect = intersectR(c(grep("allcancer", as.vector(betaDataSets$setNames), ignore.case = TRUE), grep("allhealthy", as.vector(betaDataSets$setNames), ignore.case = TRUE)), grep("%>=0.6", as.vector(betaDataSets$columnNames)))

betaDataVar = betaData[order(apply(betaData[,columnSelect], 1, var), decreasing = TRUE),]

probePlot(rownames(betaDataVar)[1:4], betaData, betaDataSets, 6, plotType = "lineplot")

require(manipulate)

manipulate({
  probePlot(rownames(betaDataVar)[1:9], betaData, betaDataSets, betaValue)
}, betaValue = slider(0,10))




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
    # Multiple plot function
    #
    # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
    # - cols:   Number of columns in layout
    # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
    #
    # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
    # then plot 1 will go in the upper left, 2 will go in the upper right, and
    # 3 will go all the way across the bottom.
    #
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      require(grid)
      
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



intersectR = function(...) {
  arguments = list(...)
  
  if(length(arguments) == 1) { #unlist from recursive calls
    arguments = arguments[[1]]
  }
  
  if(class(arguments) != "list" | length(arguments) == 1) { #something is not right here
    stop("something is wrong")
  }
  
  if(length(arguments) > 2) {
    arg1 = arguments[[length(arguments)]]
    arg2 = intersectR(arguments[1:(length(arguments)- 1)])
    
    return(intersect(arg1, arg2))
  } else {
    res = intersect(arguments[[1]], arguments[[2]])
    
    return(res)
  }
}