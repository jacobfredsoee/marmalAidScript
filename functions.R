saveResults = function(dataNames, countdata, fileName, dataColNames, probeInfo = NULL) {
  
  print(paste(Sys.time(), "|", "preparing data"))
  setNames = matrix(sapply(dataNames, function(dataName) {
    rep(dataName, ncol(countdata) / length(dataNames))
  }), nrow = 1)
  columnNames = matrix(rep(dataColNames, length(dataNames)), nrow = 1)
  
  if(!is.null(probeInfo)) {
    setNames = matrix(c(rep("probeInfo", ncol(probeInfo)), setNames), nrow = 1)
    
    columnNames = matrix(c(colnames(probeInfo), columnNames), nrow = 1)
    
    countdata = cbind(probeInfo, countdata)
  }
  seperator = ";"
  write.table(setNames, fileName, sep = seperator, col.names = FALSE, row.names = "setNames")
  write.table(columnNames, fileName, sep = seperator, col.names = FALSE, row.names = "columnNames", append = TRUE)
  
  cutpoint = 5000
  counter = 1
  limit = nrow(countdata)
  #limit = 15023
  
  print(paste(Sys.time(), "|", "Writing data to", fileName))
  startTime = Sys.time()
  while(counter * cutpoint < limit) {  
    
    startLine = ((counter - 1) * cutpoint) + 1
    stopLine = counter * cutpoint
    
    write.table(countdata[startLine:stopLine,], fileName, sep = seperator, col.names = FALSE, row.names = TRUE, append = TRUE)
    
    timeLeft = (limit/cutpoint - counter) * (difftime(Sys.time(), startTime, units = "mins") / counter)
    
    print(paste(Sys.time(), " | ", round(cutpoint * counter / limit * 100, 2), "% done", " | estimated time left: ", round(timeLeft, 2), " minutes", sep = ""))
    
    counter = counter + 1
  }
  startLine = ((counter - 1) * cutpoint) + 1
  write.table(countdata[startLine:limit,], fileName, sep = seperator, col.names = FALSE, row.names = TRUE, append = TRUE)
  
  print(paste(Sys.time(), "|", "writing rds file to", paste(scriptDir, "output.rds", sep = "/")))
  saveRDS(countdata, paste(scriptDir, "output.rds", sep = "/"))
  saveRDS(list(setNames = setNames, columnNames = columnNames), paste(scriptDir, "outputSets.rds", sep = "/"))
  
  print(paste(Sys.time(), " | ", "All done", sep = ""))
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