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

regionPlot = function(chromosome, startCoord, stopCoords, betaData, betaDataSets) {
  
}