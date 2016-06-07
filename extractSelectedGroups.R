print(paste(Sys.time(), "|", "starting", basename(sys.frame(1)$ofile)))
scriptDir = dirname(sys.frame(1)$ofile)

PROBE_NUMBER = 485511
COLUMN_NUMBER = 11
BASEDIR = sub("(.*)/.+$","\\1",scriptDir)

source(paste(scriptDir, "functions.R", sep = "/"))

sampleSheet = read.csv(file = paste(scriptDir, "sampleGroups.csv", sep = "/"), sep = ";", stringsAsFactors = FALSE)

groups = sapply(unique(sampleSheet$Group), function(groupName) {
  
  subset(sampleSheet, Group == groupName)$Name
})

print(paste(Sys.time(), "|", length(groups[names(groups) == "0"][[1]]), "out of", length(unlist(groups)), "datasets not included"))
groups = groups[names(groups) != "0"]

resultMatrix = matrix(0, nrow = PROBE_NUMBER, ncol = length(groups) * COLUMN_NUMBER * 2)

for(i in 1:length(groups)) {
  group = groups[i]
  
  start = ((i - 1) * COLUMN_NUMBER * 2) + 1
  stop = start + COLUMN_NUMBER - 1
  
  for(datasetName in group[[1]]) {
    dataset = readRDS(paste(BASEDIR, "/ProcessedData/", datasetName, ".rds", sep = ""))

    resultMatrix[,start:stop] = resultMatrix[,start:stop] + dataset[,1:COLUMN_NUMBER]
    
    print(paste(Sys.time(), "|", "loading:", datasetName, "->", names(groups)[i]))
  }
  
  startPer = ((i - 1) * COLUMN_NUMBER * 2) + 1 + COLUMN_NUMBER
  stopPer = startPer + COLUMN_NUMBER - 1
  
  resultMatrix[,startPer:stopPer] = resultMatrix[,start:stop] / max(resultMatrix[,start:stop])
}

rownames(resultMatrix) = rownames(dataset)

print(paste(Sys.time(), "|", "loading: probe info data"))
probeInfo = read.csv(paste(scriptDir, "450k_essentialprobesinfo.csv", sep = "/"), sep = ";", header = TRUE)
rownames(probeInfo) = probeInfo$TARGETID
probeInfo = probeInfo[rownames(resultMatrix),]

saveResults(dataNames = names(groups),
            countdata = resultMatrix,
            fileName = paste(scriptDir, "output.csv", sep = "/"),
            dataColNames = colnames(dataset),
            probeInfo = probeInfo)



