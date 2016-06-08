## The functions and examples used here require the following packages: ggplot2, reshape2, grid, manipulate
## You only need to run the following once ever, unless you change version of R
if(!"ggplot2" %in% rownames(installed.packages())) {
  install.packages("ggplot2")
}
if(!"reshape2" %in% rownames(installed.packages())) {
  install.packages("reshape2")
}
if(!"grid" %in% rownames(installed.packages())) {
  install.packages("grid")
}
if(!"manipulate" %in% rownames(installed.packages())) {
  install.packages("manipulate")
}

## Run the following once every time R is restarted:
source("O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/functions.R")
source("O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/probePlot.R")
source("O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/coordPlot.R")
require(manipulate)
betaData = readRDS("O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/output.rds")
betaDataSets = readRDS("O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/outputSets.rds")


## Propeplot
## Description

# A selection of probes you want to look at
probes = c("cg14686949", "cg16618605", "cg10276549", "cg01480550")

# Pass the probe names, the betaData and the betaDataSets to the function. We want a lineplot in this case
# This type of plot displays the percentage of samples that had a betavalue above a certain threshold. It
# does this for all samplesets included in the betaData file and draws a graph for each probe 
probePlot(probes, betaData, betaDataSets, plotType = "lineplot")

# Another posibility is to get a 
probes = c("cg00018204", "cg00206332", "cg01480550", "cg04929599", 
           "cg05880897", "cg07053880", "cg08783090", "cg10276549", 
           "cg12204574", "cg14686949", "cg16410962", "cg16618605",
           "cg18748981", "cg19930575", "cg25528646", "cg27049053")

probePlot(probes, betaData, betaDataSets, betaCut = 5, plotType = "boxplot")
