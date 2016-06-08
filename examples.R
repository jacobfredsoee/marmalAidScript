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


## propePlot
## Description

# A selection of probes you want to look at
probes = c("cg14686949", "cg16618605", "cg10276549", "cg01480550")

# Pass the probe names, the betaData and the betaDataSets to the function. We want a lineplot in this case
# This type of plot displays the percentage of samples that had a betavalue above a certain threshold. It
# does this for all samplesets included in the betaData file and draws a graph for each probe 
probePlot(probes, betaData, betaDataSets, plotType = "lineplot")

# Another posibility is to get an overview over the beta values for a set of probes in all the datasets at
# a given beta value cutoff
probes = c("cg21928406", "cg13283952", "cg22507154", "cg03760839",
           "cg20699586", "cg11100804", "cg08658787", "cg12111714",
           "cg10919344", "cg19247475", "cg04102163", "cg25351606",
           "cg26345105", "cg07539798", "cg08886154", "cg25133753")

# Same function as above, but now with a betaCut set (5 means % of probes with a betaValue >= 0.5), and as
# a boxplot
probePlot(probes, betaData, betaDataSets, betaCut = 5, plotType = "boxplot")

# Using the manipulate package makes it easy to change beta value cutoff (note the little gear that appears
# in the upper left hand corner on the plot area)
manipulate({
  probePlot(probes, betaData, betaDataSets, betaValue, plotType = "boxplot")
}, betaValue = slider(min = 0, max = 10, initial = 5))


## coordPlot
## Description

# This functions allows to survey the beta values of a given genomic region for all datasets at given beta
# value cutoffs
startCoord = 1
stopCoord = 500000
chr = 1

# Pass the betaData, betaDataSets, the start and stop coordinate for the genomic region of interrest, as well
# as the chromosome. In this case we want a lineplot
coordPlot(betaData, betaDataSets, startCoord, stopCoord, chr, plotType = "lineplot")

# It is also possible to only view selected beta value cutoffs by setting the betaCut variable
coordPlot(betaData, betaDataSets, startCoord, stopCoord, chr, betaCut = c(2,5,8), plotType = "lineplot")

# Or if only a single betaCut is chosen, it will show that in one plot with all the datasets
coordPlot(betaData, betaDataSets, startCoord, stopCoord, chr, betaCut = 8, plotType = "lineplot")

# In all cases it is possible to display gene names by setting showGeneNames to TRUE
coordPlot(betaData, betaDataSets, startCoord, stopCoord, chr, betaCut = 8, plotType = "lineplot", showGeneNames = TRUE)

# Finally, by setting plotType to boxplot, it displays how many percentage of samples are above a given beta 
# value cutoff for all datasets. 
coordPlot(betaData, betaDataSets, startCoord, stopCoord, chr, betaValue, plotType = "boxplot", showGeneNames = TRUE)

# Using the manipulate package makes it easy to change beta value cutoff (note the little gear that appears
# in the upper left hand corner on the plot area)
manipulate({
  coordPlot(betaData, betaDataSets, startCoord, stopCoord, chr, betaValue, plotType = "boxplot", showGeneNames = TRUE)
}, betaValue = slider(min = 0, max = 10, initial = 5))
