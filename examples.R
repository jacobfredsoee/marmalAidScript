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
