probes = c("cg21928406", "cg13283952", "cg22507154", "cg03760839")

basedir = "O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/ProcessedData/betavalues/"

sampleinfo = read.csv(paste(basedir, probes[1], ".csv", sep = ""), sep = ";", header = FALSE, stringsAsFactors = FALSE)[1:2,]

betadata = sapply(probes[1:length(probes)], function(probe) {
  print(probe)
  unlist(read.csv(paste(basedir, probe, ".csv", sep = ""), sep = ";", header = FALSE, stringsAsFactors = FALSE, skip = 2)[-1])
})

completeData = data.frame(sampleName = unlist(rev(rev(sampleinfo[1,])[-1])), group = unlist(sampleinfo[2,-1]), betadata)

ExcludeSamples = c("GSM990563", "GSM1013175", "GSM1013182", "GSM1013225", "GSM1013228", "GSM1013232", "GSM1013183", "GSM1013189", "GSM1013196", "GSM1013178", "GSM1013186", "GSM1013181", "GSM1013179", "GSM1013195", "GSM1013180", "GSM1013176", "GSM1013197", "GSM1013226", "GSM1013199", "GSM1013191", "GSM1013187", "GSM1013202", "GSM1013233", "GSM1013210", "GSM1013177", "GSM1013230", "GSM1013200", "GSM1013203", "GSM1013221", "GSM1013236", "GSM1013194", "GSM1013213", "GSM1013215", "GSM1013201", "GSM1013219", "GSM1013235", "GSM1013218", "GSM1013212", "GSM1013185", "GSM1013217", "GSM1013229", "GSM1013206", "GSM1013214", "GSM1013208", "GSM1013220", "GSM1013198", "GSM1013216", "GSM1013192", "GSM1013207", "GSM1013211", "GSM1013184", "GSM1013234", "GSM1013227", "GSM1013190", "GSM1013222", "GSM1013223", "GSM1013231", "GSM1013209", "GSM1013188", "GSM1013224", "GSM1013193", "GSM1013204", "GSM1013205")

completeData = subset(completeData, !sampleName %in% ExcludeSamples)

write.table(completeData, "O:/HE_MOMA-Data/MICROARRAY/Prostata/450K_MarmalAid/Scripts/extractedBetas.csv", sep = ";", dec = ",", row.names = FALSE)
