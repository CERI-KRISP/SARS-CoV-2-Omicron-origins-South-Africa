# install.packages("devtools"); library(devtools)
# install_github("sdellicour/seraphim/unix_OS") # for Unix systems
## install_github("sdellicour/seraphim/windows") # for Windows systems
# install.packages("diagram")

library(diagram)
require(seraphim)
library(geodata)

setwd("")

loadstuff<- TRUE
if(loadstuff){

      # 1. Extracting the spatio-temporal information contained in posterior trees

      treefile<- ""

      localTreesDirectory = "Tree_extractions"
      allTrees = scan(file=treefile, what="", sep="\n", quiet=T)
      burnIn = 0
      randomSampling = FALSE
      nberOfTreesToSample = 100
      #mostRecentSamplingDatum =  #BA lineage

      coordinateAttributeName = "location"

      #treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)


      # 2. Extracting the spatio-temporal information embedded in the MCC tree

      treefile<- ""

      #source("mccExtractions.r")
      #mcc_tre = readAnnotatedNexus(treefile)
      #mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
      #write.csv(mcc_tab, "", row.names=F, quote=F)

      # 3. Estimating the HPD region for each time slice
      mcc_tab=read.csv("")
      
      nberOfExtractionFiles = nberOfTreesToSample
      prob80 = 0.80; precision = 0.025
      prob95 = 0.95; precision = 0.025
      prob99 = 0.99; precision = 0.025
      startDatum = min(mcc_tab[,"startYear"])

      polygons85 = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob80, startDatum, precision))
      polygons95 = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob95, startDatum, precision))
      polygons99 = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob99, startDatum, precision))
      

      # 4.1 spatial boundaries

      reloadborders<-TRUE
      if(reloadborders){
          borders = getData("GADM", country="ZAF", level=1)  #gadm(country = "ZAF", level = 1)
          borders2 = getData("GADM", country="ZAF", level=0)
      }
}

# 5. Generating the dispersal history plot

pdf('',width=6, height=6.3,bg="white")

ptsize<- 0.4

#Plot boundaries
plot(borders, add=F, lwd=0.5, border="darkgrey")
plot(borders2, add=T,lwd=1, border="grey20")

plot(polygons85[[1]], axes=F, col="#014636bf", add=T, border=NA)
plot(polygons95[[1]], axes=F, col="#01463680", add=T, border=NA)
plot(polygons99[[1]], axes=F, col="#01463640", add=T, border=NA)

xs<- mcc_tab[1,"startLon"]
ys<- mcc_tab[1,"startLat"]

points(xs, ys, pch=21, col="gray20", cex=.6)

dev.off()

polygons95



### BA example (95% HPD region) ###

loadstuff<- TRUE
if(loadstuff){
  
  # 1. Extracting the spatio-temporal information contained in posterior trees
  
  treefile<- ""
  
  localTreesDirectory = "Tree_extractions"
  allTrees = scan(file=treefile, what="", sep="\n", quiet=T)
  burnIn = 0
  randomSampling = FALSE
  nberOfTreesToSample = 100
  mostRecentSamplingDatum =  #BA1
  
  coordinateAttributeName = "location"
  
  treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)
  
  
  # 2. Extracting the spatio-temporal information embedded in the MCC tree
  
  treefile<- ""
  
  source("mccExtractions.r")
  mcc_tre = readAnnotatedNexus(treefile)
  mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
  write.csv(mcc_tab, "", row.names=F, quote=F)
  
  # 3. Estimating the HPD region for each time slice
  mcc_tab=read.csv("")
  
  nberOfExtractionFiles = nberOfTreesToSample
  #prob80 = 0.80; precision = 0.025
  prob95 = 0.95; precision = 0.01
  #prob99 = 0.99; precision = 0.025
  startDatum = min(mcc_tab[,"startYear"])
  
  #polygons85 = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob80, startDatum, precision))
  polygons95 = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob95, startDatum, precision))
  #polygons99 = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob99, startDatum, precision))
  
  
  # 4.1 spatial boundaries
  
  reloadborders<-TRUE
  if(reloadborders){
    borders = getData("GADM", country="ZAF", level=1)  #gadm(country = "ZAF", level = 1)
    borders2 = getData("GADM", country="ZAF", level=0)
  }
}

# 5. Generating the dispersal history plot

pdf('',width=6, height=6.3,bg="white")

#ptsize<- 0.4

#Plot boundaries
plot(borders, add=F, lwd=0.5, border="darkgrey")
plot(borders2, add=T,lwd=1, border="grey20")

#plot(polygons85[[1]], axes=F, col="#014636bf", add=T, border=NA)
plot(polygons95[[1]], axes=F, col="#01463680", add=T, border=NA)
#plot(polygons99[[1]], axes=F, col="#01463640", add=T, border=NA)

xs<- mcc_tab[1,"startLon"]
ys<- mcc_tab[1,"startLat"]

points(xs, ys, pch=21, col="gray20", cex=.6)

dev.off()
