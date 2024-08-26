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
      mostRecentSamplingDatum =  #BA lineage
      
      coordinateAttributeName = "location"

      treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)


      # 2. Extracting the spatio-temporal information embedded in the MCC tree

      treefile<- ""

      source("mccExtractions.r")
      mcc_tre = readAnnotatedNexus(treefile)
      mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
      write.csv(mcc_tab, "space_time_info_combined.csv", row.names=F, quote=F)

      # 3. Estimating the HPD region for each time slice
      mcc_tab=read.csv("space_time_info_combined.csv")
      nberOfExtractionFiles = nberOfTreesToSample
      prob = 0.95; precision = 0.025
      startDatum = min(mcc_tab[,"startYear"])

      polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


      # 4.1 spatial boundaries

      reloadborders<-TRUE
      if(reloadborders){
          borders = getData("GADM", country="ZAF", level=1)  #gadm(country = "ZAF", level = 1)
          borders2 = getData("GADM", country="ZAF", level=0)
          template_raster = crop(raster(""), extent(borders))
      }
}


# 4.2 Defining the different colour scales

minYear = min(mcc_tab[,"startYear"]); maxYear = max(mcc_tab[,"endYear"])
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1

##colors have to go to max time slice - AT LEAST
n_number_colours_needed<- max(round(endYears_indices))
n_repeats_discrete<- 10
c2<- (brewer.pal(9,"PuBuGn"))

colours<- rev(rep(c(c2), each=n_repeats_discrete))

colour_scale<- colorRampPalette(colours)(n_number_colours_needed)

endYears_colours = colour_scale[round(endYears_indices)]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
{
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"15")
}

# 5. Generating the dispersal history plot

pdf('',width=6, height=6.3,bg="white")

ptsize<- 0.4
pitjit<- 0.3
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")

# Define custom color breaks
custom_breaks <- c(10, 100, 1000, 10000, 100000)

# Create a custom color palette for these breaks
custom_palette <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(length(custom_breaks) - 1)

# Plot the raster with custom color breaks
plot(template_raster, col = custom_palette, colNA = "white", breaks = custom_breaks, box = F, axes = F, legend = F)

plot(borders, add=T, lwd=0.5, border="darkgrey")
plot(borders2, add=T,lwd=1, border="grey20")

for (i in length(polygons):1)
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
for (i in 1:dim(mcc_tab)[1])
{
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=4*1.1, lty=1, lcol="grey20", arr.col=endYears_colours[i], arr.pos=1, curve=0.3, dr=NA, endhead=F)
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=3.5*1.1, lty=1, lcol=endYears_colours[i], arr.col=endYears_colours[i], arr.pos=1, curve=0.3, dr=NA, endhead=F)
}
for (i in dim(mcc_tab)[1]:1)
{
  xs<- mcc_tab[i,"startLon"]
  ys<- mcc_tab[i,"startLat"]
  xe<- jitter(mcc_tab[i,"endLon"],pitjit)
  ye<- jitter(mcc_tab[i,"endLat"],pitjit)
  if (i == 1)
  {
    points(xs, ys, pch=21, col="gray20", bg=colour_scale[1], cex=.6)
  }
  points(xe, ye, pch=21, col="gray20", bg=endYears_colours[i], cex=.6)
}

xrange<- c(xmin(template_raster), xmax(template_raster))
yrange<- c(ymin(template_raster), ymax(template_raster))
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tab[,"startYear"]); rast[2] = max(mcc_tab[,"endYear"])
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.8, legend.shrink=0.2, smallplot=c(0.40,0.80,0.14,0.155),
     legend.args=list(text="", cex=0.8, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))

dev.off()
