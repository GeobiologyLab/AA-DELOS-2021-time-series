library(vegan)
library(iCAMP)
library(tidyverse)

set.seed(3)

#load data
metadata = read.csv("data/distance-decay-metadata.csv")
load("data/dada2_pseudopool.Rdata")

#prepare the snv table
snv = seqtab.nochim
ctrl = snv[which(rownames(snv) == "Negative-Control_S65.cropped.trimmed.joined.fastq.gzjoin"),]
snv = snv[,-which(ctrl>0)]

snv = snv[which(rowSums(snv)>30000),]



#subset data
delos_ids = metadata[which(metadata$type=="fracture"),]
delos_ids = delos_ids[year(delos_ids$date) == 2020, ]$full_name
surface_ids = metadata[which(metadata$description=="surface"),]$full_name
grimsel_ids = c("GtsSB80001-20210401_S63.cropped.trimmed.joined.fastq.gzjoin","GtsPinkel-20210401_S62.cropped.trimmed.joined.fastq.gzjoin",
                "GtsISCinj2-20210401-3L_S60.cropped.trimmed.joined.fastq.gzjoin")




grimsel = data.frame(snv[which(rownames(snv)%in%grimsel_ids),])


surface = data.frame(snv[which(rownames(snv)%in%surface_ids),])

delos = snv[which(rownames(snv)%in%delos_ids),]

### Create a dataset of Grimsel and DELOS

grimsdelos = rbind(delos,grimsel)

grimsdelDM = data.frame(as.matrix(1-vegdist(grimsdelos,binary=TRUE)))
grimsdelDM = grimsdelDM[-which(rownames(grimsdelDM)%in%grimsel_ids),which(colnames(grimsdelDM)%in%gsub('-','.',grimsel_ids))]
grimsdelDM$TM = as.numeric(str_extract(rownames(grimsdelDM), "\\d+"))

### Create a dataset of DELOS fracture and surface

surfdelos = rbind(delos,surface)

surfaceDM = data.frame(as.matrix(1-vegdist(surfdelos,binary=TRUE)))
surfaceDM = surfaceDM[-which(rownames(surfaceDM)%in%surface_ids),which(colnames(surfaceDM)%in%gsub('-','.',surface_ids))]    
surfaceDM$TM = as.numeric(str_extract(rownames(surfaceDM), "\\d+"))
  



plotSorensen <- function(df, iter){
  

    r = rrarefy(df,30000)
    sorensen = vegdist(r,method="bray", binary=TRUE)
    sorensencol =  dist.3col(sorensen)
    sorensencol$sim = 1-sorensencol$dis
    sorensencol$separation =  abs(as.numeric(str_extract(sorensencol$name1, "\\d+")) - as.numeric(str_extract(sorensencol$name2, "\\d+")))
    print(sorensencol$separation)
    plot(sorensencol$separation,sorensencol$sim,col = c(rgb(0, 0, 0, alpha = 0.35)), pch = 16, cex=1.5,
         xlab="Lateral Distance or DELOS Tunnel Meter (m)",ylab="Sørensen Similarity",xlim=c(-10,5300))

      #fit <- lm(sim ~ separation, data = sorensencol)
      #abline(h=grimsmed, col = "#808000", lwd = 2, lty = 2) 
      #rect(xleft = -10, xright = 5000, ybottom = 0.015, ytop = .1198, col = rgb(128/255, 128/255, 0/255, alpha = 0.25), border = NA)
  
  
  points(grimsdelDM$TM,grimsdelDM$GtsSB80001.20210401_S63.cropped.trimmed.joined.fastq.gzjoin,
         col = "black", bg = rgb(128/255, 128/255, 0/255, alpha = .75),pch=23,cex=1.5)

  
  points(surfaceDM$TM,surfaceDM$BdSTicinoRiverBridge.20210821_S10.cropped.trimmed.joined.fastq.gzjoin,
         col = "black", bg = rgb(23/255, 190/255, 207/255,alpha=.75),pch=23,cex=1.5)

  legend("topright", legend = c("DELOS vs DELOS Comparison", "DELOS vs Surface Comparison", "DELOS vs Grimsel Comparison"),
         col = c("black", "#17becf", "#808000"), pch = c(16, 18, 18))
  
}



plotSorensen(delos,1)
