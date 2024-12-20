library(dada2)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(dplyr)
library(reshape2)
library(vegan)

set.seed(3)


md = read.csv('data/delos-surface-grimsel-sites-labels.csv')




load(file="data/dada2_pseudopool_annotations_v138.2.Rdata")
load(file="data/dada2_pseudopool.Rdata")
#remove control ASVs
rows_to_check <- c("Negative-Control_S65.cropped.trimmed.joined.fastq.gzjoin")
# Identify columns where any of the specified rows have values 1
columns_to_keep <- which(seqtab.nochim[which(rownames(seqtab.nochim)=="Negative-Control_S65.cropped.trimmed.joined.fastq.gzjoin"),]==0)
# Subset the DataFrame to keep only these columns
cleaned.seqtab <- seqtab.nochim[, columns_to_keep]
asv_table <- as.data.frame(t(seqtab.nochim))
asv_table$ASV <- rownames(asv_table)
# Convert taxa to a data frame and add rownames as a column
taxa_table <- as.data.frame(taxa)
taxa_table$ASV <- rownames(taxa_table)
# Merge the two data frames by ASV
a <- merge(asv_table, taxa_table, by = "ASV")
a <- a[-which(a$Genus=="Ralstonia"),]
expanded_a <- a

compressed_a <- expanded_a[,-which(colnames(expanded_a)%in%c("Kingdom","Class","Order","Genus","Species","KtoC"))] %>%
  group_by(Phylum) %>%
  summarize(across(where(is.numeric), sum))

compressed_a <- as.data.frame(compressed_a)
compressed_a <- compressed_a[which(complete.cases(compressed_a)),]
rownames(compressed_a) <- compressed_a$Phylum

compressed_a <- t(compressed_a[,which(colnames(compressed_a)%in%md$full_name)])


phylum <- rrarefy(compressed_a,min(rowSums(compressed_a)))

aitch <- vegdist(phylum,method="robust.aitchison")

isomap_result <- isomap(as.matrix(aitch), ndim = 2, k = 3)  # k controls the number of nearest neighbors
coords <- isomap_result$points

md_ordered <- md[match(rownames(coords), md$full_name), ]

colors <- as.factor(md_ordered$location_id)
plot(isomap_result, col = colors, pch=19,ylab='Isomap2',xlab='Isomap1',cex=1.5)
legend("bottomright", legend = levels(as.factor(md$location_id)), col = 1:length(levels(as.factor(md$location_id))), pch = 19)

## To see the which point belongs to which sample, add
#text(coords[,1], coords[,2], labels = rownames(coords), cex = 0.7, pos = 3)



