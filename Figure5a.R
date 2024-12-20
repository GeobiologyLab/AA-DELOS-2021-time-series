### This code creates the taxonomy time series plot of Figure 5a. Excuting this code will not generate the legend for the figure or x-axis labels. To do this, remove guide = "none" and add axis.text.x = element_text(angle = 45, hjust = 1)

library(dada2)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(dplyr)
library(reshape2)

load(file="data/dada2_pseudopool_annotations_v138.2.Rdata")
load(file="data/dada2_pseudopool.Rdata")

#remove control ASVs
rows_to_check <- c("Negative-Control_S65.cropped.trimmed.joined.fastq.gzjoin")

# Identify columns where any of the specified rows have values > 1
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

expanded_a$KP <- paste(expanded_a$Kingdom, "; ", expanded_a$Phylum, sep="")
expanded_a$KtoC <- paste(expanded_a$Kingdom, "; ", expanded_a$Phylum, "; ", expanded_a$Class, sep="")

expanded_a[which(expanded_a$KP=="Bacteria; Nitrospirota"),]$KP <- expanded_a[which(expanded_a$KP=="Bacteria; Nitrospirota"),]$KtoC
expanded_a[which(expanded_a$KP=="Bacteria; Verrucomicrobiota"),]$KP <- expanded_a[which(expanded_a$KP=="Bacteria; Verrucomicrobiota"),]$KtoC

compressed_a <- expanded_a[,-which(colnames(expanded_a)%in%c("Kingdom","Phylum","Class","Order","Genus","Species","KtoC"))] %>%
  group_by(KP) %>%
  summarize(across(where(is.numeric), sum))

compressed_a <- as.data.frame(compressed_a)
rownames(compressed_a) <- compressed_a$KP

file_paths <- c("data/TM444_time_series_metadata.csv","data/TM1306C_time_series_metadata.csv", "data/TM1494_time_series_metadata.csv", "data/TM2848A_time_series_metadata.csv", "data/TM4652_time_series_metadata.csv")

soi <- c()
for (file_path in file_paths){
  ff <- read.csv(file_path)
  soi <- c(soi,ff$full_name)
}

filtered_a <- compressed_a[,which(colnames(compressed_a)%in%soi)]

mpdf <- t(filtered_a)
mpdf <- mpdf/rowSums(mpdf)

dpann <- c("Archaea; Aenigmarchaeota", "Archaea; Altiarchaeota", "Archaea; Huberarchaeota", "Archaea; Iainarchaeota",
           "Archaea; Mamarchaeota", "Archaea; Micrarchaeota", "Archaea; Nanoarchaeota", "Archaea; Nanohalarchaeota",
           "Archaea; Parvarchaeota", "Archaea; Undinarchaeota", "Archaea; Woesearchaeota") 

df <- data.frame(
  "Unclassified Bacteria" = rowSums(mpdf[,which(colnames(mpdf)%in%c("Bacteria; NA","Bacteria; Incertae Sedis"))]),
  "Nitrospirota; Leptospirillia" = mpdf[,which(colnames(mpdf)=="Bacteria; Nitrospirota; Leptospirillia")],
  "Nitrospirota; Thermodesulfovibrionia" = mpdf[,which(colnames(mpdf)=="Bacteria; Nitrospirota; Thermodesulfovibrionia")],
  "Chloroflexota" = mpdf[,which(colnames(mpdf)=="Bacteria; Chloroflexota")],
  "DPANN" = rowSums(mpdf[,which(colnames(mpdf)%in%dpann)]),
  "CPR (Patescibacteria)" = mpdf[,which(colnames(mpdf)=="Bacteria; Patescibacteria")],
  "Pseudomonadota" = mpdf[,which(colnames(mpdf)=="Bacteria; Pseudomonadota")],
  "Candidatus Kryptonia" = mpdf[,which(colnames(mpdf)=="Bacteria; Candidatus Kryptonia")],
  "Omnitrophia" = mpdf[,which(colnames(mpdf)=="Bacteria; Verrucomicrobiota; Omnitrophia")],
  "Ignavibacteriota" = mpdf[,which(colnames(mpdf)=="Bacteria; Ignavibacteriota")],
  "Spirochaetota" = mpdf[,which(colnames(mpdf)=="Bacteria; Spirochaetota")],
  "Planctomycetota" = mpdf[,which(colnames(mpdf)=="Bacteria; Planctomycetota")],
  "Elusimicrobiota" = mpdf[,which(colnames(mpdf)=="Bacteria; Elusimicrobiota")],
  check.names = FALSE
  )


df$Other <- rowSums(mpdf)-rowSums(df)
rowSums(df)

column_order <- c(
  c("Other",
  sort(setdiff(colnames(df), c("Other", "Unclassified Bacteria","CPR","DPANN")))),  c("CPR","DPANN", "Unclassified Bacteria")
)


column_order <- c(
  c("Other",
    sort(setdiff(colnames(df), c("Other", "Unclassified Bacteria")))),  c("Unclassified Bacteria")
)

# Reorder the columns of the data frame
df <- df[, column_order, drop = FALSE]

plot_list <- list()

for (file_path in file_paths) {
  ff <- read.csv(file_path)
  
  sampleIDs <- ff$full_name ## sampleID matches colnames(a)
  tmp <- df[which(rownames(df)%in%sampleIDs),]
  #rownames(tmp) <- compressed_a$Tracked_Taxa
  tmp.melt <- melt(t(tmp))
  parts <- strsplit(as.character(tmp.melt$Var2), "-")
  extracted_numbers <- sapply(parts, function(x) {
    # Check if the split has at least 2 parts to avoid errors
    if (length(x) >= 2) {
      strsplit(x[2], "_")[[1]][1]
    } else {
      NA  # Return NA if the expected structure is not found
    }
  })  
  tmp.melt$Date <- as.Date(extracted_numbers,format = "%Y%m%d")
  tmp.melt$Var1[tmp.melt$Var1 == "Other"] <- NA
  plot <- ggplot(tmp.melt[complete.cases(tmp.melt), ], aes(x = Date, y = value, fill = Var1)) +
    geom_area(alpha=0.5) +
    geom_bar(stat = "identity", position = "stack", width = 5,color="black",size=0.25) +  # Set the width parameter to 0.8 for consistent bar widths
    labs(y = "Relative Abundance", x = element_blank()) + # x = "Date", fill = "Tracked_Taxa") +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(12, "Paired"),"gray"),guide="none") +  # guide = "none" removes the legend. Use the "tab20" color palette and remove the legend
    theme_minimal() +
    theme(axis.text.x = element_blank(), 
          axis.title.y = element_blank(),#axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 12)) + ylim(0, 1) +
    xlim(c(as.Date("20201015",format="%Y%m%d"),as.Date("20220101",,format="%Y%m%d")))
  
  plot_list[[file_path]] <- plot
}

# # Add x-axis labels to the bottom of the last plot
# last_plot <- plot_list[[length(plot_list)]]
# last_plot <- last_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# # Replace the last plot in the plot list with the updated plot
# plot_list[[length(plot_list)]] <- last_plot

combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))  

combined_plot
