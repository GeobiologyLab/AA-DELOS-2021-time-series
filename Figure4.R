library(plotly)
library(RColorBrewer)
library(reshape2)
library(vegan)

set.seed(3)
dev.new(width=5, height=4)
load(file = "data/merged_snv_table_v138.1_annotations.Rdata")  ## this has control snvs removed; variable named "a"

md = read.csv("data/CCA_2_Input_MetadataTable_DecNov20.csv",header=TRUE,row.names = "X")

a = a[-which(a$Genus=="Ralstonia"),]

df = a[which(a$Var1%in%rownames(md)),]
df = data.frame(dcast(df,Var1~Var2, value.var="value",fun=sum))
rownames(df) = df$Var1

df = df[,-which(colnames(df)=="Var1")]

r = rrarefy(df,30000)

env = c("depth","sulfate_ppm","T_degC","pH_field","Flowrate_lps","d18O_ppt")

m = md[,which(colnames(md)%in%env)]

m = m[which(rownames(m)%in%rownames(df)),]

df = df[order(rownames(df)),]
m = m[order(rownames(m)),]

## gradient for points. 
tunnel_meter = gsub(".*?(\\d+).*", "\\1", rownames(m))
tunnel_meter = sort(as.integer(tunnel_meter))
color_scale <- colorRampPalette(c("lightblue", "darkblue"))
sample_colors <- color_scale(length(tunnel_meter))

ord = cca(df,m)

## To see variouls labels/alternative plots
# plot(ord)
# plot(ord, display = "sites")
# plot(ord,display="species")

plot(ord, display = "sites", type = "n",cex=3)
points(ord, display = "sites", col = sample_colors, pch = 16,cex=1.5)

env_vars <- scores(ord, display = "bp", choices = c(1, 2), scaling = 2)
arrows(0, 0, env_vars[, 1]*1.75, env_vars[, 2]*1.75, angle = 20, length = 0.1, col = "red")



