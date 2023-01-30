
setwd("/Users/banana/documents/bioinformatics/cel-counts")

#load packages: -----
install.packages("Matrix", repos="http://R-Forge.R-project.org")
install.packages("lattice")
install.packages("fdrtool")
install.packages("rpart")
install.packages("ggplot2")
install.packages ("rpart.plot") 

# Install Bioconductor packages - done 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::available("")
BiocManager::install(c("limma"))
BiocManager::install(c("Biostrings"))
BiocManager::install(c("Biobase"))
BiocManager::install(c("genefilter"))
BiocManager::install("simpleaffy")
BiocManager::install("gcrma")
BiocManager::install("pd.rta.1.0")

remove.packages(      "ff")
install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz',repos=NULL)
BiocManager::install("oligo") 
BiocManager::install("oligoClasses") 

library(oligoClasses)
library(oligo)
library(affy)
library(affyPLM)
library(Matrix)
library(lattice)
library(org.Rn.eg.db)
library(fdrtool)
library(rpart)
library(ggplot2)
library(limma)
library(Biostrings)
library(Biobase)
library(genefilter)
library(BiocManager)
library(simpleaffy)

## To read in data: -----

celpath="/Users/banana/documents/bioinformatics/cel-counts/RawData/"

# import CEL files containing raw probe-level data into an R AffyBatch object
list = list.files(celpath,full.names=TRUE)
data = oligo::read.celfiles(list)
data

expression= exprs(data)
expression


#Retrieving annotation
ph = data@phenoData
ph@data
# Renames the files in a more meaningful way -
ph@data[ ,1] = c("A1","A2","A3","A4","A5","B1","B2","B3","B4","B5","C1","C2","C3","C4","C5","D1","D2","D3","D4","D5")
ph@data

ph@data[ ,1] = c("1 CeA:Dehydrated","2 CeA:Dehydrated","3 CeA:Dehydrated","4 CeA:Dehydrated","5 CeA:Dehydrated",
                 "1 CeA:Control","2 CeA:Control","3 CeA:Control","4 CeA:Control","5 CeA:Control",
                 "1 BLA:Dehydrated","2 BLA:Dehydrated","3 BLA:Dehydrated","4 BLA:Dehydrated","5 BLA:Dehydrated",
                 "1 BLA:Control","2 BLA:Control","3 BLA:Control","4 BLA:Control","5 BLA:Control")

#Retrieving probe annotation using Oligo
feat = data@featureData
feat
feat@data 

#Retrieving experiment annotation using Oligo 
exp = data@experimentData
exp


for (i in 1:20) 
  name = paste("image",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=ph@data$sample[i])
  dev.off()


Pset = fitProbeLevelModel(data)

for (i in 1:20)
{
  name = paste("pseudoimage",i,".jpg",sep="")
  jpeg(name)
  image(Pset,which=i,type="residuals",main=ph@data$sample[i])
  dev.off()
}


for (i in 1:20)
{
  name = paste("histogram",i,".jpg",sep="")
  jpeg(name)
  hist(data[,i],lwd=2,which='both',ylab='Density',xlab='Log2 intensities',main=ph@data$sample[i])
  dev.off()
}

#Plotting histograms - Together
pmexp = oligo::pm(data) 
SampleNames= vector()
logs = vector()
for (i in 1:20)
{
  SampleNames = c(SampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,Condition=SampleNames)
dataHist2 = ggplot(logData, aes(logInt, colour = Condition)) 
dataHist2 + geom_density() + labs(fill = "SampleNames", x = "Log Intensitiy", y = "Density") + theme(text = element_text(size=20))


# Principal components analysis: PCA ----- 
data
df_pca <- prcomp(pmexp)
df_out <- as.data.frame(df_pca$rotation) 
df_out 

# Create PCA groupings: Renamed each .CEL file 
df_out$Group = c("DH-CeA","DH-CeA", "DH-CeA","DH-CeA","DH-CeA","CT-CeA","CT-CeA","CT-CeA","CT-CeA","CT-CeA",
                 "DH-BLA", "DH-BLA","DH-BLA", "DH-BLA", "DH-BLA","CT-BLA","CT-BLA","CT-BLA","CT-BLA","CT-BLA")
df_out$Group1 = c("DH","DH","DH","DH","DH","CT","CT","CT","CT","CT","DH","DH","DH","DH","DH","CT","CT","CT","CT","CT" )
df_out$Group2 = c("CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA",
                 "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA","BLA" )
# Create subsets of Groups - 
df_BLA <- subset(df_out,Group2=="BLA") 
df_CeA <- subset(df_out,Group2=="CeA") 

head(df_out) 
dim(df_out) 

# Time to create plots

# Raw PCA plots for all samples ----- 1 

Title="Experimental Group vs Amygdala Region- Raw PCA"
Group="Experimental Group vs Amygdala Region"
Output="All CRaw_PCA.pdf"

#Add the percentage of PCA 
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
#Generate the plot
p<-ggplot(df_out,aes(x=PC1,y=PC2, color=Group))+ 
  xlab(percentage[1]) + ylab(percentage[2])+geom_point(size=6)+
  ggtitle(Title)+theme(axis.text=element_text(size=12),
  axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
  axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
  legend.text = element_text(size=12),
  legend.title.align = 0.5,
  legend.title = element_text(size=12),
  plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
  plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
# Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

# Raw PCA plots for dehydrated vs control -  Basolateral Amygdala ----- 2
Title="Basolateral Amygdala of Dehydrated vs Control Groups- Raw PCA"
Group="Experimental Group"
Output="BLA-DHvsCT-Raw_PCA.pdf"

df_BLA <- subset(df_out,Group2=="BLA") 


percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

# Generate the plot - Plot below is focused on the BLA of Dehydrated and control group

p<-ggplot(df_BLA,aes(x=PC1,y=PC2, color=Group1))+ 
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  geom_point(size=6)+
  ggtitle(Title)+
  theme(axis.text=element_text(size=12),
  axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
  axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
  legend.text = element_text(size=12),
  legend.title.align = 0.5,
  legend.title = element_text(size=12),
  plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
  plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
# Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

# Raw PCA plots for dehydrated vs control -  Central Amygdala ----- 3
Title="Central Amygdala of Dehydrated vs Control Groups- Raw PCA"
Group="Experimental Group"
Output="CeA-DHvsCT-Raw_PCA.pdf"

# Add the percentage of PCA 

percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

# Generate the plot - Plot below is focused on the CeA of Dehydrated and control group
p<-ggplot(df_CeA,aes(x=PC1,y=PC2, color=Group1))+ 
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  geom_point(size=6)+
  ggtitle(Title)+
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(size=12),
        legend.title.align = 0.5,
        legend.title = element_text(size=12),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
        plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
#Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

## Normalising data: -------

# Normalisation - RMA
data.rma = oligo::rma(data)
data.matrix = oligo::exprs(data.rma)

# Boxplots ------------- 

#Boxplots - RAW
pmexp = oligo::pm(data)
sampleNames = vector()
logs = vector()
for (i in 1:20)
{
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,sampleName=sampleNames)
dataBox = ggplot(logData,aes(sampleName,logInt))
dataBox + geom_boxplot() + ylim(2,16) + ggtitle("Before Normalization") + labs(x="Samples")

#Boxplots - Normalised
sampleNames = vector()
normlogs = vector()
for (i in 1:20)
{
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(data.matrix)[1]))
  normlogs = c(normlogs,data.matrix[,i])
}
normData = data.frame(norm_logInt=normlogs,sampleName=sampleNames)
dataBox = ggplot(normData, aes(sampleName,norm_logInt))
dataBox + geom_boxplot() + ylim(2,16) + ggtitle("After Normalization") + labs(x="Samples")



# Perfect match and mismatch: NO NEED - margins too large
par(mfrow=c(2,1))
hist(log2(pm(Data_affy [,1])),breaks=100,col="steelblue",main="PM",xlim=c(4,12))
hist(log2(mm(Data_affy [,1])),breaks=100,col="steelblue",main="MM",xlim=c(4,14))

dev.off()

#########################################################################

# Normalised PCA plots: 

# Renamed CEL files for it, so extract new:
mappath="/Users/banana/documents/bioinformatics/cel-counts/HmapData2/"

# import CEL files containing raw probe-level data into an R AffyBatch object
map_list = list.files(mappath,full.names=TRUE)
map_data = oligo::read.celfiles(map_list)

# Normalise new data...
data.map = oligo::rma(map_data) 
map.matrix = oligo::exprs(map_data)

# Time to prepare PCA... 

df_pca <- prcomp(map.matrix) 
df_out <- as.data.frame(df_pca$rotation)


df_out$Group = c("CT-BLA","CT-BLA","CT-BLA","CT-BLA","CT-BLA","DH-BLA", "DH-BLA","DH-BLA", "DH-BLA", "DH-BLA",
                 "CT-CeA","CT-CeA","CT-CeA","CT-CeA","CT-CeA","DH-CeA","DH-CeA", "DH-CeA","DH-CeA","DH-CeA")
df_out$Group1 = c("CT","CT","CT","CT","CT","DH","DH","DH","DH","DH","CT","CT","CT","CT","CT","DH","DH","DH","DH","DH")
df_out$Group2 = c("BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA","BLA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA")

# Create subsets of Groups
df_BLA <- subset(df_out,Group2=="BLA") 
df_CeA <- subset(df_out,Group2=="CeA") 


# 1) all data 

Title="Experimental group vs Amygdala region- RMA Normalised PCA"
Group="Experimental Group vs Amygdala Region"
Output="All-Norm_PCA.pdf"

#Generate the plot
p<-ggplot(df_out,aes(x=PC1,y=PC2, color=Group))+ 
  xlab(percentage[1]) + ylab(percentage[2])+geom_point(size=6)+
  ggtitle(Title)+theme(axis.text=element_text(size=12),
                       axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
                       axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
                       legend.text = element_text(size=12),
                       legend.title.align = 0.5,
                       legend.title = element_text(size=12),
                       plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
                       plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
#Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
library(ggplot2)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()


## 2) Basolateral amygdala 

Title="Basolateral Amygdala of Dehydrated vs Control Groups- Normalised PCA"
Group="Experimental Group"
Output="BLA-DHvsCT-Norm_PCA.pdf"


percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(df_BLA,aes(x=PC1,y=PC2, color=Group1))+ 
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  geom_point(size=6)+
  ggtitle(Title)+
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(size=12),
        legend.title.align = 0.5,
        legend.title = element_text(size=12),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
        plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
# Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

# 3) Central Amygdala 
Title="Central Amygdala of Dehydrated vs Control Groups- Normalised PCA"
Group="Experimental Group"
Output="CeA-DHvsCT-Norm_PCA.pdf"

p<-ggplot(df_CeA,aes(x=PC1,y=PC2, color=Group1))+ 
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  geom_point(size=6)+
  ggtitle(Title)+
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(size=12),
        legend.title.align = 0.5,
        legend.title = element_text(size=12),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
        plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
#Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()


############################################################################

# Removed the outlier for Normalised CeA and BLA PCA ------------ 4


# Find out which sample has the outlier

new_out$Group3 = c("CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA",
                   "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA","BLA" )

newpath="/Users/banana/documents/bioinformatics/cel-counts/Outlier2"

# import CEL files containing raw probe-level data into an R AffyBatch object
pca_list = list.files(newpath,full.names=TRUE)
pca_data = oligo::read.celfiles(pca_list)

# Normalise new data...
data.pca = oligo::rma(pca_data)
pca.matrix = oligo::exprs(data.pca)

# PCA of new data 
new_pca <- prcomp(pca.matrix) 
new_out <- as.data.frame(new_pca$rotation)

new_out$Group = c("CT-BLA","CT-BLA","CT-BLA","CT-BLA","CT-BLA","DH-BLA", "DH-BLA","DH-BLA", "DH-BLA", 
                  "CT-CeA","CT-CeA","CT-CeA","CT-CeA","CT-CeA","DH-CeA","DH-CeA", "DH-CeA","DH-CeA")
new_out$Group1 = c("CT","CT","CT","CT","CT","DH","DH","DH","DH","CT","CT","CT","CT","CT","DH","DH","DH","DH")
new_out$Group2 = c("BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA")

# Create subsets of Groups 
new_BLA <- subset(new_out,Group2=="BLA") 
new_CeA <- subset(new_out,Group2=="CeA") 


# Basolateral Amygdala 
Title="Basolateral Amygdala of Dehydrated vs Control Groups- Normalised PCA"
Group="Experimental Group"
Output="BLA-DHvsCT-Norm-Outlier_PCA.pdf"

p<-ggplot(new_BLA,aes(x=PC1,y=PC2, color=Group1))+ 
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  geom_point(size=6)+
  ggtitle(Title)+
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(size=12),
        legend.title.align = 0.5,
        legend.title = element_text(size=12),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
        plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
# Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

# Central Amygdala
Title="Central Amygdala of Dehydrated vs Control Groups- Normalised PCA"
Group="Experimental Group"
Output="CeA-DHvsCT-Norm_PCA.pdf"

p<-ggplot(new_CeA,aes(x=PC1,y=PC2, color=Group1))+ 
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  geom_point(size=6)+
  ggtitle(Title)+
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=14,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.text = element_text(size=12),
        legend.title.align = 0.5,
        legend.title = element_text(size=12),
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "mm"),
        plot.title = element_text(size=18,face="bold",hjust=0.5, margin = margin(t = 0, r = 0, b = 30, l = 0)))
p
#Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

##############################################################################

# Measuring outliers to see if they should be extracted from data 

remotes::install_github("privefl/bigutilsr")
library(bigutilsr)
install.packages("bigstatsr")
library(bigstatsr)


pca <- prcomp(X, scale. = TRUE, rank. = 10)
U <- df_pca$rotation

theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2]) + coord_equal()


# Measuring outlierness 
apply(U, 2, function(x) which( abs(x - mean(x)) > (3 * sd(x)) ))


###################################################################################


##### Heatmap clustering analysis 



#data.matrix = oligo::exprs(data.rma)
head(Biobase::pData(data.map))

exp_map <- Biobase::exprs(data.map)
PCA <- prcomp(t(exp_map), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])


dataGG$Condition = c("CT","CT","CT","CT","CT","DH","DH","DH","DH","DH","CT","CT","CT","CT","CT","DH","DH","DH","DH","DH")
dataGG$Region = c( "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA", "BLA","BLA",
                   "CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA","CeA")


annotation_for_heatmap <- 
  data.frame(Condition = dataGG$Condition,Region = dataGG$Region)

row.names(annotation_for_heatmap) <- row.names(pData(data.map))

dists <- as.matrix(dist(t(exp_map), method = "manhattan")) 
rownames(dists) <- row.names(pData(data.map))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL

diag(dists) <- NA
ann_colors <- list(
  Condition = c(DH = "chartreuse4", CT = "burlywood3"),
  Region = c(CeA = "blue4", BLA = "cadetblue2")
)
install.packages("pheatmap")
library(pheatmap)
pheatmap(dists, col = (hmcol), 
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")

###################################################################################
#######  Volcano plot: ------  

ph@data[ ,2] = c("DH_CeA","DH_CeA", "DH_CeA","DH_CeA","DH_CeA","CT_CeA","CT_CeA","CT_CeA","CT_CeA","CT_CeA",
                 "DH_BLA", "DH_BLA","DH_BLA", "DH_BLA", "DH_BLA","CT_BLA","CT_BLA","CT_BLA","CT_BLA","CT_BLA")
colnames(ph@data)[2]="conditions_regions"
ph@data
groups = ph@data$conditions_regions
f = factor(groups,levels=c("DH_CeA","CT_CeA","DH_BLA","CT_BLA")) 
design = model.matrix(~0 + f) 
colnames(design) = c("DH_CeA","CT_CeA","DH_BLA","CT_BLA") # double check this part 
data.fit = lmFit(data.matrix,design)

## Currently only comparing the CeA of DH and CT
contrast.matrix = makeContrasts(DH_CeA-CT_CeA,levels=design) 
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

names(data.fit.eb)
data.fit.eb$t[1:10,]
data.fit.eb$p.value

write.csv(data.fit.eb, file = "DHvsCT-CeA_output.csv")

# save all into an object:
tableTop <- topTable(data.fit.eb, number = 40000)
#Extract a table of the top-ranked genes from a linear model fit


p1 <- order(abs(tableTop$logFC), decreasing= TRUE)[1:100]
length(p1)

downregulated = tableTop[order(tableTop$logFC),]
?abs()

# tabletop 

# p2: index for the 25 with smallest P.Val
p2 <- order(abs(tableTop$P.Val), decreasing = FALSE)[1:100]
length(p2)
# union of p1 and p2
p <- union(p1,p2)
length(p)

#Volcano plot
EnhancedVolcano(tableTop,
                lab = rownames(tableTop),
                x = 'logFC',
                y = 'P.Value',
                title = 'Control vs Dehydrated Central Amygdala',
                pCutoff = 0.01,
                FCcutoff = 0.75,
                pointSize = 3.0,
                labSize = 6.0,
                shape = c(1, 4, 23, 25),
                colAlpha = 1)


# Old school way of doing it below: 
name = "CT-DH-CeA,Volcano.jpg"
jpeg(name)
volcanoplot(data.fit.eb,coef=1,highlight=10,names=data.fit.eb$genes$ID)
abline(v=1.5)
abline(h=-log(0.001, base = 10)) # 0.001 cuts it off at p-value 0.001
dev.off()


## Do the same - this time BLA: 
contrast.matrix = makeContrasts(DH_BLA-CT_BLA,levels=design) 
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

names(data.fit.eb)
data.fit.eb$t[1:10,]
data.fit.eb$p.value

write.csv(data.fit.eb, file = "DHvsCT-BLA_output.csv")

tableTop <- topTable(data.fit.eb, number = 40000)
view(tableTop)

#Volcano plot
EnhancedVolcano(tableTop,
                lab = rownames(tableTop),
                x = 'logFC',
                y = 'P.Value',
                title = 'Control vs Dehydrated Basolateral Amygdala',
                pCutoff = 0.01,
                FCcutoff = 0.75,
                pointSize = 3.0,
                #labSize = 6.0,
                shape = c(1, 4, 23, 25),
                colAlpha = 1)
dev.off()
## Again: BLA vs CeA of dehydrated models 

contrast.matrix = makeContrasts(DH_CeA-DH_BLA,levels=design) 
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

names(data.fit.eb)
data.fit.eb$t[1:10,]
data.fit.eb$p.value

write.csv(data.fit.eb, file = "CeAvsBLA-DH_output.csv")
tableTop <- topTable(data.fit.eb, number = 40000)

#Volcano plot
EnhancedVolcano(tableTop,
                lab = rownames(tableTop),
                x = 'logFC',
                y = 'P.Value',
                title = 'Central vs Basolateral Amygdala of Dehydrated Models',
                pCutoff = 0.01,
                FCcutoff = 1.25,
                pointSize = 3.0,
                labSize = 6.0,
                shape = c(1, 4, 23, 25),
                colAlpha = 1)

dev.off()
# Volcano plots: subsetted version CeA --------


# How to subset: (stratify by region) ]
ph@data
ph@varMetadata
data_CeA <- ph@data[ph@data[,"conditions_regions"] %in% c("DH_CeA","CT_CeA"),]
data_CeA
data.matrix_CeA <- data.matrix[,1:10] 

groups = data_CeA$conditions_regions
f = factor(groups,levels=c("DH_CeA","CT_CeA")) 
design = model.matrix(~0 + f) 
colnames(design) = c("DH_CeA","CT_CeA") 
data.fit = lmFit(data.matrix_CeA,design)

# Currently only comparing the CeA of DH and CT
contrast.matrix = makeContrasts(DH_CeA-CT_CeA,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

names(data.fit.eb)
data.fit.eb$t[1:10,]
data.fit.eb$p.value

write.csv(data.fit.eb, file = "Unannotated_output.csv")

tableTop <- topTable(data.fit.eb, number = 100)


#Volcano plot
name = "CT-DH-CeA,Volcano.jpg"
jpeg(name)
volcanoplot(data.fit.eb,coef=1,highlight=10)
dev.off()



##############################################################################################

##### 12.5 - extracting results: ######

### Creating a CeA histogram 
# First creating results for CeA:
contrast.matrix = makeContrasts(DH_CeA-CT_CeA,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit_CeA = eBayes(data.fit.con)

table_CeA <- topTable(data.fit_CeA, number = Inf)
head(table_CeA) 

data.table::update.dev.pkg() 
BiocManager::install("EBImage")
install.packages("GiNA")
library(GiNA)
library(EBImage)


# 1. Adjusted p-value
hist(table_CeA$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     main = "Control vs Dehydrated - Central Amygdala", xlab = "adjusted p-values")
# 2. Normal p-value 
hist(table_CeA$P.Value, col = brewer.pal(3, name = "Set2")[2],
     main = "Control vs Dehydrated - Central Amygdala", xlab = "p-values")

# Histograms for BLA:
contrast.matrix_BLA = makeContrasts(DH_BLA-CT_BLA,levels=design) 
data.fit.con_BLA = contrasts.fit(data.fit,contrast.matrix_BLA)
data.fit_BLA = eBayes(data.fit.con_BLA)

table_BLA <- topTable(data.fit_BLA, number = Inf)
head(table_BLA)

#First for adjusted p-value
hist(table_BLA$adj.P.Val, col = brewer.pal(3, name = "Set2")[2],
     main = "Control vs Dehydrated - Basolateral Amygdala", xlab = "adjusted p-values")
# Second for normal p-value 
hist(table_BLA$P.Value, col = brewer.pal(3, name = "Set2")[2],
     main = "Control vs Dehydrated - Basolateral Amygdala", xlab = "p-values")


#### Plots: 

data.fit.rich = eBayes(data.fit.con)
tableRich <- topTable(data.fit.rich, number = 30)


names(data.fit.rich)
setwd("/Users/banana/documents/bioinformatics/cel-counts/VP")
d <- read.csv("Attempt.csv")

str(tableRich)
D2 <- data.frame(probeID = row.names(tableRich),logFC=tableRich[,"logFC"], 
                 AveExpr=tableRich["AveExpr"], P.Value=tableRich[,"P.Value"],
                 adj.P.Val=tableRich["adj.P.Val"])
D2

D3 <- merge(D2,d[,c( "AFFY.Rat230.2.probe","Gene.name","gene.ID")],by.x = "probeID",by.y ="AFFY.Rat230.2.probe", all.x =T )
D3 
view(D3)

# C3 is CeA, B3 is BLA.... the objects are created further below

p<-ggplot(data=B3, aes(x=Gene.name , y=logFC)) +
  geom_bar(stat="identity")
p

# Horizontal bar plot
p + coord_flip()


###### 14: Pathway enrichment analysis #####

# read out csv file: 

#### CeA #####
data.fit.rich = eBayes(data.fit.con)

CeArich <- topTable(data.fit.rich, number = 30 ) 


names(data.fit.rich)
setwd("/Users/banana/documents/bioinformatics/cel-counts/Gene names")
c <- read.csv("CeA.csv")
view(c)

str(CeArich)
C2 <- data.frame(probeID = row.names(CeArich),logFC=CeArich[,"logFC"], 
                 AveExpr=CeArich["AveExpr"], P.Value=CeArich[,"P.Value"],
                 adj.P.Val=CeArich["adj.P.Val"])
C2

# Below we're linking together both data frames via probe ID
C3 <- merge(C2,c[,c( "AFFY.Rat230.2.probe","SYMBOL","Gene.stable.ID")],by.x = "probeID",by.y ="AFFY.Rat230.2.probe", all.x =T ) # all.x = T
C3 
View(C3)


# new attempt at retrieving Entrez ID 

BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)

## Retrieving CeA Entrez ID ##
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
datasets
mart <- useDataset("rnorvegicus_gene_ensembl", 
                   useMart("ensembl"))
mart
attributes = listAttributes(mart)
attributes

# removed uniprot from attributes 
#resbm = getBM(
#filters ="affy_rat230_2",
#attributes = c("affy_rat230_2","ensembl_gene_id","uniprot_gn_id","rgd_symbol", "entrezgene_id"),
#values =C3, mart = mart)
resbm = getBM(
  filters ="affy_rat230_2",
  attributes = c("affy_rat230_2","ensembl_gene_id","rgd_symbol", "entrezgene_id"),
  values =C3$probeID, mart = mart)

view(resbm)
length(unique(resbm$affy_rat230_2))
table(resbm$affy_rat230_2)

# x represents the first data from, y the second 
C4 <- merge(C3,resbm,by.x = "probeID",by.y ="affy_rat230_2", all =F )

View (C4) 
names(C4) 
# turns the log FC into normal FC

####### Gene set enrichment analysis 

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
organism = "org.Rn.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# we want the log2 fold change 
original_gene_list <- C4$logFC

# name the vector
names(original_gene_list) <- C4$entrezgene_id

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
keytypes(org.Rn.eg.db)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
gse
View(gse)

# Multiple different plots:

require(DOSE)
# enrichment map:
emapplot(gse, showCategory = 10)
# 1) Great for this data: dot plot:
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
# 2) ridgeplots
ridgeplot(gse) + labs(x = "enrichment distribution")
# 3) Great for this data:  categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 9)
# 4) Maybe for this data: 
# To view pubmed trend and see the amount of publications regarding what you're data
terms <- gse$Description[1:10]
pmcplot(terms, 2010:2018, proportion=FALSE)


##### BLA ###### Enrichment analysis 

BLArich <- topTable(data.fit.rich, number = 30) 


names(data.fit.rich)
setwd("/Users/banana/documents/bioinformatics/cel-counts/Gene names")
b <- read.csv("BLA.csv")
view(b)

str(BLArich)
B2 <- data.frame(probeID = row.names(BLArich),logFC=BLArich[,"logFC"], 
                 AveExpr=BLArich["AveExpr"], P.Value=BLArich[,"P.Value"],
                 adj.P.Val=BLArich["adj.P.Val"])
B2

# Below we're linking together both data frames via probe ID
B3 <- merge(B2,b[,c( "AFFY.Rat230.2.probe","Gene.name","Gene.stable.ID")],by.x = "probeID",by.y ="AFFY.Rat230.2.probe", all.x =T  )
B3 
View(B3)


# new attempt at retrieving Entrez ID 

# Retrieving BLA Entrez ID ##
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
datasets
mart <- useDataset("rnorvegicus_gene_ensembl", 
                   useMart("ensembl"))
mart
attributes = listAttributes(mart)
attributes


resbLA = getBM(
  filters ="affy_rat230_2",
  attributes = c("affy_rat230_2","ensembl_gene_id","rgd_symbol", "entrezgene_id"),
  values =B3$probeID, mart = mart)

view(resbLA)

# function below shows you the amount of unique ID's you data frame has
length(unique(resbLA$affy_rat230_2))

table(resbLA$affy_rat230_2)

B4 <- merge(B3,resbLA,by.x = "probeID",by.y ="affy_rat230_2", all =F )

View (B4)
names(B4) 

####### Gene set enrichment analysis 

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
organism = "org.Rn.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# we want the log2 fold change 
original_gene_BLA <- B4$logFC

# name the vector
names(original_gene_BLA) <- B4$entrezgene_id

# omit any NA values 
BLA_list<-na.omit(original_gene_BLA)

# sort the list in decreasing order (required for clusterProfiler)
BLA_list = sort(BLA_list, decreasing = TRUE)
keytypes(org.Rn.eg.db)

gse2 <- gseGO(geneList=BLA_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
gse2

#####
gse2@result$Gene.name <- paste0(str_sub(
  gse2@result$Gene.name, 1, 20),
  "...")

#####

# Multiple different plots:

require(DOSE)
# enrichment map:
emapplot(gse2, showCategory = 10)
# 1) Great for this data: dot plot:
dotplot(gse2, showCategory=20, split=".sign") + facet_grid(.~.sign)
# 2) ridgeplots
ridgeplot(gse2) + labs(x = "enrichment distribution")
# 3) Great for this data:  categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse2, categorySize="pvalue", foldChange=BLA_list, showCategory = 9)
# 4) Maybe for this data: 
# To view pubmed trend and see the amount of publications regarding what you're data
terms <- gse2$Description[1:10]
pmcplot(terms, 2010:2018, proportion=FALSE)



#if you want to look at how the genes interact with different features:

gse2 <- gseGO(geneList=BLA_list, 
              ont ="CC", 
              keyType = "ENTREZID", 
              nPerm = 10000, 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = organism, 
              pAdjustMethod = "none")

dotplot(gse2, showCategory=20, split=".sign") + facet_grid(.~.sign)





###################################################################################

# Section below seems redundant: 
data_medians <- rowMedians(Biobase::exprs(data.rma))

hist_res <- hist(data_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")

man_threshold <- 4
hist_res <- hist(data_medians, 100, col = "cornsilk", freq = FALSE, 
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)


###### Maybe - Look up why this is necessary: and if not - delete section


## All transcripts that don't have intensities greater than threshold need to be filtered out.
# First, get a list entailing the number of samples (=arrays) in the experimental groups:
no_of_samples <- 
  table(paste0(pData(data.rma)$Factor.Value.condition., "_", 
               pData(data.rma)$Factor.Value.region.))
no_of_samples

# All transcripts that don't have intensities greater than threshold are filtered out:
samples_cutoff <- min(no_of_samples)  

idx_threshold <- apply(Biobase::exprs(data.rma), 1,
                           function(x){sum(x > idx_threshold) >= samples_cutoff})
table(idx_threshold) # summarises how many genes have been filtered out 

data_filtered <- subset(data.rma,idx_threshold)

######## above section ends here ######


