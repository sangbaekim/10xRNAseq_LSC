#--------------------------------------------------------------------------
# NAME: scRNAseq_analysis_LSC.R
#
# INPUT FILES:
# - Read count file & R source
#
# User's input parameters:
# - dirIn1,2: dir for filtered_feature_bc_matrix folders
# - dirOut  : output folder for saving processed files
#
#--------------------------------------------------------------------------

#------------------------------------------------------------
# User's input parameters
#------------------------------------------------------------
dirIn1  <- "1st/filtered_feature_bc_matrix"
dirIn2  <- "2nd/filtered_feature_bc_matrix"
dirOut  <- "output"
#------------------------------------------------------------
if(!file.exists(dirOut))
  dir.create(dirOut)
setwd(dirOut) 
dirOut <- getwd()


# Install R packages
install.packages(c("Seurat","dplyr","topGo","tidyverse","org.Hs.eg.db","org.Mm.eg.db"))

# load library
library(Seurat)
library(dplyr)
library(Matrix)
library(topGO)
library(tidyverse)
library(org.Hs.eg.db)


# load data 
objAdataR <- Read10X(data.dir =dirIn)
objBdataR <- Read10X(data.dir =dirIn)

## 1) 1st data process
# numb of cell with expressed read counts(>0) across genes 
cut_cell1 <- round(dim(objAdataR)[2] * 0.03) # filtering of cell: 3% of total cells

# creat object
objAdata <- CreateSeuratObject(objAdataR, min.cells = cut_cell1)
objAdata

#use 15% and 85% quantiles to filter cells
qt1<-as.integer(quantile(objAdata@meta.data$nGene, c(0.05, 0.95)))

# QC & normalization
objAdata <-FilterCells(object = objAdata, subset.names = "nGene", low.thresholds = qt1[1], high.thresholds = qt1[2])
objAdata <-NormalizeData(objAdata)
objAdata <-ScaleData(objAdata)
objAdata <-FindVariableGenes(objAdata)
objAdata@meta.data$subset<-"1st"


# 2) 2nd data process
# numb of cell with expressed read counts(>0) across genes 
cut_cell2 <- round(dim(objBdataR)[2] * 0.03) # filtering of cell: 3% of total cells

# creat object
objBdata <- CreateSeuratObject(objBdataR, min.cells = cut_cell2)
objBdata

#use 15% and 85% quantiles to filter cells
qt2<-as.integer(quantile(objBdata@meta.data$nGene, c(0.05, 0.95)))

# QC & normalization
objBdata <-FilterCells(object = objBdata, subset.names = "nGene", low.thresholds = qt1[1], high.thresholds = qt1[2])
objBdata <-NormalizeData(objBdata)
objBdata <-ScaleData(objBdata)
objBdata <-FindVariableGenes(objBdata)
objBdata@meta.data$subset<-"1st"


# 3) data integration
numb_varGenes <- 1000
g.1 <- head(rownames(objAdata@hvg.info), numb_varGenes)
g.2 <- head(rownames(objBdata@hvg.info), numb_varGenes)

genes.use1 <- unique(c(g.1,g.2))
genes.use <- intersect(genes.use1, rownames(objAdata@scale.data))
genes.use <- intersect(genes.use1, rownames(objBdata@scale.data))

# Run CCA
obj.integrated <- RunCCA(object=objAdata, object2=objBdata, genes.use = genes.use, num.cc = 20)
obj.integrated <- AlignSubspace(obj.integrated, reduction.type = "cca", grouping.var = "subset",dims.align = 1:20)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
setwd(dirOut)
pdf("plot_feature.pdf")
DimPlot(object = obj.integrated, reduction.use = "cca", group.by = "subset", pt.size = 0.5, do.return = TRUE)
VlnPlot(object = obj.integrated, features.plot = "CC1", group.by = "subset", do.return = TRUE)
dev.off()



# 4) Perform data processing for integrated data
# UMAP, t-SNE and Clustering
obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca", dims = 1:20)
obj.integrated <- RunTSNE(obj.integrated, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
obj.integrated <- FindClusters(obj.integrated, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

# Visualization
pdf("6_plot_t-SNE_cc20.pdf")
# umap
pdf("plot_umap_integrated.pdf")
DimPlot(object = obj.integrated, reduction = "umap",label=T,label.size = 10)
DimPlot(object = obj.integrated, reduction = "umap",label=F, group.by = 'subset')
DimPlot(object = obj.integrated, reduction = "umap",label=T,label.size = 5, split.by = 'subset')
dev.off()

# 5) DEG : Find markers for each cluster.
DefaultAssay(object = obj.integrated) <- "RNA"
integ.markers <- FindAllMarkers(object = obj.integrated, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)

# select top 20, 50 genes of each cluster
iTop5  <- integ.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
iTop10 <- integ.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
iTop20 <- integ.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
iTop50 <- integ.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)

# Heatmap
pdf("plot_heatmap_top10.pdf", height = 25)
DoHeatmap(object = subset(obj.integrated,downsample=500), features = iTop10$gene) + NoLegend()
dev.off()

# write table - top50 markers
setwd(dirOut)
write.csv(iTop50, "out_table_top50markers_eachCluster.csv",row.names= F)





#----------------------------------------------------------------------------
# cell cycle analysis
#----------------------------------------------------------------------------

# read cellcycle genes
cc.genes  <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes   <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Assign Cell-Cycle Scores
obj.integrated <- CellCycleScoring(object = obj.integrated, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = F)

# view t-sne plot
setwd(dirOut)
pdf("plot_umap_cellcycle.pdf")
DimPlot(object = obj.integrated, reduction="umap",label=T,label.size = 10)
DimPlot(object = obj.integrated, reduction="umap",group.by='Phase',label=T,label.size = 10)
dev.off()

# heatmap -summary table (%)
library(pheatmap)
prop2  <- table(Idents(objAdata), objAdata@meta.data$Phase)
prop2A <- round(prop2/rowSums(prop2)*100)

setwd(dirOut)
pdf("heatmap_cellcycle.pdf", width = 4)
print(pheatmap(prop2,display_numbers = TRUE,cluster_rows = F,cluster_cols = F,fontsize = 15))
print(pheatmap(prop2A,display_numbers = TRUE,cluster_rows = F,cluster_cols = F,fontsize = 15))
dev.off()




#==================================================
### SCenic 
#==================================================
library(SCENIC)
library(Seurat)
library(dplyr)
library(Matrix)

# subset of 9 clusters from whole data
LSC <- subset(obj.integrated, idents=c(0:6,8,10))
LSC@active.assay <-'RNA'
LSC_norm <- NormalizeData(LSC, normalization.method = "LogNormalize", scale.factor = 10000)
LSC_norm$cluster <- LSC_norm@active.ident
cidents <- c("TDC_C0","PMC_C1","TDC_C2","PMC_C3","PMC_C4",'LPC_C5','LPC_C6','TAC_C8','LSC_C10')
table(LSC_norm@active.ident)

# Generate subset matrix
for (i in 1:8) {
  #i=9
  cid <- cidents[i]
  LSC_s1 <- subset(LSC_norm, idents = cid)
  LSC_s1_matrix <- as.matrix(GetAssayData(LSC_s1,slot = "counts"))
  LSC_s1_matrix <- as.data.frame(LSC_s1_matrix)
  sets_n = floor(length(colnames(LSC_s1))/20)
  LSC_s1_matrix<-t(LSC_s1_matrix)
  LSC_s1_matrix <- as.data.frame(LSC_s1_matrix)
  V<-rep(1:sets_n, each=20)
  set.seed(001) # just to make it reproducible
  V<-sample(V)
  LSC_s1_matrix_split<-split(LSC_s1_matrix,V)
  round(0.11,0)
  List<-list()
  for (j in 1:sets_n){
    #      normF<-colMeans(LSC_s1_matrix_split[[j]])
    normF<-round(colMeans(LSC_s1_matrix_split[[j]]),0)
    List[[j]] <- normF
  }
  
  LSC_s1_mean <- do.call(rbind, List)
  LSC_s1_mean <- t(LSC_s1_mean)
  LSC_s1_mean <- as.data.frame(LSC_s1_mean)
  colnames(LSC_s1_mean) <- paste0(cid,'_',colnames(LSC_s1_mean))
  head(colnames(LSC_s1_mean))
  fout <- paste0("LSC_",cid,".csv")
  write.csv(LSC_s1_mean, fout)
  print(dim(LSC_s1_mean))
  
  rm(LSC_s1)
  rm(LSC_s1_mean)
}

# use for LSC10
i=9
cid <- cidents[i]
LSC_s1 <- subset(LSC_norm, idents = cid)
#LSC_s1_matrix <- as.matrix(LSC_s1[["RNA"]]@data)
LSC_s1_matrix <- as.matrix(GetAssayData(LSC_s1,slot = "counts"))
LSC_s1_matrix <- as.data.frame(LSC_s1_matrix)
colnames(LSC_s1_matrix) <- paste0(cid,'_',colnames(LSC_s1_matrix))
head(colnames(LSC_s1_matrix))
fout <- paste0("LSC_",cid,".csv")
write.csv(LSC_s1_matrix, fout)
print(dim(LSC_s1_matrix))

















