library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(tibble)
library(openxlsx)
setwd("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ")
    sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_joinlayer.rds")
    #Save the R data file so we can use it later
    saveRDS(seurat_c, file='seurat_c_filter_before.rds')
    saveRDS(sc_subset_c,file="sc_subset_c_cluster_before.rds")
    saveRDS(sc_subset_c,file="sc_subset_c_cluster_after.rds")
    saveRDS(sc_subset_c,file="sc_subset_c_join_before.rds")
    #load 
    sc_subset_c <- readRDS("sc_subset_c_cluster_before.rds")
     #save.image(file="sox9_scRNAseq_QZ.Rdata")
    # load("sox9_scRNAseq_QZ_3.Rdata")
    #savehistory(file="SOX9_scRNAseq_QZ_3.Rhistory")
    #loadhistory(file="SOX9scRNAseq_QZ_3.Rhistory")

# Paths to CellRanger outputs
path1 <- "/Users/jackzhou/Desktop/Project_Sox9/Row data_CellRanger/Cornea1/filtered_feature_bc_matrix"
path2 <- "/Users/jackzhou/Desktop/Project_Sox9/Row data_CellRanger/Cornea2/filtered_feature_bc_matrix"

# Load data
data1 <- Read10X(data.dir = path1)
data2 <- Read10X(data.dir = path2)

# Create Seurat objects
seurat1 <- CreateSeuratObject(counts = data1, project = "Sample1")
seurat2 <- CreateSeuratObject(counts = data2, project = "Sample2")

#merge 2 datasets
seurat_c <- merge(seurat1, y = seurat2, project = "MergedDataset", add.cell.ids=c("s1","s2"))
seurat_c
head(seurat_c$orig.ident)

# count how many cells are in each sample using table()
table(seurat_c$orig.ident)

###################################################################################

# make a list of the mitochondrial genes
gene_names_c <- rownames(seurat_c)
mt_genes_c <- grep("^mt-", gene_names_c, value = TRUE)
length(mt_genes_c)
head(mt_genes_c)

# calculate mt percentage
seurat_c[["percent.MT"]] <- PercentageFeatureSet(seurat_c, pattern = "^mt-")
#Review mc%
summary(seurat_c[["percent.MT"]])

# make plots of the distribution of counts, genes, and mitochondrial expression
VlnPlot(seurat_c, features = c("nFeature_RNA","nCount_RNA","percent.MT"),pt.size=0.2)
# make scatterplots of one value vs the other
# number of UMI counts vs number of genes
FeatureScatter(seurat_c, feature1="nCount_RNA",feature2="nFeature_RNA")
# number of UMI counts vs % mitochondrial expression
FeatureScatter(seurat_c, feature1="nCount_RNA",feature2="percent.MT")
# number of genes vs % mitochondrial expression
FeatureScatter(seurat_c, feature1="nFeature_RNA",feature2="percent.MT")

#save
saveRDS(seurat_c, file = 'c_seurat.rds')
seurat_c <- readRDS('c_seurat.rds')
####################################################################################
# Filter cells for downstream analysis
#For today, our passing cells must have:
# ??? Number of genes (nFeature) per cell greater than 2000 and less than 5000
#??? UMI counts (nCount) per cell greater than 2000
#??? Less than 10% mitochondrial content
#seurat <- subset(seurat, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA>2000)

#works with 50% EGFP
sc_subset_c  <- subset(seurat_c, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & nCount_RNA>2000 & percent.MT < 10)

#filter criteria of CORE
#sc_subset_c <- subset(seurat_c, subset = nFeature_RNA < 8000  &nFeature_RNA > 200 & percent.MT < 20)

#test
#sc_subset_c <- subset(seurat_c, subset = nFeature_RNA < 8000  &nFeature_RNA > 200 & percent.MT < 10)

# check what fraction of cells we kept from each sample, using the table function
orig.counts <- table(seurat_c$orig.ident)
subset.counts <- table(sc_subset_c$orig.ident)
cell_stats <- cbind(orig.counts, subset.counts, subset.counts/orig.counts)
colnames(cell_stats) <- c("Starting Cells","Retained Cells","Fraction")
cell_stats

#save
saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_filted.rds")
#sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_filted.rds")

###########################################################################################
#core criteria: Log normalized; scale.factor = 10000
# normalization
sc_subset_c <- NormalizeData(sc_subset_c, normalization.method = "LogNormalize", scale.factor = 10000)#core
#sc_subset_c <- NormalizeData(sc_subset_c, scale.factor = 10000)

#save
saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_Normalized.rds")
#sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_Normalized.rds")
###########################################################################################
# Find variable genes. In practice, you may want to vary the number of features

#core criteria: FindVariableFeatures: nfeatures = 2000, selection.method = 'mvp'

# you select here based on the plot.
#core
#sc_subset_c <- FindVariableFeatures(sc_subset_c,nfeatures=2000, selection.method = 'mvp')
#50%
sc_subset_c <- FindVariableFeatures(sc_subset_c,nfeatures=9000, selection.method = 'mvp')#the top 9,000 variable genes 
VariableFeaturePlot(sc_subset_c)

#save
saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_selected.rds")
#sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_selected.rds")

###########################################################################################
# z-score
all.genes_c <- rownames(sc_subset_c)
sc_subset_c <- ScaleData(sc_subset_c)

###########################################################################################
#Centering and scaling data matrix;
sc_subset_c <- RunPCA(sc_subset_c,features=VariableFeatures(object=sc_subset_c),
                    npcs = 50,verbose=F)
# plots
ElbowPlot(sc_subset_c,ndims=24)
DimHeatmap(sc_subset_c,dims=1:24,cells=300,balanced=T)
###########################################################################################

#save
saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_PCA.rds")
#sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_PCA.rds")

###########################################################################################
#Dimensional reduction
#UMAP method


###########################################################################################

# Clustering
# Run clustering on the top PCs to filter out noise
#pca.dims1 = 1:10#core
#pca.dims1 = 1:16#50%
pca.dims1 = 1:12#13 cluster

sc_subset_c <- FindNeighbors(sc_subset_c, dims=pca.dims1)
#save
saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_PCA2.rds")
#sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_PCA2.rds")

##########################################################################################
#To optimize resolution by silhouette

##########################################################################################
# run clutering, with resolution 
sc_subset_c <- FindClusters(sc_subset_c, algorithm=2,  resolution=0.5)


##########################################################################################
# check the cluster silhouette needed

#save
saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_cluster.rds")
#sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_cluster.rds")

##########################################################################################
# Look at cluster IDs of the first 5 cells
head(Idents(sc_subset_c ), 2)
###########################################################################################
# tSNE plot
sc_subset_c <- RunTSNE(sc_subset_c, dims=pca.dims1)
DimPlot(sc_subset_c, reduction='tsne')
DimPlot(sc_subset_c, reduction='tsne',label=TRUE, label.size = 5)

# we can also make the tSNE plot with respect to the orignal identity
DimPlot(sc_subset_c, reduction='tsne', group.by="orig.ident")

# size of clusters
table(sc_subset_c$seurat_clusters)

# number of samples in each cluster
cluster_sample_c <- data.frame(cluster=sc_subset_c$seurat_clusters,
                               group=sc_subset_c$orig.ident)
cluster_abundance_c <- table(cluster_sample_c)
cluster_abundance_c

# see counts as percent of total
cluster_percent_c <- t(t(cluster_abundance_c)/colSums(cluster_abundance_c))
cluster_percent_c

# UMAP plot
sc_subset_c <- RunUMAP(sc_subset_c, dims = pca.dims1)

png(filename="c_umap2025.png")
DimPlot(sc_subset_c, reduction = "umap", label=TRUE, label.size = 5)
dev.off()  # Close the device to save the file
######################################################################################################################################################################################
  #save
#saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_UMAP.rds")
sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_UMAP.rds")
######################################################################################################################################################################################


##################################################################################################################
###########Differential expression analysis(DEG) on EGFP-bGhpoly expressed cells############
# Joining layers if they're not integrated
#important for SeuratObject 5.0.0. 
sc_subset_c<- JoinLayers(sc_subset_c, features = rownames(sc_subset_ccluster), assay = 'RNA') # Or another relevant assay
#saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_joinlayer.rds")
sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_joinlayer.rds")
##################################################################################################################


# Add expression information for "EGFP-bGhpolyA"
sc_subset_c$EGFP_bGhpolyA_expr <- GetAssayData(object = sc_subset_c, assay = "RNA", slot = "counts")["EGFP-bGhpolyA", ] > 0
#sc_subset_c_join$EGFP_bGhpolyA_expr <- GetAssayData(object = sc_subset_c_join, assay = "RNA", slot = "counts")["EGFP-bGhpolyA", ] > 0
sc_subset_c$EGFP_bGhpolyA_expr <- GetAssayData(object = sc_subset_c, assay = "RNA", layer = "counts")["EGFP-bGhpolyA", ] > 0

###########################################################################################
# check cells expressed with "EGFP-bGhpolyA"
# Subset cells from a specific sample group,
Idents(sc_subset_c) = "orig.ident"
sc_subset_c1 <- subset(sc_subset_c, subset = orig.ident == "Sample1")
sc_subset_c2 <- subset(sc_subset_c, subset = orig.ident == "Sample2")

counts_c<-GetAssayData(object = sc_subset_c ,  assay = "RNA", slot = "counts")
sum(counts_c["EGFP-bGhpolyA",]>0)

counts_c1<-GetAssayData(object = sc_subset_c1, assay = "RNA", slot = "counts")
counts_c2<-GetAssayData(object = sc_subset_c2, assay = "RNA", slot = "counts")

# write cluster matrix and gene expression matrix 

#write.csv(counts_c, file = "counts_c.csv")
#write.csv(counts_c1, file = "counts_c1.csv")
#write.csv(counts_c2, file = "counts_c2.csv")

sum(counts_c1["EGFP-bGhpolyA",]>0)
sum(counts_c2["EGFP-bGhpolyA",]>0)

# Count cells expressing "EGFP-bGhpolyA" and total cells per cluster
cluster_counts_c_EGFP <- sc_subset_c@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    EGFP_bGhpolyA_count = sum(EGFP_bGhpolyA_expr),
    total_count = n(),
    EGFP_bGhpolyA_percentage = 100 * EGFP_bGhpolyA_count / total_count,
    no_EGFP_bGhpolyA_percentage = 100 - EGFP_bGhpolyA_percentage  # Cells without EGFP-bGhpolyA expression
  )
cluster_counts_c_EGFP 

#save
saveRDS(sc_subset_c, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_join.rds")
#sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_join.rds")

# View the plot
barplot(cluster_counts_c_EGFP$EGFP_bGhpolyA_percentage )
barplot(cluster_counts_c_EGFP$EGFP_bGhpolyA_count )
barplot(cluster_counts_c_EGFP$total_count )

# View the result
write.xlsx(cluster_counts_c_EGFP,  "c_cluster_counts_EGFP.xlsx")
write.csv(cluster_counts_c_EGFP, file = "c_cluster_counts_EGFP.csv")

# Count cells expressing "EGFP-bGhpolyA", "orig.ident" and total cells per cluster
cluster_counts_c_EGFP2 <- sc_subset_c@meta.data %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(
    EGFP_bGhpolyA_count = sum(EGFP_bGhpolyA_expr),
    total_count = n(),
    EGFP_bGhpolyA_percentage = 100 * EGFP_bGhpolyA_count / total_count,
    no_EGFP_bGhpolyA_percentage = 100 - EGFP_bGhpolyA_percentage  # Cells without EGFP-bGhpolyA expression
  )
# View the plot
barplot(cluster_counts_c_EGFP2$EGFP_bGhpolyA_percentage )
barplot(cluster_counts_c_EGFP2$EGFP_bGhpolyA_count )
barplot(cluster_counts_c_EGFP2$total_count )

write.xlsx(cluster_counts_c_EGFP2,  "c_cluster_counts_EGFP2.xlsx")



###########################################################################################

#Differential expression with respect to clusters
# Joining layers if they're not integrated
sc_subset_c<- JoinLayers(sc_subset_c, features = rownames(sc_subset_ccluster), assay = 'RNA') # Or another relevant assay

#Differential expression with respect to clusters
markers_c <- FindAllMarkers(sc_subset_c, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_c, file = "markers_c2025.csv")

##########################################################################################
#the average gene expression in each clusters
##########################################################################################

avg_exp <- AverageExpression(sc_subset_c) #AggregateExpression
head(avg_exp$RNA)  # Check the first few rows of the average expression matrix for RNA assay
sc_subset_c <- subset(sc_subset_c, features = rownames(sc_subset_c)[rowSums(GetAssayData(sc_subset_c, slot = "counts")) > 0])
avg_exp <- AverageExpression(sc_subset_c)
head(avg_exp$RNA)
write.csv(avg_exp, file = "c_avg_exp.csv")





################################Manipulate clusters##########################################
#################Differential expression analysis(DEG) on selected cluster####################
# if you wanted to test between one cluster vs the combination of two other clusters:
##########################################################################
#wilcox
######################################################################

Idents(sc_subset_c) <- "seurat_clusters" 
#run stats with Wilcox instead of ROC
wilcox_stats_c <- FindAllMarkers(sc_subset_c, test.use="wilcox",group.by="seurat_clusters",
                                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )  
write.csv(wilcox_stats_c, file = "c_wilcox_stats.csv")

# find all markers of cluster 12
cluster12_markers_c <- FindMarkers(sc_subset_c, test.use="wilcox", ident.1 = 12, 
                                   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
write.csv(cluster12_markers_c , file = "c_cluster12_markers.csv")

#LSCs vs Diff & TA
stem_cluster <- "12"  # replace with the actual stem cell cluster ID
diff_cluster <-  c(2,3,4,5,7,10) # replace with differentiated cluster IDs
Idents(sc_subset_c) <- "seurat_clusters" 
c_clusterStemVsDif <- FindMarkers(sc_subset_c, ident.1 = stem_cluster , ident.2 = diff_cluster , test.use="wilcox")
write.csv(c_clusterStemVsDif, file = "c_clusterStemVsDif.csv")


#LSCs vs TA
stem_cluster <- "12"  # replace with the actual stem cell cluster ID
TA_cluster <-  c(5,7) # replace with TA cluster IDs
Idents(sc_subset_c) <- "seurat_clusters" 
c_clusterStemVsTA_wil<- FindMarkers(sc_subset_c, ident.1 = stem_cluster , ident.2 = TA_cluster , test.use="wilcox",
                                  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(c_clusterStemVsTA_wil, file = "c_clusterStemVsTA_wil.xlsx", rowNames=T)
c_clusterStemVsTA_wil <- read_excel("c_clusterStemVsTA_wil.xlsx",  col_names = TRUE)

c_clusterStemVsTA_roc<- FindMarkers(sc_subset_c, ident.1 = stem_cluster , ident.2 = TA_cluster , test.use="roc",
                                only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(c_clusterStemVsTA_roc, file = "c_clusterStemVsTA_roc.xlsx", rowNames=T)
c_clusterStemVsTA_roc <- read_excel("c_clusterStemVsTA_roc.xlsx",  col_names = TRUE)

c_cluster12vsTA_wil_roc_join <- inner_join(c_clusterStemVsTA_wil, c_clusterStemVsTA_roc, by="Name")
write.xlsx(c_cluster12vsTA_wil_roc_join, file = "c_cluster12vsTA_wil_roc_join.xlsx",  rowNames=T)

head()
c_12vsTA_list <- read.table(file = 'c_12vsTA_list.txt')
head(c_12vsTA_list)
c_12vsTA_list <- toupper(c_12vsTA_list$V1)
head(c_12vsTA_list)
write.table(c_12vsTA_list, file = "c_12vsTA_list_capital.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)








#########################################################################################
# Set the new cluster identities based on EGFP-bGhpolyA expression
Idents(sc_subset_c) <- sc_subset_c$EGFP_bGhpolyA_expr
# Perform DEA between cells expressing "EGFP-bGhpolyA" and those that do not
deg_results_EGFP_c <- FindMarkers(sc_subset_c, ident.1 = TRUE, ident.2 = FALSE,
                                  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                  test.use="wilcox")
# Save DEG results
write.csv(deg_results_EGFP_c, file = "c_DEG_EGFP_bGhpolyA_vs_Other.csv", row.names = TRUE)


##########################################################################
#AUC
######################################################################
#pct (Percentage of expressing cells in group )
# The p-value adjusted for multiple testing using the Benjamini-Hochberg (BH) correction to control the false discovery rate (FDR).
#AUC > 0.8: Strong marker gene for distinguishing between clusters.
#Statistical power : A measure of the statistical power of the test for each gene. 


#Receiver Operating Characteristic (ROC) Test (`test.use = “roc”)
Idents(sc_subset_c) <- "seurat_clusters" 
roc_c <- FindAllMarkers(sc_subset_c, test.use="roc", group.by="seurat_clusters",
                        only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
write.csv(roc_c, file = "c_roc.csv")
#roc_c <- read.csv(file = "roc_c2025.csv")
# number of genes with AUC>0.7 for each cluster
roc8_c <- table(roc_c[roc_c$myAUC>0.8,"cluster"])
write.csv(roc8_c, file = "c_roc8.csv")

# get top 5 genes for cluster 5
cluster12_c = roc_c[roc_c$cluster==12,]
cluster12_top5_c <- head(cluster12_c[order(cluster12_c[,1],decreasing=T),],5)
# getting top 10 genes for all clusters
library(dplyr)
top10_c <- roc_c %>% group_by(cluster) %>% top_n(n=10, wt=myAUC)
top10_c


# find all markers of cluster 12
cluster12_roc_c <- FindMarkers(sc_subset_c, test.use="roc", ident.1 = 12, 
                                   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
write.csv(cluster12_roc_c , file = "c_cluster12_markers_roc.csv")

#LSCs vs Diff & TA
c_roc_12vs2345710 <- FindMarkers(sc_subset_c, ident.1 = 12, ident.2 = c(2,3,4,5,7,10) , test.use="roc",
                                    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
write.csv(c_roc_12vs2345710, file = "c_cluster12vs2345710_roc.csv")

# Set the new cluster identities based on EGFP-bGhpolyA expression
Idents(sc_subset_c) <- sc_subset_c$EGFP_bGhpolyA_expr
# Perform DEA between cells expressing "EGFP-bGhpolyA" and those that do not
deg_roc_EGFP_c <- FindMarkers(sc_subset_c, ident.1 = TRUE, ident.2 = FALSE,
                                  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                                  test.use="roc")
# Save DEG results
write.csv(deg_roc_EGFP_c, file = "c_DEG_EGFP_vs_Other_roc.csv", row.names = TRUE)



##########################################################################################
#cell cycle
##########################################################################################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#2. Force Lowercase/Uppercase Matching
rownames(sc_subset_c) <- tolower(rownames(sc_subset_c))
s.genes <- tolower(cc.genes$s.genes)
g2m.genes <- tolower(cc.genes$g2m.genes)
sc_subset_c <- CellCycleScoring(sc_subset_c, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#2. Remove Missing Genes from Lists
s.genes <- s.genes[s.genes %in% rownames(sc_subset_c)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(sc_subset_c)]

#Assign Cell-Cycle Scores
sc_subset_c <- CellCycleScoring(sc_subset_c, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
sc_subset_c <- RunPCA(sc_subset_c, features = c(s.genes, g2m.genes))

# Create the table for cell cycle phases across seurat clusters
cycle_table_c <- table(sc_subset_c@meta.data$Phase,sc_subset_c@meta.data$seurat_clusters )
write.csv(cycle_table_c, file="c_cycle_table2025.csv", row.names = TRUE)

# Convert the table to a data frame for easier manipulation
cycle_table_c_df <- as.data.frame.matrix(cycle_table_c)

# Calculate the percentage of cells in each phase within each cluster
cycle_percentage_c <- sweep(cycle_table_c_df, 2, colSums(cycle_table_c_df), FUN = "/") * 100

# Write the percentage table to a CSV file
write.csv(cycle_percentage_c, file = "c_cycle_percentage_table2025.csv", row.names = TRUE)

# Install and load the pheatmap package (if not already installed)
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)
# Create a custom annotation for the x-axis (clusters)
cluster_annotation <- data.frame(
  Cluster = colnames(cycle_percentage_c) # Cluster names
)
rownames(cluster_annotation) <- colnames(cycle_percentage_c) # Set rownames to match heatmap column names

# Generate the heatmap
png(filename ="cycle_percentage_c_heatmap.png" )
cycle_percentage_c_heatmap <- pheatmap(
  mat = cycle_percentage_c,                # Data matrix
  color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient from blue (low) to red (high)
  cluster_rows = FALSE,                    # Disable clustering for rows (phases)
  cluster_cols = FALSE,                    # Disable clustering for columns (clusters)
  display_numbers = TRUE,                  # Show the percentage values on the heatmap
  fontsize_number = 10,                    # Font size for displayed numbers
  main = "Cell Cycle Phases Across Clusters", # Title of the heatmap
  fontsize_row = 12,                       # Font size for row labels
  fontsize_col = 12,                        # Font size for column labels
  annotation_col = cluster_annotation,      # Add the cluster annotation to the x-axis
  angle_col = 45                           # Rotate cluster numbers by 45 degrees
)
dev.off()



##########################################################################################
#Pseudotime
##########################################################################################
#Monocle, Velocyto , CytoTRACE, SlingShot.

library(monocle3)
# source the R script
source("https://uic-ric.github.io/workshop-data/scrna/importCDS2.R")
# use the importCDS2 function (our function) instead of importCDS (monocle function)
monocle_c <- importCDS2(sc_subset_c)
monocle_c

#Estimate size factors and pick which genes to use
monocle_c <- estimateSizeFactors(monocle_c)
monocle_c <- setOrderingFilter(monocle_c, sc_subset_c@assays$RNA@var.features)
#Reduce dimensions - this step may take awhile
monocle_c <- reduceDimension(monocle_c, max_components=2, method='DDRTree') #spend 2 days
#Infer the cell ordering (pseudotime). 
monocle_c <- orderCells(monocle_c)

#Try different visualizations
png(filename="monocle_c.png")
plot_cell_trajectory(monocle_c)
dev.off()
saveRDS(monocle_c, file="monocle_c.rds")

plot_cell_trajectory(monocle_data, color_by = "Pseudotime")
#Color by seurat cluster
png(filename="wh_monocle_RNA_snn_res.png")
plot_cell_trajectory(monocle_data, color_by = "RNA_snn_res.0.6")
#RNA_snn_res.0.5 is present in sc_subset_wh_join@meta.data, it will be transferred to monocle_data.
dev.off()

head(pData(monocle_c))
# Get the top marker gene for cluster 12
roc_c <- readRDS("roc_c.rds")
topgene = rownames(roc_c[roc_c$cluster==12,])[1]
# plot using gene expression as size of dot
png(filename="monocle_c12_g1_c.png")
plot_cell_trajectory(monocle_c, color_by = "RNA_snn_res.0.5", markers = topgene)
dev.off()
#Can plot pseudotime vs gene expression using color
png(filename="monocle_c12_g1_expression_c.png")
plot_cell_trajectory(monocle_c, markers = topgene, use_color_gradient=T)
dev.off
head(t(monocle_c@reducedDimS))



##########################################################################################
#Harmony Integration #https://satijalab.org/seurat/reference/harmonyintegration
#https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/quickstart.html
##########################################################################################




