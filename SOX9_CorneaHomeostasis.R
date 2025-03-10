library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(tibble)
library(openxlsx)
library( readxl)
setwd("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ")
sc_subset_c <- readRDS("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/c_sc_subset_joinlayer.rds")

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
# Check the expression information for "EGFP-bGhpolyA" in original data.
###############################################################################
# Joining layers if they're not integrated
seurat_c_EGFP<- JoinLayers(seurat_c, features = rownames(seurat_ccluster), assay = 'RNA') # Or another relevant assay

# Add expression information for "EGFP-bGhpolyA"
seurat_c_EGFP$EGFP_bGhpolyA_expr <- GetAssayData(object = seurat_c_EGFP, assay = "RNA", slot = "counts")["EGFP-bGhpolyA", ] > 0

# Subset cells from a specific sample group,
Idents(seurat_c_EGFP) = "orig.ident"
seurat_c1 <- subset(seurat_c_EGFP, subset = orig.ident == "Sample1")
seurat_c2 <- subset(seurat_c_EGFP, subset = orig.ident == "Sample2")

counts_c<-GetAssayData(object = seurat_c_EGFP ,  assay = "RNA", slot = "counts")
counts_c1<-GetAssayData(object = seurat_c1, assay = "RNA", slot = "counts")
counts_c2<-GetAssayData(object = seurat_c2, assay = "RNA", slot = "counts")

sum(counts_c["EGFP-bGhpolyA",]>0)
sum(counts_c1["EGFP-bGhpolyA",]>0)
sum(counts_c2["EGFP-bGhpolyA",]>0)

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

####################################################################################
# Filter cells for downstream analysis
sc_subset_c  <- subset(seurat_c, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & nCount_RNA>2000 & percent.MT < 10)

# check what fraction of cells we kept from each sample, using the table function
orig.counts <- table(seurat_c$orig.ident)
subset.counts <- table(sc_subset_c$orig.ident)
cell_stats <- cbind(orig.counts, subset.counts, subset.counts/orig.counts)
colnames(cell_stats) <- c("Starting Cells","Retained Cells","Fraction")
cell_stats

###########################################################################################
# normalization
sc_subset_c <- NormalizeData(sc_subset_c, normalization.method = "LogNormalize", scale.factor = 10000)

###########################################################################################
# Find variable genes. choose the top 9,000 variable genes (or 2000)
sc_subset_c <- FindVariableFeatures(sc_subset_c,nfeatures=9000, selection.method = 'mvp') 
VariableFeaturePlot(sc_subset_c)

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
#Dimensional reduction
###########################################################################################
# Clustering
# Run clustering on the top PCs to filter out noise
pca.dims1 = 1:12#13 cluster
sc_subset_c <- FindNeighbors(sc_subset_c, dims=pca.dims1)

##########################################################################################
#To optimize resolution by silhouette and clusteree      
##########################################################################################
# Load required libraries
library(Seurat)
library(clustree)

#To find the optimal number of clusters using the clustree package (Zappia and Oshlack 2018)
#Step 1: Perform Clustering at Multiple Resolutions
resolutions <- seq(0.2, 1.0, by = 0.1) # Define the range of resolutions
for (res in resolutions) {
  sc_subset_c <- FindClusters(sc_subset_c, resolution = res)
}
clustree(sc_subset_c, prefix = "RNA_snn_res.")

#Step 2: Create a Distance Matrix for Silhouette Analysis
# Extract PCA coordinates
pca_coords <- Embeddings(sc_subset_c, reduction = "pca")

# Calculate the distance matrix (e.g., Euclidean distance)
distance_matrix <- dist(pca_coords)

#Step 3: Calculate Silhouette Widths for Each Resolution
library(cluster)  # For silhouette function

# Initialize a vector to store average silhouette widths
asw_values <- c()

# Loop through resolutions and calculate silhouette widths
resolutions <- seq(0.2, 1.0, by = 0.1)
for (res in resolutions) {
  # Get cluster assignments
  cluster_assignments <- as.numeric(as.character(Idents(sc_subset_c)))
  
  # Calculate silhouette widths
  sil <- silhouette(cluster_assignments, dist = distance_matrix)
  
  # Store the average silhouette width for this resolution
  asw_values <- c(asw_values, mean(sil[, 3]))
}

# Plot silhouette widths across resolutions
plot(resolutions, asw_values, type = "b", xlab = "Resolution", ylab = "Average Silhouette Width",
     main = "Silhouette Width vs. Resolution")
asw_values

##########################################################################################
# run clutering, with resolution 
sc_subset_c <- FindClusters(sc_subset_c, algorithm=2,  resolution=0.5)

##########################################################################################
##########################################################################################
#to use Library of factoextra to assess average silhouette width of each cluste
# Install and load libraries
#install.packages("factoextra")
#install.packages("cluster")
library(factoextra)
library(cluster)

# Get cluster assignments as numeric values
cluster_assignments <- as.numeric(as.character(Idents(sc_subset_c)))

# Compute PCA distance matrix
pca_coords <- Embeddings(sc_subset_c, reduction = "pca")
distance_matrix <- dist(pca_coords)

# Perform silhouette analysis
sil <- silhouette(cluster_assignments, dist = distance_matrix)

# Visualize silhouette plot
fviz_silhouette(sil)

# Extract average silhouette widths
silhouette_summary <- summary(sil)
cluster_avg_widths <- silhouette_summary$clus.avg.widths  # Average silhouette width per cluster
overall_avg_width <- silhouette_summary$avg.width         # Overall average silhouette width

# Print the results
print(cluster_avg_widths)  # Average silhouette width for each cluster
print(overall_avg_width)   # Overall average silhouette width

##########################################################################################
# Look at cluster IDs of the first 5 cells
head(Idents(sc_subset_c ), 5)
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

#########################################################################################
# UMAP plot
sc_subset_c <- RunUMAP(sc_subset_c, dims = pca.dims1)

png(filename="c_umap2025.png")
DimPlot(sc_subset_c, reduction = "umap", label=TRUE, label.size = 5)
dev.off()  # Close the device to save the file
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

#in order to save it to excel with gene name.
df <- as.data.frame(avg_exp$RNA)
df$Gene <- rownames(df)  # Add gene names as a column
write.xlsx(df, file = "c_avg_exp.xlsx", rowNames = TRUE)  

################################Manipulate clusters##########################################
#################Differential expression analysis(DEG) on selected cluster####################
# if you wanted to test between one cluster vs the combination of two other clusters:
##########################################################################
#wilcox
######################################################################
#AUC
######################################################################
#pct (Percentage of expressing cells in group )
# The p-value adjusted for multiple testing using the Benjamini-Hochberg (BH) correction to control the false discovery rate (FDR).
# marker gene for distinguishing between clusters: AUC > 0.7: Acceptable; AUC > 0.8: Strong
#Statistical power : A measure of the statistical power of the test for each gene. 
#Receiver Operating Characteristic (ROC) Test (`test.use = “roc”)
##########################################################################
Idents(sc_subset_c) <- "seurat_clusters" 
#run stats with Wilcox and ROC

###################################################################################################
# find the DEG markers of LSCs-EGFP vs TA
###################################################################################################
stem_cluster <- "12"  # replace with the actual stem cell cluster ID
TA_cluster <-  c(5,7) # replace with TA cluster IDs
Idents(sc_subset_c) <- "seurat_clusters" 
#Wilcox
c_clusterStemVsTA_wil<- FindMarkers(sc_subset_c, ident.1 = stem_cluster , ident.2 = TA_cluster , test.use="wilcox",
                                    min.pct = 0.25, logfc.threshold = 0.25) #only.pos = TRUE 
head(c_clusterStemVsTA_wil)
write.xlsx(c_clusterStemVsTA_wil, file = "c_clusterStemVsTA_wil.xlsx", rowNames=T)
#Roc
c_clusterStemVsTA_roc<- FindMarkers(sc_subset_c, ident.1 = stem_cluster , ident.2 = TA_cluster , test.use="roc",
                                    min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(c_clusterStemVsTA_roc, file = "c_clusterStemVsTA_roc.xlsx", rowNames=T)

#adding the column name with 'Gene', then read the file 
c_clusterStemVsTA_wil <- read_excel("c_clusterStemVsTA_wil.xlsx",  col_names = TRUE)
c_clusterStemVsTA_roc <- read_excel("c_clusterStemVsTA_roc.xlsx",  col_names = TRUE)
avg_exp <- read_excel("c_avg_exp.xlsx",  col_names = TRUE)

#inner_join the two statistics test results from wilcoxin and roc 
c_cluster12vsTA_wil_roc_join <- inner_join(c_clusterStemVsTA_wil, c_clusterStemVsTA_roc, by="Gene")

#inner_join the two statistics test results with the original average expression data
c_cluster12vsTA_wil_roc_Aveexp <- avg_exp  %>%
  inner_join(c_clusterStemVsTA_wil, by = "Gene")%>%
  inner_join(c_clusterStemVsTA_roc, by = "Gene") 
write.xlsx(c_cluster12vsTA_wil_roc_Aveexp, "c_cluster12vsTA_wil_roc_Aveexp.xlsx", rowNames=F)

##########################################################
#pathway analysis: Metascape website
##########################################################
#Edit the genes name to capital for Metascape website using
head()
c_12vsTA_list <- read.table(file = 'c_12vsTA_list.txt')
head(c_12vsTA_list)
c_12vsTA_list <- toupper(c_12vsTA_list$V1)
head(c_12vsTA_list)
write.table(c_12vsTA_list, file = "c_12vsTA_list_capital.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

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
#Pseudotime: #Monocle, Velocyto , CytoTRACE, SlingShot.
##########################################################################################
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

#################################################################################
#Gene Set Enrichment Analysis (GSEA) on the GSEA-MSigDB website, 
#################################################################################

#Step 1: Load Your DEG Results and Load necessary libraries
library(dplyr)
library(readr)
# Read the DEG file
deg_data <- read_excel("c_clusterStemVsTA_wil.xlsx",  col_names = TRUE)
# View the first few rows
head(deg_data)

#Step 2: Prepare a Ranked Gene List for GSEA

# Rank genes by log2 fold-change
ranked_genes <- deg_data %>%
  arrange(desc(avg_log2FC)) %>% # Sort in descending order
  dplyr::select(Gene, avg_log2FC) #Force dplyr::select() Explicitly

# Convert to named vector for GSEA
gene_list <- setNames(ranked_genes$avg_log2FC, ranked_genes$Gene)

# Save ranked gene list to a file for GSEA input
write.table(ranked_genes, file = "c_ranked_genes.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Step 3: Run GSEA in R (fgsea)

# Install fgsea if not installed
if (!requireNamespace("fgsea", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("fgsea")
}

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####Step 3 pathway selection: Load MSigDB Gene Sets associated with stem cell pathway
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load necessary libraries
library(msigdbr)
library(fgsea)
library(ggplot2)
library(dplyr)
library(pheatmap)

# 1: Load MSigDB Gene Sets (Stem Cell Pathways)
msigdb <- msigdbr(species = "Mus musculus")
stem_cell_pathways <- msigdb %>%
  filter(grepl("stem|pluripotent|progenitor|self-renew|embryonic", gs_name, ignore.case = TRUE))

# 2: Prepare Gene Sets for GSEA
gene_sets <- split(stem_cell_pathways$gene_symbol, stem_cell_pathways$gs_name)
# 3: Run GSEA
fgsea_results <- fgsea(pathways = gene_sets,
                       stats = gene_list,
                       minSize = 10,
                       maxSize = 500,
                       nperm = 10000)

# Sort and show top enriched pathways
fgsea_results <- fgsea_results %>% arrange(padj)
head(fgsea_results)

# Convert list column to character for saving
fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ";"))
head(fgsea_results)

# Save results 
write.xlsx(fgsea_results, "c_GSEA_StemCell_results.xlsx")
#write.csv(fgsea_results, "c_GSEA_StemCell_results.csv", row.names = FALSE)


#select GSEA
significant_pathways <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(padj)  # Sort by significance

top_upregulated <- significant_pathways %>%
  filter(NES > 0) %>%
  arrange(desc(NES))  # Strongest positive enrichment

top_downregulated <- significant_pathways %>%
  filter(NES < 0) %>%
  arrange(NES)  # Strongest negative enrichment

moderate_size_pathways <- significant_pathways %>%
  filter(size > 15, size < 500)

library(ggplot2)

ggplot(significant_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Enriched Pathways (FGSEA)", x = "Pathway", y = "Normalized Enrichment Score") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Downregulated", "Upregulated"))

#==================================================
### SCENIC 
#==================================================
library(SCENIC)
library(dplyr)
library(Matrix)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(data.table)


head(sc_subset_c@meta.data)
sc_subset_c <- SetIdent(sc_subset_c, value = "seurat_clusters")
levels(sc_subset_c)
table(Idents(sc_subset_c))

# subset of 9 clusters from whole data
LSC <- subset(sc_subset_c, idents=c(5,7,12))
LSC@active.assay <-'RNA'
LSC_norm <- NormalizeData(LSC, normalization.method = "LogNormalize", scale.factor = 10000)
LSC_norm$cluster <- LSC_norm@active.ident

# Assign cluster names
cidents <- c("TA_C5","TA_C7",'LSC_C12')
table(LSC_norm@active.ident)
new.cluster.ids <- c("TA_C5","TA_C7",'LSC_C12')
names(new.cluster.ids) <- levels(LSC_norm)
LSC_norm <- RenameIdents(LSC_norm, new.cluster.ids)

# use for LSC_TA_C5 and C7
# Generate subset matrix
for (i in 1:2) {
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

# use for LSC12
i=3
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

################################################################################
#To identify the Key Regulators by subset matrix "LSC_LSC_C12.csv".
################################################################################
# Load necessary libraries
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(data.table)
library(dplyr)
library(Matrix)

# Step 1: Load the subset matrix (LSC_LSC_C10.csv)
exprMat <- fread("LSC_LSC_C12.csv", header = TRUE, data.table = FALSE)  # Read file as a data frame
# Convert first column (gene names) to row names
rownames(exprMat) <- exprMat[,1]  # Assign first column as row names
exprMat <- exprMat[,-1]  # Remove the first column as it's now row names

# Convert to matrix (SCENIC requires a matrix format)
exprMat <- as.matrix(exprMat)

# Check matrix format
print(dim(exprMat))
print(head(exprMat[, 1:5]))  # View first 5 columns

# Step 2: Ensure genes are in rows and cells in columns
if (ncol(exprMat) > nrow(exprMat)) {
  exprMat <- t(exprMat)
}

# Step 3: Gene Filtering (keep only highly expressed genes)
genesKept <- rownames(exprMat)[rowMeans(exprMat > 0.01) > 0.01]
exprMat_filtered <- exprMat[genesKept, ]

# Step 4: Gene Regulatory Network (GRN) Inference using GENIE3
set.seed(123)
weightMat <- GENIE3(exprMat_filtered)

# Step 5: Identify Regulons using RcisTarget

#the code doesn't works
list.files("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/SCENIC_databases/")
load("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/SCENIC_databases/motifAnnotations_mgi.RData")
#https://github.com/aertslab/RcisTarget/blob/master/data/motifAnnotations_mgi.RData
#https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/
motifAnnotations_mgi <- motifAnnotations

data(motifAnnotations_mgi)

scenicOptions <- SCENIC::initializeScenic(org = "mgi", dbDir = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/SCENIC_databases/", nCores = 10)
library(Matrix)
# Update SCENIC and dependencies
install.packages("BiocManager")
BiocManager::install("SCENIC", force = TRUE)
BiocManager::install("RcisTarget")

weightMat <- as(weightMat, "CsparseMatrix")  # Convert to sparse matrix
exprMat_filtered <- as(exprMat_filtered, "CsparseMatrix")  # Convert to sparse matrix
class(weightMat)  # Should return "dgCMatrix"
class(exprMat_filtered)  # Should return "dgCMatrix"

#issue with codes origin?
regulons <- SCENIC::runSCENIC_1_coexNetwork2modules(weightMat, exprMat_filtered)






# Run RcisTpath = # Run RcisTarget to Identify Regulons
regulons <- SCENIC::runSCENIC_2_createRegulons(weightMat, exprMat_filtered, org="mmu")

















# Step 6: Compute Regulon Activity with AUCell
aucellRankings <- AUCell_buildRankings(exprMat_filtered, nCores = 1)
regulonAUC <- AUCell_run(aucellRankings, regulons, nCores = 1)

# Step 7: Identify Key Regulators
keyRegulators <- AUCell_exploreThresholds(regulonAUC, plotHist = TRUE)
highConfRegulons <- names(which(sapply(keyRegulators$thresholds, function(x) x$aucThr > 0.5)))

# Step 8: Save Results
write.csv(regulonAUC@assays$data$AUC, "LSC_C12_RegulonAUC.csv")
write.csv(highConfRegulons, "LSC_C12_Key_Regulators.csv")

# Step 9: Visualize top key regulators
library(ggplot2)
library(ComplexHeatmap)
topTFs <- names(sort(rowMeans(regulonAUC@assays$data$AUC), decreasing = TRUE)[1:10])
Heatmap(regulonAUC@assays$data$AUC[topTFs, ], cluster_rows = TRUE, cluster_columns = TRUE)

print("SCENIC analysis completed for LSC_C12. Key Regulators Identified.")

