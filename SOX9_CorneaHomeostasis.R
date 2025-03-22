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

counts_c<-GetAssayData(object = sc_subset_c ,  assay = "RNA", layer = "counts")
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

#Trpm8
seurat_c$Trpm8_expr <- GetAssayData(object = sc_subset_c, assay = "RNA", slot = "counts")["Trpm8", ] > 0
sum(counts_c["Trpm8",]>0)


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
stem_cluster <- 12  # replace with the actual stem cell cluster ID
TA_cluster <-  c(5,7) # replace with TA cluster IDs
Idents(sc_subset_c) <- "seurat_clusters" 
#Wilcox
c_clusterStemVsTA_wil<- FindMarkers(sc_subset_c, ident.1 = stem_cluster , ident.2 = TA_cluster , test.use="wilcox",
                                    min.pct = 0.25, logfc.threshold = 0.25) #only.pos = TRUE 
head(c_clusterStemVsTA_wil)
write.xlsx(c_clusterStemVsTA_wil, file = "c_clusterStemVsTA_wil.xlsx", rowNames=T)

#filter by p_val_adj<0.05
head(c_clusterStemVsTA_wil)
c_clusterStemVsTA_wil_filter <- filter (c_clusterStemVsTA_wil, p_val_adj<0.05)%>%
    arrange(avg_log2FC)
head(c_clusterStemVsTA_wil_filter) 

#lost genes name by saving, so
#in order to save it to excel with gene name.
df <- as.data.frame.matrix(c_clusterStemVsTA_wil_filter)
head(df)

#The gene names are actually stored as row names.generate new column with Gene name.
df <- df %>% rownames_to_column(var = "Gene") # Move row names into a new column
#df <- df[, -1]
write.xlsx(df, file = "c_clusterStemVsTA_wil_filter.xlsx", rowNames = TRUE) 

#Roc
c_clusterStemVsTA_roc<- FindMarkers(sc_subset_c, ident.1 = stem_cluster , ident.2 = TA_cluster , test.use="roc",
                                    min.pct = 0.25, logfc.threshold = 0.25)
head(c_clusterStemVsTA_roc)
write.xlsx(c_clusterStemVsTA_roc, file = "c_clusterStemVsTA_roc.xlsx", rowNames=T)

#filter by myAUC
head(c_clusterStemVsTA_roc)
c_clusterStemVsTA_roc_filter <- filter (c_clusterStemVsTA_roc, myAUC>0.7)%>%
  arrange(avg_log2FC)
head(c_clusterStemVsTA_roc_filter) 
write.xlsx(c_clusterStemVsTA_roc_filter, file="c_clusterStemVsTA_roc_filter.xlsx", rownames=T)

#Overlap Wilcon and Roc
c_clusterStemVsTA_overlap_filter <- inner_join(c_clusterStemVsTA_wil_filter, c_clusterStemVsTA_roc_filter, by )

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

#==================================================
### pathway 
#==================================================

##########################################################
#pathway analysis 1: Metascape website
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

##################################################
#Using Monocle For Pseudotime Trajectory (Time permits)
##################################################
#https://monashbioinformaticsplatform.github.io/Single-Cell-Workshop/pbmc3k_tutorial.html#Finding_differentially_expressed_features_(cluster_biomarkers)

library(monocle3)
library(SeuratWrappers)

cds <- as.cell_data_set(sc_subset_c)
cds

cds <- cluster_cells(cds, k = 10, random_seed = 5)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
#install.packages("patchwork")  # Install if not already installed
library(patchwork)  # Load the package
wrap_plots(p1, p2)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE)

head(colData(cds))
head(sc_subset_c)

# Create a vector of idents to keep
selected_ids <- c("2", "3", "4","5", "7",  "9", "10",  "12")
cells_sc_subset_c <- subset(sc_subset_c, idents = selected_ids)  ## subset the PBMC seurat object to tcells
cds <- as.cell_data_set(sc_subset_c)  ## convert this to cell_data_set
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

#to order cells in pseudotime
get_earliest_principal_node <- function(cds, cell_type = "12") {
  cell_ids <- which(colData(cds)[, "ident"] == cell_type)
  
  closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)$UMAP)$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,
  ]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE)
plot_cells(cds, color_cells_by = "ident", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)












#==================================================
### Pseudotime: 
##Monocle, Velocyto , CytoTRACE, SlingShot. 
#==================================================

##########################################################################################
##Pseudotime 1:Monocle
##########################################################################################

#subset cells of LSCS, TA, Dif
head(sc_subset_c@meta.data)
sc_subset_c <- SetIdent(sc_subset_c, value = "seurat_clusters")
levels(sc_subset_c)
table(Idents(sc_subset_c))

# subset of 3 clusters from whole data
LSC <- subset(sc_subset_c, idents=c(2,3,4,5,7,9,10, 12))
LSC@active.assay <-'RNA'
LSC_norm <- NormalizeData(LSC, normalization.method = "LogNormalize", scale.factor = 10000)
LSC_norm$cluster <- LSC_norm@active.ident

# Assign cluster names
cidents <- c("Dif1_C2", "Dif2_C3", "Dif3_D4" "TA_C5","TA_C7","LSC_C9","Dif4_D10","LSC_C12")
table(LSC_norm@active.ident)
new.cluster.ids <- c("Dif1_C2", "Dif2_C3", "Dif3_D4" "TA_C5","TA_C7","LSC_C9","Dif4_D10","LSC_C12")
names(new.cluster.ids) <- levels(LSC_norm)
LSC_norm <- RenameIdents(LSC_norm, new.cluster.ids)

# use for "Dif1_C2", "Dif2_C3", "Dif3_D4" "TA_C5","TA_C7","LSC_C9","Dif4_D10"
# Generate subset matrix
for (i in 1:7) {
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



#=====================================================================================
### Gene Set Enrichment Analysis (GSEA) from Molecular Signatures Database (MSigDB)
#=====================================================================================

#################################################################################
#GSEA 1 : on the GSEA-MSigDB website, 
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
#####Step 3 pathway selection: Load MSigDB Gene Sets associated with Hallmark pathway
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(fgsea)
library(msigdbr)

category = "H" #"H": hallmark gene set; # C5: GO gene sets

# Load MSigDB Gene Sets (Hallmark Pathway )
hallmark_gene_sets <- msigdbr(species = "Mus musculus", category = category)#Homo sapiens
# Convert to list format for fgsea
msigdb_list <- split(hallmark_gene_sets$gene_symbol, hallmark_gene_sets$gs_name)

# Run GSEA
fgsea_results <- fgsea(pathways = msigdb_list, 
                       stats = gene_list, 
                       minSize = 15, 
                       maxSize = 500)


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####Step 3 alternate pathway selection 2: Load MSigDB Gene Sets associated with stem cell pathway
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load necessary libraries
library(msigdbr)
library(fgsea)
library(ggplot2)
library(dplyr)
library(pheatmap)

# 1: Load MSigDB Gene Sets (Stem Cell Pathways)
#msigdb <- msigdbr(species = "Mus musculus")
#stem_cell_pathways <- msigdb %>%
#  filter(grepl("stem|pluripotent|progenitor|self-renew|embryonic", gs_name, ignore.case = TRUE))

# 2: Prepare Gene Sets for GSEA
#gene_sets <- split(stem_cell_pathways$gene_symbol, stem_cell_pathways$gs_name)
# 3: Run GSEA
#fgsea_results <- fgsea(pathways = gene_sets,
 #                      stats = gene_list,
#                       minSize = 10,
 #                      maxSize = 500,
  #                     nperm = 10000)

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####Step 4
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sort and show top enriched pathways
fgsea_results <- fgsea_results %>% arrange(padj)
head(fgsea_results)

# Convert list column to character for saving
fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ";"))
head(fgsea_results)

# Save results # change the name based on GO or Hallmark.
write.xlsx(fgsea_results, "c_GSEA_Hallmark_results.xlsx")
#write.xlsx(fgsea_results, "c_GSEA_StemCell_results.xlsx")
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
#change the significant_pathways to top_upregulated for plot in upregulate genes only.

## Generate gene File (Pathway and Genes list)
library(tidyverse)  # Load required package

colnames(significant_pathways)
# Convert 'significant_pathways' to long format
GSEA_list <- significant_pathways %>%
  select(pathway = pathway, leadingEdge, pval,padj, log2err, ES, NES,size ) %>%  # Select relevant columns
  separate_rows(leadingEdge, sep = ";") %>%  # Split 'leadingEdge' into multiple rows
  rename(Gene = leadingEdge)  # Rename column to 'Gene'

# View the first few rows
head(GSEA_list)

# View sorted pathway names
significant_pathways <- significant_pathways %>% arrange(desc(NES))
significant_pathways$pathway

HALLMARK1<- significant_pathways[significant_pathways$pathway=="HALLMARK_KRAS_SIGNALING_UP", ]
HALLMARK2<- significant_pathways[significant_pathways$pathway=="HALLMARK_IL2_STAT5_SIGNALING", ]
HALLMARK3<- significant_pathways[significant_pathways$pathway=="HALLMARK_G2M_CHECKPOINT", ]

HALLMARK4 <- significant_pathways[significant_pathways$pathway=="HALLMARK_ESTROGEN_RESPONSE_EARLY", ]
HALLMARK5<- significant_pathways[significant_pathways$pathway=="HALLMARK_INTERFERON_GAMMA_RESPONSE", ]
HALLMARK6<- significant_pathways[significant_pathways$pathway=="HALLMARK_ESTROGEN_RESPONSE_LATE", ]

# Select only the pathway and leadingEdge columns
HALLMARK1 <- HALLMARK1 %>%
  select(Pathway = pathway,leadingEdge, pval,padj, log2err, ES, NES,size) %>%
  separate_rows(leadingEdge, sep = ";") %>%
  rename(Gene = leadingEdge)

HALLMARK2 <- HALLMARK2 %>%
  select(Pathway = pathway,leadingEdge, pval,padj, log2err, ES, NES,size) %>%
  separate_rows(leadingEdge, sep = ";") %>%
  rename(Gene = leadingEdge)

HALLMARK3 <- HALLMARK3 %>%
  select(Pathway = pathway, leadingEdge, pval,padj, log2err, ES, NES,size) %>%
  separate_rows(leadingEdge, sep = ";") %>%
  rename(Gene = leadingEdge)

GSEA_overlap <- inner_join(HALLMARK2 ,HALLMARK3, by ="Gene" )
head(GSEA_overlap) 
#Overlap: HALLMARK1n2n3: 0. 1n2: Car2, Rabgap1l.1n3: 0. 2n3: 0

c_LSCSvsTA_filter <- read_excel("c_clusterStemVsTA_wil_filter.xlsx", sheet =1, col_names = TRUE)
GSEA_overlap <- inner_join(HALLMARK2 ,c_LSCSvsTA_filter, by ="Gene" )
head(GSEA_overlap$Gene)
write.xlsx(GSEA_overlap, "c_LSCSvsTA_wil_filter_overlapGSEA_HALLMARK_IL2_STAT5_SIGNALING.xlsx")

setwd("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/GSEA")  
write.xlsx(GSEA_overlap, "c_GSEA_LSCSvsTA_HALLMARK_G2M_CHECKPOINT.xlsx")

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

# subset of 3 clusters from whole data
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
#SCENIC continue: To identify the Key Regulators by subset matrix "LSC_LSC_C12.csv".
################################################################################
# Load necessary libraries
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(data.table)
library(dplyr)
library(Matrix)

# Step 1: Load the subset matrix (LSC_LSC_C12.csv)
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
#################
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

#the code doesn't works
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


#=================================================================================
#Pathway Enrichment in DAVID
#=================================================================================

#https://davidbioinformatics.nih.gov/conversion.jsp
#c_clusterStemVsTA_wil_genelist.txt
#official_gene_symbol
#Mus musculus
#gene list
#gene ontology: goterm_bp_fat, 
#pathway: kegg
#Functional Annotation Chart
#2009 chart records

#convert the filtered gene symbol to ENSEMBL in c_clusterStemVsTA_wil_genelist by DAVID. Then read the list. 
  #1123 genes.
    
    genelist <- read.table("c_clusterStemVsTA_wil_genelist_convert.txt", sep = "\t", header = TRUE)
    head(genelist, 5)
    colnames(genelist)[1] <- "Gene_Symbol"  
    colnames(genelist)[2] <- "ENSEMBL_ID"  
    
    #prepare the gene list for DAVID
    deg_data <- read_excel("c_clusterStemVsTA_wil_filter.xlsx", col_names = TRUE)
    head(deg_data)
    
    #Extract the Upregulated Genes
    upregulated_genes <- deg_data %>%
      filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
      pull(Gene)
    
    # Write to a text file
    write.table(upregulated_genes, "C_clusterStemVsTA_wil_Upregulated_Genes.txt", 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    #Extract the Downregulated Genes
    downregulated_genes <- deg_data %>%
      filter(avg_log2FC< 0 & p_val_adj < 0.05) %>%
      pull(Gene)
    
    # Write to a text file
    write.table(downregulated_genes, "C_clusterStemVsTA_wil_Downregulated_Genes.txt", 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    #Go to DAVID: https://david.ncifcrf.gov/tools.jsp, 
    #Functional Annotation, generate Pathway Enrichment Analysis, generate csv file.
    #Current Gene List: C_clusterStemVsTA_wil_Upregulated_Genes
    #Current Background: Mus musculus
    #786 DAVID IDs
    
    # Read the CSV file
    
    Upregulated_DAVID <- read.csv("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/DAVID/C_clusterStemVsTA_wil_Upregulated_DAVID.csv", sep = "\t")
    write.xlsx(Upregulated_DAVID,"/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/DAVID/C_clusterStemVsTA_wil_Upregulated_DAVID.xlsx", rowNames=F)
    Upregulated_DAVID_pfilted <- Upregulated_DAVID  %>%
      filter(PValue<0.05 &  Bonferroni < 0.05 & Benjamini < 0.05& FDR< 0.05) 
    head(Upregulated_DAVID_pfilted)
    write.xlsx(Upregulated_DAVID_pfilted,"/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/DAVID/Upregulated_DAVID_pfilted.xlsx", rowNames=F)
    #162 Functional Annotation chart records left
    
    #david websit maintaining.
    downregulated_DAVID_pfilted <- downregulated_DAVID  %>%
      filter(PValue<0.05 &  Bonferroni < 0.05 & Benjamini < 0.05& FDR< 0.05) 
    head(downregulated_DAVID_pfilted)
    write.xlsx(Upregulated_DAVID_pfilted,"/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/DAVID/Upregulated_DAVID_pfilted.xlsx", rowNames=F)
    
    filtered_results <- Upregulated_DAVID_pfilted %>%
        filter(grepl("stem cell maintenance|differentiation|pluripotency|proliferation|stem cell|epithelial cell", Term, ignore.case = TRUE))#
    print(filtered_results)
    write.xlsx(filtered_results,"/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/DAVID/C_clusterStemVsTA_wil_Upregulated_DAVID_stem.xlsx", rowNames=F)
    head(filtered_results)
    #6 Functional Annotation chart records left.
    
    #list the gene in pathway
    filtered_results$Term
    #[1] "GO:0042127~regulation of cell population proliferation"         
    #[2] "GO:0008285~negative regulation of cell population proliferation"
    #[3] "GO:0008284~positive regulation of cell population proliferation"
    #[4] "GO:0045597~positive regulation of cell differentiation"   # did not match       
    #[5] "GO:0050678~regulation of epithelial cell proliferation"         
    #[6] "GO:0045595~regulation of cell differentiation" 
    
    David1 <-filtered_results[ filtered_results$Term=="GO:0042127~regulation of cell population proliferation"  ,]
    head(David1)
    # Select only the pathway and genes
    David1_genes <- David1 %>%
      select(pathway=Term, Genes) %>%
      separate_rows(Genes, sep = ",") %>%
      rename(Gene = Genes)
    
    #capitalize only the first letter of each gene symbol while keeping the rest in lowercase
    library(dplyr)
    library(stringr)
    # Trim spaces and capitalize first letter only
    David1_genes <- David1_genes %>%
      mutate(Gene = str_to_title(str_trim(Gene)))
    head(David1_genes)
    
    David2 <-filtered_results[ filtered_results$Term=="GO:0008285~negative regulation of cell population proliferation"  ,]
    head(David1)
    # Select only the pathway and genes
    David2_genes <- David2 %>%
      select(pathway=Term, Genes) %>%
      separate_rows(Genes, sep = ",") %>%
      rename(Gene = Genes)
    # Trim spaces and capitalize first letter only
    David2_genes <- David2_genes %>%
      mutate(Gene = str_to_title(str_trim(Gene)))
    head(David2_genes)
    
    David3 <-filtered_results[ filtered_results$Term=="GO:0008284~positive regulation of cell population proliferation"  ,]
    head(David1)
    # Select only the pathway and genes
    David3_genes <- David3 %>%
      select(pathway=Term, Genes) %>%
      separate_rows(Genes, sep = ",") %>%
      rename(Gene = Genes)
    # Trim spaces and capitalize first letter only
    David3_genes <- David3_genes %>%
      mutate(Gene = str_to_title(str_trim(Gene)))
    head(David3_genes)
    
    David4 <-filtered_results[ filtered_results$Term=="GO:0050678~regulation of epithelial cell proliferation"  ,]
    head(David1)
    #Select only the pathway and genes
    David4_genes <- David4 %>%
      select(pathway=Term, Genes) %>%
      separate_rows(Genes, sep = ",") %>%
      rename(Gene = Genes)
    # Trim spaces and capitalize first letter only
    David4_genes <- David4_genes %>%
      mutate(Gene = str_to_title(str_trim(Gene)))
    head(David4_genes)
    
    David5 <-filtered_results[ filtered_results$Term=="GO:0045595~regulation of cell differentiation"  ,]
    head(David1)
    #Select only the pathway and genes
    David5_genes <- David5 %>%
      select(pathway=Term, Genes) %>%
      separate_rows(Genes, sep = ",") %>%
      rename(Gene = Genes)
    head(David5_genes)
    # Trim spaces and capitalize first letter only
    David5_genes <- David5_genes %>%
      mutate(Gene = str_to_title(str_trim(Gene)))

#overlap genes between David1-5
    David_overlap_gene <- David1_genes%>% inner_join( David2_genes, by="Gene")%>%
      inner_join(David3_genes, by="Gene" )%>%inner_join(David4_genes, by="Gene" )%>%inner_join(David5_genes, by="Gene" )
    print(David_overlap_gene)
    write.xlsx(David_overlap_gene, "overlap_gene_David1&2&3&4&5.xlsx")
    #overlap genes total by David1&2:  60 
    #overlap genes total by David1,2,3,4,5: 8
    
    
#search the genes expression 
    c_LSCSvsTA_filter <- read_excel("c_clusterStemVsTA_wil_filter.xlsx", sheet =1, col_names = TRUE)
    David_overlap <- inner_join(David1_genes,c_LSCSvsTA_filter, by="Gene")
    head(David_overlap )
    print(David_overlap )

#=================================================================================
# to generate a Cytoscape regulatory network from DAVID Pathway Results  
#=================================================================================
# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

# View the columns in your DAVID file (optional)
david_results <- read_excel("/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/DAVID/C_clusterStemVsTA_wil_Upregulated_DAVID_stem.xlsx", col_names = TRUE)

colnames(david_results)

# Extract Pathway Names and Gene List from DAVID results
pathway_data <- david_results %>%
  select(Term, Genes, PValue) %>%
  mutate(Pathway = str_extract(Term, "^[^\\(]+"), # Extract pathway name
         GeneList = strsplit(Genes, ", ")) %>% 
  unnest(GeneList)

# View the processed data (Pathway - Gene Mapping)
head(pathway_data)

# Save the Nodes File for Cytoscape
write.csv(pathway_data,"/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/upDAVID_Pathway_Gene_Network.csv")

# Generate Cytoscape Network Files (Nodes/Edges)
#Nodes File: Lists each gene and pathway with p-values.
#Edges File: Defines connections between genes and pathways.
# Generate Nodes File (Genes + Pathways)
nodes <- data.frame(
  ID = unique(c(pathway_data$Pathway, pathway_data$GeneList)),
  Label = unique(c(pathway_data$Pathway, pathway_data$GeneList)),
  Type = ifelse(unique(c(pathway_data$Pathway, pathway_data$GeneList)) %in% pathway_data$Pathway, "Pathway", "Gene")
)

# Save Nodes File
write_csv(nodes, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/upDAVID_Nodes.csv")
write.xlsx(nodes, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/upDAVID_Nodes.xlsx")

# Generate Edges File (Pathway ↔ Gene)
edges <- data.frame(
  Source = pathway_data$Pathway,
  Target = pathway_data$GeneList,
  Interaction = "associated_with"
)

# Save Edges File
write_csv(edges, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/upDAVID_Edges.csv")
write.xlsx(edges, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/upDAVID_Edges.xlsx")

#Download Cytoscape : https://cytoscape.org/download.html

#=================================================================================
#  scNetViz in Cytoscape   
#=================================================================================
#
#generate filtered expression matrix from the 12 cluster

# Subset Seurat object for Cluster 12
cluster12_cells <- subset(sc_subset_c, idents = "12")  # Change "12" to any cluster of interest
#Raw counts (recommended for further processing):
cluster12_matrix <- as.data.frame(GetAssayData(cluster1_cells, layer = "counts"))
#Log-normalized expression (for visualization):
#cluster12_matrix <- as.data.frame(GetAssayData(cluster1_cells, slot = "data"))
write.csv(cluster12_matrix, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/c_expression_matrix_cluster12.csv", row.names = TRUE)
#If extracting multiple clusters, adjust the filename accordingly:
#write.csv(as.data.frame(GetAssayData(selected_clusters, slot = "counts")), x`x`"expression_matrix_clusters_1_3_5.csv", row.names = TRUE)
  

#####filtered expression matrix into Cytoscape using scNetViz

# Remove genes expressed in fewer than 10 cells
filtered_cluster12_matrix <- cluster12_matrix[rowSums(cluster12_matrix > 0) >= 10, ]

# Save filtered matrix
write.csv(filtered_cluster12_matrix, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/c_filtered_expression_matrix_cluster12.csv", row.names = TRUE)

#Filter for Highly Variable Genes(option)
VariableFeatures(cluster12_cells) <- FindVariableFeatures(cluster12_cells, selection.method = "vst", nfeatures = 2000)
filtered_matrix_variable <- cluster12_matrix[VariableFeatures(cluster12_cells), ]

write.csv(filtered_matrix_variable, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/c_variable_genes_expression_matrix_cluster12.csv", row.names = TRUE)

#Export Metadata for Cytoscape
cluster12_metadata <- cluster12_cells@meta.data
write.csv(cluster12_metadata, "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/c_filtered_cell_metadata_cluster12.csv", row.names = TRUE)

###Import into Cytoscape (scNetViz)

##Convert CSV to MTX Format in R

# Load necessary library
library(Seurat)
library(Matrix)

# Convert to sparse matrix format
sparse_matrix <- as(as.matrix(filtered_matrix_variable), "dgCMatrix")

# Save as MTX format
writeMM(sparse_matrix, file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/c_filtered_expression_matrix_cluster12.mtx")

# Save genes (rows) as TSV
write.table(rownames(filtered_matrix_variable), file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/c_genes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save cells (columns) as TSV
write.table(colnames(filtered_matrix_variable), file = "/Users/jackzhou/Desktop/Project_Sox9/sox9_bioinfo_QZ/Cytoscape_Network/c_barcodes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

