############Figure 6A-B############
library(RColorBrewer)
library(slingshot)
library(SingleCellExperiment)
scobj <- readRDS("CD8.rds")
scale.data<- scobj@assays$integrated@scale.data
scale.gene <- rownames(scale.data)
counts <- scobj@assays$RNA@counts
counts <- counts[scale.gene,]
sim <- as.SingleCellExperiment(scobj)
umap = scobj@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
# Add the result of the dimensionality reduction to the SingleCellExperiment object
reducedDims(sim) = SimpleList(UMAP = umap)
meta = scobj@meta.data
# Add meta.data to the SingleCellExperiment object
colData(sim)$sampleId = scobj$name_day
colData(sim)$celltype = scobj@active.ident

sim <- slingshot(sim, clusterLabels='celltype',# Select the column name of the cell annotation in colData
                 reducedDim='UMAP',
                 start.clus="CD8_CCR7",#Choose a root
                 reweight = T)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) 
plotcol <- colors[cut(sim$slingPseudotime_2, breaks=100)] 
plotcol[is.na(plotcol)] <- "lightgrey"
plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("right", legend = paste0("lineage",1:3), col = unique(brewer.pal(6,"Set1")), inset=0.8, pch = 16)

coldata <- colData(sim)
coldata <- data.frame(celltype = coldata@listData$celltype,   sampleId = coldata@listData$sampleId,
                      plotcol = plotcol)
rownames(coldata) = sim@colData@rownames
filter_cell <- dplyr::filter(coldata, plotcol != "lightgrey")
filter_cell <- rownames(filter_cell)
head(filter_cell)
filter_counts <- counts[,filter_cell]
# Randomly selected 2000 cells
set.seed(111)
scell <- sample(colnames(filter_counts), size = 2000)
filter_counts = filter_counts[, scell]
# Re-convert the expression matrix to a SingleCellExperiment object
filter_sim <- SingleCellExperiment(assays = List(counts = filter_counts))
# colData
filter_coldata = colData(sim)[scell, 1:3]
filter_sim@colData = filter_coldata
rd = reducedDim(sim)
filter_rd <- rd[scell,]
reducedDims(filter_sim) <- SimpleList(UMAP = filter_rd)

library(tradeSeq)
# Fit negative binomial model
counts <- filter_sim@assays@data$counts
crv <- SlingshotDataSet(filter_sim)
set.seed(111)
icMat <- evaluateK(counts = counts, sds = crv,k = 3:10,                  
                   nGenes = 500,                  
                   verbose = T)
#  we  pick  nknots  =  6.
set.seed(111)
pseudotime  <-  slingPseudotime(crv,  na  =  FALSE)
cellWeights  <-  slingCurveWeights(crv)

system.time({  sce  <-  fitGAM(counts  =  counts,   pseudotime  =  pseudotime,   
                               cellWeights  =  cellWeights,   nknots  =  6,   verbose  =  FALSE)})
table(rowData(sce)$tradeSeq$converged)
assoRes  <-  associationTest(sce)
head(assoRes)
startRes  <-  startVsEndTest(sce)
head(startRes)
assoRes <- assoRes[order(startRes$waldStat, decreasing = TRUE),]
gene_list_lo <- rownames(assoRes[order(assoRes$waldStat, decreasing = TRUE),])[1:200]
yhatSmooth <- predictSmooth(sce, gene = gene_list_lo, nPoints = 50, tidy =
                              FALSE)
Heatmap_same_Genes<- ComplexHeatmap::Heatmap(t(scale(t(yhatSmooth))),
                                             cluster_columns = FALSE,
                                             show_row_names = T,
                                             show_column_names = F)
