rm(list = ls())

# Load packages
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)
library(readxl)
library(flowStats)
library(FlowSOMworkshop)
library(tidyverse)
library(data.table)
library(ggpubr)
library(flowAI)
library(PeacoQC)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# set workingDir
workingDir <- "201114_WorkingDirectory"
workingDirPath <- paste(PrimaryDirectory, workingDir, sep = "/")
setwd(workingDirPath)

sce <- readRDS("SCE_BMvsPB_DR.rds")
CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta12", 
                by = "cluster_id", scale = "last", bars = TRUE, perc = TRUE)
plotAbundances(sce, k = "meta12", by = "cluster_id", group_by = "condition")

# annotate clusters

annotation_table <- as.data.frame(cbind(c(1:12), c(1:12)))

colnames(annotation_table) <- c("meta12", "Clusters")
annotation_table$Clusters <- factor(annotation_table$Clusters)
sce <- mergeClusters(sce, k = "meta12", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")
# filtered_sce <- filterSCE(sce, cluster_id %in% c(paste0("C", c(1:12))), k = "cluster_annotation")

# store original_sce
# old_sce <- sce
# sce <- filtered_sce

FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")

# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition"))
contrast <- createContrast(c(0, 1, 1))

nrow(contrast) == ncol(design)

out_DA <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrast,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = FALSE)

da <- rowData(out_DA$res)
plotDiffHeatmap(sce, da, top_n = 12, all = TRUE, fdr = FDR_cutoff)

save(list = ls(), file = "workspaceDA.rds")










