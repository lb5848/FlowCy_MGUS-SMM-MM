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
library(harmony)

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

sce <- readRDS("SCE_BMvsPB.rds")
CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

# Run FlowSOM and ConsensusClusterPlus
seed <- 123456
set.seed(seed)
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
               verbose = TRUE, seed = seed)
delta_area(sce)
# Run dimensionality reduction
n_cells <- 5000
n_events <- min(n_cells(sce))
ifelse(n_cells > n_events, n_cells <- n_events, n_cells <- n_cells)
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
# sce <- runDR(sce, dr = "TSNE", cells = n_cells, features = "type", theta = 0.5, max_iter = 1000, 
#              distMethod = "euclidean",
#              PCA = TRUE, eta = eta, exaggeration_factor = 12.0)
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")
sce <- runDR(sce, dr = "DiffusionMap", cells = n_cells, features = "type", assay = "exprs")

saveRDS(sce, file = "SCE_BMvsPB_DR.rds")

display.brewer.all(colorblindFriendly = TRUE)
delta_area(sce)
cdx <- type_markers(sce)
plotMultiHeatmap(sce, k = "meta12",
                 hm1 = cdx, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = FALSE)

plotMultiHeatmap(sce, k = "meta10",
                 hm1 = cdx, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = FALSE)

plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id",
                fun = "mean", scale = "last", bars = TRUE, perc = TRUE)
plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta12", by = "cluster_id",
                fun = "mean", scale = "last", bars = TRUE, perc = TRUE)
plotAbundances(sce, k = "meta12", by = "cluster_id", group_by = "condition")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta6", by = "cluster_id",
                fun = "mean", scale = "last", bars = TRUE, perc = TRUE)
plotAbundances(sce, k = "meta6", by = "cluster_id", group_by = "condition")

CATALYST::plotDR(sce, dr = "UMAP", color_by = "meta12", facet_by = "condition")

CATALYST::plotDR(sce, dr = "DiffusionMap", color_by = "meta4", facet_by = "condition")


CATALYST::plotDR(sce, dr = "UMAP", color_by = c("CD127", "CD69"), facet_by = "condition")
plotDR(sce, dr = "UMAP", color_by = "meta12", facet_by = "condition") +  
  geom_density2d(binwidth = 0.006, colour = "black")
