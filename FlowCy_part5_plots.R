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

load(file = "workspaceDA.rds")

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", 
                by = "cluster_id", scale = "first", bars = TRUE, perc = TRUE)

# set outputPath
plotsPath <- paste(getwd(), "plots", sep = "/")
dir.create(plotsPath)

display.brewer.all(colorblindFriendly = FALSE)

# plotDiffHeatmap(sce, da, top_n = 12, all = TRUE, fdr = FDR_cutoff)
# stat_test <- as.tibble(da)
# p.adj.signif <- c(rep("ns", 4), "*", rep("ns", 2), "**", "***", "ns", "***", "****")
# group1 <- (rep("BM", nrow(stat_test)))
# group2 <- (rep("PB", nrow(stat_test)))
# y.position <- c(50, 13, 80, 50,
#                 35, 30, 5, 0.5,
#                 2.5, 11, 10, 3)
# stat.test <- cbind(stat_test, group1, group2, p.adj.signif, y.position)
# stat.test <- as.tibble(stat.test)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
# bxp <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5)
# bxp
# ggsave("abundances_stat.svg", plot = last_plot(), dpi = 300)

plotAbundances(sce, k = "cluster_annotation", by  = "cluster_id", group_by = "condition")
ggsave("abundances.svg", plot = last_plot(), dpi = 300)

svg("ExprHeatmap.svg", width = 7, height = 14)
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "median",
                scale = "last", bars = FALSE, perc = FALSE)
dev.off()

plotDR(sce, dr = "UMAP", color_by = "cluster_annotation", facet_by = "condition") + 
  geom_density2d(binwidth = 0.006, colour = "black")
ggsave("UMAP_wContours.svg", plot = last_plot())

svg("Multiheatmap.svg", width = 7, height = 14)
cdx <- type_markers(sce)
plotMultiHeatmap(sce, k = "meta12",
                 hm1 = cdx, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = FALSE)
dev.off()
