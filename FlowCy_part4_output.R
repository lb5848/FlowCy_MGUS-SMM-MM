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
outputPath <- paste(getwd(), "output", sep = "/")
dir.create(outputPath)

# save single fcs_files by cluster and by condition
# sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")
by_cluster <- sce2fcs(sce, split_by = "cluster_annotation", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(by_cluster, outdir = outputPath, filename = "bycluster")
merged <- sce2fcs(sce, split_by = NULL, keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.FCS(merged, filename = paste(outputPath, "merged.fcs", sep = "/"))
by_condition <- sce2fcs(sce, split_by = "condition", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(by_condition, outdir = outputPath, filename = "bycondition")
