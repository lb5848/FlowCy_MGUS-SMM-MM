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

# create folders//directories

dirFlowFiles <- "flow_files"
dirPath <- paste(PrimaryDirectory, dirFlowFiles, sep = "/")
dir.create(dirPath)

dirCSVfiles <- "csv_files"
CSVdirPath <- paste(dirPath, dirCSVfiles, sep = "/")
dir.create(CSVdirPath)

dirFCSfiles <- "fcs_files"
FCSdirPath <- paste(dirPath, dirFCSfiles, sep = "/")
dir.create(FCSdirPath)

workingDir <- "201114_WorkingDirectory"
workingDirPath <- paste(PrimaryDirectory, workingDir, sep = "/")
dir.create(workingDirPath)

# List files
CSVfiles <- list.files(CSVdirPath, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)

# create sample_metadata.csv
# 
# fwrite(list(fileName), "sample_metadata.csv")



convertCSVtoFCS <- function(CSVfiles, csvDirectory, csv2fcsDirectory, fileName){
  for(i in c(1:length(CSVfiles))){
    data <- fread(paste(csvDirectory, CSVfiles[i], sep = "/"))
    print(CSVfiles[i])
    print(fileName[i])
    cytofCore.write.FCS(as.matrix(data), 
                        filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                        what = "numeric")
  }
}

convertCSVtoFCS(CSVfiles = CSVfiles, csvDirectory = CSVdirPath, csv2fcsDirectory = FCSdirPath, fileName = fileName)


# Create flowSet from FCSfiles
FCSfiles <- list.files(FCSdirPath, pattern = ".fcs$", full = FALSE)
flowSet <- read.flowSet(files = FCSfiles, path = FCSdirPath, truncate_max_range = FALSE)
colnames(flowSet)

# prep sample_md
sample_md <- fread(file = paste(PrimaryDirectory, "sample_metadata.csv", sep = "/"), header = TRUE)
sample_md
is.data.table(sample_md)

sample_md[, sample_id := factor(sample_id)]
sample_md[, condition := factor(condition, 
                                levels = c("MGUS", "SMM", "NDMM"))]
levels(sample_md$condition)
colnames(sample_md)
fwrite(sample_md, "sample_metadata.csv")

# create panel_md.csv
# fwrite(list(colnames(flowSet)), "panel_md.csv")

panel_md <- fread(file = paste(PrimaryDirectory, "panel_md.csv", sep = "/"), header = TRUE)
all(colnames(flowSet) == panel_md$fcs_colname)

# preprocessing - QC

setwd(workingDirPath)

QC_dir <- "QC"
if(!dir.exists(QC_dir)){
  dir.create(QC_dir)
  dir.create(file.path(QC_dir, "flowAI"))
  dir.create(file.path(QC_dir, "PeacoQC"))
}
QC_dir_path <- QC_dir
PeacoQC_path <- paste(QC_dir_path, "PeacoQC", sep = "/")

fs_AI <- flowCore::fsApply(flowSet, function(ff){
  resQC_AI <- flow_auto_qc(fcsfiles = ff,
                           folder_results = file.path(QC_dir, paste("flowAI", gsub(".fcs", " Folder", ff@description$GUID), sep = "/")),
                           output = 1)
  return(resQC_AI)
})

fs_PeacoQC <- flowCore::fsApply(flowSet, function(ff){
  resQC <- PeacoQC(ff = ff,
                   determine_good_cells = "all",
                   channels = c(10:35),
                   plot = TRUE,
                   output_directory = PeacoQC_path)
  ff <- ff[resQC$GoodCells, ]
  return(ff)
})

flowSet[[1]]@description$`$CYT`
flowSet[[1]]@description$`$CYT` <- "FACS"

chs_of_interest <- colnames(flowSet)[10:35]
plot_aggregate(flowSet, channels = chs_of_interest, output_image = "FCSpreNorm.png")
plot_aggregate(fs_AI, channels = chs_of_interest, output_image = "FCSpreNormpostAI.png")
plot_aggregate(fs_PeacoQC, channels = chs_of_interest, output_image = "FCSpreNormpostPeacoQC.png")
# normFlowSet <- warpSet(flowSet, stains = c("TCF1","Eomes","Tbet"))
# plot_aggregate(normFlowSet, channels = chs_of_interest, output_image = "FCSpostNorm.png")

# try fs_AI, no norm needed

fs_AI[[1]]@description$`$CYT`
fs_AI[[1]]@description$`$CYT` <- "FACS"

# prepData -
colnames(sample_md)
all(colnames(fs_AI) == panel_md$fcs_colname)
flowSet[[1]]@description$`$CYT`

sce <- CATALYST::prepData(fs_AI, panel = panel_md, md = sample_md, transform = FALSE)
assay(sce, "exprs") <- assay(sce, "counts")
assays(sce)

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 9
p

n_events <- min(n_cells(sce))
n_events
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

plotNRS(sce, features = type_markers(sce), color_by = "condition")

saveRDS(sce, file = "SCE_BMvsPB.rds")






