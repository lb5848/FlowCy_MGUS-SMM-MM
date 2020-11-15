rm(list = ls())
# MacOS
# install xquartz using homebrew (brew cask install xquartz)
# install cairo (brew install cairo)


install_packages <- TRUE

if(!require(devtools)){
  BiocManager::install("devtools")
}

if(install_packages){
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
  BiocManager::install(c( 
    "CATALYST", 
    "FlowSOM", 
    "RColorBrewer", 
    "Rtsne", 
    "XML", 
    "cluster", 
    "data.table", 
    "diffcyt", 
    "dplyr",
    "flowAI",
    "flowCore", 
    "flowStats", 
    "flowViz", 
    "flowWorkspace", 
    "ggcyto", 
    "ggplot2", 
    "ggpubr", 
    "ggthemes", 
    "readxl", 
    "scales", 
    "scater", 
    "scran", 
    "stringr", 
    "tidyverse", 
    "uwot" 
  ))
  
  devtools::install_github("nolanlab/cytofCore")
  devtools::install_github("saeyslab/PeacoQC")
  remotes::install_github("saeyslab/FlowSOM_workshop")
}

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