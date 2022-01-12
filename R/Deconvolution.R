library(Seurat)
sc_data <- readRDS("/data/allen_cortex_dwn.rds")
sc_data <- Seurat::SCTransform(sc_data , verbose = FALSE)
sc_data <-sc_data %>% Seurat::RunPCA(verbose = FALSE) %>%Seurat::RunUMAP(dims = 1:30, verbose = FALSE) %>%Seurat::FindNeighbors(dims = 1:30, verbose = FALSE) %>%Seurat::FindClusters(verbose = FALSE)

sc_data <-readRDS("/scRNA-seq/reference.rds")
