library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
print(args)
input_file = args[1]
cluster = args[2]             
proportion = as.numeric(args[3])


##--- Load Data
data <- readRDS(paste0("./Data/",input_file))     # SeuratObject
DefaultAssay(data) <-"RNA"

##--- Add metadata
meta <- data@meta.data %>% as_tibble %>% select(ID)

pheno <- read_csv("/work22/home/redahiro/analysis/COVID-19_scRNAseq/Pheno_data/Phenotype_analysis_220702.csv")

meta <- meta %>% left_join(pheno, by="ID") %>% select(ID,Age,Sex)

data[["Age"]] <- meta$Age
data[["Sex"]] <- meta$Sex

##--- SCE conversion
sce <- as.SingleCellExperiment(data)
altExps(sce) <- NULL                           

reducedDim(sce, "PCA", withDimnames=TRUE) <- data[['harmony']]@cell.embeddings   
reducedDim(sce, "UMAP", withDimnames=TRUE) <- data[['umap']]@cell.embeddings

##--------------------------------------##
##    Differential abundance testing    ##
##--------------------------------------##

##--- Step1: Create a Milo object
milo <- Milo(sce)

##--- Step2: Construct KNN graph
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA") 

##--- Step3: Defining representative neighbourhoods on the KNN graph
milo <- makeNhoods(milo, prop = proportion, k = 30, d = 30, refined = TRUE, reduced_dims = "PCA")    #  ALL PBMC,T/NK cells: 0.05, Mono/DC,B cells: 0.1

# check average neighbourhood size: over 5 x N_samples is required 
p.NhSize = plotNhoodSizeHist(milo)
ggsave(file = paste0("./Milo/",cluster,"/Distribution_Neighbourhood.size_ALL.png"), plot = p.NhSize, width = 10, height = 7, dpi = 100)

##--- Step4: Counting cells in neighbourhoods
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="ID")

##--- Step5: Defining experimental design
design <- data.frame(colData(milo))[,c("ID", "Status", "Age", "Sex")]
design <- distinct(design)
rownames(design) <- design$ID

##--- Step6: Computing neighbourhood connectivity---##
milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")

##----- Step7: Testing -----##
# The last component of the formula or last column of the model matrix are by default the test variable.
da_results <- testNhoods(milo, design = ~ Age+Sex+Status, design.df = design)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

milo <- buildNhoodGraph(milo)
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1)     # alpha: FDR

ggsave(file = paste0("./Milo/",cluster,"/Milo_UMAP_Case_vs_HC.png"), plot = nh_graph_pl, width = 8.5, height = 7.5, dpi = 300)
saveRDS(milo, paste0("./Milo/",cluster,"/Milo_Case_vs_Ctrl.rds"))

##-----------------------------------------##
##         downstream analysis             ##
##-----------------------------------------##

#-- Assign a cell type label to each neighbourhood
da_results <- annotateNhoods(milo, da_results, coldata_col = "l3")    # l3: the most finest annotation

p.DA_celltype = plotDAbeeswarm(da_results, group.by = "l3")
ggsave(file = paste0("./Milo/",cluster,"/Celltype_fraction_DA.analysis_Case_vs_HC.png"), plot = p.DA_celltype, width = 7.5, height = 7.5, dpi = 300)
ggsave(file = paste0("./Milo/",cluster,"/Celltype_fraction_DA.analysis_Case_vs_HC.pdf"), plot = p.DA_celltype, width = 7.5, height = 7.5, dpi = 300)




