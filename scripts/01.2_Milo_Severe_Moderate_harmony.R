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



##--- Load Data & pick up COVID-19 samples
data <- readRDS(paste0("./Data/",input_file))            # SeuratObject

DefaultAssay(data) <-"RNA"
data <- subset(data, subset=Status==1)



##--- Add metadata
meta <- data@meta.data %>% as_tibble %>% select(ID)

pheno <- read_csv("/work22/home/redahiro/analysis/COVID-19_scRNAseq/Pheno_data/Phenotype_analysis_220702.csv")

meta <- meta %>% left_join(pheno, by="ID") %>% select(ID,days_after_onset,days_after_steroids,Age,Sex)

data[["Onset"]] <- meta$days_after_onset
data[["Steroids"]] <- meta$days_after_steroids
data[["Age"]] <- meta$Age
data[["Sex"]] <- meta$Sex


##--- SCE conversion
sce <- as.SingleCellExperiment(data)
altExps(sce) <- NULL                        

reducedDim(sce, "PCA", withDimnames=TRUE) <- data[['harmony']]@cell.embeddings      
reducedDim(sce, "UMAP", withDimnames=TRUE) <- data[['umap']]@cell.embeddings

rm(data)

##--------------------------------------##
##    Differential abundance testing    ##
##--------------------------------------##

##--- Step1: Create a Milo object
milo <- Milo(sce)


##--- Step2: Construct KNN graph
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")


##--- Step3: Defining representative neighbourhoods on the KNN graph
milo <- makeNhoods(milo, prop = proportion, k = 30, d = 30, refined = TRUE, reduced_dims = "PCA")   #  ALL PBMC,T/NK cells: 0.05, Mono/DC,B cells: 0.1

# check average neighbourhood size: over 5 x N_samples is required 
p.NhSize = plotNhoodSizeHist(milo)
ggsave(file = paste0("./Milo/",cluster,"/Distribution_Neighbourhood.size_Case.png"), plot = p.NhSize, width = 10, height = 7, dpi = 100)


##--- Step4: Counting cells in neighbourhoods
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="ID")

##--- Step5: Defining experimental design
design <- data.frame(colData(milo))[,c("ID", "Severity", "Age", "Sex", "Onset","Steroids")]
design <- distinct(design)
rownames(design) <- design$ID


##--- Step6: Computing neighbourhood connectivity---##
milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")


##----- Step7: Testing -----##
# The last component of the formula or last column of the model matrix are by default the test variable.
da_results <- testNhoods(milo, design = ~ Age+Sex+Onset+Steroids+Severity, design.df = design)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 


milo <- buildNhoodGraph(milo)
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1)     # alpha: FDR

ggsave(file = paste0("./Milo/",cluster,"/Milo_UMAP_Severe_Moderate.png"), plot = nh_graph_pl, width = 8.5, height = 7.5, dpi = 300)

nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=1)       # All nodes are colored by their log2 fold change
ggsave(file = paste0("./Milo/",cluster,"/Milo_UMAP_Severe_Moderate_FDR1.png"), plot = nh_graph_pl, width = 8.5, height = 7.5, dpi = 300)

saveRDS(milo, paste0("./Milo/",cluster,"/Milo_Severe_vs_Moderate.rds"))


##-----------------------------------------##
##         downstream analysis             ##
##-----------------------------------------##

# Beeswarm plot
p.DA_celltype = plotDAbeeswarm(da_results, group.by = "l3")                   # l3: the most finest annotation
ggsave(file = paste0("./Milo/",cluster,"/Celltype_fraction_DA.analysis_Severe_vs_Moderate.png"), plot = p.DA_celltype, width = 7.5, height = 7.5, dpi = 300)
ggsave(file = paste0("./Milo/",cluster,"/Celltype_fraction_DA.analysis_Severe_vs_Moderate.pdf"), plot = p.DA_celltype, width = 7.5, height = 7.5, dpi = 300)

# BoxPlot
da_results_t <- as_tibble(da_results) 

      bp <- ggplot(da_results_t, aes(x = l3, y = logFC, fill = l3)) +
               geom_hline(yintercept = 0, size=0.75, linetype = "dashed") +
               stat_boxplot(aes(fill=l3), position=position_dodge(.9), coef = 1.5) +
               # scale_fill_manual(values = c("#8549BA","#14D4F4","#F7AFFF","#00B22D","#99E2A6","#F36F21","#9BDEE0")) +
               # scale_y_continuous(breaks = seq(-10,10,by=1), limits=c(-3.25,3.25)) +
               labs(x = "", y = "log2FC") + 
               theme_classic() +
               coord_flip()
     
      bp = bp + theme(axis.text.x = element_text(size = 12),
           axis.text.y = element_text(size = 12),
           axis.title.x = element_text(size = 12),
           axis.title.y = element_text(size = 12),
           legend.position="none",
           axis.ticks.length = unit(0.3, "cm"))
    
      ggsave(file = paste0("./Milo/",cluster,"/milo_Severe.Moderate_boxplot.png"), plot = bp, width = 7, height = 5, dpi = 300)
      ggsave(file = paste0("./Milo/",cluster,"/milo_Severe.Moderate_boxplot.pdf"), plot = bp, width = 7, height = 5, dpi = 300)




